# run_optimizer.R
#
# Background entry point. Loops: choose a configuration -> sample a fresh trial -> evaluate
# under each CV scheme -> store -> checkpoint a report. Checkpointed after every evaluation, so
# it is fully resumable: kill it and re-launch to continue from the store.
#
#   nohup Rscript run_optimizer.R > logs/run.out 2>&1 &   # launch
#   tail -f logs/run.out                                  # watch
#   touch state/STOP                                      # stop cleanly

library(tidyverse)
here::i_am("run_optimizer.R")

source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

# One optimization step: pick config, pick trial, evaluate under each scheme, store. Returns a
# status record for logging.
#
# Failure handling (R/conditions.R): a `fatal` is not caught here and halts the loop; a
# `sample_failed` returns sampling_failed = TRUE with nothing stored, since it is not a
# config's fault; an `infeasible` is recorded as a failed eval and the step completes.
optimizer_step <- function(con, settings, conn = NULL) {
  choice <- choose_config(con, settings)
  # choose_trial, not sample_trial: a started trial must reach settings$trial_replication
  # distinct configurations before fresh trials are drawn, and the chosen config's own trials
  # are excluded so a replication step adds a trial rather than repeating one (R/optimizer.R).
  trial  <- tryCatch(choose_trial(con, settings, conn, cfg_hash = config_hash(choice$cfg)),
                     optimizer_sample_failed = function(e) {
                       message("  trial sampling: ", conditionMessage(e)); NULL })
  if (is.null(trial)) return(list(source = choice$source, sampling_failed = TRUE))

  # The optimizer targets ONE scheme per run (settings$optimize_scheme); CV0 and
  # CV00 are distinct tasks, optimized in separate runs against the shared store.
  scheme <- settings$optimize_scheme
  ev <- evaluate_config_on_trial(choice$cfg, trial, scheme, settings, conn)
  # Record the trial's target-domain attributes alongside the score so a later run
  # can train the surrogate on only its own domain (and scheme) slice of this store.
  # (Simulated trials have no program/location/year -> NA, which is correct: sim
  # runs use no target_domain, so nothing is filtered out.)
  store_eval(con, choice$cfg, trial$id, scheme, ev$score, ev$n_test,
             ev$status, ev$reason, ev$detail %||% NA_character_, ev$seconds,
             study_name    = trial$study_name %||% NA_character_,
             program_name  = trial$program    %||% NA_character_,
             location_name = trial$location   %||% NA_character_,
             year          = trial$year       %||% NA_integer_,
             # What this evaluation cost in memory, and under which marker-density budget --
             # the two numbers that decide how many workers this machine can carry.
             peak_rss_mb   = ev$peak_rss_mb %||% NA_real_,
             peak_r_mb     = ev$peak_r_mb   %||% NA_real_,
             rss_mb        = ev$rss_mb      %||% NA_real_,
             worker        = settings$worker_id %||% NA_character_,
             dosage_budget = settings$dosage_budget_bytes %||% NA_real_,
             # How em_combine derived each partial's EM weight. Like dosage_budget, this is
             # not a config parameter, so without it rows from before and after the
             # 2026-07-31 switch would be averaged together (EM_COMBINE_COMPARISON.md item 1).
             em_df_method  = "effective_n",
             build         = settings$build %||% OPTIMIZER_BUILD)
  # scores/statuses kept as length-1 vectors so the main-loop logging is unchanged.
  list(source = choice$source, trial = trial$id,
       scores = ev$score, statuses = ev$status, ei = choice$ei %||% NA_real_,
       # Prefer the true RSS peak for the log; the heap peak understates it (R/memory.R).
       # Must test is.finite, not %||%: off Linux peak_rss_mb is NA rather than NULL, and
       # %||% only falls back on NULL -- which logged "peak=NA" instead of the heap figure.
       peak_mb = local({ p <- ev$peak_rss_mb %||% NA_real_
                         if (is.finite(p)) p else ev$peak_r_mb %||% NA_real_ }))
}

# The main loop. Reusable from tests with a small max_iters.
run_optimizer <- function(settings = optimizer_settings(), conn = NULL) {
  message(sprintf("optimizer build %s | scheme %s | %s mode",
                  settings$build %||% OPTIMIZER_BUILD, settings$optimize_scheme,
                  if (isTRUE(settings$simulate)) "SIMULATE" else "real"))
  dir.create(dirname(settings$db_path), showWarnings = FALSE, recursive = TRUE)
  dir.create(settings$log_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(settings$cache_dir, showWarnings = FALSE, recursive = TRUE)

  # Warm the work cache from its durable backup. RESTORING stays leader-only: N workers
  # pulling one tree at startup fight for the disk while the others read it. The others wait
  # on cache_ready_file, or they would re-download what the rsync is about to deliver.
  # Flushing back is every worker's job -- LESSONS #25.
  leader <- isTRUE(settings$is_leader %||% TRUE)
  ready_file <- settings$cache_ready_file
  if (leader) {
    if (!is.null(ready_file)) unlink(ready_file)       # stale flag from a previous run
    restore_cache_from_backup(settings)
    if (!is.null(ready_file)) {
      dir.create(dirname(ready_file), showWarnings = FALSE, recursive = TRUE)
      file.create(ready_file)
    }
  } else if (!is.null(ready_file) && !is.null(settings$cache_backup_dir) &&
             nzchar(settings$cache_backup_dir)) {
    # Bounded wait: if the leader is absent or its restore failed, start anyway rather than
    # stall the whole run. (No backup configured means there is nothing to wait for.)
    waited <- 0; limit <- 60 * 30
    while (!file.exists(ready_file) && waited < limit) { Sys.sleep(5); waited <- waited + 5 }
    if (file.exists(ready_file))
      message(sprintf("worker %s: cache restored by the leader after %.0f s -- starting",
                      settings$worker_id %||% "?", waited))
    else
      message(sprintf("worker %s: no cache-ready signal after %.0f min -- starting anyway",
                      settings$worker_id %||% "?", limit / 60))
  }

  # Flush the cache on the way out -- EVERY worker, not just the leader. A leader-only exit
  # hook does nothing when worker 1 is mid-evaluation at the wall clock (SIGKILL: on.exit never
  # runs) or was OOM-killed earlier, which is precisely when a final flush is wanted. The short
  # floor guarantees a last capture even if a periodic sync just ran, while leaving the stamp
  # to cap the end-of-job stampede at about one extra rsync -- and at that moment the scheduler
  # is about to kill everything, so a queue of full tree-walks is the wrong thing to start.
  on.exit(sync_cache_to_backup(settings, min_age_minutes = 2), add = TRUE)

  # Real mode needs a live BrAPI connection (made once, reused), logged in from the
  # T3_USERNAME/T3_PASSWORD environment credentials. Simulate mode never touches the network.
  if (!settings$simulate && is.null(conn)) {
    conn <- t3_connect(settings)
  }

  con <- open_store(settings$db_path)
  on.exit(close_store(con), add = TRUE)
  # A build that invalidated earlier rows starts with less history than the store suggests.
  # Say so once, loudly, rather than let a silently-empty surrogate look like a fresh start.
  local({
    all_n  <- nrow(read_evals(con))
    kept_n <- nrow(filter_evals_to_build(read_evals(con), settings$build %||% OPTIMIZER_BUILD))
    if (all_n > 0)
      message(sprintf("store: %d rows, %d usable by build %s (%d retired by BUILD_CHANGES)",
                      all_n, kept_n, settings$build %||% OPTIMIZER_BUILD, all_n - kept_n))
  })

  start <- Sys.time()
  start_n <- n_evals(con)
  iter <- 0
  consec_sample_fail <- 0L
  fatal_hit <- FALSE
  message(sprintf("[%s] optimizer start (simulate=%s, worker=%s%s); %d evaluations already in store",
                  format(start), settings$simulate, settings$worker_id %||% "1",
                  if (leader) ", leader" else "", start_n))

  repeat {
    if (file.exists(settings$stop_file)) {
      message("STOP file present -> halting cleanly."); break
    }
    if (iter >= settings$max_iters) { message("iteration budget reached."); break }

    # A fatal condition halts the run; a sample_failed bubbles up as
    # step$sampling_failed; an infeasible (trial, config) was already recorded by
    # evaluate.R and the step returns normally. Unexpected errors are logged and
    # the loop continues (a stray bug should not kill a multi-day run).
    step <- tryCatch(optimizer_step(con, settings, conn),
                     optimizer_fatal = function(e) {
                       message("FATAL: ", conditionMessage(e), " -> halting."); fatal_hit <<- TRUE; NULL },
                     error = function(e) { message("step error: ", conditionMessage(e)); NULL })
    if (fatal_hit) break
    iter <- iter + 1

    if (!is.null(step) && isTRUE(step$sampling_failed)) {
      consec_sample_fail <- consec_sample_fail + 1L
      message(sprintf("[%s] iter %d  trial sampling failed (%d consecutive)",
                      format(Sys.time()), iter, consec_sample_fail))
      if (consec_sample_fail >= settings$max_sample_fail) {
        message("too many consecutive trial-sampling failures (", consec_sample_fail,
                ") -> halting; check target_domain / min_trial_acc / network."); break
      }
    } else {
      consec_sample_fail <- 0L
      if (!is.null(step)) {
        message(sprintf("[%s] iter %d  src=%-16s trial=%s  scores=%s  status=%s  peak=%s",
                        format(Sys.time()), iter, step$source, step$trial,
                        paste(sprintf("%.3f", step$scores), collapse = "/"),
                        paste(step$statuses, collapse = "/"),
                        fmt_mb(step$peak_mb)))
      }
    }
    # EVERY worker writes the report (write_report renders to a temp file and renames, so
    # concurrent writes cannot interleave). Deliberately NOT leader-only: the report
    # summarizes the whole store, and a worker can only write it between evaluations, so
    # restricting it to worker 1 meant the file went stale for as long as worker 1's current
    # evaluation ran -- hours, while the other seven kept adding rows nobody could see.
    if (iter %% settings$checkpoint_every == 0) {
      tryCatch(write_report(con, settings),
               error = function(e) message("report error: ", conditionMessage(e)))
    }
    # Cache backup (additive rsync). Bounds what an abrupt kill loses to roughly one interval
    # of freshly-downloaded cache. Called unconditionally: sync_cache_to_backup() decides
    # whether it is due, from a stamp file every worker can see, so no worker's long evaluation
    # can hold up everyone else's backup. Not leader-gated -- see LESSONS #25.
    sync_cache_to_backup(settings)
    # Copy the store to durable storage: db_path is on local disk, and the store is the one
    # file whose loss costs real work. Any worker may do it, throttled on the backup's own
    # mtime -- see should_backup_now() in R/store.R.
    if (should_backup_now(settings)) {
      # backup_store reports its own failure; note the consequence here so a run whose backups
      # are all failing says so in the log rather than only at the moment /workdir is wiped.
      if (!backup_store(con, settings$db_backup_path))
        message("  the store is NOT being backed up -- a loss of ", dirname(settings$db_path),
                " would lose this run")
    }
  }

  write_report(con, settings)          # every worker, as in the loop
  # The final backup is the one that matters most -- say plainly whether it happened. Every
  # worker, not just the leader: if worker 1 is the one that gets OOM-killed, a leader-only
  # final backup means no worker takes one at all.
  if (!is.null(settings$db_backup_path)) {
    if (backup_store(con, settings$db_backup_path))
      message("store backed up to ", settings$db_backup_path)
    else
      message("FINAL STORE BACKUP FAILED -- copy ", settings$db_path, " off this disk by hand")
  }
  message(sprintf("[%s] optimizer stop; %d evaluations in store (this run: %d)",
                  format(Sys.time()), n_evals(con), n_evals(con) - start_n))
  invisible(con)
}

# Auto-run only when invoked directly (Rscript run_optimizer.R), not when this
# file is source()d by a test or another script. At the top level of a script
# run, sys.nframe() == 0; inside a source() call it is > 0.
if (sys.nframe() == 0L && !interactive()) {
  run_optimizer(optimizer_settings())
}
