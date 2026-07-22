# run_optimizer.R
#
# Background entry point. Sources the subsystem, then loops:
#   choose a configuration -> sample a fresh random trial -> evaluate under each
#   CV scheme -> store -> checkpoint a report.
# It checkpoints after every evaluation, so the process is fully resumable: kill
# it any time and re-launch to continue from the SQLite store. Stop it cleanly by
# creating the stop-file (state/STOP) or hitting the iteration / wall-clock
# budget in settings.R.
#
# Launch in the background from this directory:
#   nohup Rscript run_optimizer.R > logs/run.out 2>&1 &
# Watch it:
#   tail -f logs/run.out        # live log
#   cat state/report.md         # current best pipeline + learning curve
# Stop it:
#   touch state/STOP

library(tidyverse)
here::i_am("run_optimizer.R")

source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

# One optimization step: pick config, pick trial, evaluate under each scheme,
# store. Returns a short status record for logging.
#
# Failure handling (see R/conditions.R):
#   * a fatal condition (bad settings / unimplemented method) is NOT caught here
#     -- it propagates to run_optimizer()'s loop, which halts.
#   * a sample_failed condition (no usable trial) is caught and reported back as
#     sampling_failed = TRUE; no eval is stored (it is not a config's fault). The
#     loop tolerates a few in a row, then stops.
#   * an infeasible (trial, config, scheme) is recorded by evaluate_config_on_trial
#     as a failed eval (the failure log) and the step completes normally; the NEXT
#     step samples a new trial and chooses a new configuration.
optimizer_step <- function(con, settings, conn = NULL) {
  choice <- choose_config(con, settings)
  trial  <- tryCatch(sample_trial(settings, conn),
                     optimizer_sample_failed = function(e) {
                       message("  trial sampling: ", conditionMessage(e)); NULL })
  if (is.null(trial)) return(list(source = choice$source, sampling_failed = TRUE))

  recs <- purrr::map(settings$schemes, function(scheme) {
    ev <- evaluate_config_on_trial(choice$cfg, trial, scheme, settings, conn)
    store_eval(con, choice$cfg, trial$id, scheme, ev$score, ev$n_test,
               ev$status, ev$reason, ev$detail %||% NA_character_, ev$seconds)
    ev
  },
  .progress = "Evaluate config on trials")
  scores   <- vapply(recs, function(r) r$score, numeric(1))
  statuses <- vapply(recs, function(r) r$status, character(1))
  list(source = choice$source, trial = trial$id,
       scores = scores, statuses = statuses, ei = choice$ei %||% NA_real_)
}

# The main loop. Reusable from tests with a small max_iters.
run_optimizer <- function(settings = optimizer_settings(), conn = NULL) {
  dir.create(dirname(settings$db_path), showWarnings = FALSE, recursive = TRUE)
  dir.create(settings$log_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(settings$cache_dir, showWarnings = FALSE, recursive = TRUE)

  # Real mode needs a live BrAPI connection (made once, reused), logged in from the
  # T3_USERNAME/T3_PASSWORD environment credentials. Simulate mode never touches the network.
  if (!settings$simulate && is.null(conn)) {
    conn <- t3_connect(settings)
  }

  # Bug oracle: before spending hours, confirm the pipeline can still predict
  # trials we KNOW are feasible. A canary failure means the code is hiding real
  # data (a name/parse/column bug), so the "infeasible" verdicts cannot be
  # trusted. We warn loudly but do not halt -- the choice to proceed is yours.
  if (!settings$simulate && !is.null(settings$canary_trials)) {
    tryCatch(check_canaries(settings, conn),
             error = function(e) message("canary check error: ", conditionMessage(e)))
  }

  con <- open_store(settings$db_path)
  on.exit(close_store(con), add = TRUE)

  start <- Sys.time()
  start_n <- n_evals(con)
  iter <- 0
  consec_sample_fail <- 0L
  fatal_hit <- FALSE
  message(sprintf("[%s] optimizer start (simulate=%s); %d evaluations already in store",
                  format(start), settings$simulate, start_n))

  repeat {
    if (file.exists(settings$stop_file)) {
      message("STOP file present -> halting cleanly."); break
    }
    if (iter >= settings$max_iters) { message("iteration budget reached."); break }
    if (as.numeric(difftime(Sys.time(), start, units = "hours")) >= settings$max_hours) {
      message("wall-clock budget reached."); break
    }

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
        message(sprintf("[%s] iter %d  src=%-16s trial=%s  scores=%s  status=%s",
                        format(Sys.time()), iter, step$source, step$trial,
                        paste(sprintf("%.3f", step$scores), collapse = "/"),
                        paste(step$statuses, collapse = "/")))
      }
    }
    if (iter %% settings$checkpoint_every == 0) {
      tryCatch(write_report(con, settings),
               error = function(e) message("report error: ", conditionMessage(e)))
    }
  }

  write_report(con, settings)
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
