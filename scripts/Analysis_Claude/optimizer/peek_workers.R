# peek_workers.R
#
# Read-only answer to "what is every worker doing right now, and has anything died?".
#
# Three signals, because no one of them is sufficient:
#
#   claims + pid liveness   catches a worker killed MID-EVALUATION -- its claim outlives it,
#                           and a claim whose owner is gone is the record of what it died on.
#                           Misses a worker that died BETWEEN evaluations, holding no claim.
#   process count           how many are alive now. Says nothing about which ids.
#   worker-log mtimes       catches a worker gone quiet. Silence is not proof of death: a
#                           single evaluation has run for 33 hours (LESSONS #28).
#
# WORKER 1 IS NOT SPECIAL ANY MORE. Its only exclusive jobs are restoring the cache and
# clearing stale claims at startup; backups are every worker's job (LESSONS #25). After
# startup it is an ordinary worker, so the question is which workers are alive, not whether
# that one is.
#
# WHERE THIS WORKS. Claims live in the LIVE store only, which on a cluster is node-local
# scratch belonging to the job that wrote it. So this must run on the node the workers are on
# -- from a separate allocation there is nothing to see, and the backup's copy of the table is
# a stale VACUUM INTO artifact rather than the truth. It says so rather than reporting an
# empty run.
#
# Run from THIS directory (R reads .Renviron and here() from the working directory):
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript peek_workers.R

if (!file.exists("peek_workers.R"))
  stop("run this FROM the optimizer directory:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript peek_workers.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("peek_workers.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

if (.on_login_node()) stop(.login_node_message("peek_workers.R"), call. = FALSE)
s <- optimizer_settings()

# NOT resolve_read_store(): that falls back to the backup, whose claims table is whatever
# happened to be in flight when the last VACUUM INTO ran. Reporting those as live work would
# be worse than reporting nothing.
if (!file.exists(s$db_path)) {
  stop("no live store at ", s$db_path, "\n",
       "  In-flight work is recorded only there, so this has to run ON THE NODE the workers\n",
       "  are on. Under SLURM db_path is node-local scratch belonging to that job, and a\n",
       "  separate allocation cannot see it. Attach to the running job instead:\n",
       "    squeue -u $USER                              # find the job id\n",
       "    srun --overlap --jobid=<jid> --pty bash      # a shell ON that node\n",
       "    cd ", here::here(), " && ./container/run_in_container.sh exec peek_workers.R\n",
       "  (For finished work, which IS in the backup, use peek_backup.R or report_timing.R.)",
       call. = FALSE)
}

con <- open_store(s$db_path)
on.exit(close_store(con), add = TRUE)
cl <- claims_snapshot(con)
n_ev <- n_evals(con)

cat("store  : ", s$db_path, " (", n_ev, " evals)\n", sep = "")
cat("host   : ", Sys.info()[["nodename"]], "\n", sep = "")
cat("scheme : ", s$optimize_scheme, "\n\n", sep = "")

# --- 1. in-flight work, from the claims table -------------------------------
if (!nrow(cl)) {
  cat("IN FLIGHT: nothing.\n")
  cat("  Every worker is between evaluations, or none is running. The process count below\n")
  cat("  distinguishes those two.\n\n")
} else {
  cat("IN FLIGHT (", nrow(cl), " evaluation(s)):\n\n", sep = "")
  state <- ifelse(is.na(cl$alive), "other-host", ifelse(cl$alive, "alive", "DEAD"))
  out <- data.frame(
    worker  = cl$worker,
    pid     = cl$pid,
    trial   = substr(as.character(cl$trial_id), 1, 18),
    config  = substr(cl$config_hash, 1, 8),
    held_h  = sprintf("%.1f", cl$hours_held),
    state   = state,
    note    = ifelse(state == "DEAD", "<- died holding this", ""))
  print(out, row.names = FALSE)
  cat("\n")
  dead <- sum(state == "DEAD")
  if (dead > 0) {
    cat("** ", dead, " claim(s) held by a process that is gone.**\n", sep = "")
    cat("  That worker was killed mid-evaluation -- OOM is the usual cause. The row above is\n")
    cat("  what it was working on. For how much memory it was using at the time, find its pid\n")
    cat("  in the pids_rss column of ", file.path(s$log_dir, "memory_*.tsv"), ".\n", sep = "")
    cat("  Nothing needs cleaning up: the next claim on that pair reclaims it (dead owner).\n\n")
  }
}

# --- 1b. workers that died earlier ------------------------------------------
# The claims table is self-healing: choose_trial purges a dead owner's claim within seconds,
# so the live table almost never catches a death. claims_reaped is the copy taken before that
# purge, and is what makes "has anything broken since I launched it" answerable at all.
rp <- reaped_claims(con)
if (nrow(rp)) {
  cat("\nWORKERS LOST MID-EVALUATION (", nrow(rp), " since this store was created):\n\n", sep = "")
  print(data.frame(
    worker = rp$worker, pid = rp$pid,
    trial  = substr(as.character(rp$trial_id), 1, 18),
    config = substr(rp$config_hash, 1, 8),
    started = rp$ts, lost_at = rp$reaped_ts),
    row.names = FALSE)
  cat("\n  Each row is a worker that held a claim and then vanished -- OOM is the usual cause.\n")
  cat("  The evaluation was lost, not corrupted: no row was written, and the pair is free to\n")
  cat("  be re-run. For memory at the time, find the pid in ",
      file.path(s$log_dir, "memory_*.tsv"), ".\n", sep = "")
  cat("  Some of these are ordinary shutdowns: a job killed at its wall clock reaps whatever\n")
  cat("  was in flight when the next leader starts. Compare lost_at against your job history.\n")
} else {
  cat("\nWORKERS LOST MID-EVALUATION: none since this store was created.\n")
}

# --- 2. how many worker processes are actually running ----------------------
# Same idea as run_workers.sh:99, but matched in R rather than piped through a shell: the
# `[r]un_optimizer` bracket trick has to survive two levels of quoting to reach grep, and does
# not reliably do so.
n_proc <- tryCatch(
  sum(grepl("run_optimizer[.]R", system2("ps", c("-e", "-o", "args="), stdout = TRUE))),
  error = function(e) NA_integer_)
cat("PROCESSES: ", if (is.na(n_proc)) "?" else n_proc,
    " run_optimizer.R running on this host\n", sep = "")

# --- 3. per-worker log activity ---------------------------------------------
# Catches the death the claims table cannot see: a worker that exited between evaluations.
logs <- list.files(s$log_dir, pattern = "^run_w[0-9]+[.]out$", full.names = TRUE)
if (!length(logs)) {
  cat("LOGS     : none in ", s$log_dir, "\n", sep = "")
} else {
  age <- as.numeric(difftime(Sys.time(), file.mtime(logs), units = "mins"))
  ord <- order(as.integer(sub(".*run_w([0-9]+)[.]out$", "\\1", logs)))
  cat("\nLOGS in ", s$log_dir, ":\n\n", sep = "")
  print(data.frame(
    worker    = sub(".*run_w([0-9]+)[.]out$", "\\1", logs[ord]),
    quiet_min = sprintf("%.0f", age[ord]),
    holding   = ifelse(sub(".*run_w([0-9]+)[.]out$", "\\1", logs[ord]) %in% as.character(cl$worker),
                       "yes", "no")), row.names = FALSE)
  cat("\n")
  # Quiet AND holding nothing is the suspicious combination: a worker mid-evaluation is quiet
  # for legitimate hours, but one holding no claim should be logging an iteration line.
  susp <- age[ord] > 30 &
    !(sub(".*run_w([0-9]+)[.]out$", "\\1", logs[ord]) %in% as.character(cl$worker))
  if (any(susp))
    cat("** worker(s) ", paste(sub(".*run_w([0-9]+)[.]out$", "\\1", logs[ord])[susp], collapse = ", "),
        " have logged nothing for >30 min and hold no claim -- look for an exit at the end of\n",
        "   their log. A worker mid-evaluation is quiet legitimately; one that is not, is not.\n",
        sep = "")
}
