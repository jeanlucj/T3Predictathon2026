# test_concurrency.R
#
# Several workers may run against ONE evals.sqlite and ONE cache (see run_workers.sh). The
# properties that makes safe are not visible from single-process testing, so they are
# exercised here with REAL concurrent processes:
#
#   1. WAL + busy_timeout: N processes inserting at once lose no rows and raise no errors.
#   2. The migration guard: N processes opening a PRE-MIGRATION store at the same moment all
#      survive the ALTER TABLE race (only one can win it).
#   3. The cache lock: N processes wanting one uncached project produce exactly ONE download.
#   4. Atomic cache writes: a reader never sees a partial file.
#   5. The store backup happens even when no worker is the leader, and its interval is shared
#      across processes rather than honoured N times over.
#
# Run: Rscript tests/test_concurrency.R

library(tidyverse)
here::i_am("tests/test_concurrency.R")
source(here::here("settings.R"))
invisible(lapply(list.files(here::here("R"), "[.]R$", full.names = TRUE), source))

ok <- 0L; fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) { ok <<- ok + 1L }
  else { fail <<- fail + 1L; cat("  FAIL:", msg, "\n") }
}

tmp <- file.path(tempdir(), paste0("conc_", Sys.getpid()))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

# Launch `n` background Rscript processes running `body`, and wait for all of them. Each gets
# its worker number as a trailing argument. Returns their exit statuses.
spawn <- function(n, body, ...) {
  f <- file.path(tmp, "worker_body.R")
  writeLines(c(
    'suppressMessages(library(tidyverse))',
    sprintf('here::i_am("%s")', "run_optimizer.R"),
    sprintf('setwd("%s")', here::here()),
    'source(here::here("settings.R"))',
    'invisible(lapply(list.files(here::here("R"), "[.]R$", full.names = TRUE), source))',
    'ARG <- commandArgs(trailingOnly = TRUE)',
    sprintf('TMP <- "%s"', tmp),
    body), f)
  pids <- lapply(seq_len(n), function(i)
    system2("Rscript", c(f, as.character(i), ...), wait = FALSE,
            stdout = file.path(tmp, paste0("w", i, ".out")),
            stderr = file.path(tmp, paste0("w", i, ".err"))))
  # No portable way to wait on a detached pid from R, so poll for the sentinel each worker
  # writes on exit. The budget is deliberately generous: four concurrent R startups each
  # sourcing tidyverse and the whole R/ tree can take a while on a loaded machine, and a
  # timeout that fires on slowness rather than on a defect trains you to ignore this test.
  deadline <- Sys.time() + 600
  repeat {
    done <- length(list.files(tmp, pattern = "^done_"))
    if (done >= n || Sys.time() > deadline) break
    Sys.sleep(0.5)
  }
  done <- length(list.files(tmp, pattern = "^done_"))
  # On a shortfall, say WHICH workers never finished and show their stderr, so "too slow" is
  # distinguishable from "crashed" without re-running.
  if (done < n) {
    missing <- setdiff(seq_len(n), as.integer(sub("^done_", "",
                                                 list.files(tmp, pattern = "^done_"))))
    cat("  worker(s) did not finish within", round(as.numeric(difftime(Sys.time(), deadline,
        units = "secs")) + 600), "s:", paste(missing, collapse = ", "), "\n")
    for (i in missing) {
      ef <- file.path(tmp, paste0("w", i, ".err"))
      if (file.exists(ef)) cat("   w", i, " stderr: ",
                               paste(utils::tail(readLines(ef, warn = FALSE), 3),
                                     collapse = " | "), "\n", sep = "")
    }
  }
  done
}

cat("1. concurrent writers -----------------------------------------------\n")
db <- file.path(tmp, "evals.sqlite")
N_W <- 4L; N_ROWS <- 50L
n_done <- spawn(N_W, c(
  'con <- open_store(file.path(TMP, "evals.sqlite"))',
  'set.seed(as.integer(ARG[1]))',
  sprintf('for (i in seq_len(%d)) {', N_ROWS),
  '  store_eval(con, sample_config(), paste0("t", i), "CV0", runif(1), 100L, "ok",',
  '             worker = ARG[1], peak_r_mb = 100 * as.integer(ARG[1]),',
  '             peak_rss_mb = 300 * as.integer(ARG[1]), dosage_budget = 2e9)',
  '}',
  'close_store(con)',
  'file.create(file.path(TMP, paste0("done_", ARG[1])))'))
check(n_done == N_W, sprintf("all %d writers finished (got %d)", N_W, n_done))

errs <- unlist(lapply(list.files(tmp, "^w[0-9]+[.]err$", full.names = TRUE),
                      function(f) grep("locked|BUSY|Error", readLines(f, warn = FALSE), value = TRUE)))
check(length(errs) == 0, paste("no lock/errors in worker stderr:", paste(head(errs, 3), collapse = " | ")))

con <- open_store(db)
n <- n_evals(con)
check(n == N_W * N_ROWS, sprintf("all rows present: expected %d, got %d", N_W * N_ROWS, n))
check(identical(tolower(DBI::dbGetQuery(con, "PRAGMA journal_mode")$journal_mode), "wal"),
      "store is in WAL mode")
e <- read_evals(con)
check(length(unique(e$worker)) == N_W, "every worker's rows are attributed to it")
check(all(is.finite(e$peak_r_mb)) && all(is.finite(e$peak_rss_mb)),
      "both memory columns round-trip (heap peak and true RSS peak)")
close_store(con)

cat("2. concurrent migration of a pre-migration store ---------------------\n")
# A store with the ORIGINAL schema: no detail/domain/memory columns. Every worker will want
# to add all of them at once.
old_db <- file.path(tmp, "old.sqlite")
oc <- DBI::dbConnect(RSQLite::SQLite(), old_db)
invisible(DBI::dbExecute(oc, "CREATE TABLE evals (id INTEGER PRIMARY KEY AUTOINCREMENT,
  config_hash TEXT, config_json TEXT, trial_id TEXT, scheme TEXT, score REAL,
  n_test INTEGER, status TEXT, reason TEXT, seconds REAL, ts TEXT)"))
DBI::dbDisconnect(oc)
unlink(list.files(tmp, "^done_", full.names = TRUE))

n_done <- spawn(4L, c(
  # No stagger: the point is to have them collide inside open_store.
  'con <- open_store(file.path(TMP, "old.sqlite"))',
  'store_eval(con, sample_config(), "t1", "CV0", 0.5, 10L, "ok", worker = ARG[1])',
  'close_store(con)',
  'file.create(file.path(TMP, paste0("done_", ARG[1])))'))
check(n_done == 4L, sprintf("all 4 migrating workers finished (got %d)", n_done))
merrs <- unlist(lapply(list.files(tmp, "^w[0-9]+[.]err$", full.names = TRUE),
                       function(f) grep("Error", readLines(f, warn = FALSE), value = TRUE)))
check(length(merrs) == 0, paste("no worker died in the migration race:",
                                paste(head(merrs, 3), collapse = " | ")))
oc <- open_store(old_db)
have <- DBI::dbListFields(oc, "evals")
check(all(c("detail", "study_name", "peak_r_mb", "rss_mb", "worker", "dosage_budget",
            "peak_rss_mb") %in% have),
      "every migrated column exists exactly once")
check(DBI::dbGetQuery(oc, "SELECT COUNT(*) n FROM evals")$n == 4, "all 4 rows landed")
DBI::dbDisconnect(oc)

cat("3. cache lock: one download for N workers ---------------------------\n")
# get_project_dosage's expensive branch is stubbed with a marker file, so the test measures
# the LOCK rather than the network: each worker appends a line when it does "the work". With
# the lock held, only the first should get there; the rest wait and read the result.
unlink(list.files(tmp, "^done_", full.names = TRUE))
cache <- file.path(tmp, "cache"); dir.create(cache, showWarnings = FALSE)
n_done <- spawn(4L, c(
  'st <- list(cache_dir = TMP_CACHE <- file.path(TMP, "cache"), worker_id = ARG[1],',
  '           lock_wait_minutes = 2, lock_stale_minutes = 90)',
  'result_file <- file.path(TMP_CACHE, "result.rds")',
  'v <- .with_cache_lock(st, "proj_999",',
  '  ready    = function() file.exists(result_file),',
  '  on_ready = function() readRDS(result_file),',
  '  expr     = {',
  '    cat(ARG[1], "\n", file = file.path(TMP, "downloads.log"), append = TRUE)',
  '    Sys.sleep(6)                     # a "download" long enough for the others to queue',
  '    saveRDS("payload", result_file)',
  '    "payload"',
  '  })',
  'stopifnot(identical(v, "payload"))',
  'file.create(file.path(TMP, paste0("done_", ARG[1])))'))
check(n_done == 4L, sprintf("all 4 lock workers finished (got %d)", n_done))
dl <- if (file.exists(file.path(tmp, "downloads.log")))
  length(readLines(file.path(tmp, "downloads.log"))) else 0L
check(dl == 1L, sprintf("exactly one worker did the download (got %d)", dl))
check(!dir.exists(file.path(cache, "locks", "proj_999.lock")), "the lock was released")

cat("4. atomic cache writes ----------------------------------------------\n")
st <- list(cache_dir = file.path(tmp, "cache2"))
p <- .cache_save(st, "dosage", "1_sz9", matrix(1:20, 4))
check(file.exists(p), "cache file written")
check(length(list.files(dirname(p), pattern = "^[.]tmp")) == 0, "no temporary file left behind")
check(identical(readRDS(p), matrix(1:20, 4)), "value round-trips")
# A partially-written file must never be visible to the glob .find_densest_dosage uses.
check(length(list.files(dirname(p), pattern = "^dosage_")) == 1,
      "exactly one dosage file matches the glob")

cat("5. store backup does not depend on the leader ------------------------\n")
# Both properties are checked with NO leader among the workers: a backup still happens, and the
# mtime throttle is shared across processes rather than honoured once per process. LESSONS #25.
unlink(list.files(tmp, "^done_", full.names = TRUE))
bk_db <- file.path(tmp, "bk.sqlite")
bk_dest <- file.path(tmp, "bk_backup.sqlite")
n_done <- spawn(3L, c(
  'con <- open_store(file.path(TMP, "bk.sqlite"))',
  'dest <- file.path(TMP, "bk_backup.sqlite")',
  'cfg <- sample_config()',
  'made <- 0L',
  'for (i in 1:4) {',
  '  store_eval(con, cfg, trial_id = paste0(ARG[1], "_", i), scheme = "CV0",',
  '             score = runif(1), n_test = 30L, status = "ok")',
  # db_backup_minutes = 0 here would disable it; 1/60 min = 1 s keeps the test quick while
  # still exercising a throttle that is genuinely shorter than the loop.
  # is_leader = FALSE on EVERY worker: the gate must not consult it. 1/60 min = 1 s keeps the
  # test quick while still exercising a throttle genuinely shorter than the loop.
  '  st <- list(db_backup_path = dest, db_backup_minutes = 1/60, is_leader = FALSE)',
  '  if (should_backup_now(st) && backup_store(con, dest)) made <- made + 1L',
  '  Sys.sleep(0.3)',
  '}',
  'close_store(con)',
  # One file per worker, not a shared append: three processes appending to one log interleave
  # and tear a line, which made this test intermittently read NA.
  'writeLines(as.character(made), file.path(TMP, paste0("backups_", ARG[1], ".txt")))',
  'file.create(file.path(TMP, paste0("done_", ARG[1])))'))
check(n_done == 3L, sprintf("all 3 backup workers finished (got %d)", n_done))
check(file.exists(bk_dest), "a backup exists though NO worker was the leader")
# Shared throttle: 12 iterations over ~3.6 s at a 1 s interval. A per-process timestamp would
# let all three fire every second (up to ~12); a shared one caps it near the elapsed seconds.
bf <- list.files(tmp, pattern = "^backups_", full.names = TRUE)
nb <- if (length(bf)) sum(vapply(bf, function(p) as.integer(readLines(p)[1]), integer(1))) else 0L
check(nb >= 1L, sprintf("at least one backup was taken (got %d)", nb))
check(nb <= 6L, sprintf("the throttle is shared, not per-worker (got %d backups, want <= 6)", nb))
# Whatever the race, every renamed candidate must be a complete, readable database.
bc <- DBI::dbConnect(RSQLite::SQLite(), bk_dest)
nrows <- tryCatch(DBI::dbGetQuery(bc, "SELECT COUNT(*) n FROM evals")$n, error = function(e) -1L)
DBI::dbDisconnect(bc)
check(nrows > 0, sprintf("the backup is a valid store with rows (got %d)", nrows))
check(!length(list.files(tmp, pattern = "^bk_backup[.]sqlite[.]tmp")),
      "no half-written backup temp file was left behind")
# The gate itself, on its own. is_leader is varied to pin down that it is NOT consulted --
# this is the specific regression that let the backup go 24 h stale.
gate <- function(dest, mins, leader)
  should_backup_now(list(db_backup_path = dest, db_backup_minutes = mins, is_leader = leader))
missing_bk <- file.path(tmp, "no_such_backup.sqlite")
check(is.infinite(backup_age_minutes(missing_bk)), "a missing backup is infinitely old")
check(is.infinite(backup_age_minutes(NULL)), "an unconfigured backup path is infinitely old")
check(gate(missing_bk, 30, leader = FALSE), "a NON-leader backs up when none exists")
check(gate(missing_bk, 30, leader = TRUE),  "so does the leader")
check(!gate(bk_dest, 30, leader = TRUE),  "a fresh backup is not redone, even by the leader")
check(!gate(bk_dest, 30, leader = FALSE), "nor by a non-leader")
check(!gate(NULL, 30, FALSE), "no backup path configured (local mode) means no backup")
check(!gate(missing_bk, 0, FALSE), "db_backup_minutes = 0 disables backups")
Sys.setFileTime(bk_dest, Sys.time() - 3600)     # one hour old, interval 30 min
check(gate(bk_dest, 30, leader = FALSE), "a stale backup is refreshed by a non-leader")

cat(sprintf("\n%d checks passed, %d failed\n", ok, fail))
quit(status = if (fail > 0) 1L else 0L)
