# test_subtasks.R  --  Tier 1: deterministic core of the pipeline.
#
# Each test states its ORACLE (the expected answer, knowable without reading the
# implementation) in a comment, so reviewing the test is faster than auditing the
# code. See README_DEVELOP.md (test suite). All offline. Run: Rscript tests/test_subtasks.R

suppressMessages(library(tidyverse))
here::i_am("tests/test_subtasks.R")
source(here::here("settings.R"))
invisible(lapply(list.files(here::here("R"), "[.]R$", full.names = TRUE), source))

ok <- 0L; fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) ok <<- ok + 1L
  else { fail <<- fail + 1L; cat("  FAIL:", msg, "\n") }
}
approx <- function(a, b, tol = 1e-6) all(abs(a - b) <= tol)

# Minimal config builders (only the fields each subtask reads). Overrides in `...`
# are given short names and auto-prefixed with the subtask, so e.g.
# kcfg("vanRaden_single", ridge = 0.01) sets kernel.ridge.
.prefixed <- function(prefix, ...) {
  ov <- list(...); if (length(ov)) names(ov) <- paste0(prefix, names(ov)); ov
}
kcfg <- function(method, ...) modifyList(list(
  kernel.method = method, kernel.ridge = 0, kernel.maf = 0.05,
  kernel.max_missing = 0.5, kernel.impute = "mean_round", kernel.rkhs_theta = 1),
  .prefixed("kernel.", ...))
pcfg <- function(method, ...) modifyList(list(
  pheno_prep.method = method, pheno_prep.z_thr = Inf, pheno_prep.standardize = "none",
  pheno_prep.ge_weighting = "none", pheno_prep.ge_bandwidth = 1),
  .prefixed("pheno_prep.", ...))
prcfg <- function(method, ...) modifyList(list(
  predict_post.method = method, predict_post.blend_obs_w = 0, predict_post.min_overlap = 2),
  .prefixed("predict_post.", ...))
no_trial <- list(lat = NA_real_, long = NA_real_, year = NA_integer_)

# ===========================================================================
cat(".brapi_try (retry through a flaky server)\n")
# Oracle: retries a transient error and returns the eventual success; a persistent error
# raises after exactly `tries` attempts; tries=1 disables retry. base_delay=0 keeps it fast.
n <- 0L; th <- function() { n <<- n + 1L; if (n < 3L) stop("boom"); "ok" }
check(identical(suppressMessages(.brapi_try(th, tries = 4, base_delay = 0)), "ok") && n == 3L,
      "retries a transient error, then succeeds (3rd attempt)")
n2 <- 0L; th2 <- function() { n2 <<- n2 + 1L; stop("always") }
check(inherits(suppressMessages(try(.brapi_try(th2, tries = 3, base_delay = 0), silent = TRUE)),
               "try-error") && n2 == 3L, "persistent error raises after exactly `tries` attempts")
n3 <- 0L; th3 <- function() { n3 <<- n3 + 1L; stop("x") }
suppressMessages(try(.brapi_try(th3, tries = 1, base_delay = 0), silent = TRUE))
check(n3 == 1L, "tries=1 disables retry (single attempt)")

# ===========================================================================
cat(".brapi_tries accessor + auth detection + re-login\n")
# Oracle: the tries default is single-sourced; NULL/absent settings falls back to 4.
check(.brapi_tries(list(brapi_tries = 7)) == 7L, ".brapi_tries reads settings$brapi_tries")
check(.brapi_tries(NULL) == 4L && .brapi_tries(list()) == 4L, ".brapi_tries default is 4")

# Oracle: an unauthorized 401 surfaces as a WARNING and as a response $status reason, and
# both are recognized (exact strings observed live from wheat.triticeaetoolbox.org).
check(.is_auth_warning("Unauthorized (HTTP 401).") &&
      .is_auth_warning("You must login and have permission to access this BrAPI call.") &&
      !.is_auth_warning("Internal Server Error (HTTP 500)."),
      ".is_auth_warning matches 401 wording, not 500")
# search shape: $status is a per-page LIST of http_status lists.
auth_search <- list(status = list(page0 = list(category = "Client error",
              reason = "Unauthorized", message = "Client error: (401) Unauthorized")),
              combined_data = list())
# wizard shape: $status is a FLAT http_status list (category/reason/message atomic). This is
# the shape that crashed the first implementation ("$ operator is invalid for atomic vectors").
auth_wizard <- list(status = list(category = "Client error", reason = "Unauthorized",
              message = "Client error: (401) Unauthorized"), data = list(ids = character()))
ok_wizard   <- list(status = list(category = "Success", reason = "OK",
              message = "Success: (200) OK"), data = list(ids = c("1", "2")))
check(.response_auth_failed(auth_search) && .response_auth_failed(auth_wizard),
      ".response_auth_failed detects 401 in both search (nested) and wizard (flat) $status")
check(!.response_auth_failed(ok_wizard) &&
      !.response_auth_failed(list(status = list(page0 = list(reason = "OK")))) &&
      !.response_auth_failed("not a response") &&
      !.response_auth_failed(list(data = 1)),
      ".response_auth_failed is FALSE (no crash) on success/other shapes")

# Oracle: on a 401 (emitted as a warning), .brapi_try re-logs in via conn/settings ONCE and
# retries, returning the eventual success. Uses a fake conn + throwaway env credentials.
old_u <- Sys.getenv("T3_USERNAME"); old_p <- Sys.getenv("T3_PASSWORD")
Sys.setenv(T3_USERNAME = "u", T3_PASSWORD = "p")
logged_in <- FALSE
# The fake conn must behave like a real BrAPI connection: login() mutates it IN PLACE,
# storing the token t3_login() then verifies. An environment, not a list, so it can.
fake_conn <- new.env(parent = emptyenv())
fake_conn$auth_token <- NULL
fake_conn$login <- function(username, password) {
  logged_in <<- TRUE; fake_conn$auth_token <- "tok"; invisible(NULL) }
tries_seen <- 0L
th_auth <- function() {
  tries_seen <<- tries_seen + 1L
  if (!logged_in) { warning("Unauthorized (HTTP 401)."); return(list(status = list(
    page0 = list(reason = "Unauthorized")), combined_data = list())) }
  "DATA"
}
res <- suppressMessages(.brapi_try(th_auth, conn = fake_conn, settings = list(brapi_tries = 3),
                                   base_delay = 0))
check(identical(res, "DATA") && logged_in && tries_seen == 2L,
      "401 -> t3_login() once -> retry succeeds")

# Oracle: BrAPI's login() assigns resp$content$access_token UNCONDITIONALLY, so a rejected
# password leaves auth_token NULL and returns NORMALLY -- no error, no exception. Unchecked,
# every later call goes out anonymous and the run reports data-shaped failures (empty
# searches, "no descriptor") for what is a credentials problem. t3_login must catch it.
mute_conn <- new.env(parent = emptyenv())
mute_conn$auth_token <- NULL
mute_conn$login <- function(username, password) invisible(NULL)   # no token, no error
bc <- tryCatch(t3_login(mute_conn), error = function(e) e)
check(inherits(bc, "t3_bad_credentials"), "login issuing no token -> t3_bad_credentials")
check(grepl("REJECTED", conditionMessage(bc)), "the message says the login was REJECTED")

# Oracle: .brapi_try fails FAST on bad creds, as it does on missing ones -- retrying just
# re-sends the same rejected password.
n_bad <- 0L
th_bad <- function() { n_bad <<- n_bad + 1L; warning("Unauthorized (HTTP 401).")
                       list(status = list(reason = "Unauthorized")) }
res_bad <- suppressMessages(tryCatch(
  .brapi_try(th_bad, conn = mute_conn, settings = list(brapi_tries = 4), base_delay = 0),
  error = function(e) e))
check(inherits(res_bad, "t3_bad_credentials") && n_bad == 1L,
      "bad creds: .brapi_try fails fast (1 attempt, no retry storm)")

# Oracle: absent credentials -> a t3_missing_credentials condition (not a generic error), so
# .brapi_try can fail fast on it rather than retry.
Sys.unsetenv("T3_USERNAME"); Sys.unsetenv("T3_PASSWORD")
mc <- tryCatch(t3_login(fake_conn), error = function(e) e)
check(inherits(mc, "t3_missing_credentials"), "missing creds -> t3_missing_credentials condition")
check(grepl("RESTART R", conditionMessage(mc)), "the message tells the user to RESTART R")

# Oracle: on missing creds, .brapi_try does NOT burn its retry budget -- it fails FAST (one
# thunk call), after LOUDLY reporting the re-login failure (no more silent "unauthorized").
n_thunk <- 0L
th401 <- function() { n_thunk <<- n_thunk + 1L; warning("Unauthorized (HTTP 401).")
                      list(status = list(reason = "Unauthorized")) }
msgs <- character()
res <- withCallingHandlers(
  tryCatch(.brapi_try(th401, conn = fake_conn, settings = list(brapi_tries = 4), base_delay = 0),
           error = function(e) e),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
check(inherits(res, "t3_missing_credentials") && n_thunk == 1L,
      "missing creds: .brapi_try fails fast (1 attempt, no retry storm)")
check(any(grepl("re-login FAILED", msgs)), "re-login failure is reported LOUDLY, not swallowed")

# Oracle: a TRANSIENT login failure (not missing creds) does not fail fast -- login is retried
# across attempts; here it succeeds on the 2nd try and the call then returns.
Sys.setenv(T3_USERNAME = "u", T3_PASSWORD = "p")
login_calls <- 0L; authed <- FALSE
flaky_conn <- new.env(parent = emptyenv())
flaky_conn$auth_token <- NULL
flaky_conn$login <- function(username, password) {
  login_calls <<- login_calls + 1L
  if (login_calls < 2L) stop("Timeout was reached")   # first login attempt fails transiently
  authed <<- TRUE; flaky_conn$auth_token <- "tok"; invisible(NULL) }
th_recover <- function() { if (!authed) { warning("Unauthorized (HTTP 401).")
                           list(status = list(reason = "Unauthorized")) } else "DATA" }
check(identical(suppressMessages(.brapi_try(th_recover, conn = flaky_conn,
                 settings = list(brapi_tries = 5), base_delay = 0)), "DATA") && login_calls >= 2L,
      "transient login failure is retried across attempts, then recovers")
Sys.unsetenv("T3_USERNAME"); Sys.unsetenv("T3_PASSWORD")
if (nzchar(old_u)) Sys.setenv(T3_USERNAME = old_u)
if (nzchar(old_p)) Sys.setenv(T3_PASSWORD = old_p)

# ===========================================================================
cat("VCF-download retry budget (.vcf_download_plan + .ensure_project_vcf give-up)\n")
# Oracle: effort decreases with each prior failure and stops at the cap.
rm(list = ls(envir = .vcf_download_fails), envir = .vcf_download_fails)   # clean slate
bs <- list(brapi_tries = 4, vcf_max_download_attempts = 3)
p0 <- .vcf_download_plan("P", bs)
check(!p0$skip && p0$tries == 4L, "no prior failures: full effort (tries = brapi_tries)")
.note_vcf_download_fail("P"); p1 <- .vcf_download_plan("P", bs)
check(!p1$skip && p1$tries == 3L && p1$base_delay < p0$base_delay,
      "one prior failure: fewer tries and shorter backoff")
.note_vcf_download_fail("P"); check(.vcf_download_plan("P", bs)$tries == 2L, "two prior failures: tries = 2")
.note_vcf_download_fail("P"); check(.vcf_download_plan("P", bs)$skip, "at the cap (3 failures): skip")
check(.vcf_download_plan("Q", bs)$tries == 4L, "the budget is per-project (Q unaffected)")
.clear_vcf_download_fail("P"); check(!.vcf_download_plan("P", bs)$skip, "clearing (a success) resets the project")

# Integration: a conn whose vcf_archived always errors (simulated timeout). .ensure_project_vcf
# must stop invoking it once the cap is reached, and never storm indefinitely.
etmp <- tempfile("cache_vcf_"); dir.create(etmp)
es <- modifyList(optimizer_settings(),
                 list(cache_dir = etmp, brapi_tries = 2, vcf_max_download_attempts = 3))
rm(list = ls(envir = .vcf_download_fails), envir = .vcf_download_fails)
calls <- 0L
bad_conn <- list(vcf_archived = function(output, genotyping_project_id) {
  calls <<- calls + 1L; stop("Timeout was reached") })
for (i in 1:6)
  suppressMessages(try(.ensure_project_vcf("8217", bad_conn, es), silent = TRUE))
# 3 failing attempts, each doing (decreasing) in-call retries, then attempts 4-6 skip with 0 calls.
check(calls > 0L && .vcf_download_fails[["8217"]] == 3L,
      "download failures are recorded and capped at vcf_max_download_attempts")
calls_at_cap <- calls
suppressMessages(try(.ensure_project_vcf("8217", bad_conn, es), silent = TRUE))
check(calls == calls_at_cap, "once capped, .ensure_project_vcf skips WITHOUT calling vcf_archived")
unlink(etmp, recursive = TRUE)
rm(list = ls(envir = .vcf_download_fails), envir = .vcf_download_fails)

# ===========================================================================
cat("sweep_rich_trials oracle (coverage + valid configs + degeneracy rule)\n")
ov <- .oracle_variants()
# Oracle: every subtask METHOD is exercised by at least one sweep variant (else the code
# oracle has a blind spot). SUBTASKS is the single source of truth for the method list.
meth_seen <- function(st) unique(vapply(ov, function(v) as.character(v$cfg[[paste0(st, ".method")]]), character(1)))
for (st in names(SUBTASKS)) {
  miss <- setdiff(SUBTASKS[[st]]$methods, meth_seen(st))
  check(length(miss) == 0, sprintf("sweep covers every %s method (missing: %s)", st, paste(miss, collapse = ",")))
}
# Oracle: each variant is a well-formed config (repair_config left no applicable *.method NA).
check(all(vapply(ov, function(v) all(vapply(names(SUBTASKS),
        function(st) !is.na(v$cfg[[paste0(st, ".method")]]), logical(1))), logical(1))),
      "every sweep variant is a complete config")
# Oracle: the degeneracy rule is exactly direct_blup under CV00 (and nothing else).
# Oracle: NOTHING is excused as degenerate. direct_blup x CV00 looks like it should be -- a fit
# on training lines alone has no random effect for the masked focal lines -- but predict_test()
# predicts those through the kernel, so it is a real test and excusing it would hide a
# regression in that guard.
check(!.oracle_degenerate(list(predict_post.method = "direct_blup"), "CV00") &&
      !.oracle_degenerate(list(predict_post.method = "direct_blup"), "CV0") &&
      !.oracle_degenerate(list(predict_post.method = "cond_expectation"), "CV00"),
      ".oracle_degenerate: no cell is excused (direct_blup x CV00 now predicts via the kernel)")
# Oracle: `only` selects variants by label substring; empty -> all; no match -> error.
check(length(.select_variants(ov, NULL)) == length(ov), ".select_variants: NULL keeps all")
sel <- .select_variants(ov, "em_combine")
check(length(sel) == 1 && grepl("em_combine", sel[[1]]$label), ".select_variants: 'em_combine' picks the one kernel variant")
check(length(.select_variants(ov, "model=")) >= 3, ".select_variants: 'model=' picks all model variants")
check(inherits(try(.select_variants(ov, "nonesuch"), silent = TRUE), "try-error"), ".select_variants: no match -> error")

# Oracle: .note_geno_once messages once per key per session, then stays silent for that key.
rm(list = ls(envir = .geno_note_seen), envir = .geno_note_seen)
n_msg <- 0L
withCallingHandlers({ .note_geno_once("k1", "hi"); .note_geno_once("k1", "hi"); .note_geno_once("k2", "yo") },
                    message = function(m) { n_msg <<- n_msg + 1L; invokeRestart("muffleMessage") })
check(n_msg == 2L, ".note_geno_once: fires once per key (2 distinct keys, repeat suppressed)")
rm(list = ls(envir = .geno_note_seen), envir = .geno_note_seen)

# ===========================================================================
cat("remote_server path resolution + cache backup/restore\n")
# Oracle: local mode puts state under the project dir and disables cache backup; remote mode
# puts state under OPTIMIZER_HOME and points the backup there; remote with OPTIMIZER_HOME unset
# fails LOUDLY (not silently at the filesystem root). `remote_server` is the settings.R global.
# HERMETIC: this block forces both `remote_server` and OPTIMIZER_HOME to known values and
# restores BOTH afterward, so it does not depend on -- or clobber -- the user's real settings
# (a run with remote_server = TRUE and OPTIMIZER_HOME set from .Renviron must be left intact).
old_rs <- remote_server
old_op <- Sys.getenv("OPTIMIZER_HOME", unset = NA_character_)
old_or <- Sys.getenv("OPTIMIZER_REMOTE", unset = NA_character_)

# Oracle: remote_server is ENV-DRIVEN (so the tracked settings.R is unedited on every machine):
# OPTIMIZER_REMOTE (truthy/falsy) wins; else it auto-detects from OPTIMIZER_HOME being set.
Sys.unsetenv("OPTIMIZER_HOME"); Sys.unsetenv("OPTIMIZER_REMOTE")
check(isFALSE(.detect_remote_server()), "no OPTIMIZER_HOME/REMOTE -> local (FALSE)")
Sys.setenv(OPTIMIZER_HOME = "/tmp/opt_perm_test")
check(isTRUE(.detect_remote_server()), "OPTIMIZER_HOME set -> remote (auto-detect)")
Sys.setenv(OPTIMIZER_REMOTE = "false")
check(isFALSE(.detect_remote_server()), "OPTIMIZER_REMOTE=false overrides the auto-detect")
Sys.unsetenv("OPTIMIZER_HOME"); Sys.setenv(OPTIMIZER_REMOTE = "yes")
check(isTRUE(.detect_remote_server()), "OPTIMIZER_REMOTE=yes -> remote even without OPTIMIZER_HOME")
Sys.unsetenv("OPTIMIZER_HOME"); Sys.unsetenv("OPTIMIZER_REMOTE")

# From here, force the global into known LOCAL state (the source-time value may be remote if
# the env had OPTIMIZER_HOME); optimizer_settings() reads this global, not the env.
remote_server <<- FALSE

# Oracle: the startup canary check is a boolean opt-out flag (value is the user's choice).
check({ v <- optimizer_settings()$run_startup_canary; is.logical(v) && length(v) == 1 && !is.na(v) },
      "run_startup_canary is a single logical (TRUE/FALSE opt-out flag)")

# local_overrides = FALSE throughout this block: these assertions are about how paths are
# DERIVED from remote_server/OPTIMIZER_HOME, and a settings.local.R that legitimately moves
# db_path (e.g. onto /workdir so WAL works for parallel workers) would otherwise fail them.
# The override MECHANISM is tested separately, just below.
sl <- optimizer_settings(local_overrides = FALSE)
check(sl$db_path == file.path(here::here(), "state", "evals.sqlite") && is.null(sl$cache_backup_dir),
      "local mode: state under project dir, cache backup disabled")

remote_server <<- TRUE
Sys.setenv(OPTIMIZER_HOME = "/tmp/opt_perm_test")
sr <- optimizer_settings(local_overrides = FALSE)
check(sr$db_path == "/tmp/opt_perm_test/state/evals.sqlite" &&
      sr$log_dir == "/tmp/opt_perm_test/logs" &&
      sr$cache_backup_dir == "/tmp/opt_perm_test/cache",
      "remote mode: state + backup under OPTIMIZER_HOME")
check(sr$cache_dir == here::here("cache"), "remote mode: cache still on the work disk")
# The store backup that makes a local-disk db_path safe to lose (parallel workers require it).
check(sr$db_backup_path == "/tmp/opt_perm_test/state/evals_backup.sqlite",
      "remote mode: db_backup_path under OPTIMIZER_HOME")
Sys.unsetenv("OPTIMIZER_HOME")
check(inherits(try(optimizer_settings(local_overrides = FALSE), silent = TRUE), "try-error"),
      "remote_server = TRUE with OPTIMIZER_HOME unset -> loud error (not root paths)")

# (The override MECHANISM itself -- .apply_overrides layering and .local_overrides parsing --
# is covered by the "settings.local.R overrides" block below, using temp files.)

# Oracle: *_path / *_file name a FILE, *_dir names a DIRECTORY, and giving the wrong KIND is
# rejected at settings time. The mistake that motivated this is db_backup_path set to a
# directory: file.rename() then returns FALSE (it does NOT error), so the backup silently
# never happened. Shape only -- a file whose parent does not exist yet is normal.
local({
  base <- list(db_path = "/a/evals.sqlite", db_backup_path = NULL,
               report_path = "/a/report.md", stop_file = "/a/STOP",
               cache_ready_file = "/a/.cache_ready", cache_dir = "/a/cache",
               log_dir = "/a/logs", cache_backup_dir = NULL)
  bad <- function(...) inherits(try(.validate_paths(modifyList(base, list(...))), silent = TRUE),
                                "try-error")
  adir <- tempfile("vdir_"); dir.create(adir)
  afile <- tempfile("vfile_"); writeLines("x", afile)

  check(!bad(),                            ".validate_paths: file paths whose parents do not exist are fine")
  check(bad(db_path = adir),               ".validate_paths: db_path as a directory errors")
  check(bad(db_backup_path = adir),        ".validate_paths: db_backup_path as a directory errors")
  check(bad(db_path = "/a/b/"),            ".validate_paths: a trailing separator errors (basename() hides it)")
  check(bad(cache_dir = afile),            ".validate_paths: cache_dir as a file errors (the mirror mistake)")
  check(!bad(db_backup_path = NULL, cache_backup_dir = NULL),
                                           ".validate_paths: NULL backup paths are fine (local mode)")
  check(!bad(cache_dir = adir, log_dir = adir),
                                           ".validate_paths: existing directories for *_dir are fine")
  unlink(c(adir, afile), recursive = TRUE)
})

# Oracle: the per-evaluation RSS peak is only trusted when the kernel PROVES it reset the
# high-water mark. VmHWM is monotonic, so a kernel that ignores the clear_refs write would
# otherwise have every row record the whole worker's high-water mark -- a plausible-looking
# number that is not a per-evaluation peak. Verified by shadowing the /proc readers, so this
# runs on any platform (the real reset only exists on Linux).
local({
  keep <- list(r = .rss_peak_reset, h = proc_peak_rss_mb, s = proc_rss_mb)
  on.exit({ .rss_peak_reset <<- keep$r; proc_peak_rss_mb <<- keep$h; proc_rss_mb <<- keep$s
            rm(list = ls(.mem_env), envir = .mem_env) }, add = TRUE)
  probe <- function(reset_ok, hwm, rss) {
    rm(list = ls(.mem_env), envir = .mem_env)          # clear the cached answer
    .rss_peak_reset  <<- function() reset_ok
    proc_peak_rss_mb <<- function() hwm
    proc_rss_mb      <<- function() rss
    suppressMessages(.rss_peak_resettable())
  }
  check(isTRUE(probe(TRUE, 2000, 1990)),        "rss peak: reset honoured (VmHWM fell to VmRSS) -> usable")
  check(isFALSE(probe(TRUE, 9000, 1990)),       "rss peak: VmHWM stayed high -> NOT usable (the dangerous case)")
  check(isFALSE(probe(FALSE, 2000, 1990)),      "rss peak: clear_refs write failed -> not usable")
  check(isFALSE(probe(TRUE, NA_real_, 1990)),   "rss peak: /proc unreadable -> not usable")
  # And the NA propagates rather than a stale number reaching the store.
  invisible(probe(TRUE, 9000, 1990))
  check(is.na(mem_peak_rss_mb()), "rss peak: mem_peak_rss_mb() is NA when the reset is unusable")
})

# Oracle: backup_store REPORTS a failed backup. file.rename returns FALSE (with a warning)
# rather than erroring when dest is a directory, so a tryCatch on `error` alone saw nothing and
# the run reported a healthy backup that never happened.
local({
  dbf <- tempfile(fileext = ".sqlite")
  con <- open_store(dbf)
  store_eval(con, sample_config(), "t1", "CV0", 0.5, 10L, "ok")
  good <- tempfile(fileext = ".sqlite")
  check(isTRUE(backup_store(con, good)) && file.exists(good),
        "backup_store: writes the backup file and returns TRUE")
  adir <- tempfile("bdir_"); dir.create(adir)
  msgs <- character()
  res <- withCallingHandlers(backup_store(con, adir),
           message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                                   invokeRestart("muffleMessage") })
  check(isFALSE(res), "backup_store: returns FALSE when dest is a directory")
  check(any(grepl("FAILED", msgs)) && any(grepl("DIRECTORY", msgs)),
        "backup_store: SAYS SO -- a failed backup is never silent")
  check(length(list.files(dirname(adir), pattern = "[.]tmp[0-9]+$")) == 0,
        "backup_store: cleans up its temp file on failure")
  close_store(con); unlink(c(dbf, good)); unlink(adir, recursive = TRUE)
})

remote_server <<- old_rs                                          # restore the global and env vars
if (is.na(old_op)) Sys.unsetenv("OPTIMIZER_HOME")   else Sys.setenv(OPTIMIZER_HOME   = old_op)
if (is.na(old_or)) Sys.unsetenv("OPTIMIZER_REMOTE") else Sys.setenv(OPTIMIZER_REMOTE = old_or)

# Oracle: cache backup is additive and excludes raw_project/; restore ADDITIVELY reconciles
# the work cache from the backup -- filling gaps in a PARTIALLY-populated cache (the case that
# regressed) without deleting work-only files. (rsync-dependent; skipped cleanly if absent.)
if (nzchar(Sys.which("rsync"))) {
  work <- tempfile("work_"); bak <- tempfile("bak_"); dir.create(work)
  dir.create(file.path(work, "dosage"), recursive = TRUE)
  dir.create(file.path(work, "raw_project"), recursive = TRUE)
  saveRDS(1L, file.path(work, "dosage", "dosage_1_sz1.rds"))
  saveRDS(2L, file.path(work, "dosage", "dosage_2_sz1.rds"))
  writeLines("x", file.path(work, "raw_project", "raw_project_1.vcf"))
  bset <- list(cache_dir = work, cache_backup_dir = bak)
  sync_cache_to_backup(bset)
  check(file.exists(file.path(bak, "dosage", "dosage_1_sz1.rds")),
        "cache backup copies dosage files to the backup dir")
  check(!file.exists(file.path(bak, "raw_project", "raw_project_1.vcf")),
        "cache backup EXCLUDES the transient raw_project/ VCFs")
  # restore into an EMPTY work cache -> full warm-up.
  empty <- tempfile("empty_"); dir.create(empty)
  restore_cache_from_backup(list(cache_dir = empty, cache_backup_dir = bak))
  check(file.exists(file.path(empty, "dosage", "dosage_2_sz1.rds")),
        "restore warms an empty work cache from the backup")
  # restore into a PARTIAL work cache -> fills the MISSING file, keeps a work-only file (the
  # regression: the old empty-guard skipped this entirely).
  part <- tempfile("part_"); dir.create(file.path(part, "dosage"), recursive = TRUE)
  saveRDS(1L, file.path(part, "dosage", "dosage_1_sz1.rds"))   # has #1, missing #2
  saveRDS(9L, file.path(part, "work_only.rds"))                # not in the backup
  restore_cache_from_backup(list(cache_dir = part, cache_backup_dir = bak))
  check(file.exists(file.path(part, "dosage", "dosage_2_sz1.rds")),
        "restore fills a MISSING file in a partially-populated work cache")
  check(file.exists(file.path(part, "work_only.rds")),
        "restore is additive: a work-only file is not deleted")
  unlink(c(work, empty, part, bak), recursive = TRUE)
} else message("  (rsync not found -- skipping cache backup/restore integration checks)")

# ===========================================================================
cat("settings.local.R overrides (untracked, machine-specific)\n")
# Oracle: .apply_overrides layers overrides on defaults and warns on an unknown (typo) key.
check(identical(.apply_overrides(list(a = 1, b = 2), list()), list(a = 1, b = 2)),
      ".apply_overrides: no overrides -> defaults unchanged")
check(identical(.apply_overrides(list(a = 1, b = 2), list(b = 9))$b, 9),
      ".apply_overrides: an override replaces a default")
check(inherits(tryCatch(.apply_overrides(list(a = 1), list(typo = 3)), warning = function(w) w),
               "warning"), ".apply_overrides: an unknown key warns (typo guard)")
# Oracle: .local_overrides reads settings_override from a file; absent -> list(); non-list -> error.
lf1 <- tempfile(fileext = ".R")
writeLines("settings_override <- list(dosage_budget_bytes = 16e9, simulate = FALSE)", lf1)
ov1 <- .local_overrides(lf1)
check(ov1$dosage_budget_bytes == 16e9 && isFALSE(ov1$simulate), ".local_overrides reads settings_override")
check(identical(.local_overrides(tempfile(fileext = ".R")), list()), ".local_overrides: absent file -> empty list")
lf2 <- tempfile(fileext = ".R"); writeLines("settings_override <- 42", lf2)
check(inherits(try(.local_overrides(lf2), silent = TRUE), "try-error"),
      ".local_overrides: non-list settings_override -> error")
unlink(c(lf1, lf2))
# End-to-end via optimizer_settings() -- only if there is no REAL settings.local.R to clobber.
lf <- here::here("settings.local.R")
if (!file.exists(lf)) {
  writeLines("settings_override <- list(dosage_budget_bytes = 7e9, run_startup_canary = FALSE)", lf)
  s_ov <- suppressWarnings(optimizer_settings())
  check(s_ov$dosage_budget_bytes == 7e9 && isFALSE(s_ov$run_startup_canary),
        "optimizer_settings() applies settings.local.R overrides")
  invisible(file.remove(lf))
} else message("  (a real settings.local.R exists -- skipping the end-to-end override check)")

# ===========================================================================
cat("cached() / category-partitioned cache paths\n")
# Oracle: a cache entry lives at cache/<category>/<category>_<identifier>.rds; reads fall
# back to the pre-migration FLAT path so an un-migrated cache still hits; writes go nested;
# valid= gates the write. Uses a throwaway cache dir so the real cache is untouched.
.ctmp <- tempfile("cache_test_"); dir.create(.ctmp)
cs <- modifyList(optimizer_settings(), list(cache_dir = .ctmp))

# path construction: nested, with the singleton (no identifier) special case.
check(.cache_path(cs, "acc", "10676") == file.path(.ctmp, "acc", "acc_10676.rds"),
      "nested path is cache/<category>/<category>_<identifier>.rds")
check(.cache_path(cs, "trial_catalog") == file.path(.ctmp, "trial_catalog", "trial_catalog.rds"),
      "singleton path (no identifier) is cache/<category>/<category>.rds")
check(.cache_path(cs, "raw_project", "8130", ext = "vcf") ==
        file.path(.ctmp, "raw_project", "raw_project_8130.vcf"),
      "ext= controls the extension (raw VCF)")

# write + read round-trip through cached().
v <- cached(cs, "acc", "A1", expr = c("x", "y"))
check(file.exists(file.path(.ctmp, "acc", "acc_A1.rds")), "cached() writes the nested file")
check(identical(cached(cs, "acc", "A1", expr = stop("must not evaluate on a hit")), c("x", "y")),
      "cached() re-read hits the nested file (expr not evaluated)")

# flat fallback: a legacy flat file is found even though nothing is nested.
saveRDS("legacy", file.path(.ctmp, "acc_L1.rds"))
check(identical(cached(cs, "acc", "L1", expr = stop("must not evaluate on a hit")), "legacy"),
      "flat fallback: cached() reads a pre-migration flat file")
check(!is.na(.cache_existing(cs, "acc", "L1")) &&
        is.na(.cache_existing(cs, "acc", "NOPE")),
      ".cache_existing finds a present key (flat) and NA for an absent one")

# nested is preferred over flat when both exist.
saveRDS("flat", file.path(.ctmp, "obs_B1.rds")); .cache_save(cs, "obs", "B1", "nested")
check(identical(readRDS(.cache_existing(cs, "obs", "B1")), "nested"),
      "nested wins when both nested and flat exist")

# valid= gates the write: an invalid (empty) result is returned but NOT cached.
cached(cs, "acc", "E1", valid = function(a) length(a) > 0, expr = character())
check(is.na(.cache_existing(cs, "acc", "E1")), "valid= gates: an invalid result is not cached")

# .find_densest_dosage globs BOTH the nested dosage folder and the legacy flat cache, and
# returns the DENSEST (smallest-thin) file when several linger.
.cache_save(cs, "dosage", "700_sz123", matrix(0L, 1, 1))                 # thin 1, nested
saveRDS(matrix(0L, 1, 1), file.path(.ctmp, "dosage_701_sz9.rds"))        # legacy flat
.cache_save(cs, "dosage", "702_thin5_sz1", matrix(0L, 1, 5))            # thin 5
.cache_save(cs, "dosage", "702_thin2_sz1", matrix(0L, 1, 5))            # thin 2 (denser)
check(!is.null(.find_densest_dosage(cs, "700")), ".find_densest_dosage finds a nested dosage")
check(!is.null(.find_densest_dosage(cs, "701")), ".find_densest_dosage finds a legacy-flat dosage")
check(is.null(.find_densest_dosage(cs, "999")), ".find_densest_dosage returns NULL when absent")
check(.find_densest_dosage(cs, "702")$thin == 2L, ".find_densest_dosage picks the densest (thin 2 over 5)")
unlink(.ctmp, recursive = TRUE)

# ===========================================================================
cat("score_predictions\n")
# Oracle: correlation of a vector with a positive affine transform of itself is 1.
obs <- setNames(c(1, 2, 3, 4, 5, 6), letters[1:6])
r <- score_predictions(setNames(2 * obs + 10, names(obs)), obs)
check(r$status == "ok" && r$n_test == 6 && approx(r$score, 1), "self-correlation = 1")
# Oracle: a negative transform gives -1.
r <- score_predictions(setNames(-obs, names(obs)), obs)
check(approx(r$score, -1), "negated = -1")
# Oracle: only shared accession names are joined (extras on each side ignored).
p <- c(a = 1, b = 2, c = 3, d = 4, e = 5, f = 6, z = 99)   # z not in obs
o <- c(obs, y = 99)                                         # y not in pred
r <- score_predictions(p, o)
check(r$n_test == 6, "join keeps only shared names")
# Oracle: fewer than 5 overlapping accessions is a defined failure.
r <- score_predictions(c(a = 1, b = 2, c = 3, d = 4), obs)
check(is.na(r$score) && r$status == "too_few_overlap", "n<5 -> too_few_overlap")
# Oracle: constant predictions have no variance -> defined failure.
r <- score_predictions(setNames(rep(5, 6), names(obs)), obs)
check(is.na(r$score) && r$status == "constant", "constant pred -> constant")

# Oracle: a JOIN THAT SUCCEEDED but whose values are all non-finite is a MODELLING failure,
# not a coverage one, and must not hide behind too_few_overlap. Both used to report
# n_test = 0 (it is counted after the finite filter), which made them indistinguishable in
# the store -- and sent a real 2026-07-31 investigation down two wrong paths.
r <- score_predictions(setNames(rep(NaN, 6), names(obs)), obs)
check(r$status == "non_finite" && r$n_test == 0,
      "names join but predictions are all NaN -> non_finite, not too_few_overlap")
check(grepl("pred non-finite 6", r$reason), "the reason names the side that went non-finite")
r <- score_predictions(setNames(2 * obs, names(obs)), setNames(rep(Inf, 6), names(obs)))
check(r$status == "non_finite" && grepl("obs non-finite 6", r$reason),
      "non-finite OBSERVATIONS are reported as such, not blamed on the predictions")

# Oracle: nothing joined at all stays too_few_overlap -- that IS a coverage problem.
r <- score_predictions(c(x = 1, y = 2, z = 3), obs)
check(r$status == "too_few_overlap" && r$n_test == 0,
      "no shared names -> still too_few_overlap (a coverage problem)")
check(grepl("0 joined", r$reason), "too_few_overlap now records how many joined vs stayed finite")

# Oracle: a partial non-finite set still scores when >=5 finite pairs survive.
p_part <- setNames(2 * obs, names(obs)); p_part["a"] <- NA
r <- score_predictions(p_part, obs)
check(r$status == "ok" && r$n_test == 5, "one NaN prediction drops that line, the rest score")

# ===========================================================================
cat(".vcf_to_dosage (orientation + encoding; guards the transposition bug)\n")
vcf_lines <- c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
        "S1", "S2", sep = "\t"),
  paste("1", "100", "m1", "A", "G", ".", ".", ".", "GT", "0/0", "1/1", sep = "\t"),
  paste("1", "200", "m2", "A", "G", ".", ".", ".", "GT", "0/1", "0/0", sep = "\t"),
  paste("1", "300", "m3", "A", "G", ".", ".", ".", "GT", "1/1", "0/1", sep = "\t"))
vpath <- tempfile(fileext = ".vcf"); writeLines(vcf_lines, vpath)
# Oracle (hand-decoded): S1 = (m1=0/0=0, m2=0/1=1, m3=1/1=2); S2 = (2, 0, 1).
d <- suppressWarnings(.vcf_to_dosage(vpath, c("S1", "S2"), thin = 1))
check(identical(dim(d), c(2L, 3L)), "dosage is accessions x markers (2x3)")
check(identical(rownames(d), c("S1", "S2")), "rownames are samples")
check(identical(colnames(d), c("m1", "m2", "m3")), "colnames are marker IDs")
check(approx(d["S1", ], c(0, 1, 2)) && approx(d["S2", ], c(2, 0, 1)), "encoding + orientation")
# Oracle: subsetting to one sample yields one row.
d1 <- suppressWarnings(.vcf_to_dosage(vpath, "S2", thin = 1))
check(nrow(d1) == 1 && approx(d1["S2", ], c(2, 0, 1)), "single-sample subset")
# Oracle: thin=2 keeps markers 1 and 3 (every 2nd, genome-wide).
dt <- suppressWarnings(.vcf_to_dosage(vpath, c("S1", "S2"), thin = 2))
check(identical(colnames(dt), c("m1", "m3")), "marker thinning keeps every 2nd")

# Oracle: a malformed variant line (wrong field count) is SKIPPED, not fatal -- the good
# markers still come through. (This is the real 8114 failure: one 3-field variant record.)
bad_lines <- c(vcf_lines[1:5],
               paste("1", "250", "BAD", sep = "\t"),            # 3 fields, not 11
               vcf_lines[6])
bpath <- tempfile(fileext = ".vcf"); writeLines(bad_lines, bpath)
db <- .vcf_to_dosage(bpath, c("S1", "S2"), thin = 1)
check(identical(colnames(db), c("m1", "m2", "m3")), "malformed variant line skipped, good markers kept")

# Oracle: a transposed / non-VCF header (markers as columns) is REJECTED, not misread.
# (The real 14511: the "#CHROM" row holds chromosome labels; field 2 is not POS.)
tr_lines <- c("##fileformat=VCFv4.2",
              paste("#CHROM", "1A", "1A", "1B", "2A", sep = "\t"),
              paste("POS", "100", "200", "300", "400", sep = "\t"))
tpath <- tempfile(fileext = ".vcf"); writeLines(tr_lines, tpath)
check(inherits(try(.vcf_to_dosage(tpath, NULL, 1L), silent = TRUE), "try-error"),
      "transposed / non-VCF header rejected")

# Oracle: .vcf_stat validates + counts in one pass.
st <- .vcf_stat(vpath)
check(st$n_samples == 2 && st$n_markers == 3 && identical(st$samples, c("S1", "S2")),
      ".vcf_stat: 2 samples, 3 markers")
check(inherits(try(.vcf_stat(tpath), silent = TRUE), "try-error"), ".vcf_stat rejects non-VCF")

# Oracle: .cache_thin -- the densest thin that fits the budget; depends only on size, NOT on
# any requested thin. A project that fits caches at thin 1 (full markers).
# Integer inputs on purpose: n_samples*n_markers must not overflow 32-bit int (683L*7.47ML).
check(.cache_thin(100L, 1000L, 2e9) == 1, ".cache_thin: a project that fits caches at thin 1")
# 683 x 7.47M x 4 = 20.4 GB; /2 GB budget -> ceil 10.2 -> 11.
check(.cache_thin(683L, 7467224L, 2e9) == 11, ".cache_thin: oversized project thinned to fit (no int overflow)")
check(.cache_thin(683L, 7467224L, 4e9) == 6, ".cache_thin: a bigger budget caches denser")

# ===========================================================================
cat("get_project_dosage: derive requested thin by column-subset (no re-parse)\n")
# Oracle: a project is cached ONCE at its densest thin; a requested marker_thin is served by
# keeping every k = max(1, floor(marker_thin / e))-th marker of that cache. Uses a throwaway
# cache + conn = NULL (never reached -- always a cache hit), so it is fully offline.
.dtmp <- tempfile("cache_dose_"); dir.create(.dtmp)
ds <- modifyList(optimizer_settings(), list(cache_dir = .dtmp))
Xfull <- matrix(as.integer(1:120), nrow = 4,
                dimnames = list(paste0("s", 1:4), paste0("m", 1:30)))   # 4 samples x 30 markers
.cache_save(ds, "dosage", "P1_sz100", Xfull)                            # cached at thin 1 (full)

d1  <- get_project_dosage("P1", NULL, NULL, ds, marker_thin = 1L)
d3  <- get_project_dosage("P1", NULL, NULL, ds, marker_thin = 3L)
d10 <- get_project_dosage("P1", NULL, NULL, ds, marker_thin = 10L)
check(identical(d1, Xfull), "thin 1 from a full cache returns the whole matrix")
check(identical(d3, Xfull[, seq(1L, 30L, by = 3L), drop = FALSE]),
      "thin 3 from a full cache == every 3rd marker")
check(identical(d10, Xfull[, seq(1L, 30L, by = 10L), drop = FALSE]),
      "thin 10 from a full cache == every 10th marker")

# Oracle: from a cache already at thin e=2, request r=5 -> k=floor(5/2)=2 -> every 2nd marker
# of the (already-thinned) cache; and request r=2 < e -> k=1 -> the cache as-is (densest we have).
Xe2 <- matrix(as.integer(1:40), nrow = 4,
              dimnames = list(paste0("s", 1:4), paste0("m", 1:10)))     # a thin=2 cache
.cache_save(ds, "dosage", "P2_thin2_sz100", Xe2)
check(identical(get_project_dosage("P2", NULL, NULL, ds, marker_thin = 5L),
                Xe2[, seq(1L, 10L, by = 2L), drop = FALSE]),
      "e=2, r=5 -> k=2 (denser-multiple derivation)")
check(identical(get_project_dosage("P2", NULL, NULL, ds, marker_thin = 2L), Xe2),
      "e=2, r<e -> k=1 -> serve the densest cache as-is")

# Oracle: sample subsetting still applies after the column subset.
check(identical(get_project_dosage("P1", c("s1", "s3"), NULL, ds, marker_thin = 3L),
                Xfull[c("s1", "s3"), seq(1L, 30L, by = 3L), drop = FALSE]),
      "keep_samples subsets rows after the thin subsets columns")
unlink(.dtmp, recursive = TRUE)

# Oracle: re-densify -- a bigger budget re-downloads a cache that a smaller-budget machine
# thinned (so a big server is not pinned to a laptop-thinned warmed cache). A fake conn writes
# a small full VCF; the pre-seeded coarse thin4 cache must be replaced by a full-density parse.
rdmp <- tempfile("redense_"); dir.create(rdmp)
rds  <- modifyList(optimizer_settings(),
                   list(cache_dir = rdmp, dosage_budget_bytes = 1e12,   # huge -> ideal thin = 1
                        dosage_redensify = TRUE, brapi_tries = 1))
.cache_save(rds, "stat", "PZ", list(n_samples = 3L, n_markers = 3L))    # ideal thin here = 1
.cache_save(rds, "dosage", "PZ_thin4_sz1",                              # a COARSE thin-4 cache
            matrix(0L, 3, 1, dimnames = list(c("S1", "S2", "S3"), "m1")))
vln <- c("##fileformat=VCFv4.2",
  paste("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","S1","S2","S3", sep = "\t"),
  paste("1","10","m1","A","G",".",".",".","GT","0/0","0/1","1/1", sep = "\t"),
  paste("1","20","m2","A","G",".",".",".","GT","0/1","1/1","0/0", sep = "\t"),
  paste("1","30","m3","A","G",".",".",".","GT","1/1","0/0","0/1", sep = "\t"))
fconn <- list(vcf_archived = function(output, genotyping_project_id) writeLines(vln, output))
rmsg <- character()
dz <- withCallingHandlers(get_project_dosage("PZ", NULL, fconn, rds, marker_thin = 1L),
        message = function(m) { rmsg <<- c(rmsg, conditionMessage(m)); invokeRestart("muffleMessage") })
check(any(grepl("re-densify", rmsg)), "re-densify: a coarser-than-affordable cache triggers a re-download")
check(!is.null(dz) && ncol(dz) == 3L, "re-densify: re-parsed at full density (3 markers, was thinned to 1)")
check(length(list.files(file.path(rdmp, "dosage"), "thin4")) == 0, "re-densify: the coarse thin4 cache is superseded")
# And the opposite: with dosage_redensify = FALSE, the coarse cache is served as-is (no download).
.cache_save(rds, "dosage", "PZ2_thin4_sz1", matrix(0L, 3, 1, dimnames = list(c("S1","S2","S3"), "m1")))
.cache_save(rds, "stat", "PZ2", list(n_samples = 3L, n_markers = 3L))
rds_off <- modifyList(rds, list(dosage_redensify = FALSE))
check(ncol(get_project_dosage("PZ2", NULL, NULL, rds_off, marker_thin = 1L)) == 1L,
      "dosage_redensify = FALSE keeps the coarse cache (no re-download; conn = NULL never called)")
unlink(rdmp, recursive = TRUE)

# ===========================================================================
cat(".qc_markers\n")
set.seed(1)
good <- matrix(rep(c(0, 1, 2, 0, 1, 2), 60), nrow = 6)        # 60 polymorphic markers
colnames(good) <- paste0("g", 1:60)
Xqc <- cbind(good,
             mono = rep(0, 6),                  # MAF 0 -> dropped
             hm   = c(0, 1, NA, NA, NA, NA),     # 4/6 missing > 0.5 -> dropped
             imp  = c(0, 2, 2, 2, 2, NA))        # one NA -> imputed
rownames(Xqc) <- paste0("id", 1:6)
q <- .qc_markers(Xqc, kcfg("vanRaden_single"))
check(!("mono" %in% colnames(q)), "monomorphic marker dropped")
check(!("hm" %in% colnames(q)), "high-missing marker dropped")
check("imp" %in% colnames(q) && !anyNA(q), "kept marker imputed (no NA remains)")
# Oracle: mean of non-NA imp values is (0+2+2+2+2)/5 = 1.6; mean_round -> 2.
check(q["id6", "imp"] == 2, "mean_round imputation = round(1.6) = 2")
check(.qc_markers(Xqc, kcfg("vanRaden_single", impute = "mean"))["id6", "imp"] == 1.6,
      "mean imputation = 1.6")
# Oracle: fewer than 50 surviving markers is an error.
check(inherits(try(.qc_markers(good[, 1:10], kcfg("vanRaden_single")), silent = TRUE),
               "try-error"), "errors below 50 markers")

# ===========================================================================
cat(".merge_markers\n")
mk <- function(rows, cols) { m <- matrix(1, length(rows), length(cols),
  dimnames = list(rows, cols)); m }
d1 <- mk(c("a", "b", "c", "d"), paste0("m", 1:60))
d2 <- mk(c("c", "d", "e", "f"), paste0("m", 11:80))   # shares m11..m60 = 50 markers
mg <- .merge_markers(list(d1, d2))
# Oracle: 50 shared markers; accessions a,b,c,d,e,f with duplicates c,d dropped.
check(ncol(mg) == 50, "merge keeps the 50-marker intersection")
check(setequal(rownames(mg), c("a", "b", "c", "d", "e", "f")) && !anyDuplicated(rownames(mg)),
      "every accession once; duplicates c,d dropped")
# Oracle: when an accession sits in two projects, its row must come from the project with
# the RICHER marker build (more markers = the better-genotyped call of the same line).
rich <- matrix(2, 2, 70, dimnames = list(c("c", "d"), paste0("m", 11:80)))   # 70 markers, all 2s
poor <- matrix(0, 2, 60, dimnames = list(c("c", "d"), paste0("m", 1:60)))    # 60 markers, all 0s
check(all(.merge_markers(list(poor = poor, rich = rich))["c", ] == 2),
      "duplicated accession takes its row from the richer marker build")
# Oracle: <50 shared markers -> fall back to the project with more accessions.
small1 <- mk(c("a", "b", "c", "d"), paste0("m", 1:40))
small2 <- mk(c("a", "b", "c", "d", "e", "f"), paste0("p", 1:40))   # 0 shared markers
check(identical(.merge_markers(list(small1, small2)), small2), "no-overlap -> largest project")
check(identical(.merge_markers(list(d1)), d1), "single project returned as-is")

# ===========================================================================
cat("build_kernel\n")
set.seed(2)
Xc <- matrix(sample(0:2, 6 * 60, replace = TRUE), nrow = 6)
Xc[2, ] <- Xc[1, ]                                    # id2 is a clone of id1
rownames(Xc) <- paste0("id", 1:6); colnames(Xc) <- paste0("m", 1:60)
all6 <- rownames(Xc)
K <- build_kernel(kcfg("vanRaden_single"), list(Xc), all6)
check(isSymmetric(unname(K)), "VanRaden GRM symmetric")
# Oracle: genetically identical lines have identical relationship vectors.
check(approx(K["id1", ], K["id2", ], tol = 1e-8), "clone lines -> identical GRM rows")
# Oracle: adding ridge raises only the diagonal, by exactly ridge.
K0 <- build_kernel(kcfg("vanRaden_single", ridge = 0), list(Xc), all6)
K1 <- build_kernel(kcfg("vanRaden_single", ridge = 0.01), list(Xc), all6)
check(approx(diag(K1) - diag(K0), 0.01, 1e-8), "ridge adds exactly 0.01 to diagonal")
check(approx((K1 - K0)[upper.tri(K1)], 0, 1e-8), "ridge leaves off-diagonals unchanged")
# Oracle: em_combine with a single project is the plain VanRaden GRM.
check(approx(build_kernel(kcfg("em_combine"), list(Xc), all6), K0, 1e-8),
      "em_combine, one project -> VanRaden")
# RKHS Gaussian kernel.
Kr <- build_kernel(kcfg("rkhs_gaussian", rkhs_theta = 1, ridge = 0), list(Xc), all6)
check(isSymmetric(unname(Kr)), "RKHS kernel symmetric")
check(approx(diag(Kr), 1, 1e-8), "RKHS diagonal = 1")
check(all(Kr > 0 & Kr <= 1 + 1e-9), "RKHS entries in (0,1]")
check(approx(Kr["id1", "id2"], 1, 1e-8), "RKHS clone pair = 1")

# ===========================================================================
cat(".vanraden (population frequencies; the A.mat coding bug)\n")
set.seed(7)
npop <- 200; nmk <- 300
Xp <- matrix(sample(0:2, npop * nmk, replace = TRUE), nrow = npop,
             dimnames = list(paste0("g", 1:npop), paste0("m", 1:nmk)))
need4 <- c("g3", "g17", "g42", "g99")
# Oracle (the whole point of the change): the GRM built over `need` with population
# allele frequencies IS the [need, need] submatrix of the whole population's GRM.
# So "QC + GRM on the population, then subset" and this are the same numbers.
Kfull <- .vanraden(Xp, rownames(Xp))
Ksub  <- .vanraden(Xp, need4)
check(approx(Ksub, Kfull[need4, need4], 1e-10),
      "GRM on needed rows w/ population freqs == submatrix of the full-population GRM")

# Oracle: hand-computed VanRaden on a panel where EVERY marker is alt-major (p > 0.5).
# This is the case rrBLUP::A.mat silently destroys: it assumes {-1,0,1}, so on {0,1,2}
# its internal freq = p + 0.5 > 1, its MAF goes negative, and min.MAF drops the marker.
set.seed(9)
pmaj <- runif(nmk, 0.7, 0.95)                        # every marker alt-major
Xmaj <- sapply(pmaj, function(pk) rbinom(60, 2, pk))
dimnames(Xmaj) <- list(paste0("g", 1:60), paste0("m", 1:nmk))
phat  <- colMeans(Xmaj) / 2
Wman  <- sweep(Xmaj, 2, 2 * phat, "-")
Kman  <- tcrossprod(Wman) / (2 * sum(phat * (1 - phat)))   # VanRaden, by hand
check(approx(.vanraden(Xmaj, rownames(Xmaj)), Kman, 1e-10),
      "alt-major panel: .vanraden matches hand-computed VanRaden")
# ... and the regression itself: A.mat, fed the same {0,1,2} matrix, does NOT. It reads
# every alt-major marker as MAF < 0 and drops it, so on this panel it drops ALL of them
# and returns an entirely NaN GRM. (On real panels it drops ~half and returns a GRM that
# merely looks plausible -- which is why this went unnoticed.)
Abad <- suppressWarnings(tryCatch(rrBLUP::A.mat(Xmaj), error = function(e) NULL))
check(is.null(Abad) || anyNA(Abad) || !isTRUE(approx(Abad, Kman, 1e-6)),
      "regression: rrBLUP::A.mat on {0,1,2} dosages does not give VanRaden (the bug)")
check(is.null(Abad) || all(is.na(Abad)),
      "regression: A.mat drops every alt-major marker -> all-NaN GRM")

# Oracle: the {0,1,2} assumption now lives in .vanraden's `p`, so it must be guarded.
# Handing it rrBLUP's {-1,0,1} coding must FAIL LOUDLY, not silently return a wrong p.
check(inherits(try(.vanraden(Xp - 1, rownames(Xp)), silent = TRUE), "try-error"),
      "guard: {-1,0,1} coding is rejected, not silently mis-centred")

# Oracle: QC decisions depend on the population they are estimated from. QC on a
# 5-accession subset keeps a different marker set than QC on the population.
qc_pop <- colnames(.qc_markers(Xp, kcfg("vanRaden_single")))
qc_sub <- colnames(.qc_markers(Xp[1:5, , drop = FALSE], kcfg("vanRaden_single")))
check(!identical(qc_pop, qc_sub), "QC on 5 accessions != QC on the 200-accession population")

# ===========================================================================
cat(".bridge_accessions + em_combine stitching\n")
set.seed(11)
mkpanel <- function(rows, mcols) {
  m <- matrix(sample(0:2, length(rows) * length(mcols), replace = TRUE),
              nrow = length(rows), dimnames = list(rows, mcols))
  m
}
# Two DIFFERENT protocols (disjoint marker sets). need = {n1, n2}: n1 only in panel A,
# n2 only in panel B. b1..b3 are genotyped in BOTH -> the bridge.
bridge3 <- paste0("b", 1:3)
pA <- mkpanel(c("n1", bridge3, paste0("xa", 1:8)), paste0("mA", 1:120))   # protocol A markers
pB <- mkpanel(c("n2", bridge3, paste0("xb", 1:8)), paste0("mB", 1:120))   # protocol B markers
need2 <- c("n1", "n2")

# Oracle: bridge = accessions in >= 2 panels; disjoint panels -> none.
check(setequal(.bridge_accessions(list(pA, pB)), bridge3),
      ".bridge_accessions finds accessions in >= 2 panels")
check(length(.bridge_accessions(list(pA, mkpanel(paste0("z", 1:6), paste0("mC", 1:120))))) == 0,
      "disjoint panels -> no bridge accessions")

# Oracle (the crux): with need alone, panel A has a single needed row (n1) -> .vanraden NULL,
# uncombinable. Keeping the bridge accessions makes the panel a valid GRM.
check(is.null(.vanraden(.qc_markers(pA, kcfg("em_combine")), need2)),
      "need-only: panel with one needed accession is uncombinable (NULL)")
check(!is.null(.vanraden(.qc_markers(pA, kcfg("em_combine")), union(need2, bridge3))),
      "with bridges: the same panel yields a valid GRM")

# Oracle: em_combine end to end -- feasible BECAUSE of the bridges, and the bridges are
# eliminated so the result is exactly need x need. (Runs the real covariance_combiner.)
Kem <- build_kernel(kcfg("em_combine", ridge = 0), list(A = pA, B = pB), need2)
check(setequal(rownames(Kem), need2), "em_combine result is exactly need x need (bridges dropped)")
check(isSymmetric(unname(Kem)), "em_combine stitched GRM is symmetric")
check(all(is.finite(Kem)), "em_combine stitched GRM is finite (cross-panel block estimated)")

# Oracle (the sweep-oracle regression): a marker-POOR panel (<50 markers -> .qc_markers throws
# too_few_markers) must be DROPPED, not sink the whole combine. Rich panel survives -> valid K.
pRich <- mkpanel(c("n1", "n2", paste0("r", 1:8)), paste0("mR", 1:120))   # passes QC
pPoor <- mkpanel(c("n1", "n2", paste0("r", 1:8)), paste0("mP", 1:10))    # 10 markers -> QC fails
Kdrop <- build_kernel(kcfg("em_combine", ridge = 0), list(rich = pRich, poor = pPoor), need2)
check(setequal(rownames(Kdrop), need2) && all(is.finite(Kdrop)),
      "em_combine drops a marker-poor panel instead of failing (too_few_markers no longer sinks it)")

# Oracle: a SINGULAR partial covariance must not crash the combine. Without the pre-combine
# ridge this raises "Lapack routine dgesv: system is exactly singular"; with it, the combine
# succeeds.
#
# The duplicated rows must be ones that REACH the combiner -- inside `keep` (need + bridges),
# since .vanraden only builds over those. Duplicating a row outside `keep` tests nothing.
# Bridges are also the realistic case: lines genotyped in more than one project are exactly
# where re-called duplicates occur.
pClA <- pA; pClB <- pB
pClA["b1", ] <- pClA["b2", ]        # duplicate BRIDGE rows -> rank-deficient partial GRM
pClB["b1", ] <- pClB["b2", ]
outcome <- tryCatch(
  { Kc <- build_kernel(kcfg("em_combine", ridge = 0), list(A = pClA, B = pClB), need2)
    if (is.matrix(Kc) && all(is.finite(Kc))) "matrix" else "bad-matrix" },
  optimizer_infeasible = function(e) "infeasible",
  error = function(e) paste("UNCAUGHT:", conditionMessage(e)))
check(outcome %in% c("matrix", "infeasible"),
      paste("em_combine survives a singular partial covariance (got:", outcome, ")"))

# ===========================================================================
cat(".best_panel (focal coverage, not just max union)\n")
# The suspect configuration, reproduced: a BIG panel carrying training lines only, and a
# SMALL panel carrying the focal lines. Old rule maximized coverage of union(train, focal)
# and picked the big one -> test_in = 0 -> insufficient_geno_overlap flagged suspect.
tr_only <- mkpanel(paste0("t", 1:60), paste0("mT", 1:200))               # 60 training lines
fo_some <- mkpanel(c(paste0("t", 1:25), paste0("f", 1:10)), paste0("mF", 1:200))  # 25 train + 10 focal
focal10 <- paste0("f", 1:10)
need_bp <- union(paste0("t", 1:60), focal10)
pl <- list(train_only = tr_only, focal_covering = fo_some)

check(identical(.best_panel(pl, need_bp), pl$train_only),
      ".best_panel with focal = NULL keeps the old max-union behaviour (train_only wins)")
check(identical(.best_panel(pl, need_bp, focal10, min_test = 5, min_train = 20),
                pl$focal_covering),
      ".best_panel prefers the focal-covering panel over a bigger training-only panel")

# No behaviour change when both panels clear the guards: max union still wins.
both_a <- mkpanel(c(paste0("t", 1:50), paste0("f", 1:10)), paste0("mA", 1:200))
both_b <- mkpanel(c(paste0("t", 1:25), paste0("f", 1:10)), paste0("mB", 1:200))
check(identical(.best_panel(list(a = both_a, b = both_b), need_bp, focal10, 5, 20), both_a),
      ".best_panel still takes max union coverage when both panels are feasible")

# Nothing qualifies -> fall back to the best focal coverage (without focal lines there is
# nothing to predict, so that is the failure worth avoiding).
thin_a <- mkpanel(c(paste0("t", 1:8), paste0("f", 1:2)), paste0("mC", 1:200))
thin_b <- mkpanel(c(paste0("t", 1:5), paste0("f", 1:7)), paste0("mD", 1:200))
check(identical(.best_panel(list(a = thin_a, b = thin_b), need_bp, focal10, 5, 20), thin_b),
      ".best_panel falls back to max FOCAL coverage when no panel clears the guards")
check(identical(.best_panel(list(only = thin_a), need_bp, focal10, 5, 20), thin_a),
      ".best_panel with one panel returns it regardless")

# ===========================================================================
cat(".effective_n / .center_dfs (em_combine degrees of freedom)\n")
# Oracle: df is a partial's relative WEIGHT in the Wishart-EM likelihood. .effective_n is the
# Galwey (2009) effective number of independent samples; it must not move when the matrix is
# rescaled, because the combiner is fed unit-mean-diagonal partials.
In <- diag(20)
check(abs(.effective_n(In) - 20) < 1e-8, ".effective_n(I_n) == n (a fully independent panel)")
check(abs(.effective_n(In) - .effective_n(2 * In)) < 1e-8, ".effective_n is scale-invariant")

set.seed(101)
# Rank-deficient panel: fewer markers than accessions -> effective_n well below nrow. This is
# the case the accession count overstates, and the reason for the switch.
Mrd <- matrix(rbinom(60 * 15, 2, 0.3), 60, 15)
prd <- colMeans(Mrd) / 2; Wrd <- sweep(Mrd, 2, 2 * prd, "-")
Grd <- tcrossprod(Wrd) / (2 * sum(prd * (1 - prd))); Grd <- Grd / mean(diag(Grd))
check(.effective_n(Grd) < 0.5 * nrow(Grd),
      ".effective_n is far below nrow on a rank-deficient panel (nrow overstates it)")
check(.effective_n(Grd) >= 1, ".effective_n is bounded below by 1")

# Oracle: .center_dfs keeps the ORDERING but re-centres and caps the spread -- that cap is
# what stops a big panel outweighting a small one 10:1.
d <- .center_dfs(c(2000, 200, 600), 60, 15)
check(identical(order(d), order(c(2000, 200, 600))), ".center_dfs preserves the ordering of m_eff")
check(abs(mean(d) - 60) < 1e-8, ".center_dfs re-centres on mean_df")
check(stats::sd(d) <= 15 + 1e-8, ".center_dfs caps the spread at sd_df")
check(max(d) / min(d) < 2, ".center_dfs keeps the weight ratio well under the 10:1 nrow gives")
check(all(.center_dfs(c(5, 5, 5), 60, 15) == 60), "identical m_eff -> all dfs == mean_df")
check(all(.center_dfs(c(1, 1e6), 60, 15) >= 1), ".center_dfs floors dfs at 1")

# Oracle (the deviation from the sibling): df is measured BEFORE the ridge, so the searchable
# kernel.ridge cannot move the panel weights. Same two panels, ridge at both ends of the
# config-space range -> identical dfs.
dfs_at <- function(r) {
  gl <- list(A = .vanraden(.qc_markers(pA, kcfg("em_combine")), union(need2, bridge3)),
             B = .vanraden(.qc_markers(pB, kcfg("em_combine")), union(need2, bridge3)))
  unridged <- lapply(gl, function(g) g / mean(diag(g)))
  .center_dfs(vapply(unridged, .effective_n, numeric(1)), 60, 15)
}
check(isTRUE(all.equal(dfs_at(1e-5), dfs_at(1e-2))),
      "em_combine dfs are independent of kernel.ridge (measured pre-ridge)")

# Oracle: equal dfs are an unweighted average, so the combined result must be identical
# whatever the common value is -- only the RATIO between partials reaches the M-step.
eqA <- build_kernel(kcfg("em_combine", ridge = 0), list(A = pA, B = pB), need2)
check(all(is.finite(eqA)) && setequal(rownames(eqA), need2),
      "em_combine still produces a finite need x need kernel under effective-n dfs")

# ===========================================================================
cat(".group_by_panel / .prune_redundant\n")
gset <- optimizer_settings()
pmk <- function(rows, cols) matrix(1, length(rows), length(cols),
                                   dimnames = list(rows, cols))
# Same panel (v1 vs v2 reference build): v2 contains ALL of v1's markers.
v1 <- pmk(paste0("a", 1:10), paste0("m", 1:100))
v2 <- pmk(paste0("a", 1:10), paste0("m", 1:140))       # 100/100 = 100% containment
# A different panel: shares only 40 of the smaller 100 markers.
other <- pmk(paste0("b", 1:10), c(paste0("m", 1:40), paste0("x", 1:60)))
g <- .group_by_panel(list(v1 = v1, v2 = v2, other = other), gset)
check(any(vapply(g, function(x) setequal(x, c("v1", "v2")), logical(1))),
      "100% marker containment -> same protocol group")
check(any(vapply(g, function(x) identical(x, "other"), logical(1))),
      "40% marker containment -> its own group")

# Oracle: identical accession sets -> keep the richer marker build (the v1/v2 case).
pr <- .prune_redundant(list(v1 = v1, v2 = v2), gset)
check(identical(names(pr), "v2"), "identical accessions -> project with more markers survives")
# Oracle: >=90% shared accessions but not identical -> keep the one with MORE accessions,
# even though it has fewer markers.
big   <- pmk(paste0("a", 1:100), paste0("m", 1:100))
small <- pmk(paste0("a", 1:95),  paste0("m", 1:140))   # 95/95 = 100% of the smaller set
pr2 <- .prune_redundant(list(big = big, small = small), gset)
check(identical(names(pr2), "big"), "92% shared accessions -> project with more accessions survives")
# Oracle: distinct accession sets are NOT redundant.
indep <- pmk(paste0("z", 1:50), paste0("m", 1:100))
check(length(.prune_redundant(list(big = big, indep = indep), gset)) == 2,
      "disjoint accessions -> both projects kept")

# ===========================================================================
cat(".onehop_filter (subtask C: focal_plus_onehop's hop)\n")
# Three panels: `focal` genotypes the focal trial's lines; `near` shares 3 lines with it
# (a bridge); `far` shares none. Marker sets are irrelevant here -- grouping already
# happened -- so give each its own.
fmat <- function(rows) matrix(1, length(rows), 5, dimnames = list(rows, paste0("m", 1:5)))
focal <- fmat(c("f1", "f2", "f3", "s1", "s2", "s3"))
near  <- fmat(c("s1", "s2", "s3", "n1", "n2"))          # 3 accessions shared with focal
far   <- fmat(c("x1", "x2", "x3"))                       # none shared
panels <- list(focal = focal, near = near, far = far)
test_acc <- c("f1", "f2", "f3")

# Oracle: the focal panel is always kept; `far` bridges nothing so it is never admitted;
# `near` is admitted exactly while min_bridge <= 3 (it shares 3).
for (mb in 1:5) {
  keptn <- names(.onehop_filter(panels, test_acc, mb))
  check("focal" %in% keptn, sprintf("one-hop: the focal panel is kept (min_bridge=%d)", mb))
  check(!("far" %in% keptn), sprintf("one-hop: an unbridged panel is dropped (min_bridge=%d)", mb))
  check(("near" %in% keptn) == (mb <= 3),
        sprintf("one-hop: a 3-bridge panel is admitted iff min_bridge<=3 (min_bridge=%d)", mb))
}
# Oracle: nested -- raising min_bridge can only shrink the admitted set.
sets <- lapply(1:5, function(mb) names(.onehop_filter(panels, test_acc, mb)))
check(all(vapply(2:5, function(i) all(sets[[i]] %in% sets[[i - 1]]), logical(1))),
      "one-hop: the admitted set is nested in min_bridge")
# Oracle: this is what distinguishes focal_plus_onehop from all_projects -- at min_bridge=1
# it still drops `far`, which all_projects would keep.
check(length(sets[[1]]) < length(panels),
      "one-hop at min_bridge=1 is still narrower than all_projects")
# Oracle: no panel genotypes the focal lines -> nothing to hop from -> pass through
# unchanged (run_pipeline's overlap check judges the trial, not this filter).
check(length(.onehop_filter(panels, c("zz1", "zz2"), 2)) == length(panels),
      "one-hop: no focal coverage -> list returned unchanged")
# Oracle: a single panel is returned as-is regardless of threshold.
check(length(.onehop_filter(panels["focal"], test_acc, 5)) == 1,
      "one-hop: a lone panel is never filtered away")

# ===========================================================================
cat(".loo_lambda / .ridge_blup (subtask E: how lambda is chosen)\n")
set.seed(11)
nL <- 25
ZL <- matrix(rnorm(nL * 40), nL, 40)
KL <- tcrossprod(ZL) / 40
diag(KL) <- diag(KL) + 1e-6
yL  <- as.numeric(KL %*% rnorm(nL)) + rnorm(nL, sd = 0.5)
ycL <- yL - mean(yL)

# Oracle: the PRESS identity the implementation relies on. Brute force -- actually refit
# on n-1 points and predict the held-out one -- must agree with r_i/(1-h_ii) to numerical
# precision. This is the assumption that makes the grid search cheap; if it is wrong, the
# selected lambda is wrong in a way nothing downstream would reveal.
loo_brute <- function(K, y, lam) {
  vapply(seq_len(nrow(K)), function(i) {
    Kii <- K[-i, -i, drop = FALSE]
    y[i] - as.numeric(K[i, -i, drop = FALSE] %*% solve(Kii + lam * diag(nrow(Kii)), y[-i]))
  }, numeric(1))
}
loo_press <- function(K, y, lam) {
  H <- K %*% solve(K + lam * diag(nrow(K)))
  (y - as.numeric(H %*% y)) / (1 - diag(H))
}
check(approx(loo_brute(KL, ycL, 0.3), loo_press(KL, ycL, 0.3), tol = 1e-8),
      "PRESS identity: shortcut LOO residuals equal an explicit leave-one-out refit")

# Oracle: .loo_lambda returns the grid point minimizing brute-force LOO MSE.
gridL <- 10^seq(-4, 4, length.out = 30)
mseL  <- vapply(gridL, function(l) mean(loo_brute(KL, ycL, l)^2), numeric(1))
check(approx(.loo_lambda(KL, ycL), gridL[which.min(mseL)], tol = 1e-9),
      ".loo_lambda picks the grid lambda with the lowest leave-one-out MSE")
check(.loo_lambda(KL, ycL) %in% gridL, ".loo_lambda returns a point of its own grid")

# Oracle: .ridge_blup solves u = K(K + lambda I)^-1 (y - mean y) -- computed here the
# long way -- and reports mean(y) as the intercept.
idsL <- paste0("g", seq_len(nL))
rownames(KL) <- colnames(KL) <- idsL
fitL <- .ridge_blup(yL, KL, 0.7, idsL)
uL   <- as.numeric(KL %*% solve(KL + 0.7 * diag(nL), ycL))
check(approx(as.numeric(fitL$u), uL) && identical(names(fitL$u), idsL),
      ".ridge_blup = K(K+lambda I)^-1 (y - mean y), named by training id")
check(approx(fitL$mu, mean(yL)), ".ridge_blup intercept is the training mean")
# Oracle: shrinkage is monotone in lambda -- more ridge, smaller effects.
norms <- vapply(c(0.01, 0.1, 1, 10, 100),
                function(l) sqrt(sum(.ridge_blup(yL, KL, l, idsL)$u^2)), numeric(1))
check(all(diff(norms) < 0), ".ridge_blup shrinks monotonically as lambda grows")

# Oracle: train_model dispatches on lambda_select. "fixed" must use the given lambda
# (identical to calling .ridge_blup with it) and REML must NOT (it estimates its own),
# so the two branches are genuinely different fits rather than one code path.
mcfg <- function(rule, lam = 1) list(model.method = "gblup_rrblup",
                                     model.lambda_select = rule, model.lambda_fixed = lam,
                                     model.include_E = "no")
yn    <- stats::setNames(yL, idsL)
f_fix <- train_model(mcfg("fixed", 0.7), yn, KL, idsL, character(0), NULL, NULL)
check(approx(as.numeric(f_fix$u), uL), "train_model[fixed] fits at exactly lambda_fixed")
f_rem <- train_model(mcfg("reml"), yn, KL, idsL, character(0), NULL, NULL)
f_loo <- train_model(mcfg("loo"),  yn, KL, idsL, character(0), NULL, NULL)
check(!approx(f_rem$lambda, f_fix$lambda, tol = 1e-3) ||
      !approx(as.numeric(f_rem$u), uL, tol = 1e-6),
      "train_model[reml] estimates its own lambda rather than reusing lambda_fixed")
check(approx(f_loo$lambda, .loo_lambda(KL, ycL)), "train_model[loo] uses the LOO-selected lambda")
check(all(vapply(list(f_fix, f_rem, f_loo),
                 function(f) identical(f$kind, "gblup") && length(f$u) == nL &&
                             is.finite(f$mu) && !is.matrix(f$mu), logical(1))),
      "every lambda rule returns the same shape (kind/u/mu) for predict_test")
# Oracle: an absent or NA rule falls back to REML rather than erroring (crossover can
# hand train_model a config whose lambda_select came from a block that lacked it).
f_na <- train_model(list(model.method = "gblup_rrblup", model.lambda_select = NA,
                         model.include_E = "no"), yn, KL, idsL, character(0), NULL, NULL)
check(approx(f_na$lambda, f_rem$lambda, tol = 1e-8), "an NA lambda_select falls back to REML")

# ===========================================================================
cat(".trial_similarity: an unknown year/coordinate is UNRELATED, not NA\n")
# Oracle: 28% of the T3 catalogue has year = NA. Before 2026-07-31 such a candidate made d2
# NA -> similarity NA -> weight NA -> weighted.mean NA -> an NA TARGET for every accession
# phenotyped there -> non-finite predictions -> a scored-as-empty test set, several subtasks
# downstream. A candidate we cannot place must score 0, not poison the arithmetic.
sim_cat <- tibble::tibble(
  study_db_id = c("A", "B", "C"),
  latitude    = c(42.5, 42.6, NA_real_),
  longitude   = c(-76.5, -76.4, -76.3),
  year        = c(2020L, NA_integer_, 2021L))
sim_settings <- modifyList(optimizer_settings(), list(simulate = TRUE))
local({
  # Stub the catalogue: .trial_similarity's environmental branch reads only these columns.
  tc <- trial_catalog
  assign("trial_catalog", function(conn, settings) sim_cat, envir = globalenv())
  on.exit(assign("trial_catalog", tc, envir = globalenv()), add = TRUE)
  focal <- list(id = "F", accessions = character(), lat = 42.5, long = -76.5, year = 2020)
  s <- .trial_similarity(focal, c("A", "B", "C"), conn = NULL, sim_settings, "environmental")
  check(all(is.finite(s)), ".trial_similarity: no NA even when a candidate lacks year/coords")
  check(s[["A"]] > 0.9, "a candidate at the focal position scores ~1")
  check(s[["B"]] == 0 && s[["C"]] == 0, "candidates with unknown year / coordinates score 0")
  # The consequence that mattered: weights derived from this are usable.
  w <- exp(log(pmax(s, 1e-6)) / 1.5)
  check(all(is.finite(w)) && is.finite(stats::weighted.mean(c(1, 2, 3), w)),
        "weights from the fixed similarity give a finite weighted.mean")

  # End to end: build_targets under env_gaussian weighting, with an unplaceable trial in the
  # training set. This is the exact path that produced trial 11033's non-finite predictions.
  tobs <- tibble::tibble(
    study_id       = rep(c("A", "B", "C"), each = 4),
    germplasm_name = rep(paste0("g", 1:4), times = 3),
    value          = c(5000, 5200, 4800, 5100, 5050, 5250, 4850, 5150, 4990, 5190, 4790, 5090),
    rep = "1", block = "1", unit_id = NA_character_, col = NA_character_, row = NA_character_)
  cfg_ge <- modifyList(canary_config(),
    list(pheno_prep.method = "raw_mean", pheno_prep.ge_weighting = "env_gaussian",
         pheno_prep.ge_bandwidth = 1.5))
  tg <- build_targets(cfg_ge, tobs, focal, conn = NULL, sim_settings)
  check(length(tg) == 4 && all(is.finite(tg)),
        "build_targets: env_gaussian over an unplaceable trial yields FINITE targets")
})

# ===========================================================================
cat("build_targets + .blue_per_trial\n")
mkobs <- function(study, germ, value) tibble::tibble(
  study_id = study, germplasm_name = germ, value = value, rep = "r1", block = "b1")
# raw_mean. Oracle: mean of per-trial means. g1: trial A mean(10,12)=11, trial B=20 -> 15.5; g2=8.
obs_rm <- dplyr::bind_rows(
  mkobs("A", "g1", 10), mkobs("A", "g1", 12), mkobs("A", "g2", 8), mkobs("B", "g1", 20))
t_rm <- build_targets(pcfg("raw_mean"), obs_rm, no_trial, NULL, list())
check(approx(t_rm["g1"], 15.5) && approx(t_rm["g2"], 8), "raw_mean = mean of per-trial means")
# trial_center. Oracle: removing trial means leaves only a global shift, so adding
# a constant to one trial shifts ALL targets equally (variance of the change = 0).
t1 <- build_targets(pcfg("trial_center"), obs_rm, no_trial, NULL, list())
obs_shift <- obs_rm |> dplyr::mutate(value = value + ifelse(study_id == "B", 100, 0))
t2 <- build_targets(pcfg("trial_center"), obs_shift, no_trial, NULL, list())
common <- intersect(names(t1), names(t2))
check(approx(stats::var(t2[common] - t1[common]), 0, 1e-8), "trial_center removes trial effect")
# blue_lm. Oracle: single-level design -> germ coefficients are the group means.
obs_blue <- dplyr::bind_rows(
  mkobs("A", "g1", 10), mkobs("A", "g1", 10), mkobs("A", "g2", 20), mkobs("A", "g2", 20))
t_bl <- build_targets(pcfg("blue_lm"), obs_blue, no_trial, NULL, list())
check(approx(t_bl["g1"], 10) && approx(t_bl["g2"], 20), "blue_lm recovers group means")
# blue_lm fallback. Oracle: a one-germplasm trial -> mean. g1 = mean(10,14) = 12.
t_fb <- build_targets(pcfg("blue_lm"),
  dplyr::bind_rows(mkobs("A", "g1", 10), mkobs("A", "g1", 14)), no_trial, NULL, list())
check(approx(t_fb["g1"], 12), "blue_lm single-line fallback = mean")
# outlier removal. Oracle: with enough normal points, z_thr=3 flags the 500; g1 -> 10.
obs_out <- dplyr::bind_rows(
  mkobs("A", "g1", rep(10, 20)), mkobs("A", "g1", 500), mkobs("A", "g2", rep(20, 5)))
t_o3  <- build_targets(pcfg("blue_lm", z_thr = 3),   obs_out, no_trial, NULL, list())
t_oInf <- build_targets(pcfg("blue_lm", z_thr = Inf), obs_out, no_trial, NULL, list())
check(approx(t_o3["g1"], 10, 1e-6), "z_thr=3 removes the outlier (g1 -> 10)")
check(t_oInf["g1"] > 30, "z_thr=Inf keeps the outlier (g1 inflated)")
# two_stage_blup. Oracle: BLUP shrinks the spread below the raw germ effects.
if (requireNamespace("lme4", quietly = TRUE)) {
  obs_ts <- dplyr::bind_rows(
    mkobs("A", "g1", 4), mkobs("A", "g2", 11), mkobs("A", "g3", 14),
    mkobs("B", "g1", 6), mkobs("B", "g2", 9),  mkobs("B", "g3", 16))
  t_ts <- suppressWarnings(build_targets(pcfg("two_stage_blup"), obs_ts, no_trial, NULL, list()))
  raw  <- c(g1 = 5, g2 = 10, g3 = 15)             # per-germ means across the two trials
  check(stats::sd(t_ts - mean(t_ts)) < stats::sd(raw - mean(raw)),
        "two_stage_blup shrinks spread below raw effects")
} else cat("  (skipped two_stage_blup: lme4 not installed)\n")
# standardize. Oracle: per-trial z-scores have mean 0 and sd 1.
obs_sd <- mkobs("A", paste0("g", 1:5), 1:5)
t_sd <- build_targets(pcfg("raw_mean", standardize = "per_trial_z"), obs_sd, no_trial, NULL, list())
check(approx(mean(t_sd), 0, 1e-9) && approx(stats::sd(t_sd), 1, 1e-3), "per_trial_z -> mean 0, sd 1")

# ===========================================================================
cat("predict_test\n")
# Conditional expectation. Oracle: a test line that is a clone of training line t1
# (its K row to the training set equals t1's) gets prediction mu + u[t1].
Kp <- matrix(0, 3, 3, dimnames = list(c("t1", "t2", "x"), c("t1", "t2", "x")))
diag(Kp) <- 2; Kp["t1", "t2"] <- Kp["t2", "t1"] <- 0.5
Kp["x", "t1"] <- Kp["t1", "x"] <- 2; Kp["x", "t2"] <- Kp["t2", "x"] <- 0.5   # x clone of t1
fit <- list(kind = "gblup", u = c(t1 = 1, t2 = -0.5), mu = 10)
pc <- predict_test(prcfg("cond_expectation"), fit, Kp, c("t1", "t2"), "x",
                   targets = c(t1 = 0, t2 = 0), scheme = "CV00")
check(approx(pc["x"], 10 + 1), "cond_expectation clone -> mu + u[t1]")
# Blending. Oracle: under CV0, blend = w*BLUE + (1-w)*genomic; under CV00 it's a no-op.
Kd <- matrix(0, 4, 4, dimnames = list(c("a", "t1", "t2", "t3"), c("a", "t1", "t2", "t3")))
fitd <- list(kind = "gblup", u = c(a = 0, t1 = 0, t2 = 0, t3 = 0), mu = 10)
tg <- c(a = 100, t1 = 0, t2 = 0, t3 = 0)
p0 <- predict_test(prcfg("direct_blup", blend_obs_w = 0.5), fitd, Kd,
                   c("t1", "t2", "t3"), "a", tg, "CV0")
p00 <- predict_test(prcfg("direct_blup", blend_obs_w = 0.5), fitd, Kd,
                    c("t1", "t2", "t3"), "a", tg, "CV00")
check(approx(p0["a"], 0.5 * 100 + 0.5 * 10), "CV0 blends observed BLUE")
check(approx(p00["a"], 10), "CV00 blend is a no-op")
# Fallback. Oracle: too little genotyped overlap -> every prediction is mean(targets).
pf <- predict_test(prcfg("cond_expectation", min_overlap = 3), fit, Kp,
                   c("t1", "t2"), "x", targets = c(t1 = 4, t2 = 8), scheme = "CV0")
check(approx(pf["x"], 6), "fallback below min_overlap -> mean(targets)")

# direct_blup when the fit has NO random effect for the focal lines -- the rrBLUP backbone
# fitted on training lines only, which is the norm under CV00. Oracle: it must predict through
# the kernel, giving exactly the conditional expectation, NOT collapse to the intercept.
# `fit` above carries u for t1/t2 only, so fit$u["x"] is NA.
pd <- suppressMessages(predict_test(prcfg("direct_blup"), fit, Kp,
                                    c("t1", "t2"), "x", targets = c(t1 = 0, t2 = 0), "CV00"))
check(approx(pd["x"], pc["x"]),
      "direct_blup with no test-line effect falls back to the conditional expectation")
check(!approx(pd["x"], fit$mu),
      "...and therefore does NOT collapse to the intercept (the CV00 `constant` bug)")
# Oracle: when the fit DOES carry the focal lines (the sommer G+E path), direct_blup uses them
# and is not redirected -- the fallback must be conditional, not unconditional.
fit_all <- list(kind = "gblup", u = c(t1 = 1, t2 = -0.5, x = 3), mu = 10)
pk <- predict_test(prcfg("direct_blup"), fit_all, Kp, c("t1", "t2"), "x",
                   targets = c(t1 = 0, t2 = 0), scheme = "CV00")
check(approx(pk["x"], 13), "direct_blup uses the model's own effect when it exists")

# ===========================================================================
cat("mask_cv\n")
mo <- dplyr::bind_rows(
  mkobs("A", "f1", 1), mkobs("A", "f2", 2), mkobs("A", "x1", 3), mkobs("B", "x2", 4))
focal <- c("f1", "f2")
m0  <- mask_cv(mo, focal, "CV0")
m00 <- mask_cv(mo, focal, "CV00")
check(nrow(m0) == 4, "CV0 keeps all rows")
check(nrow(m00) == 2 && !any(m00$germplasm_name %in% focal), "CV00 drops focal-accession rows")
check(setequal(setdiff(m0$germplasm_name, m00$germplasm_name), focal),
      "CV0 minus CV00 = exactly the focal accessions")

# ===========================================================================
cat("evals store: domain attributes persisted + surrogate trains per-domain\n")
# Oracle: the store keeps each trial's target-domain attributes so a run can learn
# only from its own domain's slice of this shared archive.
dbp <- tempfile(fileext = ".sqlite"); con <- open_store(dbp)
flds <- DBI::dbListFields(con, "evals")
check(all(c("study_name", "program_name", "location_name", "year") %in% flds),
      "open_store: evals table has the four domain columns")
cfgA <- sample_config()
store_eval(con, cfgA, "101", "CV0", 0.40, 50L, "ok",
           study_name = "T_Cornell_24", program_name = "Cornell", location_name = "Ithaca", year = 2024)
store_eval(con, cfgA, "102", "CV0", 0.20, 50L, "ok",
           study_name = "T_OSU_24", program_name = "OSU", location_name = "Columbus", year = 2024)
ev <- read_evals(con)
check(nrow(ev) == 2 && setequal(ev$program_name, c("Cornell", "OSU")),
      "store_eval persists program_name for each row")
# Oracle: filter_evals_to_domain keeps only the matching program; NULL domain keeps all.
cornell <- filter_evals_to_domain(ev, list(programs = "Cornell"))
check(nrow(cornell) == 1 && cornell$program_name == "Cornell",
      "filter_evals_to_domain restricts to the target program")
check(nrow(filter_evals_to_domain(ev, NULL)) == 2,
      "filter_evals_to_domain: NULL domain is no constraint")
check(nrow(filter_evals_to_domain(ev, list(programs = "Cornell", years = 2023))) == 0,
      "filter_evals_to_domain: a year mismatch excludes the row")
# Oracle: an unknown/NA-attribute row (e.g. simulated or pre-migration) is out of
# domain once a constraint is set, but retained when the domain is NULL.
store_eval(con, sample_config(), "sim_9", "CV0", 0.5, 10L, "ok")   # all attrs NA
ev2 <- read_evals(con)
check(nrow(filter_evals_to_domain(ev2, list(programs = "Cornell"))) == 1,
      "filter_evals_to_domain: NA-attribute rows are out-of-domain under a constraint")
check(nrow(filter_evals_to_domain(ev2, NULL)) == 3,
      "filter_evals_to_domain: NA-attribute rows kept when domain is NULL")
close_store(con); unlink(dbp)
# Oracle: a config seen only OUT of the current domain still counts as untried, so
# choose_config re-offers a seed rather than treating it as done.
dbp2 <- tempfile(fileext = ".sqlite"); con2 <- open_store(dbp2)
# Pin the scheme rather than inheriting it: optimizer_settings() layers in settings.local.R,
# so a machine optimizing CV00 would filter out this CV0 row and the seed would look untried
# for the wrong reason. The seed must come from the SAME scheme for the hashes to correspond.
TEST_SCHEME <- "CV0"
seed1 <- seed_configs(TEST_SCHEME)[[1]]
# Stamp the current build: these oracles are about the DOMAIN/SCHEME filters, and an
# unstamped row would additionally be retired by filter_evals_to_build, making the seed look
# untried for the wrong reason.
store_eval(con2, seed1, "900", TEST_SCHEME, 0.3, 40L, "ok",
           study_name = "elsewhere", program_name = "OSU", location_name = "X", year = 2024,
           build = OPTIMIZER_BUILD)
# simulate = FALSE so the (real-data) target-domain filter actually engages.
st <- modifyList(optimizer_settings(),
                 list(simulate = FALSE, optimize_scheme = TEST_SCHEME,
                      target_domain = list(programs = "Cornell")))
pick <- choose_config(con2, st)
check(identical(config_hash(pick$cfg), config_hash(seed1)) && grepl("^seed:", pick$source),
      "choose_config: an out-of-domain seed eval does not count as done in this domain")
# Control: with NO domain restriction, that same seed eval DOES count as done, so
# choose_config moves past seed1 to a later unevaluated seed.
st0  <- modifyList(optimizer_settings(),
                   list(simulate = FALSE, optimize_scheme = TEST_SCHEME, target_domain = NULL))
pick0 <- choose_config(con2, st0)
check(!identical(config_hash(pick0$cfg), config_hash(seed1)),
      "choose_config: with no domain, the recorded seed is not re-offered")
close_store(con2); unlink(dbp2)

# ===========================================================================
cat("aggregate_scores: Fisher-z pooling weighted by n_test\n")
# Oracle: a score is a correlation, so a 5-accession one is nearly pure noise while a
# 200-accession one is informative. Pooling must reflect that -- an unweighted mean of
# r = 0.9 (n=5) and r = 0.3 (n=200) is 0.6, which is badly wrong.
mkev <- function(scores, ns) tibble::tibble(
  config_hash = "H", config_json = "{}", score = scores, n_test = as.integer(ns))
a1 <- aggregate_scores(mkev(c(0.9, 0.3), c(5, 200)), min_n_test = 0L)
check(abs(a1$unweighted - 0.6) < 1e-8, "aggregate_scores keeps the old unweighted mean alongside")
check(a1$mean_score < 0.35 && a1$mean_score > 0.29,
      "Fisher-z weighting: n=5 r=0.9 barely moves an n=200 r=0.3 (pools near 0.3)")

# Equal n must reproduce the unweighted mean (in z space), so the change is a no-op when
# every evaluation is equally precise.
a2 <- aggregate_scores(mkev(c(0.4, 0.2), c(50, 50)), min_n_test = 0L)
check(abs(a2$mean_score - tanh(mean(atanh(c(0.4, 0.2))))) < 1e-8,
      "equal n_test -> plain mean of the z-transformed scores")

# atanh(+-1) is Inf; the clamp must keep a perfect correlation finite.
a3 <- aggregate_scores(mkev(c(1, -1), c(40, 40)), min_n_test = 0L)
check(is.finite(a3$mean_score), "r = +-1 does not produce Inf through atanh")

# The floor drops scores too small to inform, without dropping the ROW from n.
a4 <- aggregate_scores(mkev(c(0.9, 0.3), c(5, 200)), min_n_test = 10L)
check(abs(a4$mean_score - 0.3) < 1e-8 && a4$n == 2 && a4$n_ok == 1,
      "min_n_test drops the low-n score from the pool but still counts the row")
check(is.na(aggregate_scores(mkev(c(0.9), c(5)), min_n_test = 10L)$mean_score),
      "a config with only sub-floor scores gets mean_score = NA")

# ===========================================================================
cat("filter_evals_to_build (retire rows a later build invalidated)\n")
bev <- tibble::tibble(
  config_json = c(config_to_json(modifyList(sample_config(), list(kernel.method = "em_combine"))),
                  config_to_json(modifyList(sample_config(), list(kernel.method = "em_combine"))),
                  config_to_json(modifyList(sample_config(), list(kernel.method = "two_stage_blup")))),
  build = c(NA_character_, "0.7.1", NA_character_))
ch <- list(list(build = "0.7.1", what = "em_combine df",
                affects = function(cfg) identical(cfg$kernel.method, "em_combine")))
kept <- filter_evals_to_build(bev, "0.7.1", ch)
check(nrow(kept) == 2, "filter_evals_to_build drops the pre-build row whose config matches")
check(all(!is.na(kept$build) | !grepl("em_combine", kept$config_json)),
      "the dropped row is the NA-build em_combine one, not the same-build one")
check(nrow(filter_evals_to_build(bev, "0.7.1", list())) == 3,
      "no change rules -> nothing is retired")
check(nrow(filter_evals_to_build(bev[0, ], "0.7.1", ch)) == 0,
      "filter_evals_to_build handles an empty store")
# A row stamped with a LATER build survives a rule from an earlier one.
bev2 <- tibble::tibble(config_json = bev$config_json[1], build = "0.7.2")
check(nrow(filter_evals_to_build(bev2, "0.7.2", ch)) == 1,
      "a row from a later build is not retired by an earlier build's rule")

# ===========================================================================
cat("evals store: optimizer targets ONE scheme (per-scheme slice of the archive)\n")
# Oracle: filter_evals_to_scheme keeps only the target scheme; NULL keeps all; and
# it composes with filter_evals_to_domain.
dbp3 <- tempfile(fileext = ".sqlite"); con3 <- open_store(dbp3)
cfgS <- sample_config()
store_eval(con3, cfgS, "201", "CV0",  0.4, 50L, "ok", program_name = "Cornell")
store_eval(con3, cfgS, "202", "CV00", 0.3, 50L, "ok", program_name = "Cornell")
store_eval(con3, cfgS, "203", "CV0",  0.2, 50L, "ok", program_name = "OSU")
evS <- read_evals(con3)
check(nrow(filter_evals_to_scheme(evS, "CV0")) == 2 &&
      all(filter_evals_to_scheme(evS, "CV0")$scheme == "CV0"),
      "filter_evals_to_scheme keeps only the target scheme")
check(nrow(filter_evals_to_scheme(evS, NULL)) == 3,
      "filter_evals_to_scheme: NULL scheme is no constraint")
check(nrow(filter_evals_to_domain(evS, list(programs = "Cornell")) |>
             filter_evals_to_scheme("CV0")) == 1,
      "domain + scheme filters compose (Cornell & CV0 -> 1 row)")
close_store(con3); unlink(dbp3)
# Oracle: a config seen only under the OTHER scheme counts as untried here, so
# choose_config re-offers a seed; once it exists under this scheme it is done.
dbp4 <- tempfile(fileext = ".sqlite"); con4 <- open_store(dbp4)
seedA <- seed_configs("CV00")[[1]]   # same scheme the settings below target
store_eval(con4, seedA, "900", "CV0", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)  # CV0 only
stCV00 <- modifyList(optimizer_settings(),
                     list(simulate = FALSE, target_domain = NULL, optimize_scheme = "CV00"))
pickB <- choose_config(con4, stCV00)
check(identical(config_hash(pickB$cfg), config_hash(seedA)) && grepl("^seed:", pickB$source),
      "choose_config: a seed evaluated only under CV0 is untried when optimizing CV00")
store_eval(con4, seedA, "901", "CV00", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)  # now under CV00
pickC <- choose_config(con4, stCV00)
check(!identical(config_hash(pickC$cfg), config_hash(seedA)),
      "choose_config: once the seed has a CV00 eval it is not re-offered")
close_store(con4); unlink(dbp4)
# Oracle: optimize_scheme must be a single scheme within `schemes`; a bad override
# is rejected by optimizer_settings(). Guard so a real settings.local.R is untouched.
lfs <- here::here("settings.local.R")
if (!file.exists(lfs)) {
  writeLines('settings_override <- list(optimize_scheme = "CVxx")', lfs)
  check(inherits(try(optimizer_settings(), silent = TRUE), "try-error"),
        "optimizer_settings(): an optimize_scheme outside `schemes` errors")
  invisible(file.remove(lfs))
} else message("  (a real settings.local.R exists -- skipping the optimize_scheme validation check)")

# ===========================================================================
cat(sprintf("\nTier 1 subtask tests: %d passed, %d failed\n", ok, fail))
if (fail > 0) quit(status = 1)
