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
  writeLines("settings_override <- list(dosage_budget_bytes = 7e9, max_iters = 3L)", lf)
  s_ov <- suppressWarnings(optimizer_settings())
  check(s_ov$dosage_budget_bytes == 7e9 && identical(s_ov$max_iters, 3L),
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

# no_cache() is the OTHER veto -- for when the value cannot show what went wrong. Oracle: the
# caller gets the value, unmarked, and nothing is written, so the next call retries.
nv <- cached(cs, "acc", "N1", expr = no_cache(c("p", "q")))
check(identical(nv, c("p", "q")), "no_cache(): the caller still gets the value")
check(is.null(attr(nv, ".no_cache")), "no_cache(): the marker is stripped before returning")
check(is.na(.cache_existing(cs, "acc", "N1")), "no_cache(): nothing is written")
check(identical(cached(cs, "acc", "N1", expr = c("second", "try")), c("second", "try")),
      "no_cache(): the next call re-evaluates rather than serving a cached failure")
check(!is.na(.cache_existing(cs, "acc", "N1")), "... and THAT result is cached")
# The two vetoes compose: either one alone suppresses the write.
cached(cs, "acc", "N2", valid = function(a) TRUE, expr = no_cache("x"))
check(is.na(.cache_existing(cs, "acc", "N2")), "no_cache() vetoes even when valid() says yes")
check(inherits(tryCatch(no_cache(NULL), error = function(e) e), "error"),
      "no_cache(NULL) errors rather than silently caching NULL")

# The negative cache: a structural failure recorded once, so it is not retried every run.
check(!.neg_cache_hit(cs, "unparseable", "P9"), "negative cache: absent key is not a hit")
.neg_cache_mark(cs, "unparseable", "P9", "no usable genotypes")
check(.neg_cache_hit(cs, "unparseable", "P9"), "negative cache: a marked key is a hit")
check(identical(readRDS(.cache_existing(cs, "unparseable", "P9"))$reason, "no usable genotypes"),
      "negative cache: the reason is recorded for the reader")

# .project_stat is the sanctioned read of the stat cache from outside R/genotypes.R.
check(is.null(.project_stat(cs, "P404")), ".project_stat is NULL for a project never parsed")
.cache_save(cs, "stat", "P8", list(n_samples = 10L, n_markers = 2000L))
check(identical(.project_stat(cs, "P8")$n_markers, 2000L), ".project_stat reads a stored stat")

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

# get_observations through the veto. Oracle: a FAILED fetch must leave no obs_<sid>.rds --
# storing the empty tibble it degrades to means "this trial has no phenotypes" for 30 days
# (LESSONS #7) -- while a successful but genuinely empty response is a real answer and is kept.
local({
  resp <- function(recs) list(combined_data = recs)
  stub <- function(fail) list(search = function(path, body) {
    if (fail && grepl("observations$", path)) stop("simulated HTTP 500")
    resp(list())
  })
  ob <- modifyList(cs, list(brapi_tries = 1))
  invisible(get_observations("S_FAIL", stub(TRUE), ob))
  check(is.na(.cache_existing(ob, "obs", "S_FAIL")),
        "get_observations: a failed fetch caches NOTHING (no 30-day empty)")
  invisible(get_observations("S_OK", stub(FALSE), ob))
  check(!is.na(.cache_existing(ob, "obs", "S_OK")),
        "get_observations: a successful empty response IS cached")
})
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
cat("geno_select = best_single_project (same criterion as .best_panel)\n")
# Oracle: "best single project" must mean best AT PREDICTING THE FOCAL TRIAL. A project
# holding 60 training lines and none of the focal ones cannot predict it at any size, so the
# 35-line project that carries the focal lines is the right pick even though it covers less
# of union(train, focal). Selecting on union coverage alone handed the kernel a GRM with no
# test rows and every such eval failed test_in = 0.
#
# choose_geno_sources needs the network, so its two data calls are stubbed in the global
# environment (which is where R/ was sourced, so the real function finds the stubs). This
# exercises the real selection path rather than a copy of its rule.
.real_pfa <- projects_for_accessions
.real_gpd <- get_project_dosage
gs_settings <- list(dosage_total_budget_bytes = Inf, merge_containment = 0.95,
                    redundant_acc_overlap = 0.90, min_test_acc = 5L, min_train_acc = 20L)
gs_cfg <- function(method) list(geno_select.method = method, geno_select.min_bridge = 1)
gs_train <- paste0("t", 1:60)
gs_focal <- paste0("f", 1:10)

# Disjoint marker sets, so .group_by_panel leaves each project its own panel and the name of
# the surviving group is the project id that was chosen.
gs_stub <- function(panels) {
  projects_for_accessions <<- function(need, conn, settings) names(panels)
  get_project_dosage      <<- function(pid, keep, conn, settings, marker_thin = 1L)
    panels[[as.character(pid)]]
}
gs_run <- function(method, panels) {
  gs_stub(panels)
  choose_geno_sources(gs_cfg(method), gs_train, gs_focal, NULL, gs_settings)
}

gs_panels <- list(
  big_train_only = mkpanel(gs_train, paste0("mBIG", 1:200)),            # focal 0,  train 60
  small_focal    = mkpanel(c(paste0("t", 1:25), gs_focal), paste0("mFOC", 1:200)))  # focal 10, train 25
out_bsp <- gs_run("best_single_project", gs_panels)
check(length(out_bsp) == 1 && identical(names(out_bsp), "small_focal"),
      "best_single_project takes the focal-covering project over a bigger training-only one")
check(length(intersect(rownames(out_bsp[[1]]), gs_focal)) == length(gs_focal),
      "the panel best_single_project returns carries all the focal lines")

# No project covers the focal lines: the fallback must still hand back a panel rather than
# an empty list, so run_pipeline's own overlap check reports the trial rather than a crash here.
gs_none <- list(a = mkpanel(gs_train, paste0("mA", 1:200)),
                b = mkpanel(paste0("t", 1:30), paste0("mB", 1:200)))
out_none <- gs_run("best_single_project", gs_none)
check(length(out_none) == 1 && nrow(out_none[[1]]) > 0,
      "best_single_project falls back to a panel when no project covers the focal lines")

# Anti-drift: the project-level selection and .best_panel must agree on the same input. This
# is the check that keeps the two from separating again -- they already did once.
check(identical(out_bsp[[1]],
                .best_panel(gs_panels, union(gs_train, gs_focal), gs_focal, 5L, 20L)),
      "best_single_project and .best_panel choose the same panel from the same input")

projects_for_accessions <- .real_pfa
get_project_dosage      <- .real_gpd

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
# Oracle: the grouping is a property of the PANELS, not the config, so it is memoized -- but
# the key must include marker COUNTS, since the same project served at a different density
# (.dosage_thin_plan) is a different panel and can group differently.
local({
  g1 <- .group_by_panel(list(v1 = v1, v2 = v2, other = other), gset)
  g2 <- .group_by_panel(list(v1 = v1, v2 = v2, other = other), gset)
  check(identical(g1, g2), ".group_by_panel is stable across calls (memo returns the same grouping)")
  # Same ids, a DIFFERENT marker set on v2 -- 60 markers of which only 20 are in v1, so
  # containment fails (20 < 0.95*60) and v1/v2 must now separate. If the memo keyed on ids
  # alone it would wrongly return the earlier "v1+v2 together" grouping.
  v2alt <- v2[, c(1:20, 101:140), drop = FALSE]
  g3 <- .group_by_panel(list(v1 = v1, v2 = v2alt, other = other), gset)
  check(length(g3) == 3,
        "a differently-thinned panel is re-evaluated, not served from the memo")
  # Same panels, looser threshold -> everything merges. Must miss the memo keyed at 0.95.
  loose <- modifyList(gset, list(merge_containment = 0.30))
  g4 <- .group_by_panel(list(v1 = v1, v2 = v2alt, other = other), loose)
  check(length(g4) == 1, "changing merge_containment misses the memo and regroups")
})

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
# Oracle: filter_evals_to_universe keeps only the pinned trials; NULL is no constraint.
c101 <- filter_evals_to_universe(ev, "101")
check(nrow(c101) == 1 && c101$trial_id == "101",
      "filter_evals_to_universe restricts to the pinned trials")
check(nrow(filter_evals_to_universe(ev, NULL)) == 2,
      "filter_evals_to_universe: NULL universe is no constraint")
check(nrow(filter_evals_to_universe(ev, character())) == 2,
      "filter_evals_to_universe: an empty universe is no constraint")
check(nrow(filter_evals_to_universe(ev, "999")) == 0,
      "filter_evals_to_universe: a trial outside the universe is excluded")
# Oracle: membership is decided on the id alone, so a row whose stored study_name no longer
# matches the catalogue's is still in the universe its trial belongs to.
store_eval(con, sample_config(), "101", "CV0", 0.5, 10L, "ok", study_name = "renamed_later")
ev2 <- read_evals(con)
check(nrow(filter_evals_to_universe(ev2, "101")) == 2,
      "filter_evals_to_universe: a renamed trial's rows stay in the universe")
close_store(con); unlink(dbp)

# ===========================================================================
cat("run provenance: the universe is pinned once and read back
")
# Oracle: the run row round-trips, and the same settings resolve to the same run_id -- which is
# what lets a restart resume a run rather than start a second one beside it.
dbpR <- tempfile(fileext = ".sqlite"); conR <- open_store(dbpR)
sR   <- modifyList(optimizer_settings(), list(simulate = TRUE, optimize_scheme = "CV00"))
ridR <- run_id_for(sR, "0.0.0")
check(identical(ridR, run_id_for(sR, "0.0.0")), "run_id_for is stable for the same settings")
check(!identical(ridR, run_id_for(modifyList(sR, list(trial_replication = 9)), "0.0.0")),
      "run_id_for changes when a search-defining setting changes")
check(identical(ridR, run_id_for(modifyList(sR, list(worker_id = "7")), "0.0.0")),
      "run_id_for ignores per-worker settings")
check(is.null(read_run(conR, ridR)), "read_run: no row before the run is recorded")
record_run(conR, ridR, sR, "0.0.0", c("E1", "E2", "E3"), c("n1", "n2", "n3"))
check(setequal(run_universe(read_run(conR, ridR)), c("E1", "E2", "E3")),
      "record_run/run_universe round-trip the pinned trial ids")
record_run(conR, ridR, sR, "0.0.0", c("Z9"))       # a second writer must not overwrite
check(setequal(run_universe(read_run(conR, ridR)), c("E1", "E2", "E3")),
      "record_run: the first writer wins, so all workers read one universe")
check(is.null(run_universe(NULL)), "run_universe: no row means no constraint")
close_store(conR); unlink(dbpR)
# Oracle: a config seen only OUTSIDE this run's universe still counts as untried, so
# choose_config re-offers a seed rather than treating it as done.
dbp2 <- tempfile(fileext = ".sqlite"); con2 <- open_store(dbp2)
# Pin the scheme rather than inheriting it: optimizer_settings() layers in settings.local.R,
# so a machine optimizing CV00 would filter out this CV0 row and the seed would look untried
# for the wrong reason. The seed must come from the SAME scheme for the hashes to correspond.
TEST_SCHEME <- "CV0"
seed1 <- seed_configs(TEST_SCHEME)[[1]]
# Stamp the current build: these oracles are about the UNIVERSE/SCHEME filters, and an
# unstamped row would additionally be retired by filter_evals_to_build, making the seed look
# untried for the wrong reason.
store_eval(con2, seed1, "900", TEST_SCHEME, 0.3, 40L, "ok",
           study_name = "elsewhere", build = OPTIMIZER_BUILD)
st <- modifyList(optimizer_settings(),
                 list(simulate = FALSE, optimize_scheme = TEST_SCHEME))
pick <- choose_config(con2, st, universe = c("901", "902"))
check(identical(config_hash(pick$cfg), config_hash(seed1)) && grepl("^seed:", pick$source),
      "choose_config: a seed eval outside the universe does not count as done in this run")
check(!identical(config_hash(choose_config(con2, st, universe = c("900", "901"))$cfg),
                config_hash(seed1)),
      "choose_config: the same row inside the universe does count as done")
# Control: with NO domain restriction, that same seed eval DOES count as done, so
# choose_config moves past seed1 to a later unevaluated seed.
st0  <- modifyList(optimizer_settings(),
                   list(simulate = FALSE, optimize_scheme = TEST_SCHEME, target_domain = NULL))
pick0 <- choose_config(con2, st0)
check(!identical(config_hash(pick0$cfg), config_hash(seed1)),
      "choose_config: with no domain, the recorded seed is not re-offered")
close_store(con2); unlink(dbp2)

# ===========================================================================
cat("aggregate_scores: trial-adjusted config estimates (incumbent + elites)\n")
# Oracle (the point of the change): a config that drew EASY trials must not outrank a better
# config that drew HARD ones. Trials differ ~2x more than configs on the real store, so the
# raw mean is confounded with which trials a config happened to get.
#   cA is genuinely better than cB, but cB was run on the two easy trials.
blup_ev <- tibble::tibble(
  config_hash = c("cA","cA","cB","cB","cC","cC","cD","cD"),
  trial_id    = c("hard1","hard2","easy1","easy2","hard1","easy1","hard2","easy2"),
  score       = c(0.30,0.32, 0.38,0.40, 0.22,0.44, 0.20,0.42),
  n_test      = 200L, config_json = "{}")
agg_raw  <- aggregate_scores(blup_ev, adjust_trial = FALSE)
agg_adj  <- aggregate_scores(blup_ev, adjust_trial = TRUE)
check(identical(attr(agg_raw, "estimator"), "pooled"), "adjust_trial = FALSE keeps the pooled estimator")
check(identical(attr(agg_adj, "estimator"), "blup"),   "adjust_trial = TRUE uses the BLUP estimator")
rawA <- agg_raw$mean_score[agg_raw$config_hash == "cA"]
rawB <- agg_raw$mean_score[agg_raw$config_hash == "cB"]
adjA <- agg_adj$mean_score[agg_adj$config_hash == "cA"]
adjB <- agg_adj$mean_score[agg_adj$config_hash == "cB"]
check(rawB > rawA, "raw mean ranks the easy-trial config ABOVE the better one (the bug)")
check((adjB - adjA) < (rawB - rawA),
      "the trial adjustment SHRINKS that spurious advantage")
vc <- attr(agg_adj, "var_comps")
check(!is.null(vc) && all(is.finite(vc)) && all(c("sd_trial","sd_config","sd_resid") %in% names(vc)),
      "variance components are returned for the report")
check(vc[["sd_trial"]] > 0, "a real trial effect is detected in data that has one")
check(is.null(attr(agg_adj, "estimator_note")), "a clean fit records no note")
check(is.finite(attr(agg_adj, "median_weight")),
      "the median fit weight is published for the report's sd_resid rescaling")

# Oracle: a fit that WARNS is kept, not discarded. tryCatch(warning = ...) threw away good
# estimates -- and with them every `se`, which .contenders() needs, so contender replication
# silently stopped. The note records the warning; the components still arrive.
local({
  real <- .blup_scores
  # Shadow it in the global env, where aggregate_scores resolves it, so the fit does exactly
  # what it does now AND warns on the way out.
  .blup_scores <<- function(evals) { r <- real(evals); warning("simulated convergence grumble"); r }
  on.exit(.blup_scores <<- real, add = TRUE)
  a <- aggregate_scores(blup_ev)
  check(identical(attr(a, "estimator"), "blup"), "a warned fit is still used")
  check(all(is.finite(attr(a, "var_comps"))), "and its variance components survive")
  check(all(is.finite(a$se)), "and `se` arrives, without which .contenders() is empty")
  check(grepl("grumble", attr(a, "estimator_note") %||% ""),
        "and the warning is recorded for the report")
})

# Oracle: when the design cannot support the fit, the REASON is recorded rather than nothing.
local({
  one_each <- tibble::tibble(               # one config per trial -> trial effect unidentifiable
    config_hash = c("a", "b"), trial_id = c("t1", "t2"), score = c(0.3, 0.4),
    n_test = c(50L, 50L), config_json = c("{}", "{}"))
  a <- aggregate_scores(one_each)
  check(identical(attr(a, "estimator"), "pooled"), "it falls back to pooled")
  check(grepl("not separable", attr(a, "estimator_note") %||% ""),
        "and names the reason -- the report prints this instead of a bare 'pooled'")
})

# Oracle: BLUPs shrink by REPLICATION -- same raw mean, fewer reps -> pulled further to the mean.
shrink_ev <- tibble::tibble(
  config_hash = c(rep("many", 8), "few", rep("filler", 8)),
  trial_id    = c(paste0("t", 1:8), "t1", paste0("t", 1:8)),
  score       = c(rep(0.45, 8), 0.45, rep(0.05, 8)),
  n_test      = 200L, config_json = "{}")
sa <- aggregate_scores(shrink_ev, adjust_trial = TRUE)
if (identical(attr(sa, "estimator"), "blup")) {
  gm <- mean(sa$mean_score, na.rm = TRUE)
  d_many <- abs(sa$mean_score[sa$config_hash == "many"] - gm)
  d_few  <- abs(sa$mean_score[sa$config_hash == "few"]  - gm)
  check(d_few < d_many,
        "same raw score, 1 rep vs 8: the 1-rep config is shrunk further toward the mean")
} else message("  (BLUP unavailable -- skipping the replication-shrinkage oracle)")

# Oracle: every fallback path returns the pooled estimate rather than erroring.
one_cfg_per_trial <- tibble::tibble(config_hash = c("a","b","c"), trial_id = c("t1","t2","t3"),
                                    score = c(0.2,0.3,0.4), n_test = 100L, config_json = "{}")
fb <- aggregate_scores(one_cfg_per_trial, adjust_trial = TRUE)
check(identical(attr(fb, "estimator"), "pooled"),
      "one config per trial is unidentifiable -> falls back to pooled")
check(all(is.finite(fb$mean_score)), "the fallback still produces finite scores")
single <- tibble::tibble(config_hash = "a", trial_id = "t1", score = 0.3,
                         n_test = 100L, config_json = "{}")
check(identical(attr(aggregate_scores(single), "estimator"), "pooled"),
      "a single row falls back to pooled without error")
check(nrow(aggregate_scores(single[0, ])) == 0, "an empty slice returns an empty tibble")

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
cat("trial discovery: coverage gating and the attempted set\n")
# The wizard call this replaces was 440 s of a 795 s evaluation -- 99.9% of .find_related.
# Discovery may only go local when the index covers every trial the CATALOGUE can offer as a
# candidate, because .germ_overlap filters candidates to the catalogue.
local({
  cdir <- tempfile("cov_"); dir.create(file.path(cdir, "acc"), recursive = TRUE)
  st <- list(cache_dir = cdir)
  put <- function(id, acc) saveRDS(acc, file.path(cdir, "acc", paste0("acc_", id, ".rds")))
  put("c1", c("a", "b")); put("c2", c("b", "c"))
  idx <- .trial_index(st)
  # A stub connection is enough: .index_covers_catalog only needs trial_catalog(), which the
  # cache serves, so no BrAPI call happens.
  saveRDS(tibble::tibble(study_db_id = c("c1", "c2", "c3"), location_name = "L",
                         latitude = 1, longitude = 1, year = 2024L),
          {dir.create(file.path(cdir, "trial_catalog")); file.path(cdir, "trial_catalog", "trial_catalog.rds")})
  check(!.index_covers(idx, c("c1","c2","c3"), st, "acc"),
        "coverage is FALSE while a catalogue trial is neither indexed nor attempted")
  # c3 has no accessions at all -- get_trial_accessions would not cache it, so only the
  # ATTEMPTED record can lift coverage. Without it, coverage could never reach 100%.
  .note_attempted(st, "acc", "c3")
  check(setequal(.attempted(st, "acc"), "c3"), "attempted trials are recorded")
  check(.index_covers(idx, c("c1","c2","c3"), st, "acc"),
        "an ATTEMPTED but empty trial counts as covered (else coverage never completes)")
  .note_attempted(st, "acc", "c9")
  check(setequal(.attempted(st, "acc"), c("c3", "c9")), "attempted records accumulate, not overwrite")
  unlink(cdir, recursive = TRUE)
})

# Oracle: local discovery yields the SAME candidates and counts as counting from an explicit
# candidate list. This is what licenses replacing the wizard.
local({
  cdir <- tempfile("disc_"); dir.create(file.path(cdir, "acc"), recursive = TRUE)
  st <- list(cache_dir = cdir)
  set.seed(303)
  pool <- paste0("g", 1:120); ids <- paste0("d", 1:40)
  for (i in ids) saveRDS(sample(pool, sample(5:40, 1)),
                         file.path(cdir, "acc", paste0("acc_", i, ".rds")))
  idx <- .trial_index(st)
  bad <- 0L
  for (r in 1:10) {
    acc <- sample(pool, 30)
    hits <- unlist(idx[intersect(acc, names(idx))], use.names = FALSE)
    tab  <- table(hits)
    disc <- names(tab)                                   # local discovery
    # ground truth: every trial sharing >= 1 accession, by brute force
    truth <- Filter(function(sid) length(intersect(
      as.character(readRDS(file.path(cdir, "acc", paste0("acc_", sid, ".rds")))), acc)) > 0, ids)
    if (!setequal(disc, truth)) { bad <- bad + 1L; next }
    cnt <- as.integer(tab[truth])
    ref <- vapply(truth, function(sid) length(intersect(
      as.character(readRDS(file.path(cdir, "acc", paste0("acc_", sid, ".rds")))), acc)), integer(1))
    if (!identical(cnt, as.integer(ref))) bad <- bad + 1L
  }
  check(bad == 0L,
        "local discovery finds EXACTLY the trials sharing an accession, with exact counts")
  unlink(cdir, recursive = TRUE)
})

# ===========================================================================
cat("local wizards: the projects index and memo isolation\n")
# The wizard answers `accessions` BY genotyping project, and there are only ~110 projects in
# the crop -- so accession -> projects is an inversion of a cheap primary map, exactly like
# accession -> trials. Verified against the live wizard on real accession sets: 0 mismatches.
local({
  cdir <- tempfile("pi_"); dir.create(file.path(cdir, "proj_acc"), recursive = TRUE)
  st <- list(cache_dir = cdir)
  putp <- function(pid, acc) saveRDS(acc, file.path(cdir, "proj_acc", paste0("proj_acc_", pid, ".rds")))
  putp("p1", c("a", "b")); putp("p2", c("b", "c"))
  pidx <- .project_index(st)
  check(setequal(pidx[["b"]], c("p1", "p2")), ".project_index inverts project -> accessions")
  check(setequal(attr(pidx, "keys"), c("p1", "p2")), ".project_index records the projects it saw")
  check(!.index_covers(pidx, c("p1","p2","p3"), st, "proj_acc"),
        "coverage is FALSE while a project in the universe is unseen")
  .note_attempted(st, "proj_acc", "p3")
  check(.index_covers(pidx, c("p1","p2","p3"), st, "proj_acc"),
        "an attempted-but-empty project counts as covered")
  # The trials and projects indices must not collide: same helper, different categories.
  dir.create(file.path(cdir, "acc"))
  saveRDS(c("a", "zzz"), file.path(cdir, "acc", "acc_t1.rds"))
  check(setequal(.trial_index(st)[["a"]], "t1") && setequal(.project_index(st)[["a"]], c("p1")),
        "the trial and project indices are independent")
  unlink(cdir, recursive = TRUE)
})

# Oracle: the RAM memo must key on the CACHE DIRECTORY too. Two caches with the same file
# count -- and, in a fast test, the same newest mtime -- would otherwise share an entry and
# one would be served the other's index.
local({
  d1 <- tempfile("m1_"); d2 <- tempfile("m2_")
  for (d in c(d1, d2)) dir.create(file.path(d, "acc"), recursive = TRUE)
  saveRDS(c("x"), file.path(d1, "acc", "acc_A.rds"))
  saveRDS(c("y"), file.path(d2, "acc", "acc_B.rds"))
  # FORCE the collision: identical file count AND identical mtime, so the signature alone
  # cannot tell the two caches apart. Without this the test passes either way, because temp
  # files created milliseconds apart already differ.
  tt <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  Sys.setFileTime(file.path(d1, "acc", "acc_A.rds"), tt)
  Sys.setFileTime(file.path(d2, "acc", "acc_B.rds"), tt)
  i1 <- .trial_index(list(cache_dir = d1))
  i2 <- .trial_index(list(cache_dir = d2))
  check(identical(names(i1), "x") && identical(names(i2), "y"),
        "two cache dirs with identical file counts get their OWN index, not a shared memo")
  unlink(c(d1, d2), recursive = TRUE)
})

# ===========================================================================
cat(".trial_index (inverted accession -> trials lookup)\n")
# .germ_overlap decides which trials become the TRAINING SET, so an index that is fast but
# subtly different silently changes every score and would read as a modelling result rather
# than a bug. These oracles demand exact equivalence, not agreement.
local({
  cdir <- tempfile("cache_"); dir.create(file.path(cdir, "acc"), recursive = TRUE)
  st <- list(cache_dir = cdir)
  put <- function(id, acc) saveRDS(acc, file.path(cdir, "acc", paste0("acc_", id, ".rds")))
  put("t1", c("a", "b", "c"))
  put("t2", c("b", "c", "d"))
  put("t3", c("z"))
  idx <- .trial_index(st)
  check(setequal(names(idx), c("a","b","c","d","z")), ".trial_index keys on every accession")
  check(setequal(idx[["b"]], c("t1","t2")), "an accession maps to every trial containing it")
  check(identical(idx[["z"]], "t3"), "a single-trial accession maps to just that trial")
  check(is.null(.trial_index(list(cache_dir = tempfile()))),
        ".trial_index returns NULL when there are no accession caches")

  # LEGACY FLAT LAYOUT. .cache_existing() reads cache/acc_<id>.rds as well as the nested
  # cache/acc/acc_<id>.rds, so the index must too. Globbing only the nested path left a
  # server with an older cache reporting "covered 0 of 432 candidates" -- the index was not
  # empty, it was never built.
  fdir <- tempfile("flat_"); dir.create(fdir)
  saveRDS(c("p", "q"), file.path(fdir, "acc_f1.rds"))       # flat, no acc/ subdirectory
  saveRDS(c("q", "r"), file.path(fdir, "acc_f2.rds"))
  fidx <- .trial_index(list(cache_dir = fdir))
  check(!is.null(fidx) && setequal(fidx[["q"]], c("f1", "f2")),
        ".trial_index reads the LEGACY FLAT cache layout, not just the nested one")
  check(setequal(attr(fidx, "keys"), c("f1", "f2")), "flat-layout trials are recorded as indexed")
  # Both layouts at once: the nested copy wins, as in .cache_existing.
  dir.create(file.path(fdir, "acc"))
  saveRDS(c("p", "NESTED"), file.path(fdir, "acc", "acc_f1.rds"))
  # Force every mtime equal before re-reading. Dedup by basename already leaves the file
  # COUNT unchanged, so this makes the memo signature's other component useless too -- which
  # is what a filesystem with coarse timestamp granularity does for free. This test passed on
  # APFS by luck and failed on SciNet; pinning the mtimes makes it deterministic everywhere.
  stamp <- as.POSIXct("2026-08-06 12:00:00", tz = "UTC")
  for (f in list.files(fdir, recursive = TRUE, full.names = TRUE)) Sys.setFileTime(f, stamp)
  bidx <- .trial_index(list(cache_dir = fdir))
  check("NESTED" %in% names(bidx) && setequal(bidx[["p"]], c("f1")),
        "with both layouts present the nested file takes precedence")
  unlink(fdir, recursive = TRUE)

  # An accession repeated within one trial must not double-count that trial's overlap.
  put("t4", c("q", "q", "r"))
  idx <- .trial_index(st)
  check(identical(idx[["q"]], "t4"), "a duplicated accession yields the trial once, not twice")

  # Adding a trial must invalidate the cached index rather than serve a stale one.
  put("t5", c("a", "new"))
  idx2 <- .trial_index(st)
  check("new" %in% names(idx2) && setequal(idx2[["a"]], c("t1","t5")),
        "a new accession cache file forces a rebuild")

  # THE ORACLE THAT LICENSES THE CHANGE: tabulation == the per-candidate intersect loop.
  set.seed(77)
  ids <- paste0("s", 1:60)
  pool <- paste0("g", 1:200)
  for (i in ids) put(i, sample(pool, sample(10:60, 1)))
  idx3 <- .trial_index(st)
  bad <- 0L
  for (rep in 1:10) {
    acc  <- sample(pool, 40)
    cand <- sample(ids, 25)
    loop <- vapply(cand, function(sid)
      length(intersect(as.character(readRDS(file.path(cdir,"acc",paste0("acc_",sid,".rds")))), acc)),
      integer(1))
    hits <- unlist(idx3[intersect(acc, names(idx3))], use.names = FALSE)
    tb <- table(hits); v <- as.integer(tb); names(v) <- names(tb)
    fromidx <- ifelse(cand %in% names(v), v[cand], 0L); fromidx[is.na(fromidx)] <- 0L
    if (!identical(as.integer(loop), as.integer(fromidx))) bad <- bad + 1L
  }
  check(bad == 0L, "index tabulation EQUALS the per-candidate intersect loop, element for element")

  # The index must record WHICH trials it saw, so a candidate absent from the tabulation can
  # be read as "zero overlap" rather than "unknown". Without it the index is all-or-nothing:
  # one unindexed candidate sent every candidate back to the loop, which is exactly why the
  # first cut of this showed no improvement on the server (803 s -> 736 s).
  seen <- attr(idx3, "keys")
  check(!is.null(seen) && all(ids %in% seen),
        ".trial_index records the trial ids it indexed")
  check(!("never_indexed" %in% seen), "a trial with no accession cache is not claimed as indexed")

  # A trial that IS indexed but shares nothing with the query must score 0, not fall back.
  put("iso", c("zzz_unique_1", "zzz_unique_2"))
  idx4 <- .trial_index(st)
  hits <- unlist(idx4[intersect(c("g1","g2"), names(idx4))], use.names = FALSE)
  tb <- table(hits); v <- as.integer(tb); names(v) <- names(tb)
  check("iso" %in% attr(idx4, "keys") && !("iso" %in% names(v)),
        "an indexed trial with zero overlap is absent from the tabulation (reads as 0)")
  unlink(cdir, recursive = TRUE)
})

# ===========================================================================
cat("choose_trial (trial_replication) + trial_id as a surrogate feature\n")
# Oracle: .sim_trial must be a FUNCTION of the id. Revisiting a simulated trial has to return
# the SAME trial, or trial_replication can never be satisfied offline and trial_id is
# degenerate in every simulate-mode test.
check(identical(.sim_trial("simtrial_42"), .sim_trial("simtrial_42")),
      ".sim_trial is deterministic in the id")
check(!identical(.sim_trial("simtrial_42")$heritability, .sim_trial("simtrial_43")$heritability),
      ".sim_trial gives different trials different attributes")
local({                                    # and it must not disturb the caller's RNG stream
  set.seed(99); a <- runif(3)
  set.seed(99); invisible(.sim_trial("simtrial_7")); b <- runif(3)
  check(isTRUE(all.equal(a, b)), ".sim_trial restores the global RNG state")
})

dbp6 <- tempfile(fileext = ".sqlite"); con6 <- open_store(dbp6)
tset <- function(w = 1, r = 2) modifyList(optimizer_settings(),
  list(simulate = TRUE, optimize_scheme = "CV00", trial_replication = r,
       worker_id = as.character(w)))
# Empty store -> nothing to revisit, fall through to a fresh sample.
check(grepl("^simtrial_", choose_trial(con6, tset())$id),
      "empty store: choose_trial falls through to sample_trial")
c1 <- sample_config(); c2 <- sample_config()
store_eval(con6, c1, "simtrial_111", "CV00", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)
check(identical(choose_trial(con6, tset())$id, "simtrial_111"),
      "a trial with 1 config and trial_replication=2 is REVISITED")
store_eval(con6, c2, "simtrial_111", "CV00", 0.4, 40L, "ok", build = OPTIMIZER_BUILD)
check(!identical(choose_trial(con6, tset())$id, "simtrial_111"),
      "once it has 2 distinct configs it leaves the backlog")
# trial_replication = 1 disables revisiting outright.
store_eval(con6, c1, "simtrial_222", "CV00", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)
check(!identical(choose_trial(con6, tset(r = 1))$id, "simtrial_222"),
      "trial_replication = 1 disables revisiting")
# Workers must not pile onto the same backlog trial (the 0.7.4 lockstep bug, again).
store_eval(con6, c1, "simtrial_333", "CV00", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)
store_eval(con6, c1, "simtrial_444", "CV00", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)
wpicks <- vapply(1:3, function(w) choose_trial(con6, tset(w = w))$id, character(1))
check(dplyr::n_distinct(wpicks) == 3,
      "three workers revisit three DIFFERENT backlog trials")
# A row under the other scheme must not count toward a trial's config tally.
store_eval(con6, c2, "simtrial_222", "CV0", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)
check(identical(choose_trial(con6, tset())$id, "simtrial_222"),
      "an eval under the OTHER scheme does not satisfy the replication constraint")
close_store(con6); unlink(dbp6)

# ===========================================================================
cat(".trial_backlog: what may hold the gate open\n")
# Oracle: a non-empty backlog stops fresh trials being drawn, so anything that can sit in it
# unsatisfiably stalls the whole search. The two views the caller passes are what keep that
# from happening -- evidence decides who is owed replication, `seen` (built from every computed
# row, exactly as claim_eval keys them) decides what has already been done.
bl_ev <- function(...) tibble::tibble(trial_id = c(...),
                                      config_hash = paste0("c", seq_along(c(...))),
                                      score = 0.3)
check(setequal(.trial_backlog(bl_ev("A", "B"), r = 2L), c("A", "B")),
      ".trial_backlog: under-replicated trials are owed a revisit")
check(!length(.trial_backlog(bl_ev("A", "A"), r = 2L)),
      ".trial_backlog: a trial at the replication target leaves the backlog")
check(!length(.trial_backlog(bl_ev("A"), seen = "A", r = 2L)),
      ".trial_backlog: a trial this config already ran is not re-proposed")
check(identical(.trial_backlog(bl_ev("A", "B"), universe = "B", r = 2L), "B"),
      ".trial_backlog: a trial outside the universe cannot hold the gate open")
check(setequal(.trial_backlog(bl_ev("A", "B"), universe = NULL, r = 2L), c("A", "B")),
      ".trial_backlog: a NULL universe is no constraint")
check(!length(.trial_backlog(bl_ev()[0, ], r = 2L)),
      ".trial_backlog: an empty store has no backlog")
# The regression this all exists for: rows that were COMPUTED but are not evidence (their trial
# is no longer in the universe) must not be able to keep proposing that trial, or selection and
# claim_eval disagree forever -- proposed every iteration, claimable never.
local({
  computed <- bl_ev("A", "OUT")
  evidence <- filter_evals_to_universe(computed, "A")
  seen     <- unique(as.character(computed$trial_id))
  check(!length(.trial_backlog(evidence, seen, "A", r = 2L)),
        ".trial_backlog: a computed-but-not-evidence trial is never re-proposed")
})

# ===========================================================================
cat("choose_config (config_replication) and the trial it is paired with\n")
# Oracle: a configuration short of config_replication evaluations must come back, and must come
# back paired with a trial it has NOT already been run on -- the real pipeline is deterministic
# in (config, trial, scheme), so repeating a pair recomputes a known score and double-counts it.
dbp7 <- tempfile(fileext = ".sqlite"); con7 <- open_store(dbp7)
cset <- function(w = 1, cr = 2, tr = 1) modifyList(optimizer_settings(),
  list(simulate = TRUE, optimize_scheme = "CV00", config_replication = cr,
       trial_replication = tr, worker_id = as.character(w)))
put <- function(cfg, trial, scheme = "CV00", score = 0.3)
  store_eval(con7, cfg, trial, scheme, score, 40L, "ok", build = OPTIMIZER_BUILD)

# Retire the seed phase first: every seed already replicated, so it is out of both backlogs.
for (s in seed_configs("CV00")) for (t in c("simtrial_900", "simtrial_901")) put(s, t, score = 0.25)

cA <- sample_config(); hA <- config_hash(cA)
put(cA, "simtrial_910", score = 0.5)
bl <- function(st = cset()) {
  e <- filter_evals_to_scheme(read_evals(con7), st$optimize_scheme)
  unlist(.replication_backlog(e, aggregate_scores(e), st), use.names = FALSE)
}
check(hA %in% bl(), "a config with 1 of 2 evaluations is in the replication backlog")
check(identical(choose_config(con7, cset())$source, "replicate"),
      "and the worker whose slot is due takes from it, source = replicate")
# The pairing: simtrial_910 is the only trial in the TRIAL backlog, so without cfg_hash it is
# what choose_trial returns -- and with it, it must not be.
check(identical(choose_trial(con7, cset(tr = 2))$id, "simtrial_910"),
      "the trial backlog offers the under-replicated trial")
check(!identical(choose_trial(con7, cset(tr = 2), cfg_hash = hA)$id, "simtrial_910"),
      "cfg_hash excludes a trial that configuration has already been run on")

put(cA, "simtrial_911", score = 0.4)
check(!(hA %in% unlist(.replication_backlog(
          filter_evals_to_scheme(read_evals(con7), "CV00"),
          local({ e <- filter_evals_to_scheme(read_evals(con7), "CV00")
                  a <- aggregate_scores(e); a$se <- NA_real_; a }),   # no SEs -> no contenders
          cset()))),
      "at config_replication evaluations and not a contender, it leaves the backlog")

cB <- sample_config()
put(cB, "simtrial_912")
check(!(config_hash(cB) %in% unlist(.replication_backlog(
          filter_evals_to_scheme(read_evals(con7), "CV00"),
          local({ e <- filter_evals_to_scheme(read_evals(con7), "CV00")
                  a <- aggregate_scores(e); a$se <- NA_real_; a }),
          cset(cr = 1)))),
      "config_replication = 1 leaves a once-evaluated config out of the backlog")
# A row under the OTHER scheme must not count toward the tally.
put(cB, "simtrial_913", scheme = "CV0")
check(config_hash(cB) %in% bl(),
      "an eval under the OTHER scheme does not count toward config_replication")

# Workers must spread over the backlog, and those past its end must explore rather than pile on.
cC <- sample_config()
put(cC, "simtrial_914")
check(length(bl()) >= 2, "the backlog holds more than one config")

# Oracle: the backlog counts EVALUATIONS, so it drains even where only one trial exists (the
# sim_fixed_trial case). Counting distinct trials would leave cE stuck below the target forever.
cE <- sample_config()
put(cE, "simtrial_fixed"); put(cE, "simtrial_fixed")
# --- trial exhaustion -------------------------------------------------------
# Oracle: a config already run on every eligible trial has nothing left to learn from. The
# universe comes from the catalogue restricted to target_domain, so this runs offline against a
# pre-written trial_catalog cache -- no network, and no descriptor is ever built.
local({
  ctmp <- tempfile("exh_"); dir.create(ctmp)
  catl <- tibble::tibble(study_db_id = c("E1", "E2", "E3"),
                         study_name = c("e1", "e2", "e3"),
                         program_name = "P", location_name = "L", year = 2025L)
  es <- modifyList(optimizer_settings(),
    list(simulate = FALSE, cache_dir = ctmp, optimize_scheme = "CV00",
         target_domain = list(programs = NULL, years = NULL, locations = NULL,
                              trials = c("e1", "e2", "e3"))))
  .cache_save(es, "trial_catalog", NULL, catl)
  check(setequal(eligible_trial_ids(NULL, es), c("E1", "E2", "E3")),
        "eligible_trial_ids returns the domain's trials, offline from the cached catalogue")

  dbx <- tempfile(fileext = ".sqlite"); cx <- open_store(dbx)
  cX <- sample_config(); hX <- config_hash(cX)
  for (t in c("E1", "E2")) store_eval(cx, cX, t, "CV00", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)
  ex <- read_evals(cx)
  ax <- local({ a <- aggregate_scores(ex); a$se <- 0.5; a })      # force it to be a contender
  check(hX %in% unlist(.replication_backlog(ex, ax, es, universe = c("E1","E2","E3"))),
        "a contender with 2 of 3 eligible trials is still in the backlog")
  store_eval(cx, cX, "E3", "CV00", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)
  ex <- read_evals(cx); ax <- local({ a <- aggregate_scores(ex); a$se <- 0.5; a })
  check(!(hX %in% unlist(.replication_backlog(ex, ax, es, universe = c("E1","E2","E3")))),
        "once it covers the whole domain it is dropped, contender or not")
  check(hX %in% unlist(.replication_backlog(ex, ax, es, universe = NULL)),
        "and with no universe (simulate mode) the domain test does not apply")
  close_store(cx); unlink(dbx)

  # .sample_unseen_trial signals rather than handing back a duplicate. The stub sampler always
  # yields the one id the config has already seen -- the case the setdiff cannot catch, where
  # unseen trials exist but none is usable.
  sim1 <- modifyList(optimizer_settings(), list(simulate = TRUE, sim_fixed_trial = TRUE))
  hit <- tryCatch(.sample_unseen_trial(sim1, NULL, seen = "simtrial_fixed", tries = 3L),
                  optimizer_trials_exhausted = function(e) e)
  check(inherits(hit, "optimizer_trials_exhausted"),
        ".sample_unseen_trial raises trials_exhausted instead of repeating a pair")
  check(!inherits(tryCatch(.sample_unseen_trial(sim1, NULL, seen = character()),
                           optimizer_trials_exhausted = function(e) e),
                  "optimizer_trials_exhausted"),
        "... and returns normally when the trial has not been seen")
  unlink(ctmp, recursive = TRUE)
})

# Oracle: a SATURATED optimizer must not look like broken trial sampling. sim_fixed_trial gives
# exactly one trial, so every replication pick is instantly exhausted -- the end state a bounded
# target domain reaches. run_optimizer.R does not auto-run when sourced (it checks sys.nframe).
local({
  source(here::here("run_optimizer.R"))
  st <- modifyList(optimizer_settings(), list(
    simulate = TRUE, sim_fixed_trial = TRUE, sim_noise_sd = 0, optimize_scheme = "CV0",
    db_path = tempfile(fileext = ".sqlite"), stop_file = tempfile(),
    report_path = tempfile(fileext = ".md"), log_dir = tempdir(), cache_dir = tempdir(),
    db_backup_path = NULL, cache_backup_dir = NULL, cache_ready_file = NULL,
    max_iters = 20, checkpoint_every = 1000))
  set.seed(11)
  out <- capture.output(run_optimizer(st), type = "message")
  cn <- open_store(st$db_path); ev <- read_evals(cn); close_store(cn)
  check(!any(grepl("too many consecutive trial-sampling failures", out)),
        "exhaustion never trips the max_sample_fail halt")
  check(nrow(ev) == 20L,
        "no iteration is idled: an exhausted pick falls through to exploration in the same step")
  check(!nrow(dplyr::filter(dplyr::count(ev, config_hash, trial_id, scheme), n > 1)),
        "and no (config, trial, scheme) pair is ever evaluated twice")
  unlink(st$db_path)
})

# --- the replication share --------------------------------------------------
# Oracle: one worker in `replicate_every` replicates at any instant, and a single worker
# alternates over time -- keyed on the row count so N workers never fall into lockstep.
local({
  due <- function(n, w, every = 3L) (n + w) %% every == 0L
  check(sum(vapply(1:9, function(w) due(100L, w), logical(1))) == 3L,
        "share: 3 of 9 workers are due at any one row count")
  check(sum(vapply(1:30, function(n) due(n, 1L), logical(1))) == 10L,
        "share: a single worker is due on 1 iteration in 3 as the store grows")
  check(!identical(which(vapply(1:9, function(w) due(100L, w), logical(1))),
                   which(vapply(1:9, function(w) due(101L, w), logical(1)))),
        "share: which workers are due shifts as rows land, so it is not always the same ones")
})

# Oracle: a stop-file left over from a PREVIOUS job must not kill this one. It lives on durable
# storage, so it outlives the job that consumed it, and the loop tests it before iteration 1.
# The complement matters just as much: a STOP that appears once the loop is running still stops
# it, which is what makes clearing the stale one safe.
local({
  source(here::here("run_optimizer.R"))
  base <- function(stop_file) modifyList(optimizer_settings(), list(
    simulate = TRUE, optimize_scheme = "CV0", max_iters = 3, n_random_init = 5, ntree = 20,
    db_path = tempfile(fileext = ".sqlite"), stop_file = stop_file,
    report_path = tempfile(fileext = ".md"), log_dir = tempdir(), cache_dir = tempdir(),
    db_backup_path = NULL, cache_backup_dir = NULL, cache_ready_file = NULL,
    checkpoint_every = 1000))

  sf <- tempfile(); file.create(sf)                    # a STOP left behind by the last job
  st <- base(sf)
  invisible(capture.output(run_optimizer(st), type = "message"))
  cn <- open_store(st$db_path); n <- n_evals(cn); close_store(cn); unlink(st$db_path)
  check(n == 3L, "a stale stop-file does not halt the next run")
  check(!file.exists(sf), "and the leader clears it")

  # Still stops when the file appears mid-run: written from the report hook, which the loop
  # calls after every iteration.
  sf2 <- tempfile(); st2 <- base(sf2)
  st2$checkpoint_every <- 1
  st2$report_path <- tempfile(fileext = ".md")
  st2$max_iters <- 20
  local({
    orig <- write_report
    assign("write_report", function(con, settings) { file.create(sf2); orig(con, settings) },
           envir = globalenv())
    on.exit(assign("write_report", orig, envir = globalenv()), add = TRUE)
    invisible(capture.output(run_optimizer(st2), type = "message"))
  })
  cn2 <- open_store(st2$db_path); n2 <- n_evals(cn2); close_store(cn2); unlink(st2$db_path)
  check(n2 < 20L, "a stop-file created mid-run still halts the loop")
  unlink(c(sf, sf2))
})

# Oracle: the BASE floor is not rationed. Only the contender tier is, so config_replication
# stays a guarantee rather than a suggestion when replicate_every > 1.
local({
  source(here::here("run_optimizer.R"))
  st <- modifyList(optimizer_settings(), list(
    simulate = TRUE, optimize_scheme = "CV0", max_iters = 40, n_random_init = 10, ntree = 40,
    db_path = tempfile(fileext = ".sqlite"), stop_file = tempfile(),
    report_path = tempfile(fileext = ".md"), log_dir = tempdir(), cache_dir = tempdir(),
    db_backup_path = NULL, cache_backup_dir = NULL, cache_ready_file = NULL,
    checkpoint_every = 1000, replicate_every = 3L))
  set.seed(9)
  invisible(capture.output(run_optimizer(st), type = "message"))
  cn <- open_store(st$db_path); ev <- read_evals(cn); close_store(cn)
  reps <- table(ev$config_hash)
  # The config the run was on when max_iters cut it off gets a pass: one first evaluated on the
  # final iteration had no remaining step in which to reach its floor. Everything else must.
  last <- ev$config_hash[which.max(ev$id)]
  reps <- reps[!(names(reps) == last & as.integer(reps) < st$config_replication)]
  check(min(as.integer(reps)) >= st$config_replication,
        "every config reaches config_replication even though only 1 worker in 3 replicates")
  check(max(as.integer(reps)) > st$config_replication,
        "and a contender goes past it -- the ramp is live")
  unlink(st$db_path)
})

# --- contenders and the trial schedule --------------------------------------
local({
  mk <- function(score, se) tibble::tibble(config_hash = paste0("h", seq_along(score)),
                                           mean_score = score, se = se)
  a <- mk(c(0.50, 0.48, 0.20), c(0.02, 0.02, 0.02))
  cn <- .contenders(a, z = 1)
  check("h1" %in% cn, "the leader is always a contender")
  check("h2" %in% cn, "one whose bound exactly reaches the leader is a contender (>=, not >)")
  check(!("h3" %in% cn), "one whose optimistic bound falls short is not")
  check(length(.contenders(mk(seq(0.5, 0.4, length.out = 20), rep(0.2, 20)), z = 1)) == 8L,
        "the contender set is capped at 8")
  check(length(.contenders(mk(c(0.5, 0.2), c(NA_real_, NA_real_)), z = 1)) == 0L,
        "no standard errors (pooled fallback) means no contenders")
  ts <- function(n) .trial_target(modifyList(optimizer_settings(), list(trial_replication = 2)), n)
  check(identical(c(ts(0), ts(25), ts(100), ts(400)), c(2L, 3L, 4L, 6L)),
        "trial_replication schedule: 2 / 3 / 4 / 6 at 0 / 25 / 100 / 400 scored configs")
  check(identical(.trial_target(modifyList(optimizer_settings(),
                                           list(trial_replication = 1)), 400), 1L),
        "trial_replication = 1 disables the schedule outright")
})

close_store(con7); unlink(dbp7)

# Oracle: trial_id is only safe as a feature once every trial has >= 2 configs.
mk_ev <- function(n_cfg) tibble::tibble(
  trial_id = rep(c("A","B"), each = n_cfg),
  config_hash = paste0("c", seq_len(2 * n_cfg)), score = 0.3)
check(!trial_feature_usable(mk_ev(1)), "trial_id NOT usable at 1 config per trial")
check(trial_feature_usable(mk_ev(2)),  "trial_id usable at 2 configs per trial")
check(!trial_feature_usable(tibble::tibble(trial_id = character(), config_hash = character())),
      "trial_feature_usable is FALSE on an empty slice")

# Oracle: the marginal prediction is the mean over trials of the per-trial predictions, and
# scoring never needs a level the forest has not seen.
local({
  set.seed(5); n <- 60
  tr <- factor(sample(c("t1","t2","t3"), n, TRUE))
  X  <- data.frame(a = runif(n), b = runif(n), trial_id = tr)
  y  <- X$a * 2 + as.numeric(tr) * 0.5 + rnorm(n, 0, 0.05)
  m  <- fit_surrogate(X, y, ntree = 40, min_obs = 10)
  newX <- data.frame(a = c(0.2, 0.8), b = c(0.5, 0.5))
  pm <- predict_surrogate_marginal(m, newX, levels(tr))
  # Recompute the marginal by hand: predict on each trial, average across trials.
  byhand <- vapply(levels(tr), function(lv) {
    g <- newX; g$trial_id <- factor(lv, levels = levels(tr))
    predict_surrogate(m, g)$mean }, numeric(2))
  check(isTRUE(all.equal(pm$mean, rowMeans(byhand), tolerance = 1e-8)),
        "marginal mean == mean over trials of the per-trial predictions")
  check(length(pm$sd) == 2 && all(is.finite(pm$sd)), "marginal sd is finite, one per candidate")
  check(pm$mean[2] > pm$mean[1], "the marginal still ranks on the config feature (a)")
})

# ===========================================================================
cat("fit_surrogate_merf (trial as a random effect outside the tree)\n")
# Oracle: the trial effect is fitted OUTSIDE the forest, so the returned model predicts the
# trial-marginal value with no trial_id feature and no marginalisation step -- and it still
# ranks candidates on the config feature with a large trial effect present.
local({
  set.seed(11); n <- 150
  tr <- factor(sample(c("t1", "t2", "t3", "t4"), n, TRUE))
  X  <- data.frame(a = runif(n), b = runif(n))
  y  <- X$a * 2 + as.numeric(tr) * 3 + rnorm(n, 0, 0.05)      # trial effect dwarfs the signal
  m  <- fit_surrogate_merf(X, y, tr, ntree = 40, min_obs = 10)
  if (is.null(m)) {
    check(!requireNamespace("lme4", quietly = TRUE), "merf returns NULL only when lme4 is absent")
  } else {
    check(inherits(m, "rf_surrogate"), "merf returns an rf_surrogate")
    check(!("trial_id" %in% m$feat_names), "merf keeps trial_id OUT of the features")
    pr <- predict_surrogate(m, data.frame(a = c(0.2, 0.8), b = c(0.5, 0.5)))
    check(all(is.finite(pr$mean)) && all(is.finite(pr$sd)), "merf predictions are finite")
    check(pr$mean[2] > pr$mean[1], "merf ranks on the config feature despite the trial effect")
    # The EI floor is 0.25 * sd(y_obs); the fit ran on y - b, so keeping the residualised
    # response would shrink the floor and explore less.
    check(isTRUE(all.equal(sort(m$y_obs), sort(y))), "merf stores the UNRESIDUALISED y_obs")
  }
  check(is.null(fit_surrogate_merf(X, y, factor(rep("only", n)), ntree = 40, min_obs = 10)),
        "merf returns NULL with a single trial, so the caller falls back")
  check(is.null(fit_surrogate_merf(X[1:5, ], y[1:5], tr[1:5], ntree = 40, min_obs = 10)),
        "merf returns NULL below min_obs")
})

# Oracle: surrogate_method selects the fit, and an un-fittable choice falls through rather than
# failing. One trial cannot support merf or blocked, so both must land on pooled.
local({
  sset <- function(m) modifyList(optimizer_settings(),
    list(simulate = TRUE, optimize_scheme = "CV00", surrogate_method = m,
         n_random_init = 8, trial_replication = 1, config_replication = 1))
  dbm <- tempfile(fileext = ".sqlite"); cm <- open_store(dbm)
  for (i in 1:14) store_eval(cm, sample_config(), "one_trial", "CV00", 0.2 + i / 100,
                             40L, "ok", build = OPTIMIZER_BUILD)
  for (m in c("merf", "blocked", "pooled")) {
    got <- choose_config(cm, sset(m), replicate = FALSE)
    check(identical(got$source, "acquisition") || grepl("^seed:|^random", got$source),
          sprintf("surrogate_method = %s falls back to pooled on a single-trial store", m))
  }
  close_store(cm); unlink(dbm)
})

# ===========================================================================
cat("restore_store_from_backup (merging the durable backup into the work store)\n")
# Oracle: the merge is ADDITIVE, like the cache rsync it mirrors. Rows already on the work disk
# survive; a (config, trial, scheme) cell already present is not duplicated.
local({
  mkstore <- function(rows) {
    p <- tempfile(fileext = ".sqlite"); con <- open_store(p)
    for (r in rows) store_eval(con, r$cfg, r$trial, r$scheme, r$score, 40L, "ok",
                               build = r$build %||% OPTIMIZER_BUILD)
    close_store(con); p
  }
  cA <- sample_config(); cB <- sample_config()
  bakp <- mkstore(list(list(cfg = cA, trial = "T1", scheme = "CV0",  score = 0.3),
                       list(cfg = cA, trial = "T2", scheme = "CV0",  score = 0.4),
                       list(cfg = cB, trial = "T1", scheme = "CV0",  score = 0.5),
                       list(cfg = cA, trial = "T1", scheme = "CV00", score = 0.6)))
  sset <- function(db) modifyList(optimizer_settings(),
                                  list(db_path = db, db_backup_path = bakp))

  # 1. empty work store -> everything arrives.
  w1 <- tempfile(fileext = ".sqlite")
  r1 <- restore_store_from_backup(sset(w1))
  check(identical(r1$inserted, 4L) && identical(r1$skipped, 0L),
        "empty work store: every backup row is inserted")
  c1 <- open_store(w1)
  check(nrow(read_evals(c1)) == 4L, "and they are readable afterwards")
  # A CV0 and a CV00 row for ONE (config, trial) are different tasks and must both survive.
  check(sum(read_evals(c1)$trial_id == "T1" &
            read_evals(c1)$config_hash == config_hash(cA)) == 2L,
        "scheme is part of the key: the CV0 and CV00 rows both survive")
  close_store(c1)

  # 2. idempotent -- the property that makes it safe for the leader to run every startup.
  r2 <- restore_store_from_backup(sset(w1))
  check(identical(r2$inserted, 0L) && identical(r2$skipped, 4L),
        "running it again inserts nothing")

  # 3. partly-populated work store: only the missing cells arrive, local rows are untouched.
  w3 <- mkstore(list(list(cfg = cA, trial = "T1", scheme = "CV0", score = 0.99),
                     list(cfg = cB, trial = "T9", scheme = "CV0", score = 0.11)))
  r3 <- restore_store_from_backup(sset(w3))
  check(identical(r3$inserted, 3L), "a partly-populated store takes only the cells it lacks")
  c3 <- open_store(w3); e3 <- read_evals(c3); close_store(c3)
  check(any(e3$trial_id == "T9"), "the work store's own row is not lost")
  check(sum(e3$trial_id == "T1" & e3$scheme == "CV0" &
            e3$config_hash == config_hash(cA)) == 1L, "and the collided cell is not duplicated")
  check(isTRUE(e3$score[e3$trial_id == "T1" & e3$scheme == "CV0" &
                        e3$config_hash == config_hash(cA)] == 0.99),
        "local wins a collision -- it is the row this build computed")

  # 4. duplicates ALREADY in the backup collapse to the newest (MAX(id)).
  cdup <- open_store(bakp)
  store_eval(cdup, cB, "T1", "CV0", 0.77, 40L, "ok", build = OPTIMIZER_BUILD)  # 2nd cB/T1/CV0
  close_store(cdup)
  w4 <- tempfile(fileext = ".sqlite")
  r4 <- restore_store_from_backup(sset(w4))
  check(identical(r4$inserted, 4L), "a duplicated cell in the backup is inserted once")
  c4 <- open_store(w4); e4 <- read_evals(c4); close_store(c4)
  check(isTRUE(e4$score[e4$config_hash == config_hash(cB) & e4$trial_id == "T1"] == 0.77),
        "and the surviving row is the most recent one")

  # 5. a backup from an older build lacks the migrated columns; those are the ones worth
  #    rescuing, so a fixed column list would fail exactly where it must not.
  oldp <- tempfile(fileext = ".sqlite")
  co <- DBI::dbConnect(RSQLite::SQLite(), oldp)
  DBI::dbExecute(co, "CREATE TABLE evals (id INTEGER PRIMARY KEY AUTOINCREMENT,
     config_hash TEXT, config_json TEXT, trial_id TEXT, scheme TEXT, score REAL,
     n_test INTEGER, status TEXT, reason TEXT, seconds REAL, ts TEXT)")
  DBI::dbExecute(co, "INSERT INTO evals (config_hash, config_json, trial_id, scheme, score,
     n_test, status, ts) VALUES (?,?,?,?,?,?,?,?)",
     params = list(config_hash(cA), config_to_json(cA), "T7", "CV0", 0.42, 40L, "ok", "2026-01-01"))
  DBI::dbDisconnect(co)
  w5 <- tempfile(fileext = ".sqlite")
  r5 <- restore_store_from_backup(modifyList(optimizer_settings(),
                                             list(db_path = w5, db_backup_path = oldp)))
  check(identical(r5$inserted, 1L), "a backup missing the migrated columns still merges")

  # 6. no backup configured, and a configured path that does not exist: both are no-ops.
  check(is.null(restore_store_from_backup(modifyList(optimizer_settings(),
          list(db_path = tempfile(), db_backup_path = NULL)))),
        "no backup configured is a no-op")
  check(is.null(restore_store_from_backup(modifyList(optimizer_settings(),
          list(db_path = tempfile(), db_backup_path = tempfile())))),
        "a missing backup file is a no-op, not an error")

  # store_summary reads any store file, and NULL rather than erroring on a missing one.
  ss <- store_summary(bakp)
  check(!is.null(ss) && ss$rows == 5L && ss$n_config == 2L,
        "store_summary counts rows and distinct configs")
  check(is.null(store_summary(tempfile())), "store_summary is NULL on a missing file")
  unlink(c(bakp, w1, w3, w4, w5, oldp))
})

# ===========================================================================
cat("evals store: optimizer targets ONE scheme (per-scheme slice of the archive)\n")
# Oracle: filter_evals_to_scheme keeps only the target scheme; NULL keeps all; and
# it composes with filter_evals_to_universe.
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
check(nrow(filter_evals_to_universe(evS, c("201", "202")) |>
             filter_evals_to_scheme("CV0")) == 1,
      "universe + scheme filters compose (201/202 & CV0 -> 1 row)")
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

# ---------------------------------------------------------------------------
# Oracle: CONCURRENT WORKERS MUST NOT ALL PICK THE SAME SEED. choose_config knows only what
# has FINISHED, never what is in flight, so before the first worker stores anything every
# worker sees an identical store. Six workers once spent six evaluations on Prediction1
# before any other seed was touched (observed on the server, 2026-07-31).
dbp5 <- tempfile(fileext = ".sqlite"); con5 <- open_store(dbp5)
wsettings <- function(w) modifyList(optimizer_settings(),
  list(simulate = FALSE, target_domain = NULL, optimize_scheme = "CV00", worker_id = as.character(w)))
n_seed <- length(seed_configs("CV00"))
picks  <- vapply(seq_len(n_seed), function(w) choose_config(con5, wsettings(w))$source,
                 character(1))
check(dplyr::n_distinct(picks) == n_seed,
      "empty store: workers 1..5 pick five DIFFERENT seeds (not all the first)")
check(all(grepl("^seed:", picks)), "all five picks are seeds")
# More workers than seeds -> wrap around rather than error or return NULL.
check(identical(choose_config(con5, wsettings(n_seed + 1))$source, picks[1]),
      "worker 6 wraps around to worker 1's seed")
# Unchanged for a single worker: still the first un-done seed.
check(identical(choose_config(con5, wsettings(1))$source,
                paste0("seed:", names(seed_configs("CV00"))[1])),
      "worker 1 alone still takes the FIRST un-done seed (no change unstaggered)")
# The offset indexes the UN-DONE list, not the full one: with seed 1 stored, worker 1 gets #2.
store_eval(con5, seed_configs("CV00")[[1]], "950", "CV00", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)
check(identical(choose_config(con5, wsettings(1))$source,
                paste0("seed:", names(seed_configs("CV00"))[2])),
      "with seed 1 done, worker 1 takes seed 2 (offset is into the un-done list)")
# With every seed stored, phase 1 must fall through rather than return a seed.
for (k in 2:n_seed)
  store_eval(con5, seed_configs("CV00")[[k]], paste0("95", k), "CV00", 0.3, 40L, "ok",
             build = OPTIMIZER_BUILD)
check(!grepl("^seed:", choose_config(con5, wsettings(3))$source),
      "all seeds done -> phase 1 falls through to random/acquisition")
close_store(con5); unlink(dbp5)

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
cat("claims: (config, trial) is exclusive, the same config on another trial is not\n")
# Oracle: the pipeline is deterministic in (config, trial, scheme), so two workers running one
# triple recompute a known number and count it twice. Selection cannot see in-flight work --
# it reads committed rows -- so the store has to arbitrate.
dbp6 <- tempfile(fileext = ".sqlite"); con6 <- open_store(dbp6)
check(claim_eval(con6, "cfgA", "t1", "CV0", "1"), "the first worker takes (cfgA, t1)")
check(!claim_eval(con6, "cfgA", "t1", "CV0", "2"), "a second worker is refused the SAME pair")
check(claim_eval(con6, "cfgA", "t2", "CV0", "2"),
      "but the same config on a DIFFERENT trial is allowed -- replication wants that")
check(claim_eval(con6, "cfgB", "t1", "CV0", "3"), "and a different config on the same trial is fine")
check(claim_eval(con6, "cfgA", "t1", "CV00", "2"),
      "the scheme is part of the key: the other scheme is a different cell")
release_claim(con6, "cfgA", "t1", "CV0")
check(claim_eval(con6, "cfgA", "t1", "CV0", "9"), "released, it can be taken again")

# Oracle: STALENESS IS A DEAD OWNER, NEVER ELAPSED TIME. A claim held for days is a days-long
# evaluation -- expiring it would start a duplicate of work that is still running. This is the
# distinction from the removed max_eval_minutes: slowness must cost nothing.
invisible(DBI::dbExecute(con6, "UPDATE claims SET ts = '1999-01-01 00:00:00'"))
check(!claim_eval(con6, "cfgA", "t1", "CV0", "4"),
      "an ancient claim whose owner is ALIVE is still not stealable")
# A pid above every platform's pid_max, so it cannot belong to a live process.
invisible(DBI::dbExecute(con6, "UPDATE claims SET pid = 999999999 WHERE config_hash = 'cfgA'"))
check(claim_eval(con6, "cfgA", "t1", "CV0", "4"),
      "a claim whose owner is GONE is reclaimed, however recent")
check(nrow(active_claims(con6, "CV0")) > 0, "active_claims reports what is in flight")
clear_claims(con6)
check(nrow(active_claims(con6)) == 0, "clear_claims empties the table (the leader's startup)")

# Oracle: an ALREADY EVALUATED pair is refused too, not just one in flight. This is the hole a
# claims-only test misses and that a 4-worker simulate run hit repeatedly: the holder commits
# its row and releases, and a worker whose selection read the store just before that commit
# then claims a pair that is finished. The caller's read is always stale, so the store decides.
check(duplicate_cells(con6) == 0, "a fresh store has no duplicated cells")
s1 <- seed_configs("CV0")[[1]]
store_eval(con6, s1, "t9", "CV0", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)
check(!claim_eval(con6, config_hash(s1), "t9", "CV0", "1", OPTIMIZER_BUILD),
      "a pair already evaluated under this build cannot be claimed")
check(claim_eval(con6, config_hash(s1), "t10", "CV0", "1", OPTIMIZER_BUILD),
      "the same config on an unevaluated trial still can -- replication is not blocked")
check(claim_eval(con6, config_hash(s1), "t9", "CV0", "1", "0.0.1-other"),
      "and another build may recompute it: BUILD_CHANGES has to be able to retire a row")

# Oracle: duplicate_cells counts (config, trial, scheme) triples with more than one row --
# reported so the number can be watched for growth once claims are live.
store_eval(con6, s1, "t9", "CV0", 0.3, 40L, "ok", build = OPTIMIZER_BUILD)
check(duplicate_cells(con6) == 1, "and one counts a cell evaluated twice")

# Oracle: a death is RECORDED before the claim is reclaimed. .purge_dead_claims runs on every
# choose_trial, so a dead owner's claim is gone within seconds of the death -- measured on a
# 4-worker run, where killing a claim-holder left no trace at all. claims_reaped is the copy
# taken first, and is the only thing that makes "did anything die" answerable after the fact.
clear_claims(con6)
invisible(DBI::dbExecute(con6, "DELETE FROM claims_reaped"))
check(claim_eval(con6, "cfgD", "t1", "CV0", "7"), "a worker takes a claim")
invisible(DBI::dbExecute(con6, "UPDATE claims SET pid = 999999999 WHERE config_hash = 'cfgD'"))
snap <- claims_snapshot(con6)
check(nrow(snap) == 1 && isFALSE(snap$alive[1]),
      "claims_snapshot shows it as not alive")
check(nrow(claims_snapshot(con6)) == 1,
      "and does NOT purge -- looking at the table must not destroy the evidence")
check(claim_eval(con6, "cfgD", "t1", "CV0", "8"), "the pair is reclaimable by a live worker")
rp <- reaped_claims(con6)
check(nrow(rp) == 1 && identical(rp$worker[1], "7") && identical(rp$trial_id[1], "t1"),
      "and the death is kept in claims_reaped: which worker, on which trial")
# A leader clearing stale claims at startup is the wall-clock-kill case, and counts as a death.
invisible(DBI::dbExecute(con6, "DELETE FROM claims_reaped"))
clear_claims(con6)
check(nrow(reaped_claims(con6)) == 1,
      "clear_claims records what was in flight when the last job was killed")
close_store(con6); unlink(dbp6)

# ===========================================================================
cat("resolve_read_store: prefer the live store, fall back to the durable backup\n")
# Oracle: on a cluster db_path is node-local scratch belonging to the job that wrote it, so a
# diagnostic run from any other allocation finds nothing there. It must read db_backup_path
# rather than report "no store" while a good backup sits on /project.
local({
  live <- tempfile(fileext = ".sqlite"); bak <- tempfile(fileext = ".sqlite")
  mk <- function(p, n) {
    cn <- open_store(p)
    for (i in seq_len(n))
      store_eval(cn, seed_configs("CV0")[[1]], paste0("t", i), "CV0", 0.3, 40L, "ok",
                 build = OPTIMIZER_BUILD)
    close_store(cn)
  }
  st <- function(...) modifyList(optimizer_settings(), list(...))
  quiet <- function(expr) suppressMessages(expr)

  # Both present -> the live one, because it is the only one with rows since the last backup.
  mk(live, 2); mk(bak, 1)
  with_settings <- function(s, expr) {
    real <- optimizer_settings
    optimizer_settings <<- function(...) s
    on.exit(optimizer_settings <<- real, add = TRUE)
    force(expr)
  }
  s_both <- st(db_path = live, db_backup_path = bak)
  check(identical(with_settings(s_both, quiet(resolve_read_store())), live),
        "with both present it reads the LIVE store")

  # Live absent -> the backup. This is the Ceres case: $TMPDIR belongs to another job.
  unlink(live)
  check(identical(with_settings(s_both, quiet(resolve_read_store())), bak),
        "with the live store absent it falls back to the BACKUP")

  # Live present but EMPTY -> still the backup. open_store() creates a full empty schema at any
  # path it is handed, so a node that merely started a job would otherwise shadow the backup.
  cn <- open_store(live); close_store(cn)
  check(file.exists(live), "an empty live store exists on disk")
  check(identical(with_settings(s_both, quiet(resolve_read_store())), bak),
        "a live store with zero rows does not shadow a backup that has some")

  # An explicit path wins over both, and a bad one is named.
  check(identical(quiet(resolve_read_store(bak)), bak), "an explicit path is used as given")
  check(inherits(try(quiet(resolve_read_store(tempfile())), silent = TRUE), "try-error"),
        "an explicit path that does not exist stops")
  # report_*.R hand argv[1] in for their positional form, so a FLAG must not be read as a
  # filename -- both spellings have to work in the same script.
  check(identical(with_settings(s_both, quiet(resolve_read_store("--store=/x/y.sqlite"))), bak),
        "a flag passed as the positional argument is not mistaken for a path")

  # Neither -> one error naming both, not a silent empty store.
  unlink(c(live, bak))
  e <- try(with_settings(s_both, quiet(resolve_read_store())), silent = TRUE)
  check(inherits(e, "try-error") && grepl("no store to read", conditionMessage(attr(e, "condition"))),
        "with neither present it stops, naming both paths")
})

# Oracle: a login node is refused outright. Detected by SLURM's tools being on PATH while we
# hold no allocation.
local({
  stub <- file.path(tempdir(), "loginstub"); dir.create(stub, showWarnings = FALSE)
  writeLines("#!/bin/sh\nexit 0", file.path(stub, "squeue"))
  Sys.chmod(file.path(stub, "squeue"), "0755")
  old_path <- Sys.getenv("PATH"); old_job <- Sys.getenv("SLURM_JOB_ID")
  on.exit({ Sys.setenv(PATH = old_path)
            if (nzchar(old_job)) Sys.setenv(SLURM_JOB_ID = old_job) else Sys.unsetenv("SLURM_JOB_ID")
            unlink(stub, recursive = TRUE) }, add = TRUE)
  Sys.setenv(PATH = paste(stub, old_path, sep = .Platform$path.sep))
  Sys.unsetenv("SLURM_JOB_ID")
  check(.on_login_node(), "squeue on PATH with no allocation is a login node")
  e <- try(resolve_read_store(), silent = TRUE)
  check(inherits(e, "try-error") && grepl("login node", conditionMessage(attr(e, "condition"))),
        "and resolve_read_store refuses there, whatever the store situation")
  Sys.setenv(SLURM_JOB_ID = "1")
  check(!.on_login_node(), "inside an allocation it is not a login node")
})

# ===========================================================================
cat(sprintf("\nTier 1 subtask tests: %d passed, %d failed\n", ok, fail))
if (fail > 0) quit(status = 1)
