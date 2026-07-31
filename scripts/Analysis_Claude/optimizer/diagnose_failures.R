# diagnose_failures.R
#
# Investigate the stored evals that did NOT come back `ok`. The four non-ok statuses do not
# mean the same thing, so they are not treated the same way:
#
#   suspect          infeasible(suspect = TRUE) -- the funnel looks like a BUG (data that
#                    should be visible isn't) rather than genuine infeasibility. Diagnosed.
#   too_few_overlap  score_predictions() got <5 accessions in predictions n observations.
#                    The pipeline RAN; the cliff is one step downstream of
#                    insufficient_geno_overlap. Diagnosed.
#   non_finite       predictions and observations JOINED on name, but every pair had a
#                    NaN/Inf on one side. Nothing about the data is missing -- this is a
#                    modelling failure, and the reason names which side. Diagnosed.
#   error            an uncaught R error. REPORTED, not diagnosed -- an error is usually a
#                    property of the environment (a dead connection, a full disk), not of
#                    the trial or the config, and running diagnose_trial on it would
#                    attribute an environmental failure to the data.
#   infeasible       too_few_train_trials / too_few_train_acc, raised before genotypes are
#                    touched: the config's train_select was too strict. conditions.R
#                    documents this as EXPECTED, so it is off by default (--status to add).
#
# Each diagnosed trial gets diagnose_trial twice:
#   REPRODUCE   the config actually stored on the failing eval (the only run that can
#               reproduce the recorded failure)
#   CONTROL     the permissive default canary_config, to show whether the trial fails
#               for ANY config or only for the stored one
# and then the same probe on a known-good canary trial, which calibrates everything above
# it: if the canary fails too, the run says nothing about the failures.
#
# MUST BE LAUNCHED FROM THIS DIRECTORY. R reads .Renviron from the working directory
# only (no parent walk), and T3_USERNAME / T3_PASSWORD live in ./.Renviron. Launched
# from anywhere else the connection falls back to ~/.Renviron or goes anonymous, every
# catalogue fetch 401s, and the output is a page of plausible-looking cached counts
# followed by "could not build descriptor". The preflight below stops that at second 5
# rather than hour 2.
#
#   cd /path/to/optimizer
#   nohup Rscript diagnose_failures.R --no-replay > logs/diagnose_failures.out 2>&1 &
#
# Options (all optional):
#   --status=suspect,error status(es) to act on; default
#                          suspect,too_few_overlap,non_finite,error.
#                          Add `infeasible` to include the expected-failure bucket.
#   --trials=10260,10259   trials to diagnose; default: every trial in the store carrying
#                          an eval with one of the selected statuses
#   --canary=10676         known-good trial used as the calibration control
#   --scheme=CV00          restrict to one scheme; default: the scheme on the stored eval
#   --max-projects=8       per-project dosage attribution loop (each project loads a WHOLE
#                          dosage matrix -- up to ~10 GB apiece). 0 skips it; the panel
#                          table is the authoritative coverage answer either way.
#   --auth-window-min=60   minutes of padding either side of an authentication error when
#                          listing the evals that ran against the same dead connection.
#   --no-replay            skip the closing run_pipeline() replay. The replay is a COMPLETE
#                          evaluation (hours on a big trial) and only confirms an outcome
#                          the store already recorded -- but it is what prints the funnel.

# Guard the working directory FIRST, with an actionable message. here::i_am() below also
# halts from the wrong directory, but only in the cases where it cannot find the project --
# and cwd = optimizer/logs/ is not one of them (here() still resolves to optimizer, while
# .Renviron, which is read from cwd only, silently does not).
if (!file.exists("diagnose_failures.R"))
  stop("run this FROM the optimizer directory, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  nohup Rscript diagnose_failures.R --no-replay > logs/diagnose_failures.out 2>&1 &\n",
       "  (working directory was: ", getwd(), ")")

library(tidyverse)
here::i_am("diagnose_failures.R")

source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

# --- options ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
split_csv <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",")[[1]])

o_trials  <- split_csv(opt("trials"))
o_status  <- split_csv(opt("status", "suspect,too_few_overlap,non_finite,error"))
o_canary  <- opt("canary", "10676")
o_scheme  <- opt("scheme")
o_maxproj <- as.integer(opt("max-projects", "8"))
o_authpad <- as.numeric(opt("auth-window-min", "60"))
o_replay  <- !("--no-replay" %in% args)

# `error` evals are reported, never diagnosed (see the header). Split the requested
# statuses into the two groups once, here.
REPORT_ONLY  <- "error"
o_diagnose   <- setdiff(o_status, REPORT_ONLY)
o_report     <- intersect(o_status, REPORT_ONLY)

settings <- optimizer_settings()

# --- preflight: prove the network path works BEFORE spending hours ----------
# Everything diagnose_trial prints first (accession counts, observation rows, per-project
# dosage overlaps) is served from the disk cache and looks healthy whether or not the
# connection is authenticated. Only the catalogue fetch goes to the network. So check auth
# and the catalogue here, up front, where a failure is unambiguous.
cat("=== preflight ===\n")
cat("  working directory: ", getwd(), "\n", sep = "")
if (!file.exists(".Renviron"))
  stop("no .Renviron in the working directory -- cd to the optimizer folder and re-run.\n",
       "  R reads .Renviron from the working directory only; from anywhere else the T3\n",
       "  credentials come from ~/.Renviron (or nowhere) and every catalogue fetch 401s.")
cat("  T3_USERNAME:       ", Sys.getenv("T3_USERNAME"), "\n", sep = "")
cat("  T3_PASSWORD:       ", if (nzchar(Sys.getenv("T3_PASSWORD"))) "set" else "EMPTY", "\n", sep = "")

conn <- t3_connect(settings)     # raises t3_bad_credentials if the server issues no token
cat("  login:             OK (token ", nchar(conn$auth_token), " chars)\n", sep = "")

# The catalogue is the fetch that build_trial_descriptor needs and the disk cache does not
# reliably cover (7-day max_age). Warming it here both validates the connection against a
# real query and removes it as a variable from every diagnose_trial below.
cat_tbl <- trial_catalog(conn, settings)
cat("  trial catalogue:   ", nrow(cat_tbl), " trials",
    if ("latitude" %in% names(cat_tbl))
      sprintf(", %d with coordinates", sum(!is.na(cat_tbl$latitude))) else "", "\n", sep = "")
if (!nrow(cat_tbl)) stop("trial catalogue came back EMPTY -- fix that before reading anything below")
cat("  store:             ", settings$db_path, "\n", sep = "")
cat("  statuses:          ", paste(o_status, collapse = ", "),
    if (length(o_report)) sprintf("   (%s: reported, not diagnosed)", paste(o_report, collapse = ", ")) else "",
    "\n", sep = "")
cat("  replay:            ", if (o_replay) "on (run_pipeline -- SLOW)" else "off", "\n", sep = "")

con <- open_store(settings$db_path)
evals <- read_evals(con)

# --- store-wide summary ----------------------------------------------------
# Mirrors state/report.md so this log stands on its own: what was in scope, and what was
# left out. Reading a diagnosis without knowing the denominator invites over-reading it.
cat("\n=== store summary (", nrow(evals), " evals) ===\n", sep = "")
print(evals |> dplyr::count(status, sort = TRUE), n = Inf)
non_ok <- evals |> dplyr::filter(status != "ok")
if (nrow(non_ok)) {
  cat("\n  by status x reason:\n")
  print(non_ok |> dplyr::count(status, reason, sort = TRUE), n = Inf)
}

# --- report-only statuses --------------------------------------------------
# An `error` is a property of the environment far more often than of the data. Print it
# verbatim so it can be read as what it is, and flag the auth signature specifically --
# that one means the connection was dead, so the eval records nothing about the trial.
AUTH_SIG <- "unauthor|401|login did not resolve|must login"
if (length(o_report)) {
  rep <- evals |> dplyr::filter(status %in% o_report)
  cat("\n\n########## REPORTED (not diagnosed): ", paste(o_report, collapse = ", "),
      " -- ", nrow(rep), " eval(s) ##########\n", sep = "")
  if (!nrow(rep)) {
    cat("  none in the store.\n")
  } else {
    print(rep |> dplyr::select(id, ts, trial_id, scheme, status, reason), n = Inf)
    auth <- rep |> dplyr::filter(grepl(AUTH_SIG, reason, ignore.case = TRUE))
    if (nrow(auth)) {
      cat("\n  ", nrow(auth), " of these are AUTHENTICATION failures. The connection was dead,\n",
          "  so these evals say nothing about their trial or config -- delete and re-run them.\n", sep = "")
      # Anything stored in the same window ran against that same connection. Cached reads
      # still succeeded, so a config could have been judged on partial data and recorded a
      # plausible-looking infeasible. Worth re-running; not proof of anything.
      #
      # PAD the range. A single auth error makes range() zero-width, which would report an
      # empty window and imply nothing else was affected -- the opposite of the truth, since
      # one error is the commonest case. ts is character (store_eval writes UTC via format()),
      # so parse before doing arithmetic on it.
      pt   <- function(x) as.POSIXct(x, tz = "UTC")
      pad  <- o_authpad * 60
      w    <- range(pt(auth$ts), na.rm = TRUE) + c(-pad, pad)
      window <- evals |> dplyr::filter(pt(ts) >= w[1], pt(ts) <= w[2], !(id %in% auth$id))
      cat("  window ", format(w[1]), " .. ", format(w[2]), " UTC (+/-", o_authpad,
          " min) holds ", nrow(window), " other eval(s):\n", sep = "")
      if (nrow(window)) print(window |> dplyr::count(status, reason), n = Inf)
      cat("  Those ran against the same connection. Cached reads still worked, so a verdict\n",
          "  there may reflect partial data rather than the trial. Re-run to be sure.\n", sep = "")
    }
  }
}

# --- which trials to diagnose ----------------------------------------------
if (is.null(o_trials)) {
  o_trials <- evals |>
    dplyr::filter(status %in% o_diagnose) |>
    dplyr::distinct(trial_id) |>
    dplyr::pull(trial_id) |>
    as.character()
}
cat("\n  trials to diagnose: ",
    if (length(o_trials)) paste(o_trials, collapse = ", ") else "(none)", "\n", sep = "")

# --- per-trial: reproduce the stored failure, then control against default -----
for (tid in o_trials) {
  cat("\n\n########## TRIAL ", tid, " ##########\n", sep = "")
  bad <- failed_configs(con, tid, status = o_diagnose)
  if (is.null(bad) || !nrow(bad)) {
    cat("  no matching evals stored for this trial -- skipping\n"); next
  }
  if (!is.null(o_scheme)) bad <- bad[bad$scheme == o_scheme, , drop = FALSE]
  if (!nrow(bad)) { cat("  nothing under scheme ", o_scheme, "\n", sep = ""); next }
  print(bad |> dplyr::select(id, status, reason, scheme))

  for (i in seq_len(nrow(bad))) {
    cat("\n--- REPRODUCE: stored config from eval id ", bad$id[i],
        " (", bad$status[i], " / ", bad$reason[i] %||% "?", ") ---\n", sep = "")
    diagnose_trial(tid, settings, conn, cfg = bad$cfg[[i]], scheme = bad$scheme[i],
                   max_projects = o_maxproj, replay = o_replay)
  }

  # Same trial, permissive default config. If this succeeds where the stored config failed,
  # the fault is in the config's methods; if it fails identically, the fault is the trial's
  # data (or, as in the 2026-07-31 run, the environment).
  cat("\n--- CONTROL: permissive default config ---\n")
  diagnose_trial(tid, settings, conn, cfg = canary_config(),
                 scheme = if (!is.null(o_scheme)) o_scheme else bad$scheme[1],
                 max_projects = o_maxproj, replay = o_replay)
}

# --- calibration: a trial we KNOW is feasible ------------------------------
# The control for the whole run. A canary that fails the same way as the trials above means
# the failure is environmental and every verdict above it is unreadable.
if (nzchar(o_canary)) {
  cat("\n\n########## CANARY ", o_canary, " ##########\n", sep = "")
  diagnose_trial(o_canary, settings, conn, cfg = canary_config(),
                 scheme = if (!is.null(o_scheme)) o_scheme else settings$schemes[1],
                 max_projects = o_maxproj, replay = o_replay)
}

close_store(con)
cat("\n=== done ===\n")
