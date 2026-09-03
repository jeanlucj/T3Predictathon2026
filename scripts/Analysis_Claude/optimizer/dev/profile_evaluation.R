# dev/profile_evaluation.R
#
# Where does a 29-minute evaluation actually go? Runs ONE real evaluation with per-subtask
# timing and reports the split, plus the single number that decides whether BLAS threading
# could ever help: CPU time against wall time.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript dev/profile_evaluation.R --trial=11033
#
# Options:
#   --trial=<id>     which focal trial. Default: the most-evaluated trial in the store,
#                    which is the one most likely to have a warm cache.
#   --scheme=CV00    CV scheme. Default CV00.
#   --kernel=        vanRaden_single (default) | em_combine | rkhs_gaussian
#
# READING THE RESULT
#
#   CPU / wall ~= 1.0   The process is compute-bound on ONE core. BLAS threading is either
#                       unavailable or unused -- run dev/check_blas.R next.
#   CPU / wall  < 1.0   The process spends that fraction of its life WAITING -- network,
#                       disk. No BLAS setting will change this; more workers or fewer
#                       downloads would.
#   CPU / wall  > 1.0   Threads are being used; the value is roughly how many.
#
# The subtask split then says which stage to attack. Note the wrappers NEST: build_kernel's
# time includes .vanraden and .qc_markers, which are reported separately as "of which".
#
# This does NOT modify the pipeline -- it shadows the subtask functions with timing wrappers
# in the global environment for the duration of the run.

# The optimizer ROOT, not this script's directory: `.Renviron` -- and so the T3
# credentials -- is read from the WORKING DIRECTORY only, with no parent walk. `settings.R`
# is the marker for that root. here::i_am() below also halts from the wrong place, but only
# when it cannot find the project at all, and a cwd inside the project is not one of those
# cases: here() still resolves while .Renviron silently does not.
if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript dev/profile_evaluation.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("dev/profile_evaluation.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", n, "="), "", h[1]) }
o_scheme <- opt("scheme", "CV00")
o_kernel <- opt("kernel", "vanRaden_single")
o_trial  <- opt("trial")

s <- optimizer_settings()

# Default to the trial with the most stored evaluations: most likely to have a warm cache,
# so the profile reflects compute rather than a cold download.
if (is.null(o_trial)) {
  # Only needed to pick a default: with --trial=<id> the store is never opened, so this runs
  # fine on a machine that has none.
  con <- DBI::dbConnect(RSQLite::SQLite(), resolve_read_store(what = "dev/profile_evaluation.R"))
  e <- tryCatch(DBI::dbReadTable(con, "evals"), error = function(err) NULL)
  DBI::dbDisconnect(con)
  if (!is.null(e) && nrow(e)) {
    tt <- sort(table(e$trial_id[e$status == "ok"]), decreasing = TRUE)
    o_trial <- names(tt)[1]
  }
  if (is.null(o_trial)) stop("no trial in the store; pass --trial=<id>")
}

# --- timing wrappers --------------------------------------------------------
TIMES <- new.env(parent = emptyenv())
wrap <- function(name) {
  if (!exists(name, envir = globalenv())) return(invisible(FALSE))
  orig <- get(name, envir = globalenv())
  assign(name, function(...) {
    t0 <- Sys.time()
    on.exit({
      d <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      assign(name, (get0(name, envir = TIMES, ifnotfound = 0)) + d, envir = TIMES)
    }, add = TRUE)
    orig(...)
  }, envir = globalenv())
  invisible(TRUE)
}
TOP <- c("select_training_trials", "get_observations", "build_targets",
         "choose_geno_sources", "build_kernel", "train_model", "predict_test")
# Inner functions, so the big stages break down rather than showing as one opaque block.
# A first run on real data attributed 79% of the wall clock to choose_geno_sources but only
# 17% to the inner functions then wrapped -- these close that gap.
INNER <- c("projects_for_accessions", "get_project_dosage", ".dosage_thin_plan",
           ".group_by_panel", ".prune_redundant", ".merge_markers",
           ".find_related", ".trial_similarity", ".qc_markers", ".vanraden",
           # Inside .find_related. Two hypotheses about where its 435 s went -- the projects
           # wizard, then the overlap loop -- were both wrong, so instrument every piece
           # rather than guess again. .germ_overlap is the whole body; the rest partition it.
           ".germ_overlap", ".trial_index", "trial_catalog", "get_trial_accessions",
           # Every BrAPI call goes through .brapi_try, so its total separates "talking to T3"
           # from local compute. It spans ALL stages, so it will overlap the others.
           ".brapi_try")
for (f in c(TOP, INNER)) wrap(f)

# --- run one evaluation -----------------------------------------------------
cat("BLAS        :", extSoftVersion()[["BLAS"]], "\n")
for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"))
  cat(sprintf("%-13s %s\n", v, Sys.getenv(v, "<unset>")))
cat("trial       :", o_trial, " scheme:", o_scheme, " kernel:", o_kernel, "\n\n")

conn  <- t3_connect(s)
trial <- build_trial_descriptor(o_trial, conn, s)
cfg   <- modifyList(canary_config(), list(kernel.method = o_kernel))

cat("running one evaluation -- this takes as long as a real one...\n\n")
t_all <- system.time({
  res <- tryCatch(run_pipeline(cfg, trial, o_scheme, s, conn),
                  error = function(e) { cat("  pipeline failed: ", conditionMessage(e),
                                            "\n  (timings below are still valid up to that point)\n",
                                            sep = ""); NULL })
})

# --- report -----------------------------------------------------------------
wall <- t_all[["elapsed"]]; cpu <- t_all[["user.self"]] + t_all[["sys.self"]]
get_t <- function(n) get0(n, envir = TIMES, ifnotfound = 0)
cat(sprintf("\n%-26s %9s %7s\n", "stage", "seconds", "% wall"))
for (f in TOP) {
  v <- get_t(f)
  if (v > 0) cat(sprintf("%-26s %9.1f %6.1f%%\n", f, v, 100 * v / wall))
}
cat(sprintf("%-26s %9s %7s\n", "  -- of which (nested) --", "", ""))
for (f in INNER) {
  v <- get_t(f)
  if (v > 0) cat(sprintf("  %-24s %9.1f %6.1f%%\n", f, v, 100 * v / wall))
}
acc <- sum(vapply(TOP, get_t, numeric(1)))
cat(sprintf("\n%-26s %9.1f\n", "accounted for", acc))
cat(sprintf("%-26s %9.1f\n", "TOTAL wall", wall))
cat(sprintf("%-26s %9.1f  (unattributed: setup, masking, scoring)\n", "  unaccounted", wall - acc))

cat(sprintf("\nCPU %.1fs / wall %.1fs = %.2f\n", cpu, wall, cpu / max(wall, 1e-9)))
verdict <- if (cpu / wall < 0.7) {
  "  I/O-BOUND -- most of the time is spent WAITING (network, disk). BLAS threading cannot help."
} else if (cpu / wall < 1.4) {
  "  COMPUTE-BOUND ON ONE CORE -- threads are not being used. Run dev/check_blas.R next."
} else {
  sprintf("  THREADED, roughly %.0f-way.", cpu / wall)
}
cat(verdict, "\n")
