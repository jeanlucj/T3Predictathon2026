# prewarm_indices.R
#
# Fill the two PRIMARY maps so the pipeline can answer its relationship questions locally,
# with no BrAPI wizard call at all:
#
#   project -> accessions    110 projects   ~1 min     cache/proj_acc/proj_acc_<pid>.rds
#   trial   -> accessions    7,645 trials   ~25 min    cache/acc/acc_<sid>.rds
#
# The two "by accession" directions are INVERSIONS of these, derived in RAM, never fetched.
#
# WHY. Profiled on the server, one `trials` wizard call inside .germ_overlap was 440 s of a
# 795 s evaluation -- 55% of the run. The projects wizard is cheap today only because its
# cache happens to hit; on a fresh install it is not.
#
# SAFE TO INTERRUPT, and safe to run while workers are going. Every fetch is cached as it
# completes, so re-running skips what is done; workers contribute to the same files.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   nohup Rscript prewarm_indices.R > logs/prewarm.out 2>&1 &
#
# Options:
#   --only=projects|trials   do one map only (projects first is the best value: ~1 min)
#   --limit=500              stop after this many fetches per map
#   --sleep=0.05             seconds between calls. T3 is shared; do not set 0 casually.
#   --force                  re-fetch keys already present

if (!file.exists("prewarm_indices.R"))
  stop("run this FROM the optimizer directory:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript prewarm_indices.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("prewarm_indices.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", n, "="), "", h[1]) }
o_only  <- opt("only", "both")
o_limit <- suppressWarnings(as.integer(opt("limit", NA)))
o_sleep <- suppressWarnings(as.numeric(opt("sleep", "0.05")))
o_force <- "--force" %in% args

s    <- optimizer_settings()
conn <- t3_connect(s)

# One map. `fetch` pulls and caches a single key; `universe` is the coverage denominator.
fill <- function(label, category, universe, fetch, per_call) {
  idx   <- tryCatch(.inverted_index(s, category), error = function(e) NULL)
  known <- if (is.null(idx)) character() else (attr(idx, "keys") %||% character())
  tried <- .attempted(s, category)
  todo  <- if (o_force) universe else setdiff(universe, union(known, tried))
  cat(sprintf("\n== %s ==\n  universe %d | indexed %d | attempted %d | to fetch %d\n",
              label, length(universe), sum(universe %in% known),
              sum(universe %in% tried), length(todo)))
  if (is.finite(o_limit)) todo <- head(todo, o_limit)
  if (!length(todo)) { cat("  nothing to do\n"); return(invisible(NULL)) }
  cat(sprintf("  estimated %.0f min at %.2f s/call\n", length(todo) * per_call / 60, per_call))

  t0 <- Sys.time(); empty <- 0L; failed <- 0L; batch <- character()
  for (i in seq_along(todo)) {
    k <- todo[i]
    a <- tryCatch(fetch(k), error = function(e) NULL)
    if (is.null(a)) failed <- failed + 1L else {
      if (!length(a)) empty <- empty + 1L
      # Attempted covers empty AND non-empty: an empty result is deliberately not cached, so
      # without this record a genuinely empty key holds coverage below 100% for ever.
      batch <- c(batch, k)
    }
    # Flushed in batches: the attempted file is rewritten whole, so per-key writes would be
    # O(n^2) in IO. An interrupted run simply re-tries the unflushed tail, which is cheap
    # because each fetch is itself cached.
    if (length(batch) >= 100L) { .note_attempted(s, category, batch); batch <- character() }
    if (i %% 100L == 0L || i == length(todo)) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("  [%5d/%5d] %3.0f%%  %.2f s/call  elapsed %.0f min  ETA %.0f min  (empty %d, failed %d)\n",
                  i, length(todo), 100 * i / length(todo), el / i, el / 60,
                  el / i * (length(todo) - i) / 60, empty, failed))
      utils::flush.console()
    }
    if (o_sleep > 0) Sys.sleep(o_sleep)
  }
  if (length(batch)) .note_attempted(s, category, batch)
}

# Projects first: 110 calls for complete coverage is the best value in the run.
if (o_only %in% c("both", "projects"))
  fill("projects -> accessions", "proj_acc", .all_project_ids(conn, s),
       function(p) get_project_accessions(p, conn, s), 0.54)

if (o_only %in% c("both", "trials"))
  fill("trials -> accessions", "acc", unique(as.character(trial_catalog(conn, s)$study_db_id)),
       function(t) get_trial_accessions(t, conn, s), 0.39)

cat("\n== coverage ==\n")
pi <- tryCatch(.project_index(s), error = function(e) NULL)
ti <- tryCatch(.trial_index(s),   error = function(e) NULL)
pu <- tryCatch(.all_project_ids(conn, s), error = function(e) character())
tu <- tryCatch(unique(as.character(trial_catalog(conn, s)$study_db_id)), error = function(e) character())
cat(sprintf("  projects: %d indexed, covers universe: %s\n",
            length(attr(pi, "keys") %||% character()), .index_covers(pi, pu, s, "proj_acc")))
cat(sprintf("  trials  : %d indexed, covers universe: %s\n",
            length(attr(ti, "keys") %||% character()), .index_covers(ti, tu, s, "acc")))
cat("\nBoth TRUE means the pipeline makes no relationship wizard calls at all.\n")
