# peek_backup.R
#
# Pre-flight: what is in the durable backup, what THIS run will see of it, and which phase
# iteration 1 will land in. Read-only, no network, seconds. Run it before submitting -- an empty
# backup and a full one produce completely different first hours.
#
# The starting phase is obtained by merging the backup into a TEMP COPY of the working store and
# calling the real choose_config() on it, so the answer cannot drift from what the optimizer
# actually does. The same run reports what the leader's restore will insert.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript peek_backup.R
#
# SAFETY: opens copies, never the live store (peek_failures.R explains why that matters).

if (!file.exists("peek_backup.R"))
  stop("run this FROM the optimizer directory:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript peek_backup.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("peek_backup.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

# Not resolve_read_store(): comparing the live store against the backup IS this script's job,
# so it needs both paths, not a choice between them. Only the login-node refusal is shared.
if (.on_login_node()) stop(.login_node_message("peek_backup.R"), call. = FALSE)
s <- optimizer_settings()
hr <- function(t) cat("\n", t, "\n", strrep("-", nchar(t)), "\n", sep = "")
mb <- function(b) if (is.na(b)) "-" else sprintf("%.1f MB", b / 1e6)
age <- function(t) if (is.na(t)) "-" else {
  h <- as.numeric(difftime(Sys.time(), t, units = "hours"))
  if (h < 48) sprintf("%.1f h ago", h) else sprintf("%.1f days ago", h / 24)
}
tab <- function(x) if (!length(x)) "  (none)\n" else
  paste0("  ", format(names(x), width = max(nchar(names(x)))), "  ", as.integer(x), "\n")

# --- files -----------------------------------------------------------------
hr("files")
bak  <- store_summary(s$db_backup_path)
work <- store_summary(s$db_path)
cat(sprintf("  backup : %s\n", s$db_backup_path %||% "(none configured)"))
cat(sprintf("           %s\n",
    if (!is.null(bak)) sprintf("%d rows, %s, written %s", bak$rows, mb(bak$bytes), age(bak$mtime))
    else if (is.null(s$db_backup_path) || !nzchar(s$db_backup_path))
      "nothing is copied off this machine -- expected for a local run"
    else "MISSING -- there is nothing to restore from"))
cat(sprintf("  store  : %s\n", s$db_path))
cat(sprintf("           %s\n", if (is.null(work)) "absent -- normal on a fresh node"
                               else sprintf("%d rows, %s", work$rows, mb(work$bytes))))

# Cache counts by category: the index maps are what a pre-warm builds and are tiny; dosage is
# the bulk. Reported separately because they are restored and used independently.
cdir <- s$cache_backup_dir
if (!is.null(cdir) && nzchar(cdir) && dir.exists(cdir)) {
  n <- vapply(list.dirs(cdir, recursive = FALSE),
              function(d) length(list.files(d)), integer(1))
  names(n) <- basename(names(n))
  cat(sprintf("  cache  : %s  (%d files)\n", cdir, sum(n)))
  cat(tab(n[order(-n)][seq_len(min(8, length(n)))]), sep = "")
} else {
  cat(sprintf("  cache  : %s\n", if (is.null(cdir) || !nzchar(cdir)) "(none configured)"
                                 else paste(cdir, "-- MISSING")))
}

if (is.null(bak) && is.null(work)) {
  cat("\nNothing to report: no store and no backup. A run would start at iteration 1,",
      "\nphase seed, with the five submissions as its first evaluations.\n")
  quit(status = 0)
}

# --- what is in the backup -------------------------------------------------
if (!is.null(bak)) {
  hr("backup contents")
  cat(sprintf("  %d rows over %d configs x %d trials\n", bak$rows, bak$n_config, bak$n_trial))
  cat(sprintf("  written between %s and %s\n", bak$ts_range[1], bak$ts_range[2]))
  cat("  by scheme:\n"); cat(tab(bak$by_scheme), sep = "")
  cat("  by status:\n"); cat(tab(bak$by_status), sep = "")
  cat("  by build:\n");  cat(tab(bak$by_build),  sep = "")
  # dosage_budget changes marker density and therefore scores, and is not a config parameter:
  # more than one value here means the rows are not all comparable.
  if ("dosage_budget" %in% names(bak$evals)) {
    db <- table(bak$evals$dosage_budget, useNA = "ifany")
    if (length(db) > 1) {
      cat("  by dosage_budget (rows under different budgets are NOT comparable):\n")
      cat(tab(db), sep = "")
    }
  }
}

# --- merge the backup into a temp copy of the store, then ask the real code ---
hr("what this run will see")
tmp <- file.path(tempdir(), "peek_backup_merge.sqlite")
unlink(paste0(tmp, c("", "-wal", "-shm")))
if (!is.null(work)) invisible(.copy_store_with_sidecars(s$db_path, tmp))
merged <- restore_store_from_backup(modifyList(s, list(db_path = tmp)))
con <- open_store(tmp)

evals <- read_evals(con)
slice <- evals |>
  filter_evals_to_domain(if (isTRUE(s$simulate)) NULL else s$target_domain) |>
  filter_evals_to_scheme(s$optimize_scheme) |>
  filter_evals_to_build(s$build %||% OPTIMIZER_BUILD)
cat(sprintf("  after the leader's restore: %d rows\n", nrow(evals)))
cat(sprintf("  in this domain + scheme (%s) + build (%s): %d rows over %d configs x %d trials\n",
            s$optimize_scheme, s$build %||% OPTIMIZER_BUILD, nrow(slice),
            dplyr::n_distinct(slice$config_hash), dplyr::n_distinct(slice$trial_id)))
if (nrow(evals) > nrow(slice))
  cat(sprintf("  %d row(s) are out of slice -- real work, but not evidence for THIS run\n",
              nrow(evals) - nrow(slice)))

# --- how it will start -----------------------------------------------------
hr("iteration 1")
pick <- tryCatch(choose_config(con, s), error = function(e) NULL)
if (is.null(pick)) {
  cat("  could not determine the starting phase\n")
} else {
  cat(sprintf("  phase: %s\n", pick$source))
  cat(switch(sub(":.*$", "", pick$source),
    seed        = "  -> the five submissions are not all evaluated here yet; it will baseline them first\n",
    replicate   = "  -> a config is short of config_replication evaluations; it gets another trial\n",
    random_init = sprintf("  -> still exploring: fewer than n_random_init (%d) scored configs in slice\n",
                          s$n_random_init),
    "  -> the surrogate is driving; the search resumes where it left off\n"))
}

agg <- aggregate_scores(slice)

# The variance decomposition, from the fit aggregate_scores just paid for. Worth printing here
# and not only in report.md: this is a FRESH process, so it reflects the current code even when
# the workers writing the store are running an older build and their report.md cannot.
#
# sd_resid is on the unit-weight scale (the fit weights by n_test - 3), so it is not comparable
# with the other two until divided by sqrt(median weight) -- LESSONS #19. Same rendering as
# R/report.R.
local({
  vc   <- attr(agg, "var_comps")
  note <- attr(agg, "estimator_note")
  w    <- attr(agg, "median_weight") %||% NA_real_
  if (!is.finite(w)) w <- 1
  cat(sprintf("\n  config score estimator: %s\n", attr(agg, "estimator") %||% "pooled"))
  if (!is.null(vc) && all(is.finite(vc))) {
    cat(sprintf("    sd_trial %.3f  sd_config %.3f  sd_resid %.3f per eval (%.3f at unit weight)\n",
                vc[["sd_trial"]], vc[["sd_config"]], vc[["sd_resid"]] / sqrt(max(1, w)),
                vc[["sd_resid"]]))
    # The number the replication argument rests on: trials varying more than configs is why
    # config_replication exists at all.
    cat(sprintf("    trial effect is %.1fx the config effect\n",
                vc[["sd_trial"]] / max(vc[["sd_config"]], 1e-9)))
    if (!is.null(note)) cat("    (fit warned: ", note, ")\n", sep = "")
  } else {
    cat("    no variance components: ", note %||% "reason not recorded", "\n", sep = "")
    cat("    -- and so no `se`, which means .contenders() is empty and contender\n")
    cat("       replication is idle on this slice.\n")
  }
  cand <- .contenders(agg, s$contender_z %||% 1, k = 8L)
  cat(sprintf("    contenders: %d\n", length(cand)))
})

inc <- if (nrow(agg)) incumbent_config(agg, s$incumbent_min_reps) else NULL
if (!is.null(inc)) {
  cat(sprintf("\n  incumbent so far: %.3f over %d trial(s)\n", inc$mean_score, inc$n_ok))
  cat(format_config(inc$config), sep = "\n")
} else {
  cat("\n  no incumbent yet in this slice\n")
}

close_store(con)
unlink(paste0(tmp, c("", "-wal", "-shm")))
if (!is.null(merged))
  cat(sprintf("\nThe leader will insert %d row(s) and skip %d already present.\n",
              merged$inserted, merged$skipped))
