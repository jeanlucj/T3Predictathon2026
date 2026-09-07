# dev/analyze_failures.R
#
# Where do the failures come from, and did the DOMAIN change under us?
#
# Written when a fresh run's failure rate came back at ~15% against a much lower rate in the
# previous store. Everything here is arithmetic on rows the store already holds -- no network,
# no genotypes, no re-running -- so it answers in seconds what tools/diagnose_failures.R would
# spend days on. Reach for that one only if a specific trial is still unexplained after this.
#
# TWO ARTEFACTS THIS EXISTS TO GET PAST, both in the report's own failure table:
#
#  1. failure_summary()'s per-method rate charges a failed eval to ALL SIX of its subtask
#     methods, so a subtask with no causal role shows the base rate. Five of the six columns in
#     that table are therefore noise; only the one with real spread means anything.
#  2. A marginal rate cannot separate cause from correlation. If the search co-selects one
#     method with another, both marginals move. Crossing the method with the TRIAL does
#     separate them: a method that fails everywhere is the method, one that fails on three
#     trials is those trials.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript dev/analyze_failures.R
#
# Options:
#   --store=<path>    the store to analyse   (default: resolve_read_store(), as the tools do)
#   --archive=<path>  the PREVIOUS store, for the before/after and domain diff
#                     (default: newest $OPTIMIZER_HOME/state/evals_archive_*.sqlite)
#   --no-archive      skip sections 4 and 5
#   --z=1,2           contender_z values to count at (section 6)
#   --top=20          rows to print in the longer tables
#
# SAFETY: opens a COPY of each store, with its -wal/-shm sidecars, never the live file --
# tools/inspect_failures.R explains why that matters. Read-only throughout; safe beside a run.

# The optimizer ROOT, not this script's directory: `.Renviron` -- and so the T3
# credentials -- is read from the WORKING DIRECTORY only, with no parent walk. `settings.R`
# is the marker for that root. here::i_am() below also halts from the wrong place, but only
# when it cannot find the project at all, and a cwd inside the project is not one of those
# cases: here() still resolves while .Renviron silently does not.
if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript dev/analyze_failures.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("dev/analyze_failures.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", n, "="), "", h[1]) }
o_top <- suppressWarnings(as.integer(opt("top", "20")))
o_z   <- suppressWarnings(as.numeric(trimws(strsplit(opt("z", "1,2"), ",")[[1]])))
o_noa <- "--no-archive" %in% args

s <- optimizer_settings()
hr <- function(t) cat("\n\n", t, "\n", strrep("=", nchar(t)), "\n", sep = "")

# ---- load, always from a copy ---------------------------------------------
read_store <- function(path, what) {
  tmp <- file.path(tempdir(), paste0("analyze_", what, ".sqlite"))
  con <- open_store(.copy_store_with_sidecars(path, tmp))
  on.exit(close_store(con), add = TRUE)
  list(evals = read_evals(con),
       runs  = if ("runs" %in% DBI::dbListTables(con))
                 tryCatch(DBI::dbReadTable(con, "runs"), error = function(e) NULL) else NULL)
}

# `method_of` reads one subtask's method off each row's stored config. config_from_json is the
# same parser the optimizer uses, so a config that the search can produce always parses here.
method_of <- function(evals, key)
  unname(vapply(evals$config_json, function(j)
    tryCatch(as.character(config_from_json(j)[[key]] %||% NA_character_),
             error = function(e) NA_character_), character(1)))

param_of <- function(evals, key)
  unname(vapply(evals$config_json, function(j)
    tryCatch(suppressWarnings(as.numeric(config_from_json(j)[[key]] %||% NA_real_)),
             error = function(e) NA_real_), numeric(1)))

store_path <- resolve_read_store(opt("store"), "dev/analyze_failures.R")
cur <- read_store(store_path, "cur")
cat("store   : ", store_path, "  (", nrow(cur$evals), " rows)\n", sep = "")

arch_path <- opt("archive")
if (is.null(arch_path) && !o_noa) {
  cands <- Sys.glob(file.path(s$db_backup_path %||% "state/x" |> dirname(),
                              "evals_archive_*.sqlite"))
  if (length(cands)) arch_path <- cands[order(file.info(cands)$mtime, decreasing = TRUE)][1]
}
arch <- NULL
if (!o_noa && !is.null(arch_path) && file.exists(arch_path)) {
  arch <- read_store(arch_path, "arch")
  cat("archive : ", arch_path, "  (", nrow(arch$evals), " rows)\n", sep = "")
} else if (!o_noa) {
  cat("archive : none found -- sections 4 and 5 skipped (--archive=<path> to name one)\n")
}

e <- cur$evals
e$failed <- e$status != "ok"
if (!nrow(e)) stop("the store has no rows.")

# ---- 1. failures by trial --------------------------------------------------
hr("1. Failures by trial")
by_trial <- e |>
  group_by(trial_id, study_name, program_name) |>
  summarise(n = n(), n_fail = sum(failed), rate = mean(failed), .groups = "drop") |>
  arrange(desc(rate))
print(by_trial, n = Inf)
cat(sprintf("\n  %d trial(s); overall %d/%d = %.3f\n",
            nrow(by_trial), sum(e$failed), nrow(e), mean(e$failed)))
# Concentration is the finding, or its absence is: if a handful of trials carry the failures,
# the domain is the story; if every trial fails at the overall rate, it is the config space.
top_share <- by_trial |> arrange(desc(n_fail)) |> mutate(cum = cumsum(n_fail) / sum(n_fail))
cat(sprintf("  top 3 trials hold %.0f%% of all failures; top 5 hold %.0f%%\n",
            100 * top_share$cum[min(3, nrow(top_share))],
            100 * top_share$cum[min(5, nrow(top_share))]))

cat("\nreason x trial (non-ok only):\n")
print(e |> filter(failed) |> count(study_name, status, reason) |>
        arrange(desc(n)) |> head(o_top), n = Inf)

# ---- 2. is it train_select, and where? -------------------------------------
hr("2. Cause attribution: train_select x trial")
e$train_select <- method_of(e, "train_select.method")
cat("marginal (this is what the report shows):\n")
print(e |> group_by(train_select) |>
        summarise(n = n(), n_fail = sum(failed), rate = mean(failed), .groups = "drop") |>
        arrange(desc(rate)), n = Inf)

cat("\nwithin trial -- a method that fails EVERYWHERE is the method;\n")
cat("one that fails on a few trials is those trials:\n")
print(e |> group_by(study_name, train_select) |>
        summarise(n = n(), rate = round(mean(failed), 2), .groups = "drop") |>
        pivot_wider(names_from = train_select, values_from = c(n, rate),
                    values_fill = list(n = 0L)) , n = Inf, width = Inf)

# same_program can only work where the program HAS other trials in the domain. This counts
# them, from the rows themselves rather than from the catalogue.
cat("\ndomain trials per program (same_program's training pool):\n")
print(e |> distinct(trial_id, program_name) |> count(program_name, name = "trials_in_domain") |>
        arrange(trials_in_domain), n = Inf)

# ---- 3. the `constant` bucket ----------------------------------------------
hr("3. `constant` failures: is predict_post.min_overlap the cause?")
cst <- filter(e, status == "constant")
if (!nrow(cst)) cat("  none in this store.\n") else {
  cst$min_overlap <- param_of(cst, "predict_post.min_overlap")
  cst$train_in    <- .funnel_get(cst$detail, "train_in")
  cat(sprintf("  %d constant row(s).\n", nrow(cst)))
  cat("  predict_test() returns rep(mean(targets), .) when train_in < min_overlap, which is a\n")
  cat("  zero-variance prediction and therefore a guaranteed `constant` -- after a full run.\n\n")
  print(cst |> count(min_overlap, name = "n") |> arrange(desc(n)), n = Inf)
  known <- sum(is.finite(cst$train_in) & is.finite(cst$min_overlap))
  if (known)
    cat(sprintf("\n  of %d row(s) with both recorded, %d had train_in < min_overlap\n",
                known, sum(cst$train_in < cst$min_overlap, na.rm = TRUE)))
  else
    cat("\n  no funnel on these rows (scoring statuses carry a reason, not a funnel), so the\n",
        "  min_overlap distribution above is the evidence -- compare it to the config space.\n", sep = "")
}

# ---- 4. domain diff, keyed on trial_id -------------------------------------
hr("4. Domain diff: which trials changed?")
if (is.null(arch)) cat("  no archive -- skipped.\n") else {
  # trial_id, NEVER study_name: a renamed trial keeps its id, and `study_name` on a row is the
  # name AS OF that evaluation -- which is what makes the two stores together a rename log.
  nm <- function(ev) ev |> distinct(trial_id, study_name) |>
    group_by(trial_id) |> summarise(name = paste(unique(study_name), collapse = " | "),
                                    .groups = "drop")
  a <- nm(arch$evals); b <- nm(e)
  added   <- anti_join(b, a, by = "trial_id")
  removed <- anti_join(a, b, by = "trial_id")
  common  <- inner_join(a, b, by = "trial_id", suffix = c("_old", "_new"))
  cat(sprintf("  old %d trial(s), new %d; %d in common, %d added, %d removed\n\n",
              nrow(a), nrow(b), nrow(common), nrow(added), nrow(removed)))
  cat("ADDED (never evaluated in the old store):\n");   print(added, n = Inf)
  cat("\nREMOVED (in the old store, not in the new):\n"); print(removed, n = Inf)
  renamed <- filter(common, name_old != name_new)
  if (nrow(renamed)) {
    cat("\nRENAMED (same trial id, different study_name -- the rename log):\n")
    print(renamed, n = Inf)
  }
  if (!is.null(arch$runs) || !is.null(cur$runs)) {
    cat("\nrecorded run rows (universe + target_domain as the run itself saw them):\n")
    for (nmx in c("archive", "current")) {
      r <- if (nmx == "archive") arch$runs else cur$runs
      if (is.null(r) || !nrow(r)) { cat("  ", nmx, ": no runs table\n", sep = ""); next }
      for (i in seq_len(nrow(r))) {
        u <- tryCatch(jsonlite::fromJSON(r$universe_json[i]), error = function(e) NULL)
        td <- tryCatch(jsonlite::fromJSON(r$settings_json[i])$target_domain, error = function(e) NULL)
        cat(sprintf("  %s run %s  build %s  started %s  catalog %s  universe %d trial(s)\n",
                    nmx, r$run_id[i], r$build[i] %||% "?", r$started_ts[i] %||% "?",
                    r$catalog_ts[i] %||% "?", length(u$id %||% character())))
        if (!is.null(td$trials))
          cat("    target_domain$trials: ", length(td$trials), " name(s)\n", sep = "")
      }
    }
  }
}

# ---- 5. before / after, like for like ---------------------------------------
hr("5. Failure rate: before vs after")
if (is.null(arch)) cat("  no archive -- skipped.\n") else {
  ae <- arch$evals; ae$failed <- ae$status != "ok"
  ae$train_select <- method_of(ae, "train_select.method")
  rate <- function(d) sprintf("%d/%d = %.3f", sum(d$failed), nrow(d), mean(d$failed))
  cat("  ALL trials      old ", rate(ae), "   new ", rate(e), "\n", sep = "")
  shared <- intersect(unique(ae$trial_id), unique(e$trial_id))
  ac <- filter(ae, trial_id %in% shared); bc <- filter(e, trial_id %in% shared)
  cat("  COMMON trials   old ", rate(ac), "   new ", rate(bc), "\n", sep = "")
  cat("\n  This is the decisive line. If the COMMON-trial rate is flat while the ALL-trial rate\n")
  cat("  rose, the increase is entirely the newly admitted trials and the code is exonerated.\n")
  newonly <- filter(e, !trial_id %in% shared)
  if (nrow(newonly)) cat("  NEW-only trials new ", rate(newonly), "\n", sep = "")
  cat("\nby train_select, common trials only:\n")
  print(bind_rows(mutate(ac, era = "old"), mutate(bc, era = "new")) |>
          group_by(era, train_select) |>
          summarise(n = n(), rate = round(mean(failed), 3), .groups = "drop") |>
          pivot_wider(names_from = era, values_from = c(n, rate)), n = Inf)
}

# ---- 6. how many contenders, really -----------------------------------------
hr("6. Contenders (the report caps its list at k = 8)")
# The same filters write_report() applies, so the counts are comparable to the report's.
# write_report() reads the universe off settings$trial_universe, which run_optimizer sets after
# pinning; a standalone script has to read the run row instead, as tools/check_backup.R does.
rid <- run_id_for(s, s$build %||% OPTIMIZER_BUILD)
uni <- NULL
if (!is.null(cur$runs) && rid %in% cur$runs$run_id) {
  uj  <- cur$runs$universe_json[match(rid, cur$runs$run_id)]
  uni <- tryCatch(as.character(jsonlite::fromJSON(uj)$id), error = function(e) NULL)
}
cat("  run ", rid, ": ", if (length(uni)) paste0(length(uni), "-trial universe")
    else "no run row -- not filtering to a universe", "\n", sep = "")
ev <- e |> filter_evals_to_scheme(s$optimize_scheme) |>
  filter_evals_to_build(s$build %||% OPTIMIZER_BUILD)
if (length(uni)) ev <- filter_evals_to_universe(ev, uni)
agg <- aggregate_scores(ev)
cat(sprintf("  %d distinct config(s) in the slice; %d with an se (>= 2 evals)\n",
            nrow(agg), sum(is.finite(agg$se))))
# aggregate_scores shrinks the config BLUPs in proportion to replication, so a slice with
# little replication collapses every config onto the grand mean -- correct, and uninformative.
# Say so, or "everything is a contender" reads as a finding rather than as no evidence.
cat(sprintf("  estimator: %s;  mean_score spread %.4f (%.4f to %.4f), median se %.4f\n",
            attr(agg, "estimator") %||% "?",
            diff(range(agg$mean_score, na.rm = TRUE)),
            min(agg$mean_score, na.rm = TRUE), max(agg$mean_score, na.rm = TRUE),
            stats::median(agg$se, na.rm = TRUE)))
if (isTRUE(diff(range(agg$mean_score, na.rm = TRUE)) <
           stats::median(agg$se, na.rm = TRUE)))
  cat("  ** spread < median se: the BLUPs have shrunk together, so the counts below mean\n",
      "     'nothing is ruled out yet', NOT 'these configurations are all equally good'.\n", sep = "")
for (z in o_z) {
  n_all <- length(.contenders(agg, z, k = .Machine$integer.max))
  cat(sprintf("  contender_z = %-4g -> %d configuration(s) in play  (report shows at most 8)\n",
              z, n_all))
}
ct <- .contenders(agg, o_z[1], k = o_top)
cat(sprintf("\ntop %d by UCB at z = %g:\n", length(ct), o_z[1]))
print(agg |> filter(config_hash %in% ct) |>
        mutate(ucb = mean_score + o_z[1] * se) |>
        select(config_hash, n, n_ok, mean_score, se, ucb) |>
        arrange(desc(ucb)), n = Inf)

cat("\n\ndone.\n")
