# evaluation.R
#
# Console tooling for EVALUATING the optimizer, driven by docs/EVALUATION.md from the console.
# Nothing here is called by run_optimizer():
#
#   arm_evaluation(group)   -- debug() every function in a module group, so a pipeline call
#   disarm_evaluation()        drops into the debugger there, one module at a time.
#   eval_groups()           -- print the group menu (which functions, in fast->slow order).
#   peek(x)                 -- a one-object summary tuned to the objects that flow between the
#                              six subtasks, making the silent failure signatures (all-NA
#                              rep/block, degenerate targets, zero rowname overlap, a
#                              non-symmetric K) visible. Returns its input invisibly.
#
# run_optimizer.R's source glob loads this file too. Defining the functions is harmless, but
# disarm_evaluation() before a background run or every pipeline call blocks on the debugger.

library(tidyverse)

# ---------------------------------------------------------------------------
# Module groups, ordered fast/offline -> slow/online; each is a set of functions to arm
# together. The six subtasks are separate so that arming exactly ONE lets run_pipeline() run
# straight through to it; the orchestrators sit in their own `flow` group so arming a subtask
# does not double-break.
# ---------------------------------------------------------------------------
EVAL_GROUPS <- list(
  # ---- offline, instant -------------------------------------------------
  config_space = c("sample_config", "sample_block", "crossover", "mutate_config",
                   "repair_config", "config_hash", "configs_to_features",
                   "feature_schema", "seed_configs", ".make_seed"),
  engine       = c("aggregate_scores", "get_elites", "fresh_random",
                   "propose_candidates", "incumbent_config", "choose_config",
                   "filter_evals_to_universe", "filter_evals_to_scheme",
                   "fit_surrogate", "fit_surrogate_merf", "predict_surrogate",
                   "expected_improvement"),
  store        = c("open_store", "store_eval", "read_evals", "config_to_json",
                   "config_from_json", "write_report", "method_importance",
                   "failure_summary"),
  scoring      = c("sample_trial", "evaluate_config_on_trial", "score_predictions",
                   ".sim_true", ".sim_evaluate", "mask_cv"),

  # ---- online (live BrAPI); slow ---------------------------------------
  data         = c("cached", "t3_connect", "t3_login", ".brapi_try", "trial_catalog",
                   "sample_real_trial", "build_trial_descriptor", "get_trial_accessions",
                   "get_observations", ".obs_tibble", ".obsunits_tibble",
                   "projects_for_accessions", "get_project_dosage",
                   ".ensure_project_vcf", ".vcf_complete", ".vcf_to_dosage",
                   "restore_cache_from_backup", "sync_cache_to_backup"),

  # ---- the six subtasks (arm ONE at a time) ----------------------------
  subtaskA     = c("select_training_trials", ".find_related", ".germ_overlap",
                   ".trial_similarity"),
  subtaskB     = c("mask_cv", "build_targets", ".blue_per_trial", ".per_acc_blue"),
  subtaskC     = c("choose_geno_sources", "projects_for_accessions",
                   "get_project_dosage",".ensure_project_vcf", ".vcf_to_dosage",
                   ".vcf_stat", ".group_by_panel", ".prune_redundant"),
  subtaskD     = c("build_kernel", ".merge_markers", ".qc_markers", ".vanraden",
                   ".best_panel", ".bridge_accessions"),
  subtaskE     = c("train_model", ".fit_sommer_GE", "score_predictions"),
  subtaskF     = c("predict_test"),

  # ---- orchestration / flow glue (its OWN group -> no double-break) -----
  flow         = c("run_pipeline", "optimizer_step"),

  # ---- validation -------------------------------------------------------
  diagnostics  = c("check_canaries", "diagnose_trial", "canary_anchor",
                   "calibrate_canary_trials", "print_calibration",
                   "canary_coverage", "canary_config", "canary_configs",
                   "sweep_rich_trials", ".oracle_variants", ".select_variants")
)

# The whole inner pipeline (all six subtasks, but NOT run_pipeline): arm it and call
# run_pipeline() to break once per subtask, A -> F.
EVAL_GROUPS$pipeline <- unique(unlist(
  EVAL_GROUPS[c("subtaskA", "subtaskB", "subtaskC", "subtaskD", "subtaskE", "subtaskF")],
  use.names = FALSE))

# Evaluation order: fast/offline first, so cheap bugs surface before a download is spent.
EVAL_ORDER <- c("config_space", "engine", "store", "scoring",
                "data", "subtaskA", "subtaskB", "subtaskC", "subtaskD",
                "subtaskE", "subtaskF", "flow", "diagnostics")

# ---------------------------------------------------------------------------
# Arm / disarm debug() over a group. debug(), not debugonce(): the flag persists until
# undebug()'d, and debugonce() can be neither cancelled nor seen by isdebugged(). The
# exists()/isdebugged() guards make both idempotent.
# ---------------------------------------------------------------------------
arm_evaluation <- function(groups) {
  unknown <- setdiff(groups, names(EVAL_GROUPS))
  if (length(unknown))
    stop("unknown eval group(s): ", paste(unknown, collapse = ", "),
         "\n  known groups: ", paste(names(EVAL_GROUPS), collapse = ", "))
  fns   <- unique(unlist(EVAL_GROUPS[groups], use.names = FALSE))
  armed <- 0L
  for (fn in fns)
    if (exists(fn, mode = "function") && !isdebugged(get(fn))) { debug(get(fn)); armed <- armed + 1L }
  message("Armed debug() on ", armed, " function(s) in group(s): ",
          paste(groups, collapse = ", "), ".",
          "\n  In the debugger:  c = finish this function,  n = next line,",
          "  s = step into a call,  Q = quit.",
          "\n  disarm_evaluation() when done  (ALWAYS disarm before a background run).")
  invisible(fns)
}

disarm_evaluation <- function(groups = names(EVAL_GROUPS)) {
  fns <- unique(unlist(EVAL_GROUPS[groups], use.names = FALSE))
  n <- 0L
  for (fn in fns)
    if (exists(fn, mode = "function") && isdebugged(get(fn))) { undebug(get(fn)); n <- n + 1L }
  message("Disarmed debug() on ", n, " function(s).")
  invisible(fns)
}

# Print the group menu in fast->slow order, with each group's functions.
eval_groups <- function() {
  extra <- setdiff(names(EVAL_GROUPS), EVAL_ORDER)          # e.g. `pipeline` aggregate
  for (g in c(EVAL_ORDER, extra))
    cat(sprintf("%-11s : %s\n", g, paste(EVAL_GROUPS[[g]], collapse = ", ")))
  invisible(EVAL_GROUPS)
}

# ---------------------------------------------------------------------------
# Print a compact shape/health line for one object, then return it invisibly so it can sit in a
# pipe. Covers the five object shapes that carry the pipeline's state and the failure signatures
# that otherwise pass silently. `accessions` adds a name-overlap check on a dosage matrix or
# dosage_list (the synonym-mismatch signal).
# ---------------------------------------------------------------------------
peek <- function(x, label = deparse(substitute(x)), accessions = NULL) {
  lab <- paste0("[peek] ", label, ": ")

  # named / plain numeric vector (targets, obs, pred) -----------------------
  if (is.numeric(x) && is.null(dim(x))) {
    n_na <- sum(is.na(x)); fin <- x[is.finite(x)]
    cat(lab, sprintf("numeric[%d]  distinct=%d  NA=%d  range=[%s, %s]\n",
                     length(x), dplyr::n_distinct(x[!is.na(x)]), n_na,
                     if (length(fin)) formatC(min(fin), format = "g", digits = 4) else "NA",
                     if (length(fin)) formatC(max(fin), format = "g", digits = 4) else "NA"),
        sep = "")
    if (length(fin) && dplyr::n_distinct(fin) == 1)
      cat("        !! all finite values identical -- degenerate targets/predictions\n")
    if (!is.null(names(x))) cat("        names:", paste(utils::head(names(x), 4), collapse = ", "),
                                if (length(x) > 4) "..." else "", "\n")
    return(invisible(x))
  }

  # matrix (dosage, K) ------------------------------------------------------
  if (is.matrix(x)) {
    sym <- isTRUE(nrow(x) == ncol(x) && isSymmetric(unname(x)))
    dg  <- if (nrow(x) == ncol(x)) diag(x) else NA_real_
    cat(lab, sprintf("matrix[%d x %d]  square=%s  symmetric=%s  NA=%d",
                     nrow(x), ncol(x), nrow(x) == ncol(x), sym, sum(is.na(x))), sep = "")
    if (nrow(x) == ncol(x))
      cat(sprintf("  diag=[%.3g, %.3g]", min(dg, na.rm = TRUE), max(dg, na.rm = TRUE)))
    cat("\n")
    if (!is.null(accessions) && !is.null(rownames(x))) {
      ov <- length(intersect(rownames(x), as.character(accessions)))
      cat(sprintf("        rowname overlap with accessions = %d / %d rows\n", ov, nrow(x)))
      if (ov == 0 && nrow(x) > 0)
        cat("        !! zero overlap on a non-empty matrix -- likely a synonym/name mismatch\n")
    }
    return(invisible(x))
  }

  # list of matrices (dosage_list) ------------------------------------------
  if (is.list(x) && length(x) && all(vapply(x, is.matrix, logical(1)))) {
    cat(lab, sprintf("list of %d matrices  (n_projects attr=%s)\n",
                     length(x), attr(x, "n_projects") %||% "NA"), sep = "")
    for (nm in names(x) %||% seq_along(x)) {
      d  <- x[[nm]]
      ov <- if (!is.null(accessions)) length(intersect(rownames(d), as.character(accessions))) else NA
      cat(sprintf("        %-14s [%d x %d]%s\n", nm, nrow(d), ncol(d),
                  if (!is.na(ov)) sprintf("  overlap=%d", ov) else ""))
    }
    return(invisible(x))
  }

  # data frame / tibble (train_obs) -- NA count per column ------------------
  if (is.data.frame(x)) {
    cat(lab, sprintf("tibble[%d x %d]\n", nrow(x), ncol(x)), sep = "")
    na <- vapply(x, function(col) sum(is.na(col)), integer(1))
    for (nm in names(x))
      cat(sprintf("        %-16s %s  (NA=%d%s)\n", nm, class(x[[nm]])[1], na[[nm]],
                  if (na[[nm]] == nrow(x) && nrow(x) > 0) " -- ALL NA !!" else ""))
    return(invisible(x))
  }

  # fallback ----------------------------------------------------------------
  cat(lab, sprintf("%s  length=%d\n", paste(class(x), collapse = "/"), length(x)), sep = "")
  invisible(x)
}
