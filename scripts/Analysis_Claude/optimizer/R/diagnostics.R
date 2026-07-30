# diagnostics.R
#
# Tools to answer: "is this trial REALLY infeasible, or is a bug hiding the data?"
#
# The optimizer continues past infeasible (trial, config) pairs, so a silent bug
# -- a name/key mismatch, a parsing error, a wrong column -- can masquerade as a
# trial that "can't be predicted". Three lines of defence:
#
#   1. The failure log flags funnel CLIFFS (many accessions in, ~zero genotyped
#      overlap out) as status "suspect" rather than "infeasible" (see pipeline.R).
#   2. check_canaries(): run the most permissive config on KNOWN-feasible trials.
#      A canary that comes back infeasible can only mean a bug -- it is an oracle.
#   3. diagnose_trial(): replay one trial with the full data funnel printed AND an
#      independent re-derivation of the raw counts, so a discrepancy between "what
#      the pipeline saw" and "what is actually on T3" is staring you in the face.

library(tidyverse)

# The maximally permissive configuration: broadest training set, no outlier
# trimming, every covering genotyping project, a single VanRaden GRM with loose
# QC, plain GBLUP, no blending. If real data exists for a trial, THIS should
# predict it. A failure here is a bug signal, not genuine infeasibility.
canary_config <- function() {
  .make_seed(list(
    train_select.method = "accession_overlap",
    train_select.primary_only = "yes",
    train_select.primary_min = 2,
    train_select.secondary_min = 8,
    pheno_prep.method   = "raw_mean",
    geno_select.method  = "all_projects",
    kernel.method       = "vanRaden_single",
    kernel.maf          = 0.01,
    kernel.max_missing  = 0.80,
    model.method        = "gblup_rrblup",
    model.lambda_select = "reml",
    predict_post.method = "direct_blup"
  ))
}

# ---------------------------------------------------------------------------
# FROZEN coverage canary configs (Stage 2 of the calibrate-then-freeze bootstrap).
#
# One config per focal trial, assigned from the calibration so that across the
# canaries every one of the 18 subtask methods AND every behaviour-changing
# parameter level is exercised at least once -- so a bug in ANY branch trips the
# oracle, not just the one a single permissive config happens to take.
#
# Coverage is driven by FOUR strong, data-rich trials -- Aurora 10673, Big6 10675,
# CornellMaster 10676, YT_Urb 10677 (all verified feasible and anchor-agreeing at
# calibration) -- which between them carry every demanding branch. The rest (strong OHRWW
# 10679 / TCAP 10680 and the three weak ones) get a light best_single_project filler: they
# confirm "still predictable" and are not needed for coverage. Keyed by studyDbId.
#
# Which trial exercises which branch is tabulated in EVALUATION.md sec. 9 (kept there rather
# than duplicated here, since a config edit below must be reflected in the reader's table).
# `canary_coverage()` is the machine-checkable version and is what to trust.
canary_configs <- function() {
  mk <- function(...) .make_seed(list(...))
  list(
    `10673` = mk(  # Aurora
      train_select.method = "accession_overlap", train_select.primary_only = "no",
      train_select.primary_min = 4, train_select.secondary_min = 12,
      pheno_prep.method = "two_stage_blup",
      geno_select.method = "focal_plus_onehop", geno_select.min_bridge = 2,
      kernel.method = "em_combine", kernel.impute = "mean_round",
      model.method = "gblup_sommer_GE", model.include_E = "yes",
      predict_post.method = "direct_blup"),
    `10675` = mk(  # Big6
      train_select.method = "top_k_similar", train_select.k = 15,
      train_select.similarity = "both",
      pheno_prep.method = "blue_lm", pheno_prep.z_thr = 3,
      pheno_prep.ge_weighting = "env_gaussian", pheno_prep.ge_bandwidth = 1,
      pheno_prep.standardize = "per_trial_z",
      geno_select.method = "best_single_project",
      kernel.method = "rkhs_gaussian", kernel.rkhs_theta = 1, kernel.impute = "mean",
      model.method = "rkhs", model.include_E = "no",
      predict_post.method = "cond_expectation", predict_post.blend_obs_w = 0.3),
    `10676` = mk(  # CornellMaster (fewest projects -> hosts heavy all_projects)
      train_select.method = "same_program", train_select.prog_cap = 15,
      pheno_prep.method = "trial_center", pheno_prep.z_thr = 3,
      geno_select.method = "all_projects",
      kernel.method = "vanRaden_single",
      model.method = "gblup_rrblup", model.lambda_select = "loo",
      predict_post.method = "direct_blup"),
    `10677` = mk(  # YT_Urb
      train_select.method = "accession_overlap", train_select.primary_only = "yes",
      train_select.primary_min = 2, train_select.secondary_min = 8,
      pheno_prep.method = "raw_mean",
      geno_select.method = "focal_plus_onehop", geno_select.min_bridge = 1,
      kernel.method = "vanRaden_single", kernel.impute = "mean_round",
      model.method = "gblup_rrblup", model.lambda_select = "fixed", model.lambda_fixed = 1,
      predict_post.method = "direct_blup"),
    `10679` = .canary_filler(), # OHRWW (strong filler; best_single_project keeps it light)
    `10680` = .canary_filler(), # TCAP
    `10674` = .canary_filler(), # 24Crk  (weak)
    `10678` = .canary_filler(), # AWY1   (weak)
    `10681` = .canary_filler()  # STP1   (weak)
  )
}

# A light, robust config for the filler canaries: one genotyping project only, so
# it stays cheap even on trials with many projects (e.g. OHRWW's 39).
.canary_filler <- function() {
  .make_seed(list(
    train_select.method = "accession_overlap", train_select.primary_only = "yes",
    train_select.primary_min = 2, train_select.secondary_min = 8,
    pheno_prep.method = "blue_lm", pheno_prep.z_thr = 3,
    geno_select.method = "best_single_project",
    kernel.method = "vanRaden_single",
    model.method = "gblup_rrblup", model.lambda_select = "reml",
    predict_post.method = "direct_blup"))
}

# Coverage assertion: across canary_configs(), which subtask methods and key
# branch-levels are exercised, and which (if any) are not. Run offline; a gap is a
# documented oracle blind spot, not an error (those methods are still exercised by
# seed_configs() in the main loop).
canary_coverage <- function(cfgs = canary_configs()) {
  methods_hit <- function(st) unique(vapply(cfgs, function(c) as.character(c[[paste0(st, ".method")]]), character(1)))
  level_hit   <- function(key, val) any(vapply(cfgs, function(c) identical(as.character(c[[key]]), val), logical(1)))
  num_gt0     <- function(key) any(vapply(cfgs, function(c) { v <- c[[key]]; isTRUE(is.finite(suppressWarnings(as.numeric(v))) && as.numeric(v) > 0) }, logical(1)))

  method_rows <- purrr::map_dfr(names(SUBTASKS), function(st) {
    hit <- methods_hit(st); all_m <- SUBTASKS[[st]]$methods
    tibble::tibble(target = paste0(st, " methods"),
                   covered = paste(intersect(all_m, hit), collapse = ", "),
                   missing = paste(setdiff(all_m, hit), collapse = ", "))
  })
  levels <- tibble::tibble(
    target = c("primary_only=no", "primary_only=yes", "ge_weighting=env_gaussian",
               "standardize=per_trial_z", "include_E=yes", "include_E=no",
               "lambda_select=reml", "lambda_select=loo", "lambda_select=fixed", "impute=mean",
               "impute=mean_round", "blend_obs_w>0"),
    covered = c(level_hit("train_select.primary_only","no"), level_hit("train_select.primary_only","yes"),
                level_hit("pheno_prep.ge_weighting","env_gaussian"), level_hit("pheno_prep.standardize","per_trial_z"),
                level_hit("model.include_E","yes"), level_hit("model.include_E","no"),
                level_hit("model.lambda_select","reml"), level_hit("model.lambda_select","loo"),
                level_hit("model.lambda_select","fixed"),
                level_hit("kernel.impute","mean"), level_hit("kernel.impute","mean_round"),
                num_gt0("predict_post.blend_obs_w")))
  list(methods = method_rows, levels = levels)
}

# Run each canary trial under its OWN frozen coverage config (canary_configs()),
# across every CV scheme. Together the configs exercise all 19 methods + branch-
# levels, so a failure implicates whichever method that trial's config uses.
# Two severity tiers: a failure on a STRONG trial is a hard CANARY ALARM; a
# failure on a weak trial (settings$canary_weak_trials) is a soft warning only.
# Returns a tibble of (trial, scheme, status, reason, detail, n_test, score,
# method_signature). Any non-"ok" on a strong trial -- especially "suspect" --
# means investigate the code, not the trial.
check_canaries <- function(settings, conn, configs = canary_configs()) {
  if (is.null(configs) || !length(configs)) {
    message("no canary configs; nothing to check"); return(invisible(NULL))
  }
  weak <- as.character(settings$canary_weak_trials %||% character())
  sig  <- function(cfg) paste(vapply(names(SUBTASKS),
            function(st) as.character(cfg[[paste0(st, ".method")]]), character(1)), collapse = "/")
  out <- purrr::imap_dfr(configs, function(cfg, id) {
    trial <- tryCatch(build_trial_descriptor(id, conn, settings), error = function(e) NULL)
    if (is.null(trial)) {
      return(tibble::tibble(trial = id, scheme = NA, status = "no_descriptor",
                            reason = "could not build descriptor", detail = NA,
                            n_test = NA_integer_, score = NA_real_,
                            method_signature = sig(cfg), weak = id %in% weak))
    }
    purrr::map_dfr(settings$schemes, function(scheme) {
      ev <- evaluate_config_on_trial(cfg, trial, scheme, settings, conn)
      tibble::tibble(trial = id, scheme = scheme, status = ev$status,
                     reason = ev$reason %||% NA_character_,
                     detail = ev$detail %||% NA_character_,
                     n_test = ev$n_test %||% NA_integer_, score = ev$score,
                     method_signature = sig(cfg), weak = id %in% weak)
    })
  }, .progress = "Canary trials")
  bad <- out |> dplyr::filter(status != "ok")
  hard <- bad |> dplyr::filter(!weak)
  soft <- bad |> dplyr::filter(weak)
  if (nrow(hard)) {
    message("CANARY ALARM: ", nrow(hard), " STRONG (trial, scheme) failed -- likely a BUG: ",
            paste(sprintf("%s/%s=%s", hard$trial, hard$scheme, hard$status), collapse = ", "))
  }
  if (nrow(soft)) {
    message("canary warning (weak trials, may be genuinely marginal): ",
            paste(sprintf("%s/%s=%s", soft$trial, soft$scheme, soft$status), collapse = ", "))
  }
  if (!nrow(hard)) message("canaries OK: all STRONG trials predicted under their coverage configs.")
  out
}

# ===========================================================================
# Rich-trial config sweep (a CODE-correctness oracle)
# ===========================================================================
# check_canaries runs ONE frozen config per trial across nine trials -- a COVERAGE test that
# conflates a code bug with data-adequacy (a demanding config can legitimately fail on a
# data-poor trial). This does the complement: it varies ONE method/branch at a time from a
# robust baseline and runs each variant on a couple of DATA-RICH trials (default Big6 10675 +
# YT_Urb 10677). On trials rich enough to satisfy every method's data needs, feasibility is a
# property of the CODE, not the data -- so a branch that cannot produce a prediction on EITHER
# rich trial is a genuine bug signal. Known-degenerate cells are excluded (see .oracle_degenerate).
#
# Verdict per variant:
#   * an `error` (a crash) on any (trial, scheme) is ALWAYS a bug.
#   * otherwise the variant is HEALTHY if it produced `ok` on >=1 non-degenerate (trial, scheme);
#     if it never did (all infeasible/constant/too_few_overlap across both rich trials), it is
#     SUSPECT -- the branch could not predict even on rich data.

# The robust baseline every variant perturbs by exactly one field.
.oracle_baseline <- function()
  list(train_select.method = "accession_overlap", train_select.primary_only = "no",
       train_select.primary_min = 2, train_select.secondary_min = 8,
       pheno_prep.method = "raw_mean",
       geno_select.method = "all_projects",
       kernel.method = "vanRaden_single",
       model.method = "gblup_rrblup", model.lambda_select = "reml",
       predict_post.method = "cond_expectation", predict_post.min_overlap = 3)

# One config per method / behaviour-changing branch level (single-field perturbations of the
# baseline). Each is a full .make_seed()'d config; `label` says what it exercises.
.oracle_variants <- function() {
  b <- .oracle_baseline()
  v <- function(label, ...) list(label = label, cfg = .make_seed(modifyList(b, list(...))))
  list(
    v("baseline"),
    # -- train_select
    v("train_select=top_k_similar/genomic",       train_select.method = "top_k_similar", train_select.k = 15, train_select.similarity = "genomic"),
    v("train_select=top_k_similar/environmental", train_select.method = "top_k_similar", train_select.k = 15, train_select.similarity = "environmental"),
    v("train_select=top_k_similar/both",          train_select.method = "top_k_similar", train_select.k = 15, train_select.similarity = "both"),
    v("train_select=same_program",                train_select.method = "same_program", train_select.prog_cap = 20),
    v("train_select.primary_only=yes",            train_select.primary_only = "yes"),
    # -- pheno_prep
    v("pheno_prep=blue_lm",                       pheno_prep.method = "blue_lm", pheno_prep.z_thr = 3),
    v("pheno_prep=trial_center",                  pheno_prep.method = "trial_center", pheno_prep.z_thr = 3),
    v("pheno_prep=two_stage_blup",                pheno_prep.method = "two_stage_blup"),
    v("pheno_prep.ge_weighting=env_gaussian",     pheno_prep.ge_weighting = "env_gaussian", pheno_prep.ge_bandwidth = 1),
    v("pheno_prep.standardize=per_trial_z",       pheno_prep.standardize = "per_trial_z"),
    # -- geno_select
    v("geno_select=focal_plus_onehop",            geno_select.method = "focal_plus_onehop", geno_select.min_bridge = 1),
    v("geno_select.min_bridge=4",                 geno_select.method = "focal_plus_onehop", geno_select.min_bridge = 4),
    v("geno_select=best_single_project",          geno_select.method = "best_single_project"),
    # -- kernel
    v("kernel=em_combine",                        kernel.method = "em_combine"),
    v("kernel=rkhs_gaussian",                     kernel.method = "rkhs_gaussian", kernel.rkhs_theta = 1),
    v("kernel.impute=mean",                       kernel.impute = "mean"),
    # -- model
    v("model=gblup_sommer_GE+E",                  model.method = "gblup_sommer_GE", model.include_E = "yes"),
    v("model=gblup_rrblup[lambda=loo]",           model.method = "gblup_rrblup", model.lambda_select = "loo"),
    v("model=gblup_rrblup[lambda=fixed]",         model.method = "gblup_rrblup", model.lambda_select = "fixed", model.lambda_fixed = 1),
    v("model=rkhs-E",                             model.method = "rkhs", model.include_E = "no"),
    # -- predict_post
    v("predict_post=direct_blup",                 predict_post.method = "direct_blup"),
    v("predict_post.blend_obs_w=0.3",             predict_post.blend_obs_w = 0.3))
}

# Cells EXPECTED to be degenerate (not bugs). There are none.
#
# PITFALL: do not re-add `direct_blup` x CV00 here. It looks degenerate -- the focal lines are
# masked out of training -- but predict_test() predicts them through the kernel (LESSONS.md #22),
# so it is a real test. Excusing it would hide a regression in exactly that guard.
.oracle_degenerate <- function(cfg, scheme) FALSE

# Keep the variants whose label matches any `only` pattern (fixed substring, so "em_combine"
# matches "kernel=em_combine"). NULL/empty -> keep all. Errors if a filter matches nothing.
.select_variants <- function(variants, only = NULL) {
  if (is.null(only) || !length(only)) return(variants)
  keep <- vapply(variants, function(v)
    any(vapply(only, function(p) grepl(p, v$label, fixed = TRUE), logical(1))), logical(1))
  if (!any(keep)) stop("no sweep variant matches `only` (", paste(only, collapse = ", "),
                       "). Labels: ", paste(vapply(variants, `[[`, character(1), "label"), collapse = ", "))
  variants[keep]
}

# `only`: run just the variants whose label contains any of these strings (e.g.
# only = "em_combine", or c("kernel=", "model=")). Handy for re-checking one branch after a
# fix without re-running the whole sweep. See .oracle_variants() for the labels.
sweep_rich_trials <- function(settings, conn,
                              trials = settings$oracle_trials %||% c("10675", "10677"),
                              only = NULL) {
  trials <- as.character(trials)
  descs <- stats::setNames(lapply(trials, function(id)
    tryCatch(build_trial_descriptor(id, conn, settings), error = function(e) NULL)), trials)
  variants <- .select_variants(.oracle_variants(), only)

  rows <- purrr::map_dfr(variants, function(v) {
    purrr::map_dfr(trials, function(id) {
      tr <- descs[[id]]
      purrr::map_dfr(settings$schemes, function(sc) {
        if (is.null(tr))
          return(tibble::tibble(variant = v$label, trial = id, scheme = sc, status = "no_descriptor",
                                score = NA_real_, n_test = NA_integer_, reason = NA_character_,
                                degenerate = FALSE))
        ev <- evaluate_config_on_trial(v$cfg, tr, sc, settings, conn)
        tibble::tibble(variant = v$label, trial = id, scheme = sc, status = ev$status,
                       score = ev$score, n_test = ev$n_test %||% NA_integer_,
                       reason = ev$reason %||% NA_character_, degenerate = .oracle_degenerate(v$cfg, sc))
      })
    })
  }, .progress = "Rich-trial sweep")

  verdict <- rows |>
    dplyr::group_by(variant) |>
    dplyr::summarise(
      errored = any(status == "error"),
      healthy = any(status == "ok" & !degenerate),
      # accuracy = the Pearson r that score_predictions returns; best across the ok cells.
      best_score = suppressWarnings(max(score[status == "ok"], na.rm = TRUE)),
      .groups = "drop") |>
    dplyr::mutate(best_score = ifelse(is.finite(best_score), best_score, NA_real_),
                  verdict = dplyr::case_when(errored ~ "BUG(error)",
                                             healthy ~ "ok",
                                             TRUE    ~ "SUSPECT"))
  bad <- dplyr::filter(verdict, verdict != "ok")

  disp <- rows; disp$score <- round(disp$score, 3)
  print(as.data.frame(disp[, c("variant", "trial", "scheme", "status", "score", "reason")]), row.names = FALSE)
  message("\naccuracy (best Pearson r on the ok cells) per variant:")
  vd <- verdict; vd$best_score <- round(vd$best_score, 3)
  print(as.data.frame(vd[, c("variant", "verdict", "best_score")]), row.names = FALSE)
  if (nrow(bad)) {
    message("\nSWEEP ALARM: ", nrow(bad), " variant(s) never predicted on either rich trial ",
            "(a code bug, since these trials are data-rich):")
    for (i in seq_len(nrow(bad))) message(sprintf("  %-40s %s", bad$variant[i], bad$verdict[i]))
    message("Cross-check a SUSPECT with diagnose_trial() on the same trial + config.")
  } else {
    message("\nSWEEP CLEAN: every method/branch produced a prediction on a data-rich trial.")
  }
  invisible(list(rows = rows, verdict = verdict))
}

# Replay ONE trial and print the data funnel stage by stage, alongside an
# independent re-derivation of the raw counts straight from T3. The pipeline's
# numbers and the independent numbers should agree; where they diverge -- or
# where a stage cliffs to ~zero despite plentiful input -- is the bug.
# The configs a trial actually FAILED under, recovered from the store, ready to hand back to
# diagnose_trial(). A failure is a property of (trial x config), not of the trial alone: the
# default canary_config() is deliberately permissive and uses geno_select = "all_projects", so
# diagnosing with it does not reproduce (say) a focal_plus_onehop failure and can return a false
# clean bill of health. Recover the real thing:
#
#   con <- open_store(settings$db_path)
#   bad <- failed_configs(con, "10260")
#   diagnose_trial("10260", settings, conn, cfg = bad$cfg[[1]], scheme = bad$scheme[1])
#
# `status = NULL` returns every stored eval for the trial (useful to contrast a config that
# failed against one that succeeded on the same trial).
failed_configs <- function(con, trial_id,
                           status = c("suspect", "infeasible", "error", "too_few_overlap")) {
  e <- read_evals(con)
  e <- e[as.character(e$trial_id) == as.character(trial_id), , drop = FALSE]
  if (!is.null(status)) e <- e[e$status %in% status, , drop = FALSE]
  if (!nrow(e)) { message("no stored evals for trial ", trial_id,
                          if (!is.null(status)) paste0(" with status in {", paste(status, collapse = ", "), "}")
                          else ""); return(invisible(NULL)) }
  tibble::tibble(id = e$id, status = e$status, reason = e$reason, scheme = e$scheme,
                 detail = e$detail, cfg = lapply(e$config_json, config_from_json))
}

diagnose_trial <- function(study_id, settings, conn,
                           cfg = canary_config(), scheme = settings$schemes[1]) {
  id <- as.character(study_id)
  cat("=== diagnose trial ", id, " (scheme ", scheme, ") ===\n", sep = "")
  # Print the config being diagnosed. Which one it is decides what the output means, and the
  # default is NOT the one that failed -- see failed_configs() above.
  cat("  config under test", if (identical(cfg, canary_config()))
        " (DEFAULT canary_config -- permissive, may NOT reproduce a stored failure)" else "", ":\n", sep = "")
  cat(format_config(cfg), sep = "\n")

  # --- independent ground-truth probes (a DIFFERENT path than run_pipeline) ---
  acc <- tryCatch(get_trial_accessions(id, conn, settings), error = function(e) character())
  obs <- tryCatch(get_observations(id, conn, settings), error = function(e) NULL)
  cat(sprintf("  trial accessions (wizard):      %d\n", length(acc)))
  cat(sprintf("  focal-trait observations:       %d rows over %d germplasm\n",
              if (is.null(obs)) 0L else nrow(obs),
              if (is.null(obs)) 0L else dplyr::n_distinct(obs$germplasm_name)))
  projs <- tryCatch(projects_for_accessions(acc, conn, settings), error = function(e) character())
  cat(sprintf("  genotyping projects covering:   %d  (%s)\n",
              length(projs), paste(utils::head(projs, 8), collapse = ", ")))
  # The decisive check for the synonym/name-mismatch class of bug: do the dosage
  # matrix rownames actually intersect the trial's accession names?
  #
  # PITFALL this check was blind to until 2026-07-30: it passed `keep_samples = acc`, and
  # get_project_dosage's subset_samples returns NULL when NOTHING matches. So the cliff it exists
  # to catch was swallowed and printed as "download/parse failed", and on the non-NULL branch
  # `ov == nrow(d) > 0` by construction, so the CLIFF line could not fire either way. Load the
  # project's WHOLE population and intersect here instead.
  n_cover <- 0L; best_ov <- 0L
  if (length(projs)) {
    for (pid in projs) {
      d <- tryCatch(get_project_dosage(pid, NULL, conn, settings), error = function(e) NULL)
      if (is.null(d)) { cat(sprintf("    project %-7s dosage: NULL (download/parse failed)\n", pid)); next }
      ov <- length(intersect(rownames(d), acc))
      if (ov > 0) n_cover <- n_cover + 1L
      best_ov <- max(best_ov, ov)
      cat(sprintf("    project %-7s dosage: %d samples x %d markers, overlap with accessions = %d%s\n",
                  pid, nrow(d), ncol(d), ov,
                  if (nrow(d) > 0 && ov == 0) "   <-- no name overlap in THIS project" else ""))
    }
    cat(sprintf("  projects containing >=1 focal accession: %d of %d (best single project: %d)%s\n",
                n_cover, length(projs), best_ov,
                if (n_cover == 0 && length(acc) > 0)
                  "   <-- CLIFF: the focal lines are genotyped NOWHERE under these names" else ""))
  }

  # What subtask C hands to the kernel, and what .best_panel picks out of it. This separates the
  # two causes of test_in = 0: no panel covers the focal lines (a name/data problem) versus a
  # covering panel existing but being passed over (.best_panel maximizes coverage of
  # union(train, focal), which large training sets dominate).
  # `trial` is built once here and reused by the pipeline replay below.
  trial <- tryCatch(build_trial_descriptor(id, conn, settings), error = function(e) NULL)
  if (!is.null(trial)) {
    # Replay run_pipeline's own route to train_acc (select trials -> observations -> CV mask ->
    # targets), so the panels reported are the ones the pipeline would actually see.
    panels <- tryCatch({
      focal_acc <- trial$accessions
      train_ids <- setdiff(select_training_trials(cfg, trial, conn, settings), id)
      train_obs <- mask_cv(get_observations(train_ids, conn, settings), focal_acc, scheme)
      train_acc <- names(build_targets(cfg, train_obs, trial, conn, settings))
      list(p = choose_geno_sources(cfg, train_acc, focal_acc, conn, settings),
           train_acc = train_acc, focal_acc = focal_acc)
    }, error = function(e) { cat("  (panel probe unavailable: ", conditionMessage(e), ")\n", sep = ""); NULL })

    if (!is.null(panels) && length(panels$p)) {
      pl <- panels$p; f_acc <- panels$focal_acc; t_acc <- panels$train_acc
      cat(sprintf("  panels from subtask C (%d), focal / train overlap:\n", length(pl)))
      fo <- vapply(pl, function(d) length(intersect(rownames(d), f_acc)), integer(1))
      to <- vapply(pl, function(d) length(intersect(rownames(d), t_acc)), integer(1))
      # Exactly .best_panel's criterion: coverage of union(train, focal), counted once.
      need <- union(t_acc, f_acc)
      cov  <- vapply(pl, function(d) length(intersect(rownames(d), need)), integer(1))
      for (i in seq_along(pl))
        cat(sprintf("    panel %-20s %6d rows   focal %5d   train %6d   union %6d\n",
                    names(pl)[i] %||% i, nrow(pl[[i]]), fo[i], to[i], cov[i]))
      picked <- which.max(cov)
      cat(sprintf("  .best_panel picks %s -- focal coverage %d%s\n",
                  names(pl)[picked] %||% picked, fo[picked],
                  if (fo[picked] == 0 && any(fo > 0))
                    "   <-- PASSED OVER a focal-covering panel: .best_panel, not the data"
                  else if (fo[picked] == 0)
                    "   <-- no panel covers the focal lines: a data/name problem" else ""))
    }
  }

  # --- run the pipeline and report where it lands -----------------------------
  # `trial` was built above for the panel probe.
  if (is.null(trial)) { cat("  could not build descriptor; stopping.\n"); return(invisible(NULL)) }
  res <- tryCatch(
    { po <- run_pipeline(cfg, trial, scheme, settings, conn)
      sc <- score_predictions(po$pred, po$obs)
      cat(sprintf("  PIPELINE OK: %d predictions, score=%s (status %s)\n",
                  length(po$pred), signif(sc$score, 3), sc$status)); sc },
    error = function(e) {
      if (inherits(e, "optimizer_infeasible")) {
        cat(sprintf("  PIPELINE INFEASIBLE [%s%s]: %s\n", e$code,
                    if (isTRUE(e$suspect)) " / SUSPECT-BUG" else "",
                    conditionMessage(e)))
        if (!is.null(e$funnel))
          cat("    funnel: ", funnel_string(e$funnel), "\n", sep = "")
      } else if (inherits(e, "optimizer_fatal")) {
        cat(sprintf("  PIPELINE FATAL [%s]: %s\n", e$code, conditionMessage(e)))
      } else {
        cat("  PIPELINE ERROR: ", conditionMessage(e), "\n", sep = "")
      }
      invisible(NULL)
    })
  invisible(res)
}

# ===========================================================================
# STAGE-1 calibration: validate the data path against an INDEPENDENT anchor
# BEFORE freezing any canary configs (see plan/cuddly-coalescing-teacup.md).
#
# The catch-22: calibration runs the very functions the oracle guards, so it
# cannot self-validate them. We break it with an external ground truth that
# touches NONE of our code -- the five Predictathon teams' submission files --
# and compare our pipeline's per-trial counts against it. A data-hiding bug
# shows up as our count diverging from the anchor. Freeze only once they agree;
# this is iterative (fix divergence, re-run) and human-reviewed.
# ===========================================================================

# Independent anchor: median CV0 predicted-accession count per focal trial across
# the five anonymized algorithm submissions (scripts/PredictionN/submission/).
# Pure file reads -- no BrAPI/pipeline code. `settings$canary_trials` must be a
# named vector studyName -> studyDbId (submission subfolders are named by studyName).
canary_anchor <- function(settings, preds_root = NULL) {
  # optimizer is scripts/Analysis_Claude/optimizer -> here::here("..","..") = scripts/
  if (is.null(preds_root)) preds_root <- here::here("..", "..")
  name2id <- settings$canary_trials
  if (is.null(names(name2id))) {
    warning("settings$canary_trials must be named studyName->id for the anchor")
    return(tibble::tibble())
  }
  algos <- paste0("Prediction", 1:5)
  count_pred <- function(alg, study_name) {
    f <- file.path(preds_root, alg, "submission", study_name, "CV0_Predictions.csv")
    if (!file.exists(f)) return(NA_integer_)
    tryCatch(as.integer(nrow(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))),
             error = function(e) NA_integer_)
  }
  purrr::imap_dfr(name2id, function(id, study_name) {
    counts <- vapply(algos, count_pred, integer(1), study_name = study_name)
    tibble::tibble(
      study_db_id   = as.character(id),
      study_name    = study_name,
      anchor_n_pred = suppressWarnings(as.integer(round(stats::median(counts, na.rm = TRUE)))),
      anchor_n_part = sum(is.finite(counts)))
  }, .progress = "Read submission anchors")
}

# Probe each canary trial's data shape THROUGH the real pipeline functions and
# join to the anchor. The `divergent` flag (our genotyped-focal count vs the
# anchor's predicted-accession count outside 0.5x-2x) is the validation signal:
# investigate those before trusting calibration or freezing configs. Not called
# at startup; run by hand during the calibrate-then-freeze bootstrap.
calibrate_canary_trials <- function(settings, conn, anchor = NULL, deep = FALSE) {
  if (is.null(anchor)) anchor <- canary_anchor(settings)
  ids     <- as.character(settings$canary_trials)
  id2name <- stats::setNames(names(settings$canary_trials), ids)
  weak    <- as.character(settings$canary_weak_trials %||% character())

  has_sommer <- requireNamespace("sommer", quietly = TRUE)
  has_lme4   <- requireNamespace("lme4",   quietly = TRUE)
  message(sprintf("packages: sommer=%s (sommer_GE/rkhs), lme4=%s (two_stage_blup)",
                  has_sommer, has_lme4))

  rows <- purrr::map_dfr(ids, function(id) {
    trial <- tryCatch(build_trial_descriptor(id, conn, settings), error = function(e) NULL)
    if (is.null(trial)) {
      return(tibble::tibble(study_db_id = id, study_name = id2name[id],
                            error = "no_descriptor"))
    }
    acc <- trial$accessions
    nt  <- function(block) length(setdiff(
      tryCatch(select_training_trials(block, trial, conn, settings),
               error = function(e) character()), id))

    # train_select feasibility via the REAL code paths (top_k uses the candidate
    # pool to avoid its heavy per-candidate similarity scan).
    n_ao_primary <- nt(list(train_select.method = "accession_overlap",
                            train_select.primary_only = "yes",
                            train_select.primary_min = 2, train_select.secondary_min = 8))
    n_same_prog  <- nt(list(train_select.method = "same_program",
                            train_select.prog_cap = 9999))
    n_topk_pool  <- length(tryCatch(.find_related(id, conn, settings, 1),
                                    error = function(e) character()))

    # genotyping projects covering the focal accessions.
    projects <- tryCatch(projects_for_accessions(acc, conn, settings),
                         error = function(e) character())
    # CHEAP path: project membership via the breeder wizard -- canonical accession
    # names, exercises project selection WITHOUT downloading any VCF. This is the
    # default our_n_geno_focal we compare to the anchor.
    members <- if (length(projects)) unique(unlist(lapply(projects, function(pid) {
        w <- tryCatch(.brapi_try(function() conn$wizard("accessions", list(genotyping_projects = pid)),
                                 conn = conn, settings = settings, what = "project members wizard"),
                      error = function(e) NULL)
        if (is.null(w)) character() else as.character(w$data$names)
      }), use.names = FALSE)) else character()
    geno_wizard <- length(intersect(acc, members))
    # DEEP path (opt-in): the actual dosage extraction -- downloads VCFs and matches
    # on VCF sample names. A big gap below geno_wizard localizes a VCF name-matching
    # bug (e.g. synonym sample names) rather than project selection.
    geno_dosage <- NA_integer_
    if (deep) {
      dl <- tryCatch(choose_geno_sources(
        list(geno_select.method = "all_projects"),
        acc, acc, conn, settings), error = function(e) list())
      geno_dosage <- length(intersect(acc, unique(unlist(lapply(dl, rownames),
                                                          use.names = FALSE))))
    }
    tibble::tibble(
      study_db_id       = id,
      study_name        = id2name[id],
      weak              = id %in% weak,
      n_acc             = length(acc),
      n_ao_primary      = n_ao_primary,
      n_same_prog       = n_same_prog,
      n_topk_pool       = n_topk_pool,
      n_projects        = length(projects),
      our_n_geno_focal  = geno_wizard,
      our_n_geno_dosage = geno_dosage,
      has_coords        = is.finite(trial$lat))
  }, .progress = "Calibrate canaries")

  out <- dplyr::left_join(rows, anchor, by = c("study_db_id", "study_name")) |>
    dplyr::mutate(
      ratio = ifelse(is.finite(.data$our_n_geno_focal) &
                       is.finite(.data$anchor_n_pred) & .data$anchor_n_pred > 0,
                     .data$our_n_geno_focal / .data$anchor_n_pred, NA_real_),
      divergent = is.finite(.data$ratio) & (.data$ratio < 0.5 | .data$ratio > 2))
  attr(out, "has_sommer") <- has_sommer
  attr(out, "has_lme4")   <- has_lme4
  out
}

# Pretty-print the probe-vs-anchor table (the Stage-1 deliverable).
print_calibration <- function(cal) {
  cols <- intersect(c("study_name", "weak", "n_acc", "n_ao_primary", "n_same_prog",
                      "n_topk_pool", "n_projects", "our_n_geno_focal",
                      "our_n_geno_dosage", "anchor_n_pred", "ratio", "divergent"), names(cal))
  print(as.data.frame(cal[, cols]), row.names = FALSE)
  # If the deep VCF path was run, flag where it sees far fewer genotyped focal
  # accessions than the cheap wizard path -- the VCF name-matching (synonym) signal.
  if ("our_n_geno_dosage" %in% names(cal) && any(is.finite(cal$our_n_geno_dosage))) {
    gap <- cal[is.finite(cal$our_n_geno_dosage) & is.finite(cal$our_n_geno_focal) &
                 cal$our_n_geno_focal > 0 &
                 cal$our_n_geno_dosage < 0.5 * cal$our_n_geno_focal, , drop = FALSE]
    if (nrow(gap)) {
      message("\nVCF dosage path sees far fewer genotyped accessions than project ",
              "membership -- likely a VCF sample-name (synonym) mismatch:")
      for (i in seq_len(nrow(gap)))
        message(sprintf("  %-28s wizard=%s dosage=%s",
                        gap$study_name[i], gap$our_n_geno_focal[i], gap$our_n_geno_dosage[i]))
    }
  }
  d <- cal[isTRUE(cal$divergent) | (!is.na(cal$divergent) & cal$divergent), , drop = FALSE]
  if (nrow(d)) {
    message("\nDIVERGENT (our count vs anchor off >2x) -- investigate before freezing:")
    for (i in seq_len(nrow(d)))
      message(sprintf("  %-28s our=%s anchor=%s ratio=%.2f",
                      d$study_name[i], d$our_n_geno_focal[i], d$anchor_n_pred[i], d$ratio[i]))
  } else {
    message("\nNo divergences flagged (subject to anchor availability).")
  }
  invisible(cal)
}
