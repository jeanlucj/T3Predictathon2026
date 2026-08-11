# diagnostics.R
#
# Tools to answer: "is this trial REALLY infeasible, or is a bug hiding the data?" The loop
# continues past infeasible pairs, so a name mismatch or parse error can masquerade as a trial
# that cannot be predicted. Three defences:
#
#   1. The failure log flags funnel CLIFFS as "suspect" rather than "infeasible" (pipeline.R).
#   2. check_canaries() runs a permissive config on KNOWN-feasible trials -- an oracle, since a
#      failure there can only be a bug.
#   3. diagnose_trial() replays one trial with the funnel printed beside an independent
#      re-derivation of the raw counts from T3.

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
# FROZEN coverage canary configs, one per focal trial, assigned so that across the canaries
# every subtask method and every behaviour-changing parameter level is exercised at least once
# -- a bug in ANY branch then trips the oracle.
#
# Coverage rests on four data-rich trials (Aurora 10673, Big6 10675, CornellMaster 10676,
# YT_Urb 10677), which between them carry every demanding branch; the rest get a light filler
# that only confirms "still predictable". Keyed by studyDbId.
#
# EVALUATION.md sec. 9 tabulates which trial exercises which branch; canary_coverage() is the
# machine-checkable version and is what to trust.
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

# Which subtask methods and branch-levels canary_configs() exercises, and which it misses. Run
# offline; a gap is a documented oracle blind spot, not an error -- seed_configs() still
# exercises those methods in the main loop.
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

# Run each canary trial under its OWN frozen coverage config, across every CV scheme, so a
# failure implicates whichever method that trial's config uses. Two severity tiers: a failure on
# a STRONG trial is a hard CANARY ALARM; one on a weak trial (settings$canary_weak_trials) is a
# soft warning. Returns a tibble of (trial, scheme, status, reason, detail, n_test, score,
# method_signature); any non-"ok" on a strong trial means investigate the code, not the trial.
check_canaries <- function(settings, conn, configs = canary_configs()) {
  if (is.null(configs) || !length(configs)) {
    message("no canary configs; nothing to check"); return(invisible(NULL))
  }
  weak <- as.character(settings$canary_weak_trials %||% character())
  sig  <- function(cfg) paste(vapply(names(SUBTASKS),
            function(st) as.character(cfg[[paste0(st, ".method")]]), character(1)), collapse = "/")
  out <- purrr::imap_dfr(configs, function(cfg, id) {
    derr  <- NULL
    trial <- tryCatch(build_trial_descriptor(id, conn, settings),
                      error = function(e) { derr <<- conditionMessage(e); NULL })
    if (is.null(trial)) {
      return(tibble::tibble(trial = id, scheme = NA, status = "no_descriptor",
                            reason = "could not build descriptor", detail = derr,
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
# The complement to check_canaries, which conflates a code bug with data adequacy. This varies
# ONE method or branch at a time from a robust baseline and runs each variant on DATA-RICH
# trials, where feasibility is a property of the CODE rather than the data -- so a branch that
# cannot predict on EITHER is a genuine bug signal. Known-degenerate cells are excluded.
#
# Verdict per variant: an `error` on any (trial, scheme) is always a bug; otherwise HEALTHY if
# it produced `ok` on >= 1 non-degenerate cell, and SUSPECT if it never did.

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
# PITFALL: `direct_blup` x CV00 looks degenerate but is not -- predict_test() predicts the masked
# focal lines through the kernel (LESSONS #22). Excusing it would hide a regression in that guard.
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

# `only`: run just the variants whose label contains any of these strings (e.g. "em_combine", or
# c("kernel=", "model=")) -- re-check one branch without the whole sweep. Labels: .oracle_variants().
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

# The configs a trial actually FAILED under, recovered from the store for diagnose_trial(). A
# failure is a property of (trial x config): the permissive default canary_config() will not
# reproduce, say, a focal_plus_onehop failure, and returns a false clean bill of health.
# Recover the real thing:
#
#   con <- open_store(settings$db_path)
#   bad <- failed_configs(con, "10260")
#   diagnose_trial("10260", settings, conn, cfg = bad$cfg[[1]], scheme = bad$scheme[1])
#
# `status = NULL` returns every stored eval for the trial (useful to contrast a config that
# failed against one that succeeded on the same trial).
failed_configs <- function(con, trial_id,
                           status = c("suspect", "infeasible", "error", "too_few_overlap",
                                      "non_finite")) {
  e <- read_evals(con)
  e <- e[as.character(e$trial_id) == as.character(trial_id), , drop = FALSE]
  if (!is.null(status)) e <- e[e$status %in% status, , drop = FALSE]
  if (!nrow(e)) { message("no stored evals for trial ", trial_id,
                          if (!is.null(status)) paste0(" with status in {", paste(status, collapse = ", "), "}")
                          else ""); return(invisible(NULL)) }
  tibble::tibble(id = e$id, status = e$status, reason = e$reason, scheme = e$scheme,
                 detail = e$detail, cfg = lapply(e$config_json, config_from_json))
}

# Replay ONE trial, printing the data funnel stage by stage beside an independent re-derivation
# of the raw counts from T3. Where the two diverge -- or where a stage cliffs to ~zero despite
# plentiful input -- is the bug.
#
# Two parts are expensive and optional, since the panel table is what discriminates a name
# problem from a panel-selection problem:
#   max_projects  caps the per-project attribution loop, which loads each project's whole dosage
#                 matrix (up to ~10 GB apiece). It only attributes a cliff to individual
#                 projects; the panel probe below already says whether ANY panel covers the focal
#                 lines.
#   replay        FALSE skips the final run_pipeline(), a complete evaluation costing hours on a
#                 big trial, which only confirms an outcome the store already recorded.
diagnose_trial <- function(study_id, settings, conn,
                           cfg = canary_config(), scheme = settings$schemes[1],
                           max_projects = 8L, replay = TRUE) {
  id <- as.character(study_id)
  cat("=== diagnose trial ", id, " (scheme ", scheme, ") ===\n", sep = "")
  # Which config this is decides what the output means, and the default is NOT the one that
  # failed -- see failed_configs() above.
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
  # The decisive check for the synonym/name-mismatch class of bug: do the dosage matrix rownames
  # intersect the trial's accession names? Load each project's WHOLE population and intersect
  # here -- passing keep_samples = acc returns NULL when nothing matches, swallowing the very
  # cliff this exists to catch.
  n_cover <- 0L; best_ov <- 0L
  probe <- utils::head(projs, max(0L, as.integer(max_projects)))
  if (length(probe) < length(projs))
    cat(sprintf("  (probing %d of %d projects individually -- max_projects; the panel table below\n",
                length(probe), length(projs)),
        "   is the authoritative coverage answer)\n", sep = "")
  if (length(probe)) {
    for (pid in probe) {
      d <- tryCatch(get_project_dosage(pid, NULL, conn, settings), error = function(e) NULL)
      if (is.null(d)) { cat(sprintf("    project %-7s dosage: NULL (download/parse failed)\n", pid)); next }
      ov <- length(intersect(rownames(d), acc))
      if (ov > 0) n_cover <- n_cover + 1L
      best_ov <- max(best_ov, ov)
      cat(sprintf("    project %-7s dosage: %d samples x %d markers, overlap with accessions = %d%s\n",
                  pid, nrow(d), ncol(d), ov,
                  if (nrow(d) > 0 && ov == 0) "   <-- no name overlap in THIS project" else ""))
    }
    # Denominator is what was PROBED: under max_projects this is a sample, so it can suggest a
    # cliff but not establish one. The panel table does that.
    partial <- length(probe) < length(projs)
    cat(sprintf("  projects containing >=1 focal accession: %d of %d probed (best single: %d)%s\n",
                n_cover, length(probe), best_ov,
                if (n_cover == 0 && length(acc) > 0)
                  if (partial) "   <-- no overlap in the PROBED projects; check the panel table below"
                  else "   <-- CLIFF: the focal lines are genotyped NOWHERE under these names"
                else ""))
  }

  # What subtask C hands to the kernel and what .best_panel picks out of it, separating the two
  # causes of test_in = 0: no panel covers the focal lines (a name/data problem), or a covering
  # panel exists but is passed over (.best_panel maximizes coverage of union(train, focal), which
  # large training sets dominate). `trial` is reused by the pipeline replay below.
  trial_err <- NULL
  trial <- tryCatch(build_trial_descriptor(id, conn, settings),
                    error = function(e) { trial_err <<- conditionMessage(e); NULL })
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
  # A descriptor failure is usually the catalogue fetch (auth / server), not the trial -- and the
  # counts printed above come from cache, so they look healthy regardless.
  if (is.null(trial)) {
    cat("  could not build descriptor; stopping. reason: ", trial_err %||% "unknown", "\n", sep = "")
    return(invisible(NULL))
  }
  if (!isTRUE(replay)) {
    cat("  (pipeline replay skipped -- replay = FALSE; the stored eval already records the outcome)\n")
    return(invisible(NULL))
  }
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
# STAGE-1 calibration: validate the data path against an INDEPENDENT anchor before freezing any
# canary configs.
#
# Calibration runs the very functions the oracle guards, so it cannot self-validate them. The way
# out is a ground truth that touches NONE of our code -- the five submission files -- against
# which our per-trial counts are compared. A data-hiding bug shows up as a divergence. Freeze
# only once they agree; iterative and human-reviewed.
# ===========================================================================

# Median CV0 predicted-accession count per focal trial across the five submissions
# (scripts/PredictionN/submission/). Pure file reads -- no BrAPI or pipeline code.
# `settings$canary_trials` must be named studyName -> studyDbId, since submission subfolders are
# named by studyName.
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

# Probe each canary trial's data shape THROUGH the real pipeline functions and join to the
# anchor. `divergent` -- our genotyped-focal count outside 0.5x-2x of the anchor's
# predicted-accession count -- is the validation signal; investigate those before freezing
# configs. Run by hand, not at startup.
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

    # train_select feasibility via the REAL code paths (top_k uses the candidate pool, to avoid
    # its heavy per-candidate similarity scan).
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
    # CHEAP path: project membership via the breeder wizard -- canonical accession names,
    # exercising project selection without downloading any VCF. This is the our_n_geno_focal
    # compared to the anchor.
    members <- if (length(projects)) unique(unlist(lapply(projects, function(pid) {
        w <- tryCatch(.brapi_try(function() conn$wizard("accessions", list(genotyping_projects = pid)),
                                 conn = conn, settings = settings, what = "project members wizard"),
                      error = function(e) NULL)
        if (is.null(w)) character() else as.character(w$data$names)
      }), use.names = FALSE)) else character()
    geno_wizard <- length(intersect(acc, members))
    # DEEP path (opt-in): the real dosage extraction, downloading VCFs and matching on VCF
    # sample names. A big gap below geno_wizard localizes a VCF name-matching bug (e.g. synonym
    # sample names) rather than project selection.
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
  # Where the deep VCF path sees far fewer genotyped focal accessions than the cheap wizard
  # path -- the VCF name-matching (synonym) signal.
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
