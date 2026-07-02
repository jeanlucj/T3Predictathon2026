# pipeline.R
#
# The parameterized genomic-prediction pipeline. run_pipeline() executes the six
# subtasks for one configuration on one focal trial under one CV scheme and
# returns predicted and observed per-accession yields for scoring. Each subtask
# is a dispatcher that branches on the method named in the configuration -- this
# is where "recombine subtask methods from different algorithms" actually
# happens. To add a literature method, add a branch here and an entry in
# R/config_space.R.
#
# This is the heavy, network- and compute-bound path. It is written defensively:
# any subtask that cannot be completed for a given trial/config throws, and
# evaluate.R records the failure so the optimizer learns to avoid it. The genomic
# core (GBLUP via rrBLUP) is the reliable backbone; the G+E / RKHS / EM-combine
# variants extend it.
#
# VALIDATION: this path needs live T3 access and cannot be exercised offline, so
# run it once on a single known trial (see README "First real run") and confirm
# the predicted/observed join is non-empty before launching the background loop.

library(tidyverse)

run_pipeline <- function(cfg, trial, scheme, settings, conn) {
  focal_id  <- trial$id
  focal_acc <- trial$accessions

  # --- observed truth: per-accession yield BLUE in the focal trial ---------
  focal_obs <- get_observations(focal_id, conn, settings)
  if (!nrow(focal_obs)) infeasible("no_focal_obs")
  obs <- .per_acc_blue(focal_obs)                       # named numeric over focal germ

  # --- A. select training trials -------------------------------------------
  train_ids <- select_training_trials(cfg, trial, conn, settings)
  train_ids <- setdiff(train_ids, focal_id)
  if (length(train_ids) < settings$min_train_trials)
    infeasible("too_few_train_trials", sprintf("%d < %d", length(train_ids), settings$min_train_trials))

  # --- B. preprocess phenotypes -> per-accession targets (with CV masking) -
  train_obs <- get_observations(train_ids, conn, settings)
  train_obs <- mask_cv(train_obs, focal_acc, scheme)
  if (!nrow(train_obs)) infeasible("no_train_obs_after_mask")
  targets <- build_targets(cfg, train_obs, trial, conn, settings)   # named numeric
  if (length(targets) < 10) infeasible("too_few_train_acc", length(targets))

  # --- C. select genotyping data -------------------------------------------
  train_acc <- names(targets)
  dosage_list <- choose_geno_sources(cfg, train_acc, focal_acc, conn, settings)
  n_projects  <- attr(dosage_list, "n_projects") %||% NA_integer_
  if (!length(dosage_list)) {
    # Funnel cliff: projects covering these accessions WERE found, yet not one
    # dosage matrix came back. That points at a download/parse/name bug, not at
    # the trial being genuinely ungenotyped -> flag for review, don't bury it.
    fn <- c(focal_acc = length(focal_acc), train_acc = length(train_acc),
            geno_projects = n_projects, geno_dosage = 0L)
    infeasible("no_genotype_data", sprintf("projects=%s, dosage=0", n_projects),
               funnel = fn, suspect = isTRUE(n_projects > 0))
  }

  # --- D. build relationship / kernel --------------------------------------
  K <- build_kernel(cfg, dosage_list)
  geno_acc <- rownames(K)
  train_in <- intersect(train_acc, geno_acc)
  test_in  <- intersect(focal_acc, geno_acc)
  if (length(train_in) < settings$min_train_trials || length(test_in) < 5) {
    fn <- c(focal_acc = length(focal_acc), train_acc = length(train_acc),
            geno_projects = n_projects, geno_acc = length(geno_acc),
            train_in = length(train_in), test_in = length(test_in))
    # Funnel cliff: both the phenotyped side AND the genotyped side are sizeable,
    # yet they barely intersect. Real lines usually overlap; a near-empty
    # intersection here is the classic signature of a name/key mismatch (e.g. VCF
    # samples under synonym names) rather than a genuinely unpredictable trial.
    K_ <- settings$min_trial_acc
    suspect <- length(train_acc) >= K_ && length(geno_acc) >= K_ &&
               (length(train_in) == 0 || length(test_in) == 0)
    infeasible("insufficient_geno_overlap",
               sprintf("train %d, test %d (of train_acc %d, geno_acc %d)",
                       length(train_in), length(test_in),
                       length(train_acc), length(geno_acc)),
               funnel = fn, suspect = suspect)
  }

  # --- E. train model ------------------------------------------------------
  fit <- train_model(cfg, targets[train_in], K, train_in, test_in, trial, train_obs)

  # --- F. predict ----------------------------------------------------------
  pred <- predict_test(cfg, fit, K, train_in, test_in, targets, scheme, settings)

  list(pred = pred, obs = obs)
}

# CV masking: CV0 keeps every training observation; CV00 additionally removes any
# observation on a focal-trial accession (Jarquin et al. 2017). Factored out of
# run_pipeline so the CV0/CV00 distinction -- the crux of the challenge -- is
# unit-testable in isolation.
mask_cv <- function(train_obs, focal_acc, scheme) {
  if (scheme == "CV00") dplyr::filter(train_obs, !(germplasm_name %in% focal_acc)) else train_obs
}

# ===========================================================================
# Subtask A: select training trials
# ===========================================================================
select_training_trials <- function(cfg, trial, conn, settings) {
  m <- cfg$train_select.method
  focal_id <- trial$id

  if (m == "accession_overlap") {
    # Primary tier: trials sharing >= primary_min germplasm directly with the focal trial.
    primary <- .find_related(focal_id, conn, cfg$train_select.primary_min)
    if (identical(cfg$train_select.primary_only, "yes")) return(primary)
    # Secondary tier: trials sharing germplasm with the primary pool at a higher
    # threshold (one expansion hop from a few primary trials) -- friends-of-friends.
    secondary <- character()
    for (sid in utils::head(primary, 10)) {
      secondary <- union(secondary, .find_related(sid, conn, cfg$train_select.secondary_min))
    }
    return(union(primary, secondary))
  }
  if (m == "top_k_similar") {
    cand <- .find_related(focal_id, conn, 1)
    if (!length(cand)) return(cand)
    sim <- .trial_similarity(trial, cand, conn, settings, cfg$train_select.similarity)
    return(names(utils::head(sort(sim, decreasing = TRUE), cfg$train_select.k)))
  }
  if (m == "same_program") {
    cat <- trial_catalog(conn, settings)
    same <- cat |> dplyr::filter(program_name == trial$program, study_db_id != focal_id)
    ids <- as.character(same$study_db_id)
    if (length(ids) > cfg$train_select.prog_cap) {
      same2 <- same |> dplyr::filter(location_name == trial$location | year == trial$year)
      ids <- as.character(same2$study_db_id)
      ids <- utils::head(ids, cfg$train_select.prog_cap)
    }
    return(ids)
  }
  fatal(paste("unknown train_select method", m), "bug_unknown_method")
}

.find_related <- function(study_id, conn, min_common) {
  res <- tryCatch(
    T3BrapiHelpers::find_other_studies_evaluating_same_germplasm(
      study_id, conn, min_germ_common = max(1L, as.integer(min_common))),
    error = function(e) NULL)
  if (is.null(res)) return(character())
  # find_other_studies_evaluating_same_germplasm returns a janitor tabyl with
  # columns other_study_db_id, n.
  ids <- if (is.data.frame(res)) as.character(res$other_study_db_id %||% res[[1]]) else as.character(res)
  unique(ids)
}

# Similarity of candidate trials to the focal trial: genomic = number of shared
# accessions; environmental = Gaussian on scaled lat/long/year distance.
.trial_similarity <- function(trial, cand_ids, conn, settings, kind) {
  g <- e <- stats::setNames(rep(0, length(cand_ids)), cand_ids)
  if (kind %in% c("genomic", "both")) {
    for (sid in cand_ids) {
      acc <- tryCatch(get_trial_accessions(sid, conn, settings), error = function(e) character())
      g[sid] <- length(intersect(acc, trial$accessions))
    }
    g <- g / (max(g) + 1e-9)
  }
  if (kind %in% c("environmental", "both") && is.finite(trial$lat)) {
    cat <- trial_catalog(conn, settings) |> dplyr::filter(study_db_id %in% cand_ids)
    for (sid in cand_ids) {
      r <- cat[cat$study_db_id == sid, ]
      if (!nrow(r)) next
      d2 <- (suppressWarnings(as.numeric(r$latitude)) - trial$lat)^2 +
            (suppressWarnings(as.numeric(r$longitude)) - trial$long)^2 +
            (0.1 * (suppressWarnings(as.integer(r$year)) - trial$year))^2
      e[sid] <- exp(-d2 / 10)
    }
  }
  switch(kind, genomic = g, environmental = e, both = g + e)
}

# ===========================================================================
# Subtask B: preprocess phenotypes -> one target per accession
# ===========================================================================
build_targets <- function(cfg, train_obs, trial, conn, settings) {
  m <- cfg$pheno_prep.method
  obs <- train_obs

  # Per-trial outlier removal (for the BLUE-based methods).
  if (m %in% c("blue_lm", "trial_center") && is.finite(cfg$pheno_prep.z_thr)) {
    obs <- obs |>
      dplyr::group_by(study_id) |>
      dplyr::mutate(z = (value - mean(value)) / (stats::sd(value) + 1e-9)) |>
      dplyr::filter(abs(z) <= cfg$pheno_prep.z_thr) |>
      dplyr::ungroup() |>
      dplyr::select(-z)
  }
  if (identical(cfg$pheno_prep.standardize, "per_trial_z")) {
    obs <- obs |>
      dplyr::group_by(study_id) |>
      dplyr::mutate(value = (value - mean(value)) / (stats::sd(value) + 1e-9)) |>
      dplyr::ungroup()
  }

  # Per-trial, per-accession value (BLUE or mean) feeding the cross-trial step.
  per_trial <- switch(m,
    blue_lm        = .blue_per_trial(obs),
    two_stage_blup = .blue_per_trial(obs),
    trial_center   = obs |> dplyr::group_by(study_id) |>
                       dplyr::mutate(value = value - mean(value) + mean(train_obs$value)) |>
                       dplyr::ungroup() |>
                       dplyr::group_by(study_id, germplasm_name) |>
                       dplyr::summarise(value = mean(value), .groups = "drop"),
    raw_mean       = obs |> dplyr::group_by(study_id, germplasm_name) |>
                       dplyr::summarise(value = mean(value), .groups = "drop"),
    fatal(paste("unknown pheno_prep method", m), "bug_unknown_method"))

  # Optional G×E environmental weighting of trials when aggregating.
  w <- NULL
  if (identical(cfg$pheno_prep.ge_weighting, "env_gaussian") && is.finite(trial$lat)) {
    sim <- .trial_similarity(trial, unique(per_trial$study_id), conn, settings, "environmental")
    w <- exp(log(pmax(sim, 1e-6)) / cfg$pheno_prep.ge_bandwidth)
  }

  if (m == "two_stage_blup" && requireNamespace("lme4", quietly = TRUE) &&
      dplyr::n_distinct(per_trial$study_id) > 1) {
    # Shrinkage across trials: random accession effect, fixed trial effect.
    fit <- lme4::lmer(value ~ (1 | germplasm_name) + study_id, data = per_trial)
    re <- lme4::ranef(fit)$germplasm_name[, 1]
    nm <- rownames(lme4::ranef(fit)$germplasm_name)
    return(stats::setNames(re + mean(per_trial$value), nm))
  }

  # Default cross-trial aggregation: (weighted) mean per accession.
  agg <- per_trial
  if (!is.null(w)) agg$wt <- w[agg$study_id] else agg$wt <- 1
  agg <- agg |>
    dplyr::group_by(germplasm_name) |>
    dplyr::summarise(value = stats::weighted.mean(value, wt), .groups = "drop")
  stats::setNames(agg$value, agg$germplasm_name)
}

# Per-trial BLUE: value ~ germ (+ rep + block when they vary); fallback to mean.
.blue_per_trial <- function(obs) {
  obs |>
    dplyr::group_by(study_id) |>
    dplyr::group_modify(~{
      d <- .x
      if (dplyr::n_distinct(d$germplasm_name) < 2)
        return(dplyr::summarise(dplyr::group_by(d, germplasm_name),
                                value = mean(value), .groups = "drop"))
      rhs <- "germplasm_name"
      if (dplyr::n_distinct(d$rep)   > 1) rhs <- c(rhs, "rep")
      if (dplyr::n_distinct(d$block) > 1) rhs <- c(rhs, "block")
      form <- stats::as.formula(paste("value ~ 0 +", paste(rhs, collapse = " + ")))
      fit <- tryCatch(stats::lm(form, data = d), error = function(e) NULL)
      if (is.null(fit)) return(dplyr::summarise(dplyr::group_by(d, germplasm_name),
                                                value = mean(value), .groups = "drop"))
      co <- stats::coef(fit)
      g <- co[grep("^germplasm_name", names(co))]
      names(g) <- sub("^germplasm_name", "", names(g))
      tibble::tibble(germplasm_name = names(g), value = as.numeric(g))
    }) |>
    dplyr::ungroup()
}

# Observed focal-trial BLUE used as the scoring target (simple, robust).
.per_acc_blue <- function(focal_obs) {
  b <- .blue_per_trial(focal_obs)
  stats::setNames(b$value, b$germplasm_name)
}

# ===========================================================================
# Subtask C: select genotyping data
# ===========================================================================
choose_geno_sources <- function(cfg, train_acc, test_acc, conn, settings) {
  need  <- union(train_acc, test_acc)
  thin  <- as.integer(cfg$geno_select.marker_thin)
  projs <- projects_for_accessions(need, conn, settings)
  if (!length(projs)) return(structure(list(), n_projects = 0L))

  m <- cfg$geno_select.method
  if (m == "best_single_project") {
    cover <- vapply(projs, function(pid) {
      d <- tryCatch(get_project_dosage(pid, need, conn, settings, thin), error = function(e) NULL)
      if (is.null(d)) 0L else length(intersect(rownames(d), need))
    }, integer(1))
    projs <- projs[which.max(cover)]
  }
  # focal_plus_onehop and all_projects both load every covering project; the
  # one-hop bridging is handled in the EM combine (it uses shared accessions).
  dl <- list()
  for (pid in projs) {
    d <- tryCatch(get_project_dosage(pid, need, conn, settings, thin), error = function(e) NULL)
    if (!is.null(d) && nrow(d) > 0) dl[[as.character(pid)]] <- d
  }
  # Carry how many covering projects we started from, so run_pipeline can tell a
  # genuinely ungenotyped trial (0 projects) from a "found projects but extracted
  # no dosage" bug signature (projects > 0, dosage 0).
  attr(dl, "n_projects") <- length(projs)
  dl
}

# ===========================================================================
# Subtask D: relationship matrix / kernel
# ===========================================================================
build_kernel <- function(cfg, dosage_list) {
  ridge <- cfg$kernel.ridge
  m <- cfg$kernel.method

  if (m == "em_combine" && length(dosage_list) > 1) {
    grms <- lapply(dosage_list, function(d) rrBLUP::A.mat(.qc_markers(d, cfg)))
    std  <- lapply(grms, function(g) g / mean(diag(g)))
    names_all <- unique(unlist(lapply(std, colnames)))
    idx  <- lapply(std, function(g) match(colnames(g), names_all))
    res  <- T3BrapiHelpers::covariance_combiner(
      partial_covs = std, var_indices = idx,
      degrees_freedom = vapply(dosage_list, nrow, integer(1)))
    K <- res$psi
    if (nrow(K) == length(names_all) + 1) K <- K[-1, -1]   # drop phantom var
    rownames(K) <- colnames(K) <- names_all
  } else if (m == "rkhs_gaussian") {
    X <- .merge_markers(dosage_list)
    X <- .qc_markers(X, cfg)
    D2 <- as.matrix(stats::dist(X))^2
    theta <- cfg$kernel.rkhs_theta / (stats::median(D2[D2 > 0]) + 1e-9)
    K <- exp(-theta * D2)
  } else {
    # vanRaden_single (also the fallback when em_combine has one project).
    X <- .merge_markers(dosage_list)
    K <- rrBLUP::A.mat(.qc_markers(X, cfg))
  }
  diag(K) <- diag(K) + ridge
  K
}

# Merge per-project dosage matrices on shared markers (intersection), stacking
# accessions; later projects' duplicate accessions are dropped.
.merge_markers <- function(dosage_list) {
  if (length(dosage_list) == 1) return(dosage_list[[1]])
  common <- Reduce(intersect, lapply(dosage_list, colnames))
  if (length(common) < 50) {
    # No usable shared marker set: fall back to the largest single project.
    return(dosage_list[[which.max(vapply(dosage_list, nrow, integer(1)))]])
  }
  mats <- lapply(dosage_list, function(d) d[, common, drop = FALSE])
  X <- do.call(rbind, mats)
  X[!duplicated(rownames(X)), , drop = FALSE]
}

# MAF and missingness QC + imputation, returning a centered-ready dosage matrix.
.qc_markers <- function(X, cfg) {
  miss <- colMeans(is.na(X))
  X <- X[, miss <= cfg$kernel.max_missing, drop = FALSE]
  af <- colMeans(X, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  X <- X[, maf >= cfg$kernel.maf, drop = FALSE]
  if (ncol(X) < 50) infeasible("too_few_markers", ncol(X))
  # Impute remaining missing by per-marker mean (rounded for the additive code).
  for (j in which(colMeans(is.na(X)) > 0)) {
    mu <- mean(X[, j], na.rm = TRUE)
    X[is.na(X[, j]), j] <- if (identical(cfg$kernel.impute, "mean_round")) round(mu) else mu
  }
  X
}

# ===========================================================================
# Subtask E: train model
# ===========================================================================
train_model <- function(cfg, y_train, K, train_in, test_in, trial, train_obs) {
  m <- cfg$model.method

  if (m %in% c("gblup_sommer_GE", "rkhs") && identical(cfg$model.include_E, "yes") &&
      requireNamespace("sommer", quietly = TRUE)) {
    out <- tryCatch(.fit_sommer_GE(y_train, K, train_in, train_obs),
                    error = function(e) NULL)
    if (!is.null(out)) return(out)
  }
  # Genomic-only GBLUP via rrBLUP (the reliable backbone for every GBLUP variant;
  # mixed.solve estimates the variance ratio by REML, i.e. the optimal ridge).
  Ktt <- K[train_in, train_in, drop = FALSE]
  sol <- rrBLUP::mixed.solve(y = as.numeric(y_train[train_in]), K = Ktt)
  list(kind = "gblup", u = stats::setNames(sol$u, train_in), mu = sol$beta)
}

.fit_sommer_GE <- function(y_train, K, train_in, train_obs) {
  # Multi-environment training frame: one row per (study, accession) BLUE, so the
  # study (environment) random effect can be fit alongside the genomic kernel.
  df <- train_obs |>
    dplyr::filter(germplasm_name %in% train_in) |>
    dplyr::group_by(study_id, germplasm_name) |>
    dplyr::summarise(value = mean(value), .groups = "drop") |>
    dplyr::mutate(germplasm_name = factor(germplasm_name, levels = rownames(K)),
                  study_id = factor(study_id)) |>
    as.data.frame()
  Kg <- K[levels(df$germplasm_name), levels(df$germplasm_name), drop = FALSE]
  fit <- sommer::mmer(value ~ 1,
                      random = ~ sommer::vsr(germplasm_name, Gu = Kg) + study_id,
                      data = df, verbose = FALSE)
  u <- fit$U$`u:germplasm_name`$value
  list(kind = "gblup", u = u[train_in], mu = fit$Beta$Estimate[1])
}

# ===========================================================================
# Subtask F: predict on the focal trial + post-process
# ===========================================================================
predict_test <- function(cfg, fit, K, train_in, test_in, targets, scheme, settings) {
  mu <- fit$mu
  if (length(train_in) < cfg$predict_post.min_overlap) {
    # Fallback: too little genotyped overlap to trust the genomic model.
    return(stats::setNames(rep(mean(targets), length(test_in)), test_in))
  }

  if (cfg$predict_post.method == "cond_expectation") {
    Ktt <- K[train_in, train_in, drop = FALSE]
    K21 <- K[test_in, train_in, drop = FALSE]
    Kinv <- tryCatch(solve(Ktt), error = function(e) MASS::ginv(Ktt))
    if (!is.matrix(Kinv)) { dimnames(Kinv) <- dimnames(Ktt) }
    u_test <- as.numeric(K21 %*% Kinv %*% fit$u[train_in])
    pred <- stats::setNames(mu + u_test, test_in)
  } else {
    # direct BLUP: use the model's own random effect where available.
    u_test <- fit$u[test_in]
    u_test[is.na(u_test)] <- 0
    pred <- stats::setNames(mu + u_test, test_in)
  }

  # Blend observed BLUE for test accessions that appear in training (CV0 only;
  # under CV00 those phenotypes were masked, so blending is a no-op there).
  w <- cfg$predict_post.blend_obs_w
  if (scheme == "CV0" && w > 0) {
    seen <- intersect(test_in, names(targets))
    pred[seen] <- w * targets[seen] + (1 - w) * pred[seen]
  }
  pred
}
