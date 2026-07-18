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
  # QC + allele frequencies come from each panel's full population; the GRM itself
  # covers only the accessions we need to relate.
  K <- build_kernel(cfg, dosage_list, union(train_acc, focal_acc))
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
    primary <- .find_related(focal_id, conn, settings, cfg$train_select.primary_min)
    if (identical(cfg$train_select.primary_only, "yes")) return(primary)
    # Secondary tier: pool the germplasm across ALL primary trials, then find trials
    # sharing >= secondary_min accessions with that combined set -- one lookup against
    # the primary pool's germplasm rather than a separate hop from each primary trial.
    # A candidate is judged on its overlap with the pool as a whole, so a trial that
    # shares a few accessions with each of several primary trials can now qualify.
    prim_acc <- unique(unlist(lapply(primary, function(sid)
                        get_trial_accessions(sid, conn, settings)), use.names = FALSE))
    secondary <- .find_related(prim_acc, conn, settings, cfg$train_select.secondary_min,
                               input = "accessions")
    return(union(primary, secondary))
  }
  if (m == "top_k_similar") {
    cand <- .find_related(focal_id, conn, settings, 1)
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

# Trials sharing >= min_common germplasm with a query set of accessions. `x` is either
# a trial id (input = "trial"; the query set is that trial's accessions, and the trial
# itself is excluded from the result) or an accession vector directly
# (input = "accessions"; e.g. the pooled germplasm of several trials).
#
# Derived from the acc_<sid> cache rather than from
# T3BrapiHelpers::find_other_studies_evaluating_same_germplasm(), which spends one
# wizard query PER GERMPLASM (223 queries for trial 10676) to answer for exactly one
# trial. Here one wizard query returns the candidate trials, and the overlap counts
# are set intersections of accession lists we already cache -- and those lists are
# shared across every trial and with .trial_similarity, so the cost decays to nothing
# as the cache warms.
#
# Two differences from the helper, both deliberate: overlap is counted in germplasm
# NAME space (what acc_<sid> holds, and what every downstream join uses) rather than
# germplasm_db_id space, and candidates are restricted to the focal-trait catalogue --
# a trial that never measured the trait is not a training trial, and used to pad
# run_pipeline's min_train_trials check while contributing no observations.
.overlap_memo <- new.env(parent = emptyenv())

.find_related <- function(x, conn, settings, min_common, input = c("trial", "accessions")) {
  input <- match.arg(input)
  if (input == "trial") {
    sid <- as.character(x)
    # `get_trial_accessions` is passed lazily: on a memo hit .germ_overlap never
    # forces it, so a repeat call costs nothing.
    n <- .germ_overlap(get_trial_accessions(sid, conn, settings), conn, settings,
                       exclude = sid, memo_key = sid)
  } else {
    acc <- unique(as.character(x))
    n <- .germ_overlap(acc, conn, settings,
                       memo_key = paste0("acc_", substr(rlang::hash(sort(acc)), 1, 12)))
  }
  if (!length(n)) return(character())
  names(n)[n >= max(1L, as.integer(min_common))]
}

# Named integer vector: candidate trial id -> germplasm shared with `accessions`.
# Memoized per session (in RAM, not on disk -- it is derived data), because the
# optimizer calls this on the same query set for every config and scheme; `memo_key`
# is the trial id (trial input) or a hash of the accession set (accession input).
.germ_overlap <- function(accessions, conn, settings, exclude = character(), memo_key = NULL) {
  if (!is.null(memo_key)) {
    hit <- .overlap_memo[[memo_key]]
    if (!is.null(hit)) return(hit)              # `accessions` promise left unforced
  }

  n <- tryCatch({
    acc <- unique(as.character(accessions))
    if (!length(acc)) stop("no accessions")

    # Candidate trials: any trial sharing >= 1 accession. The breeder wizard unions
    # over the accessions it is given, so batch them (as projects_for_accessions does)
    # rather than asking once per germplasm.
    batches <- split(acc, ceiling(seq_along(acc) / 500L))
    cand <- unique(unlist(purrr::map(batches, function(b) {
      w <- .brapi_try(function() conn$wizard("trials", list(accessions = b)),
                      tries = settings$brapi_tries %||% 4L, what = "trials wizard")
      as.character(w$data$ids)
    })))
    # Keep only trials that measured the focal trait, and drop the excluded trial(s).
    cand <- setdiff(intersect(cand, as.character(trial_catalog(conn, settings)$study_db_id)),
                    exclude)
    if (!length(cand)) stop("no candidates")

    stats::setNames(purrr::map_int(cand, function(sid) {
      length(intersect(get_trial_accessions(sid, conn, settings), acc))
    }, .progress = "Germplasm overlap: candidate trials"), cand)
  }, error = function(e) integer())

  # Memoize only real answers: caching an empty result would let one transient network
  # failure make a query set permanently neighbourless for the rest of the run.
  if (!is.null(memo_key) && length(n)) assign(memo_key, n, envir = .overlap_memo)
  n
}

# Similarity of candidate trials to the focal trial: genomic = number of shared
# accessions; environmental = Gaussian on scaled lat/long/year distance.
.trial_similarity <- function(trial, cand_ids, conn, settings, kind) {
  g <- e <- stats::setNames(rep(0, length(cand_ids)), cand_ids)
  if (kind %in% c("genomic", "both")) {
    # One accession fetch per candidate trial (network on a cold cache) -> progress.
    g[cand_ids] <- purrr::map_dbl(cand_ids, function(sid) {
      acc <- tryCatch(get_trial_accessions(sid, conn, settings), error = function(e) character())
      length(intersect(acc, trial$accessions))
    }, .progress = "Similarity: candidate trials")
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
# Returns one dosage matrix per PROTOCOL GROUP, each holding the group's FULL
# genotyped population (every sample in the constituent projects), not just the
# accessions this trial needs. Subtask D estimates marker QC and allele frequencies
# from that population and only then subsets the GRM -- estimating them from a
# handful of needed accessions is what made both unstable.
choose_geno_sources <- function(cfg, train_acc, test_acc, conn, settings) {
  need  <- union(train_acc, test_acc)
  thin  <- as.integer(cfg$geno_select.marker_thin)
  projs <- projects_for_accessions(need, conn, settings)
  if (!length(projs)) return(structure(list(), n_projects = 0L))

  # Load the WHOLE population of each covering project (keep_samples = NULL). This
  # is the heaviest loop in the pipeline -- one VCF download + parse per project on
  # a cold cache -- so it reports progress per project.
  dl <- purrr::map(projs, function(pid) {
    d <- tryCatch(get_project_dosage(pid, NULL, conn, settings, thin), error = function(e) NULL)
    if (!is.null(d) && nrow(d) > 0 && length(intersect(rownames(d), need))) d else NULL
  }, .progress = "Load project dosages")
  names(dl) <- as.character(projs)
  dl <- purrr::compact(dl)
  if (!length(dl)) return(structure(list(), n_projects = length(projs)))

  if (identical(cfg$geno_select.method, "best_single_project")) {
    cover <- vapply(dl, function(d) length(intersect(rownames(d), need)), integer(1))
    dl <- dl[which.max(cover)]
  }
  # focal_plus_onehop and all_projects both keep every covering project; the one-hop
  # bridging is handled in the EM combine (it uses shared accessions).

  groups <- .group_by_panel(dl, settings)
  out <- purrr::map(groups, function(pids) {
    .merge_markers(.prune_redundant(dl[pids], settings))
  })
  names(out) <- vapply(groups, function(pids) paste(pids, collapse = "+"), character(1))
  # Carry how many covering projects we started from, so run_pipeline can tell a
  # genuinely ungenotyped trial (0 projects) from a "found projects but extracted
  # no dosage" bug signature (projects > 0, dosage 0).
  attr(out, "n_projects") <- length(projs)
  out
}

# Group projects into protocols by MARKER OVERLAP, not by protocol id: a "V2"/"v2.1"
# protocol is the same protocol scored against a different reference genome, so it
# gets its own id but (near-)identical markers. Two projects join the same group when
# they share >= settings$merge_containment of the smaller panel; grouping is
# transitive (single-linkage).
.group_by_panel <- function(dl, settings) {
  ids <- names(dl)
  if (length(ids) <= 1) return(as.list(ids))
  mk  <- lapply(dl, colnames)
  thr <- settings$merge_containment %||% 0.95

  grp <- seq_along(ids)                                   # union-find by relabelling
  for (i in seq_along(ids)) for (j in seq_len(i - 1)) {
    shared <- length(intersect(mk[[i]], mk[[j]]))
    if (shared >= thr * min(length(mk[[i]]), length(mk[[j]]))) {
      grp[grp == grp[i]] <- grp[j]
    }
  }
  unname(split(ids, grp))
}

# Drop redundant projects WITHIN a protocol group -- the same lines re-called on the
# same panel. Scoped to a group on purpose: across panels, two projects with the same
# accessions are complementary (different markers), not redundant, and dropping one
# would throw away a whole marker set.
.prune_redundant <- function(dl, settings) {
  if (length(dl) <= 1) return(dl)
  thr <- settings$redundant_acc_overlap %||% 0.90
  repeat {
    ids <- names(dl)
    drop <- NULL
    for (i in seq_along(ids)) {
      for (j in seq_len(i - 1)) {
        a <- rownames(dl[[i]]); b <- rownames(dl[[j]])
        shared <- length(intersect(a, b))
        if (shared < thr * min(length(a), length(b))) next
        loser <- if (setequal(a, b)) {
          # identical lines -> keep the richer marker build (the v1/v2 case)
          if (ncol(dl[[i]]) >= ncol(dl[[j]])) j else i
        } else if (length(a) != length(b)) {
          # near-identical lines -> keep the project with more accessions
          if (length(a) > length(b)) j else i
        } else {
          if (ncol(dl[[i]]) >= ncol(dl[[j]])) j else i
        }
        drop <- ids[loser]
        break
      }
      if (!is.null(drop)) break
    }
    if (is.null(drop)) return(dl)
    dl <- dl[setdiff(names(dl), drop)]
    if (length(dl) <= 1) return(dl)
  }
}

# ===========================================================================
# Subtask D: relationship matrix / kernel
# ===========================================================================
# `dosage_list` holds one FULL-population matrix per protocol group; `need` is the
# accession set the pipeline actually has to relate (training + focal). Marker QC and
# allele frequencies come from the population; the GRM is then formed over `need`.
build_kernel <- function(cfg, dosage_list, need) {
  ridge <- cfg$kernel.ridge
  m <- cfg$kernel.method

  if (m == "em_combine" && length(dosage_list) > 1) {
    # Bridge accessions -- genotyped in >1 protocol group -- are the shared rows through
    # which covariance_combiner stitches the separately-built per-panel GRMs. Keep them in
    # each GRM alongside the needed accessions (they may not be needed themselves), then
    # drop them from the combined result. `bridge` is empty when panels are disjoint, in
    # which case this reduces to the old need-only behaviour.
    bridge <- .bridge_accessions(dosage_list)
    keep   <- union(need, bridge)
    grms <- purrr::compact(lapply(dosage_list, function(d) .vanraden(.qc_markers(d, cfg), keep)))
    if (!length(grms)) infeasible("no_geno_overlap_in_panels")
    if (length(grms) == 1) {
      K <- grms[[1]]
    } else {
      std  <- lapply(grms, function(g) g / mean(diag(g)))
      names_all <- unique(unlist(lapply(std, colnames)))
      idx  <- lapply(std, function(g) match(colnames(g), names_all))
      res  <- T3BrapiHelpers::covariance_combiner(
        partial_covs = std, var_indices = idx,
        degrees_freedom = vapply(std, nrow, integer(1)))
      K <- res$psi
      if (nrow(K) == length(names_all) + 1) K <- K[-1, -1]   # drop phantom var
      rownames(K) <- colnames(K) <- names_all
    }
    # Eliminate the bridge accessions now that stitching is done -- they were scaffolding.
    final <- intersect(rownames(K), need)
    if (length(final) < 2) infeasible("no_geno_overlap_in_panels")
    K <- K[final, final, drop = FALSE]
  } else if (m == "rkhs_gaussian") {
    X <- .qc_markers(.best_panel(dosage_list, need), cfg)
    rows <- intersect(rownames(X), need)
    if (length(rows) < 2) infeasible("no_geno_overlap_in_panels")
    D2 <- as.matrix(stats::dist(X[rows, , drop = FALSE]))^2
    # Bandwidth from the POPULATION's typical distance, so theta means the same thing
    # regardless of how many accessions this trial happens to need. Subsampled: the
    # full population distance matrix would be O(n_pop^2).
    ref  <- rownames(X)
    if (length(ref) > 300) ref <- sample(ref, 300)
    Dref <- as.matrix(stats::dist(X[ref, , drop = FALSE]))^2
    theta <- cfg$kernel.rkhs_theta / (stats::median(Dref[Dref > 0]) + 1e-9)
    K <- exp(-theta * D2)
  } else {
    # vanRaden_single (also the fallback when em_combine has one panel).
    K <- .vanraden(.qc_markers(.best_panel(dosage_list, need), cfg), need)
    if (is.null(K)) infeasible("no_geno_overlap_in_panels")
  }
  diag(K) <- diag(K) + ridge
  K
}

# Accessions genotyped in >1 protocol group -- present in the rownames of >=2 of the
# dosage matrices. These are the shared rows covariance_combiner uses to stitch the
# per-panel GRMs together (em_combine); they are kept in the GRMs sent to the combiner
# and dropped from its result. Empty when the panels have disjoint accession sets.
.bridge_accessions <- function(dosage_list) {
  counts <- table(unlist(lapply(dosage_list, rownames), use.names = FALSE))
  names(counts)[counts >= 2L]
}

# The protocol group covering the most of `need`.
.best_panel <- function(dosage_list, need) {
  if (length(dosage_list) == 1) return(dosage_list[[1]])
  cover <- vapply(dosage_list, function(d) length(intersect(rownames(d), need)), integer(1))
  dosage_list[[which.max(cover)]]
}

# VanRaden GRM over `need`, centred and scaled on the allele frequencies of the FULL
# population in X. Because K_ij depends on the other accessions only through p, this
# is exactly the [need, need] submatrix of the whole population's GRM -- at
# O(n_need^2 * m) instead of O(n_pop^2 * m).
#
# NOT rrBLUP::A.mat: A.mat documents {-1,0,1} coding and our dosages are {0,1,2}. Its
# internal freq = (colMeans(X)+1)/2 then equals p + 0.5, so its MAF = pmin(freq, 1-freq)
# goes NEGATIVE for every marker with alt frequency > ~0.5 and its min.MAF filter drops
# them silently -- about half the panel, leaving a GRM whose off-diagonals correlate
# only ~0.7 with the correct one. (Centering happened to survive; the scaling did not.)
#
# The {0,1,2} assumption has not gone away -- it has MOVED here, into `p`. So guard it:
# a negative entry means someone changed the dosage encoding (e.g. "fixed" .vcf_to_dosage
# to rrBLUP's {-1,0,1}), which would make p meaningless and reproduce the same class of
# silent bug. Fail loudly instead.
.vanraden <- function(X, need) {
  rows <- intersect(rownames(X), need)
  if (length(rows) < 2) return(NULL)
  rng <- suppressWarnings(range(X, na.rm = TRUE))
  if (is.finite(rng[1]) && rng[1] < 0)
    fatal(sprintf("dosages must be coded {0,1,2}; saw a minimum of %g", rng[1]),
          "bug_dosage_coding")
  p <- colMeans(X) / 2                                   # population allele frequency
  denom <- 2 * sum(p * (1 - p))
  if (!is.finite(denom) || denom <= 0) return(NULL)
  W <- sweep(X[rows, , drop = FALSE], 2, 2 * p, "-")
  tcrossprod(W) / denom
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
  # Richest marker build first, so that when the same accession appears in two
  # projects the row kept below comes from the better-genotyped one.
  ord  <- order(vapply(dosage_list, ncol, integer(1)), decreasing = TRUE)
  mats <- lapply(dosage_list[ord], function(d) d[, common, drop = FALSE])
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
  # mixed.solve returns u as an n x 1 matrix and beta as a 1 x 1 matrix; coerce to plain
  # numerics so downstream `mu + u_test` is scalar + vector, not array + vector (which R
  # rejects as "non-conformable arrays").
  list(kind = "gblup", u = stats::setNames(as.numeric(sol$u), train_in),
       mu = as.numeric(sol$beta))
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
  mu <- as.numeric(fit$mu)                          # never a 1x1 matrix (see train_model)
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
