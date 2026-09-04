# pipeline.R
#
# The parameterized genomic-prediction pipeline. run_pipeline() executes the six subtasks for
# one configuration on one focal trial under one CV scheme, returning predicted and observed
# per-accession yields for scoring. Each subtask is a dispatcher branching on the method named
# in the configuration -- to add a method, add a branch here and an entry in R/config_space.R.
#
# A subtask that cannot be completed for a given trial/config throws; evaluate.R records the
# failure so the optimizer learns to avoid it.
#
# Needs live T3 access and cannot be exercised offline: run it once on a known trial and
# confirm the predicted/observed join is non-empty before launching the loop.

library(tidyverse)

# Feasibility floors on the GENOTYPED overlap, in ACCESSIONS -- not TRIALS, which
# min_train_trials counts. Single-sourced so .best_panel selects against the same guard
# run_pipeline applies.
.min_train <- function(settings) as.integer(settings$min_train_acc %||% 20L)
.min_test  <- function(settings) as.integer(settings$min_test_acc  %||% 5L)

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
  # Count USABLE targets, not just present ones. A non-finite target is not a training
  # datum -- it poisons the model fit and comes back out as a non-finite prediction, which
  # scores as an empty test set far downstream. Drop them here, where the cause is legible.
  n_na    <- sum(!is.finite(targets))
  targets <- targets[is.finite(targets)]
  if (length(targets) < 10)
    infeasible("too_few_train_acc",
               sprintf("%d finite%s", length(targets),
                       if (n_na > 0) sprintf(" (%d non-finite dropped)", n_na) else ""))

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
  K <- build_kernel(cfg, dosage_list, union(train_acc, focal_acc), settings, focal_acc)
  geno_acc <- rownames(K)
  train_in <- intersect(train_acc, geno_acc)
  test_in  <- intersect(focal_acc, geno_acc)
  if (length(train_in) < .min_train(settings) || length(test_in) < .min_test(settings)) {
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

# CV0 keeps every training observation; CV00 also removes any observation on a focal-trial
# accession (Jarquin et al. 2017). Factored out so the distinction is unit-testable.
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

# Trials sharing >= min_common germplasm with a query set. `x` is a trial id (input = "trial";
# the query set is that trial's accessions and the trial itself is excluded) or an accession
# vector (input = "accessions").
#
# Derived from the acc_<sid> cache, not from find_other_studies_evaluating_same_germplasm(),
# which costs a wizard query per germplasm -- LESSONS #17.
#
# Two invariants: overlap is counted in germplasm NAME space, and candidates are restricted to
# the focal-trait catalogue -- a trial that never measured the trait is not a training trial.
.overlap_memo <- new.env(parent = emptyenv())

# Emit an informational note at most once per `key` per session (so a recurring geno-source
# observation -- e.g. a marker-poor panel -- is reported once, not on every evaluation).
.geno_note_seen <- new.env(parent = emptyenv())
.note_geno_once <- function(key, msg) {
  if (is.null(.geno_note_seen[[key]])) { message(msg); assign(key, TRUE, envir = .geno_note_seen) }
}

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

    cat_ids <- as.character(trial_catalog(conn, settings)$study_db_id)
    idx     <- tryCatch(.trial_index(settings), error = function(e) NULL)
    local_ok <- .index_covers(idx, cat_ids, settings, "acc")

    # One tabulation answers BOTH questions the wizard was used for: which trials share an
    # accession (its names) and how many they share (its values). Profiled on the server, the
    # wizard call this replaces was 440 s of a 795 s evaluation -- 99.9% of .find_related.
    tab <- NULL
    if (!is.null(idx)) {
      hits <- unlist(idx[intersect(acc, names(idx))], use.names = FALSE)
      if (length(hits)) tab <- table(hits)
    }

    if (local_ok) {
      # Local discovery: the candidates ARE the tabulated trials. No BrAPI call at all.
      .note_geno_once("trial_discovery_mode", sprintf(
        "trial discovery: LOCAL (index covers the %d-trial catalogue) -- no trials wizard call",
        length(cat_ids)))
      cand <- setdiff(intersect(names(tab) %||% character(), cat_ids), exclude)
      if (!length(cand)) stop("no candidates")
      counts <- stats::setNames(as.integer(tab[cand]), cand)
      counts[is.na(counts)] <- 0L
      counts
    } else {
      # Fallback: ask the wizard which trials to consider, exactly as before. Reached when the
      # index does not yet cover the catalogue (run tools/prepare_indices.R) or discovery is
      # pinned to "wizard".
      .note_geno_once("trial_discovery_mode", sprintf(
        "trial discovery: WIZARD (index covers %d of %d catalogue trials) -- run tools/prepare_indices.R",
        if (is.null(idx)) 0L else sum(cat_ids %in% union(attr(idx, "keys") %||% character(),
                                                         .attempted(settings, "acc"))),
        length(cat_ids)))
      batches <- split(acc, ceiling(seq_along(acc) / 500L))
      cand <- unique(unlist(purrr::map(batches, function(b) {
        w <- .brapi_try(function() conn$wizard("trials", list(accessions = b)),
                        conn = conn, settings = settings, what = "trials wizard")
        as.character(w$data$ids)
      })))
      cand <- setdiff(intersect(cand, cat_ids), exclude)
      if (!length(cand)) stop("no candidates")

      # Count from the index where it knows the candidate; only the rest need the old
      # per-trial read + intersect. Degrading candidate by candidate matters: falling back for
      # ALL of them because ONE is unindexed put the whole cost straight back.
      known  <- if (is.null(idx)) character() else (attr(idx, "keys") %||% character())
      counts <- stats::setNames(integer(length(cand)), cand)
      in_idx <- cand %in% known
      if (any(in_idx) && !is.null(tab)) {
        v <- as.integer(tab); names(v) <- names(tab)
        got <- v[cand[in_idx]]; got[is.na(got)] <- 0L
        counts[in_idx] <- as.integer(got)
      }
      miss <- cand[!in_idx]
      if (length(miss))
        counts[miss] <- purrr::map_int(miss, function(sid) {
          length(intersect(get_trial_accessions(sid, conn, settings), acc))
        }, .progress = "Germplasm overlap: unindexed candidates")
      counts
    }
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
      # A candidate with unknown coordinates or year is UNRELATED (similarity 0), never NA:
      # a quarter of the T3 catalogue has year = NA, and an NA similarity becomes an NA
      # weight in build_targets, whose weighted.mean has no na.rm.
      e[sid] <- if (is.finite(d2)) exp(-d2 / 10) else 0
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
# Extra marker thinning per project, so that all of a trial's covering projects TOGETHER fit
# dosage_total_budget_bytes. Returns thin factors keyed by project (1 = untouched). This is the
# only place the SUM is bounded -- dosage_budget_bytes caps one project, and the pipeline holds
# them all at once. Served by column-subsetting the cache, so it costs no re-download.
#
# Every project gets the SAME factor. Panels in a protocol group merge on their SHARED markers,
# and thinning near-identical panels differently keeps different markers from each, collapsing
# the intersection the merge needs. The cost is that small panels are thinned too.
#
# Sizes come from the cached stat_<pid>. A project never parsed has no stat entry and is left
# at 1; get_project_dosage applies its own budget when it downloads it.
.dosage_thin_plan <- function(projs, settings) {
  ids  <- as.character(projs)
  plan <- as.list(stats::setNames(rep(1L, length(ids)), ids))
  cap  <- settings$dosage_total_budget_bytes %||% Inf
  if (!is.finite(cap) || cap <= 0) return(plan)

  # Size each project by what is ACTUALLY ON DISK, not by dosage_budget_bytes: a cache built
  # under a larger budget is reused as-is, since get_project_dosage only re-parses to make a
  # cache denser. Sizing from the budget would under-count it and the cap would not fire.
  thin <- stats::setNames(rep(NA_real_, length(ids)), ids)
  bytes <- vapply(ids, function(pid) {
    st <- .project_stat(settings, pid)
    if (is.null(st) || is.null(st$n_markers)) return(NA_real_)
    dense <- as.numeric(st$n_samples) * as.numeric(st$n_markers) * 4
    cd <- .find_densest_dosage(settings, pid)
    # Not cached yet -> it will be parsed at THIS machine's budget when downloaded.
    e  <- if (!is.null(cd)) cd$thin else
            .cache_thin(st$n_samples, st$n_markers, settings$dosage_budget_bytes)
    thin[[pid]] <<- e
    dense / e
  }, numeric(1))

  known <- bytes[is.finite(bytes)]
  total <- sum(known)
  if (!length(known) || total <= cap) return(plan)

  # One extra factor for every project, so panels stay marker-aligned for .merge_markers
  # (see above). marker_thin is an ABSOLUTE thin, and get_project_dosage serves a cache at
  # thin e by keeping every floor(marker_thin / e)-th column -- so the target must be
  # e * factor, not the factor alone, or a project already cached coarse gets no reduction.
  factor <- max(1L, as.integer(ceiling(total / cap)))
  for (pid in names(known)) plan[[pid]] <- max(1L, as.integer(thin[[pid]] * factor))
  kept <- total / factor
  message(sprintf(
    "geno: %d covering project(s) total %.1f GB > dosage_total_budget_bytes %.1f GB -- keeping 1 marker in %d (%.1f GB)",
    length(known), total / 1e9, cap / 1e9, factor, kept / 1e9))
  plan
}

# Returns one dosage matrix per PROTOCOL GROUP, each holding the group's FULL genotyped
# population rather than only the accessions this trial needs: subtask D estimates marker QC
# and allele frequencies from the population and subsets the GRM afterwards.
choose_geno_sources <- function(cfg, train_acc, test_acc, conn, settings) {
  need  <- union(train_acc, test_acc)
  projs <- projects_for_accessions(need, conn, settings)
  if (!length(projs)) return(structure(list(), n_projects = 0L))

  # Load each covering project's WHOLE population -- QC and allele frequencies are estimated
  # on it, LESSONS #12. The heaviest loop in the pipeline on a cold cache, and where the
  # memory peak is set: every covering project is resident at once, which is what
  # .dosage_thin_plan bounds.
  plan <- .dosage_thin_plan(projs, settings)
  dl <- purrr::map(projs, function(pid) {
    d <- tryCatch(get_project_dosage(pid, NULL, conn, settings,
                                     marker_thin = plan[[as.character(pid)]] %||% 1L),
                  error = function(e) NULL)
    if (!is.null(d) && nrow(d) > 0 && length(intersect(rownames(d), need))) d else NULL
  }, .progress = "Load project dosages")
  names(dl) <- as.character(projs)
  dl <- purrr::compact(dl)
  if (!length(dl)) return(structure(list(), n_projects = length(projs)))

  # Drop a project with intrinsically < 50 markers: QC only removes markers and every kernel
  # method needs >= 50 surviving, so such a panel can never yield a GRM under any config.
  # Skipping it at the source keeps it from sinking a downstream combine; noted once per
  # session. A panel that fails only under a strict QC config is config-dependent and is left
  # to build_kernel.
  poor <- names(dl)[vapply(dl, ncol, integer(1)) < 50L]
  for (pid in poor)
    .note_geno_once(paste0("poormarker_", pid), sprintf(
      "geno: project %s has < 50 markers -- unusable for any GRM; skipping it", pid))
  dl <- dl[setdiff(names(dl), poor)]
  if (!length(dl)) return(structure(list(), n_projects = length(projs)))

  # "Best" is .best_panel_index's criterion, not max coverage of `need`: a project thick with
  # training lines and holding none of the focal ones cannot predict the trial, however large.
  # Subset rather than extract, so the project name survives for .group_by_panel.
  if (identical(cfg$geno_select.method, "best_single_project"))
    dl <- dl[.best_panel_index(dl, need, test_acc, .min_test(settings), .min_train(settings))]

  groups <- .group_by_panel(dl, settings)
  out <- purrr::map(groups, function(pids) {
    .merge_markers(.prune_redundant(dl[pids], settings))
  })
  names(out) <- vapply(groups, function(pids) paste(pids, collapse = "+"), character(1))
  # `out` holds merged copies, so the per-project originals are now dead weight -- several GB
  # Drop the reference and collect, or they stay resident through the kernel build.
  rm(dl); gc(full = TRUE)

  # focal_plus_onehop: keep only the panels reachable from the focal trial's own panel.
  # all_projects keeps everything; this is what makes the two methods different.
  if (identical(cfg$geno_select.method, "focal_plus_onehop"))
    out <- .onehop_filter(out, test_acc, cfg$geno_select.min_bridge)
  # Lets run_pipeline tell a genuinely ungenotyped trial (0 projects) from the bug signature
  # of projects found but no dosage extracted.
  attr(out, "n_projects") <- length(projs)
  out
}

# The "one hop" of `focal_plus_onehop`, applied to MERGED protocol groups: bridging is between
# panels, and two projects in one group share a panel and need no bridge.
#
# Seed = the group(s) genotyping the focal trial; admit a further group only if it shares at
# least `min_bridge` accessions with the seed. Those shared lines are what covariance_combiner
# stitches through, so a group without them costs markers and compute but says nothing about
# the test set. Deliberately ONE hop -- the transitive closure is all_projects by another name.
#
# With no group covering the focal accessions there is nothing to hop from: return the list
# unchanged and let run_pipeline's overlap check judge the trial.
.onehop_filter <- function(dl, test_acc, min_bridge) {
  if (length(dl) <= 1) return(dl)
  mb <- suppressWarnings(as.integer(min_bridge))
  if (!isTRUE(is.finite(mb))) mb <- 1L

  covers_focal <- vapply(dl, function(d) length(intersect(rownames(d), test_acc)) > 0, logical(1))
  if (!any(covers_focal)) return(dl)

  seed_acc <- unique(unlist(lapply(dl[covers_focal], rownames), use.names = FALSE))
  bridged  <- vapply(dl, function(d) length(intersect(rownames(d), seed_acc)) >= mb, logical(1))
  dl[covers_focal | bridged]
}

# Group projects into protocols by MARKER OVERLAP, not by protocol id: a "V2"/"v2.1"
# protocol is the same protocol scored against a different reference genome, so it
# gets its own id but (near-)identical markers. Two projects join the same group when
# they share >= settings$merge_containment of the smaller panel; grouping is
# transitive (single-linkage).
# Memo for .group_by_panel, per session and in RAM. The grouping is a property of the PANELS,
# not of the configuration, so it is identical across evaluations touching the same projects.
.panel_group_memo <- new.env(parent = emptyenv())

.group_by_panel <- function(dl, settings) {
  ids <- names(dl)
  if (length(ids) <= 1) return(as.list(ids))
  mk  <- lapply(dl, colnames)
  thr <- settings$merge_containment %||% 0.95

  # Key on the project ids AND their marker counts: the same project served at a different
  # marker density (.dosage_thin_plan) is a different panel and can group differently, so ids
  # alone would be an unsound key. Counts are O(1) to read and pin the density.
  o   <- order(ids)
  key <- rlang::hash(list(ids[o], vapply(mk[o], length, integer(1)), thr))
  hit <- .panel_group_memo[[key]]
  if (!is.null(hit)) return(hit)

  grp <- seq_along(ids)                                   # union-find by relabelling
  for (i in seq_along(ids)) for (j in seq_len(i - 1)) {
    shared <- length(intersect(mk[[i]], mk[[j]]))
    if (shared >= thr * min(length(mk[[i]]), length(mk[[j]]))) {
      grp[grp == grp[i]] <- grp[j]
    }
  }
  out <- unname(split(ids, grp))
  assign(key, out, envir = .panel_group_memo)
  out
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
build_kernel <- function(cfg, dosage_list, need, settings = NULL, focal = NULL) {
  ridge <- cfg$kernel.ridge
  m <- cfg$kernel.method

  if (m == "em_combine" && length(dosage_list) > 1) {
    # Bridge accessions -- genotyped in >1 protocol group -- are the shared rows through
    # which covariance_combiner stitches the separately-built per-panel GRMs. Keep them in
    # each GRM alongside the needed accessions (they may not be needed themselves), then
    # drop them from the combined result. Dropping bridges BEFORE the combine collapses
    # em_combine to block-diagonal with no error (docs/LESSONS.md #13). `bridge` is empty when
    # the panels are disjoint.
    bridge <- .bridge_accessions(dosage_list)
    keep   <- union(need, bridge)
    # Build a GRM per panel, dropping one that cannot yield it -- too few markers after THIS
    # config's QC, or too little overlap -- rather than letting one weak panel sink the
    # combine. Config-dependent, so not persisted; a marker-poor project is filtered earlier
    # in choose_geno_sources.
    grms <- purrr::compact(purrr::imap(dosage_list, function(d, nm) {
      g <- tryCatch(.vanraden(.qc_markers(d, cfg), keep), error = function(e) e)
      if (is.null(g) || inherits(g, "error")) {
        why <- if (inherits(g, "error")) conditionMessage(g) else "too little needed/bridge overlap"
        .note_geno_once(sprintf("emdrop_%s_maf%s_miss%s", nm, cfg$kernel.maf, cfg$kernel.max_missing),
          sprintf("em_combine: dropped panel %s under this QC (%s)", nm, why))
        return(NULL)
      }
      g
    }))
    if (!length(grms)) infeasible("no_geno_overlap_in_panels")
    if (length(grms) == 1) {
      K <- grms[[1]]
    } else {
      # Standardize each partial covariance to a unit mean diagonal, THEN regularize it. The
      # ridge must go on HERE, not only on the combined matrix: the EM combiner inverts these
      # partials -- LESSONS #21.
      # Diagnostic only -- the ridge above makes either case combinable. Duplicate genotype
      # rows mean two accession NAMES on one genotype: an unresolved synonym (LESSONS #6) or
      # identical lines. Rank < nrow without duplicates means too few surviving markers.
      for (nm in names(grms)) {
        g <- grms[[nm]]
        dup <- sum(duplicated(round(g, 10)))
        if (nrow(g) > 1 && qr(g)$rank < nrow(g))
          .note_geno_once(paste0("rankdef_", nm), sprintf(
            "em_combine: panel %s partial covariance is rank %d of %d (%d duplicate row(s)) -- ridged before combining",
            nm, qr(g)$rank, nrow(g), dup))
      }
      # Standardize, MEASURE, then ridge, in that order: the ridge is a searchable parameter
      # large enough to lift a rank-deficient panel's near-zero eigenvalues and inflate its
      # effective_n, so measuring first keeps df a property of the panel, not of a knob.
      unridged <- lapply(grms, function(g) g / mean(diag(g)))
      std <- lapply(unridged, function(g) { diag(g) <- diag(g) + max(ridge, 1e-6); g })
      names_all <- unique(unlist(lapply(std, colnames)))
      idx  <- lapply(std, function(g) match(colnames(g), names_all))
      # df is each partial's RELATIVE WEIGHT in the Wishart-EM likelihood, so only the ratio
      # between panels matters. Accession count would weight a 2000-line panel 10:1 over a
      # 200-line one; markers are in LD, so that overstates its independent information. The
      # effective-sample-size measure, re-centred by .center_dfs, brings it to ~1.4:1.
      dfs <- .center_dfs(vapply(unridged, .effective_n, numeric(1)),
                         settings$em_df_mean %||% 60, settings$em_df_stdev %||% 15)
      # A combine that still fails is a property of THIS trial's panels, not a bug: record it
      # as an infeasibility so the loop moves on rather than surfacing as `status = "error"`.
      res <- tryCatch(
        T3BrapiHelpers::covariance_combiner(
          partial_covs = std, var_indices = idx, degrees_freedom = dfs),
        error = function(e) e)
      if (inherits(res, "error"))
        # Carry each partial's shape and rank into the failure record, so the cause is
        # readable from the store without re-running anything.
        infeasible("em_combine_failed",
                   paste0(conditionMessage(res), " | partials: ",
                          paste(vapply(std, function(g)
                            sprintf("%dx%d rank %d", nrow(g), ncol(g), qr(g)$rank),
                            character(1)), collapse = "; ")),
                   funnel = c(panels = length(std), accessions = length(names_all)))
      K <- res$psi
      if (nrow(K) == length(names_all) + 1) K <- K[-1, -1]   # drop phantom var
      rownames(K) <- colnames(K) <- names_all
    }
    # Eliminate the bridge accessions now that stitching is done -- they were scaffolding.
    final <- intersect(rownames(K), need)
    if (length(final) < 2) infeasible("no_geno_overlap_in_panels")
    K <- K[final, final, drop = FALSE]
  } else if (m == "rkhs_gaussian") {
    X <- .qc_markers(.best_panel(dosage_list, need, focal, .min_test(settings),
                                 .min_train(settings)), cfg)
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
    K <- .vanraden(.qc_markers(.best_panel(dosage_list, need, focal, .min_test(settings),
                                           .min_train(settings)), cfg), need)
    if (is.null(K)) infeasible("no_geno_overlap_in_panels")
  }
  diag(K) <- diag(K) + ridge
  K
}

# --- em_combine degrees of freedom -----------------------------------------
# Effective number of INDEPENDENT samples behind a relationship matrix: Galwey (2009)'s
# "effective number of independent tests", (sum sqrt(lambda))^2 / sum(lambda) over the
# positive eigenvalues. Scale-invariant (so standardizing to a unit mean diagonal does not
# move it) and bounded [1, rank]. Used instead of the accession count because markers are in
# LD: a 2000-line panel does not carry 10x the independent information of a 200-line one.
.effective_n <- function(G) {
  lambda <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  lambda <- lambda[lambda > 1e-8]                  # drop ~zero / tiny-negative
  if (length(lambda) < 1) return(1)
  (sum(sqrt(lambda)))^2 / sum(lambda)
}

# Turn effective sample sizes into EM degrees of freedom: keep their relative ORDERING but
# re-centre on mean_df and cap the spread at sd_df. Only the ratios between partials matter
# to the combiner (the common scale divides out of the M-step), so this is the knob that
# decides how far apart two panels' weights may get. Floored at 1 so every df stays positive.
.center_dfs <- function(m_eff, mean_df, sd_df) {
  if (length(m_eff) <= 1) return(rep(mean_df, length(m_eff)))
  obs_sd <- stats::sd(m_eff)
  if (!is.finite(obs_sd) || obs_sd == 0) return(rep(mean_df, length(m_eff)))
  z <- (m_eff - mean(m_eff)) / obs_sd
  pmax(mean_df + z * min(sd_df, obs_sd), 1)
}

# Accessions genotyped in >1 protocol group -- present in the rownames of >=2 of the
# dosage matrices. These are the shared rows covariance_combiner uses to stitch the
# per-panel GRMs together (em_combine); they are kept in the GRMs sent to the combiner
# and dropped from its result. Empty when the panels have disjoint accession sets.
.bridge_accessions <- function(dosage_list) {
  counts <- table(unlist(lapply(dosage_list, rownames), use.names = FALSE))
  names(counts)[counts >= 2L]
}

# Which of a set of dosage matrices to predict the focal trial from -- THE selection rule,
# shared by .best_panel (protocol groups, after merging) and choose_geno_sources's
# `best_single_project` branch (raw projects, before merging). One criterion, one
# implementation: the two drifted apart once, and a project selection that ignored focal
# coverage silently starved the kernel of test rows.
#
# Two floors, then the HARMONIC MEAN of focal and training coverage.
#
# Since train = setdiff(need, focal), f_cov + t_cov = |panel n need| exactly -- so this is not a
# new pair of quantities, it is the SUM criterion replaced by a harmonic mean of the same two.
# The sum is wrong because the two coverages do different jobs and are not interchangeable:
# training lines buy prediction ACCURACY, focal lines buy the PRECISION with which that accuracy
# is measured (SE of a correlation goes as 1/sqrt(n_test - 3)). Summing lets a panel thick with
# training lines outscore one that actually covers the trial -- on canary 10676, panel 14511
# (focal 11, train 525) beat panel 11049 (focal 186, train 94) on the sum, and an 11-point
# correlation carries noise several times the config effect the search is trying to resolve.
#
# The harmonic mean is the effective sample size when precision is limited by both groups
# (1/n_eff = (1/n_f + 1/n_t)/2), and it saturates in the larger argument the way prediction
# accuracy saturates in training size: with t >> f it tends to 2f, i.e. "maximise focal, ignore
# training past the floor". The exactly-right criterion is mean GBLUP reliability over the test
# lines (c' P^-1 c, which captures relatedness rather than mere counts), but that needs a GRM
# built and inverted PER CANDIDATE PANEL -- 46 of them on trial 10253 -- which is the cost this
# function exists to avoid.
#
# If nothing clears both floors, fall back to best focal coverage: without focal lines there is
# nothing to predict. `focal` NULL gives pure max-coverage, for callers with no focal set.
.best_panel_index <- function(dosage_list, need, focal = NULL,
                              min_test = 5L, min_train = 20L) {
  if (length(dosage_list) == 1) return(1L)
  cover <- vapply(dosage_list, function(d) length(intersect(rownames(d), need)), integer(1))
  if (is.null(focal) || !length(focal)) return(which.max(cover))
  train  <- setdiff(need, focal)
  f_cov  <- vapply(dosage_list, function(d) length(intersect(rownames(d), focal)), integer(1))
  t_cov  <- vapply(dosage_list, function(d) length(intersect(rownames(d), train)), integer(1))
  ok     <- f_cov >= min_test & t_cov >= min_train
  # Guarded against 0/0: `ok` already implies both are positive, and non-ok entries are -1.
  hmean  <- ifelse(ok, 2 * f_cov * t_cov / pmax(f_cov + t_cov, 1L), -1)
  if (any(ok)) return(which.max(hmean))
  which.max(f_cov)
}

# The protocol group to build a single-panel kernel from.
.best_panel <- function(dosage_list, need, focal = NULL,
                        min_test = 5L, min_train = 20L) {
  dosage_list[[.best_panel_index(dosage_list, need, focal, min_test, min_train)]]
}

# VanRaden GRM over `need`, centred and scaled on the allele frequencies of the FULL population
# in X. K_ij depends on the others only through p, so this is exactly the [need, need]
# submatrix of the population GRM, at O(n_need^2 * m) rather than O(n_pop^2 * m).
#
# Do NOT replace with rrBLUP::A.mat -- LESSONS #11.
#
# The {0,1,2} assumption lives here now, in `p`, so guard it: a negative entry means the
# dosage encoding changed underneath us, which would make p meaningless. Fail loudly.
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

# MAF and missingness QC + imputation, returning a centred-ready dosage matrix.
#
# The pipeline's largest allocation -- X is a whole genotyped population, several GB on the big
# panels, times the number of workers -- so it makes as few copies as it can:
#
#   * ONE subset, not two: both masks are computed up front and applied together.
#   * `miss` is reused for the imputation loop rather than recomputing colMeans(is.na(X)).
#   * INTEGER IS PRESERVED on the mean_round path, whose fills are whole numbers. Assigning a
#     double into a 4-byte integer matrix coerces the WHOLE matrix to 8 bytes on the first
#     imputed column; plain `mean` needs doubles, so there the promotion is done once.
.qc_markers <- function(X, cfg) {
  miss <- colMeans(is.na(X))
  af   <- colMeans(X, na.rm = TRUE) / 2
  maf  <- pmin(af, 1 - af)
  # which() drops the NA that an all-missing column produces in `maf`, so such a column is
  # excluded rather than becoming an NA index.
  keep <- which(miss <= cfg$kernel.max_missing & maf >= cfg$kernel.maf)
  if (length(keep) < 50) infeasible("too_few_markers", length(keep))
  X    <- X[, keep, drop = FALSE]
  miss <- miss[keep]

  # Impute remaining missing by per-marker mean (rounded for the additive code).
  na_cols <- which(miss > 0)
  if (length(na_cols)) {
    round_impute <- identical(cfg$kernel.impute, "mean_round")
    keep_integer <- round_impute && is.integer(X)
    if (!keep_integer && is.integer(X)) storage.mode(X) <- "double"
    for (j in na_cols) {
      mu   <- mean(X[, j], na.rm = TRUE)
      fill <- if (round_impute) round(mu) else mu
      X[is.na(X[, j]), j] <- if (keep_integer) as.integer(fill) else fill
    }
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
  # Genomic-only GBLUP backbone. `lambda_select` chooses the ridge, the variance ratio
  # lambda = sigma2_e / sigma2_u:
  #
  #   reml   rrBLUP::mixed.solve estimates it by REML -- optimal IF the model is correct.
  #   fixed  cfg$model.lambda_fixed as given.
  #   loo    minimise leave-one-out predictive MSE over a log grid.
  #
  # REML and LOO are different estimators, not two routes to one number: LOO shrinks harder
  # when the model is wrong, which here it is -- a block-diagonal or EM-stitched GRM, and
  # phenotypes pooled over trials with heterogeneous error variances. Hence worth searching.
  Ktt  <- K[train_in, train_in, drop = FALSE]
  y    <- as.numeric(y_train[train_in])
  rule <- cfg$model.lambda_select
  if (length(rule) != 1L || is.na(rule)) rule <- "reml"

  if (identical(rule, "reml")) {
    sol <- rrBLUP::mixed.solve(y = y, K = Ktt)
    # mixed.solve returns u as an n x 1 matrix and beta as a 1 x 1 matrix; coerce to plain
    # numerics so downstream `mu + u_test` is scalar + vector, not array + vector (which R
    # rejects as "non-conformable arrays").
    return(list(kind = "gblup", u = stats::setNames(as.numeric(sol$u), train_in),
                mu = as.numeric(sol$beta),
                lambda = as.numeric(sol$Ve) / as.numeric(sol$Vu)))
  }
  lambda <- if (identical(rule, "loo")) .loo_lambda(Ktt, y - mean(y))
            else as.numeric(cfg$model.lambda_fixed)
  if (!is.finite(lambda) || lambda <= 0) lambda <- 1
  .ridge_blup(y, Ktt, lambda, train_in)
}

# GBLUP at a GIVEN lambda: u = K (K + lambda I)^-1 (y - mean(y)), the kernel-ridge
# solution. Returns the same shape as the REML branch so predict_test cannot tell them
# apart.
.ridge_blup <- function(y, Ktt, lambda, ids) {
  mu <- mean(y)
  yc <- y - mu
  A  <- Ktt + lambda * diag(nrow(Ktt))
  u  <- tryCatch(as.numeric(Ktt %*% solve(A, yc)),
                 error = function(e) as.numeric(Ktt %*% MASS::ginv(A) %*% yc))
  list(kind = "gblup", u = stats::setNames(u, ids), mu = mu, lambda = lambda)
}

# Choose lambda by leave-one-out predictive MSE over a log-spaced grid (Prediction4's method
# and grid).
#
# No refitting: the fit is the linear smoother yhat = H y, H = K (K + lambda I)^-1, so the
# leave-one-out residual is exactly r_i / (1 - h_ii) by the PRESS identity, and one solve per
# grid point gives the full LOO error. `yc` must already be centred -- re-centring inside the
# loop would break the identity.
.loo_lambda <- function(Ktt, yc, grid = 10^seq(-4, 4, length.out = 30)) {
  n <- nrow(Ktt)
  best_lambda <- 1
  best_mse    <- Inf
  for (lam in grid) {
    H <- tryCatch(Ktt %*% solve(Ktt + lam * diag(n)), error = function(e) NULL)
    if (is.null(H)) next
    resid <- yc - as.numeric(H %*% yc)
    denom <- 1 - diag(H)
    denom[abs(denom) < 1e-10] <- 1e-10      # a point the smoother interpolates
    mse <- mean((resid / denom)^2)
    if (is.finite(mse) && mse < best_mse) { best_mse <- mse; best_lambda <- lam }
  }
  best_lambda
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

  # Conditional expectation through the kernel: mu + K21 K11^-1 u.
  .cond_exp <- function() {
    Ktt <- K[train_in, train_in, drop = FALSE]
    K21 <- K[test_in, train_in, drop = FALSE]
    Kinv <- tryCatch(solve(Ktt), error = function(e) MASS::ginv(Ktt))
    if (!is.matrix(Kinv)) { dimnames(Kinv) <- dimnames(Ktt) }
    stats::setNames(mu + as.numeric(K21 %*% Kinv %*% fit$u[train_in]), test_in)
  }

  if (cfg$predict_post.method == "cond_expectation") {
    pred <- .cond_exp()
  } else {
    # direct BLUP: the model's own random effect for the focal accessions, available only when
    # the fit has one -- the sommer G+E path carries every accession in rownames(K) as a level;
    # the rrBLUP backbone is fitted on training lines alone and does not.
    #
    # Without the guard, fit$u[test_in] is all-NA -> every prediction equals `mu` and scores
    # NA. Fall back to the kernel and RECORD it -- LESSONS #22.
    u_test <- fit$u[test_in]
    if (all(is.na(u_test))) {
      .note_geno_once("direct_blup_no_test_effect", paste(
        "predict: direct_blup has no random effect for the focal accessions (the model was",
        "fitted on training lines only) -- predicting through the kernel instead"))
      pred <- .cond_exp()
    } else {
      u_test[is.na(u_test)] <- 0
      pred <- stats::setNames(mu + u_test, test_in)
    }
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
