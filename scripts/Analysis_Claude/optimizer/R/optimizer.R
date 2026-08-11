# optimizer.R
#
# The model-based + evolutionary search. choose_config() decides which configuration to
# evaluate next, given everything tried so far:
#
#   1. Any unevaluated seed, for a baseline.
#   2. Any configuration short of its config_replication evaluations.
#   3. Fresh random configurations until there are enough scores to fit a surrogate.
#   4. Otherwise fit the random-forest surrogate, generate candidates by crossover and
#      mutation of the elites plus fresh random draws, and take the highest Expected
#      Improvement.
#
# choose_trial() picks the trial, given the configuration chosen here.

library(tidyverse)

# Restrict a set of evaluations to those whose trial falls in the current target
# domain. evals.sqlite is a shared archive of EVERY evaluation ever run (across
# target domains); a specific optimization only wants evidence from its own domain,
# on the premise that the best pipeline for a domain is the one to deploy on a new
# trial FROM that domain. The domain attributes were stored on each row under the
# catalogue's own column names, so the identical predicate that restricts trial
# sampling (.apply_target_domain) restricts the eval slice here. A NULL domain (e.g.
# simulate mode) is no constraint -> all rows. Rows with unknown attributes (NA,
# e.g. pre-migration or simulated) are out-of-domain whenever a constraint is set.
filter_evals_to_domain <- function(evals, td) {
  if (is.null(td) || !nrow(evals)) return(evals)
  .apply_target_domain(evals, td)
}

# Restrict evaluations to a single CV scheme. evals.sqlite is a shared archive that
# may hold rows for BOTH schemes (e.g. from an earlier run that optimized the other
# one); a given optimization targets exactly one, since CV0 and CV00 are distinct
# tasks that can have different optimal pipelines. Composes with filter_evals_to_domain.
filter_evals_to_scheme <- function(evals, scheme) {
  if (is.null(scheme) || !nrow(evals)) return(evals)
  dplyr::filter(evals, scheme == !!scheme)
}

# What each build changed, and which configurations it invalidated. A stored score is only
# comparable with a new one if no intervening change altered what that configuration COMPUTES.
# Most changes touch one method, so a blanket "ignore everything older" throws away good rows
# -- hence a predicate per change rather than a version cutoff.
#
# `affects(cfg)` returns TRUE when a config's result would differ under the new build.
# Add an entry whenever you bump OPTIMIZER_BUILD in a way that changes results.
BUILD_CHANGES <- list(
  list(build = "0.7.1",
       what  = "em_combine df: accession count -> centred effective_n",
       affects = function(cfg) identical(cfg$kernel.method, "em_combine")),
  list(build = "0.7.1",
       what  = ".best_panel now requires focal coverage, not just max union coverage",
       affects = function(cfg) isTRUE(cfg$kernel.method %in% c("vanRaden_single", "rkhs_gaussian"))),
  list(build = "0.7.3",
       what  = paste(".trial_similarity no longer returns NA for a candidate trial with an",
                     "unknown year/coordinates; env_gaussian weighting produced NA targets",
                     "and hence non-finite predictions before this"),
       affects = function(cfg) identical(cfg$pheno_prep.ge_weighting, "env_gaussian"))
)

# Drop rows a later build invalidated. A row goes when it was produced BEFORE a change's
# build and its config matches that change's predicate. `build` NA/empty predates stamping
# and so counts as older than everything. Composes with the two filters above.
filter_evals_to_build <- function(evals, build, changes = BUILD_CHANGES) {
  if (!nrow(evals) || is.null(build) || !length(changes)) return(evals)
  row_build <- if ("build" %in% names(evals)) as.character(evals$build) else rep(NA_character_, nrow(evals))
  # numeric_version handles 0.7.2 < 0.10.1 correctly; NA sorts oldest via the explicit test.
  older <- function(a, b) is.na(a) || !nzchar(a) ||
    isTRUE(tryCatch(numeric_version(a) < numeric_version(b), error = function(e) TRUE))
  cfgs <- lapply(evals$config_json, function(j)
    tryCatch(config_from_json(j), error = function(e) list()))
  drop <- rep(FALSE, nrow(evals))
  for (ch in changes) {
    hit <- vapply(seq_len(nrow(evals)), function(i)
      older(row_build[i], ch$build) && isTRUE(ch$affects(cfgs[[i]])), logical(1))
    drop <- drop | hit
  }
  evals[!drop, , drop = FALSE]
}

# Mean score per distinct configuration (over the trials it was run on, within the
# already domain/scheme-filtered slice), with the parsed config attached. Failed
# evaluations (NA score) count toward n but not the mean; a config that only ever
# failed gets mean_score = NA.
aggregate_scores <- function(evals, min_n_test = 10L, adjust_trial = TRUE) {
  if (!nrow(evals)) {
    return(tibble::tibble(config_hash = character(), mean_score = numeric(),
                          pooled = numeric(), unweighted = numeric(),
                          n = integer(), n_ok = integer(), config_json = character()))
  }
  # A score is a CORRELATION, so its precision depends on how many accessions it was computed
  # over -- and n_test ranges from 5 (SD(r) ~ 0.7, pure noise) to several hundred. A plain
  # mean gave those equal weight, so one lucky small trial could carry a configuration into
  # get_elites. Pool the textbook way instead: Fisher z, weight n - 3, back-transform.
  # `unweighted` is kept alongside so the two can be compared on a real run.
  nt <- suppressWarnings(as.numeric(evals$n_test %||% NA_real_))
  evals <- evals |>
    dplyr::mutate(
      .n_test = ifelse(is.na(nt), 0, nt),
      # Drop scores too small to inform anything; NA n_test predates the column, so keep it
      # (its score is real, we just cannot weight it) with the minimum weight.
      .usable = is.finite(score) & (is.na(nt) | .n_test >= min_n_test),
      .z      = atanh(pmin(pmax(score, -0.999), 0.999)),
      .w      = pmax(ifelse(is.na(nt), 1, .n_test - 3), 1))
  out <- evals |>
    dplyr::group_by(config_hash) |>
    dplyr::summarise(
      pooled      = if (any(.usable)) tanh(sum(.w[.usable] * .z[.usable]) / sum(.w[.usable]))
                    else NA_real_,
      unweighted  = mean(score[is.finite(score)]),
      n           = dplyr::n(),
      n_ok        = sum(.usable),
      config_json = dplyr::first(config_json),
      .groups = "drop")
  out$mean_score <- out$pooled
  attr(out, "estimator") <- "pooled"

  # Trial-adjusted estimate. `pooled` is confounded with WHICH trials a config drew, and
  # trials vary more than configs -- LESSONS #19 -- so a config that drew easy trials
  # outranks a better one that drew hard ones. A two-way random-effects fit removes that, and
  # shrinks the config BLUPs in proportion to replication, so a single lucky draw cannot post
  # an extreme value.
  #
  # Falls back to `pooled` whenever lme4 is unavailable or the design cannot support the fit.
  # Every failure path must return a value: this cannot take a run down.
  if (isTRUE(adjust_trial)) {
    adj <- tryCatch(.blup_scores(evals), error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(adj)) {
      i <- match(out$config_hash, adj$config_hash)
      ok <- !is.na(i) & is.finite(adj$blup[i])
      out$mean_score[ok] <- adj$blup[i][ok]
      attr(out, "estimator")  <- "blup"
      attr(out, "var_comps")  <- attr(adj, "var_comps")
    }
  }
  out
}

# Config BLUPs from a two-way random-effects fit, plus the variance components. NULL when the
# design or the packages cannot support it -- the caller then keeps its pooled estimate.
.blup_scores <- function(evals) {
  if (!requireNamespace("lme4", quietly = TRUE)) return(NULL)
  d <- evals[evals$.usable, , drop = FALSE]
  if (!nrow(d) || dplyr::n_distinct(d$trial_id) < 2L) return(NULL)
  # Identifiability: with one config per trial the trial effect and the residual are the same
  # term. lme4 errors on this ("number of levels ... must be < number of observations"), so
  # check rather than rely on the error.
  per_trial <- tapply(d$config_hash, d$trial_id, function(x) length(unique(x)))
  if (max(per_trial) < 2L) return(NULL)
  if (dplyr::n_distinct(d$config_hash) < 2L) return(NULL)

  m <- suppressMessages(lme4::lmer(.z ~ 1 + (1 | trial_id) + (1 | config_hash),
                                   data = d, weights = .w))
  re <- lme4::ranef(m)$config_hash
  v  <- as.data.frame(lme4::VarCorr(m))
  sd_of <- function(g) { x <- v$sdcor[v$grp == g]; if (length(x)) x[1] else NA_real_ }
  structure(
    tibble::tibble(config_hash = rownames(re),
                   blup = tanh(as.numeric(lme4::fixef(m)[1]) + re[, 1])),
    var_comps = c(sd_trial = sd_of("trial_id"), sd_config = sd_of("config_hash"),
                  sd_resid = sd_of("Residual")))
}

# Top-k configurations by mean score (only those with at least one success).
get_elites <- function(agg, k = 8) {
  agg |>
    dplyr::filter(is.finite(mean_score)) |>
    dplyr::arrange(dplyr::desc(mean_score)) |>
    head(k) |>
    dplyr::pull(config_json) |>
    lapply(config_from_json)
}

# Draw a fresh random configuration whose hash is not already in `avoid`.
fresh_random <- function(avoid, max_tries = 200) {
  for (i in seq_len(max_tries)) {
    cfg <- sample_config()
    if (!(config_hash(cfg) %in% avoid)) return(cfg)
  }
  sample_config()  # give up on novelty after many collisions
}

# Generate a pool of candidate configurations from the elites: crossovers,
# mutations, and random draws. Repaired and de-duplicated; novel configs
# preferred (those not yet tried), but the pool is never empty.
propose_candidates <- function(agg, elites, avoid,
                               n_cross = 60, n_mut = 60, n_rand = 60) {
  cands <- list()

  if (length(elites) >= 2) {
    for (i in seq_len(n_cross)) {
      pair <- sample(length(elites), 2)
      cands[[length(cands) + 1L]] <- repair_config(
        crossover(elites[[pair[1]]], elites[[pair[2]]]))
    }
  }
  if (length(elites) >= 1) {
    for (i in seq_len(n_mut)) {
      e <- elites[[sample(length(elites), 1)]]
      cands[[length(cands) + 1L]] <- mutate_config(e, n_subtasks = sample(1:2, 1))
    }
  }
  for (i in seq_len(n_rand)) cands[[length(cands) + 1L]] <- sample_config()

  # De-duplicate by hash; prefer novel candidates but keep some retries possible.
  hashes <- vapply(cands, config_hash, character(1))
  cands  <- cands[!duplicated(hashes)]
  hashes <- hashes[!duplicated(hashes)]
  novel  <- cands[!(hashes %in% avoid)]
  if (length(novel) >= 5) novel else cands
}

# Best configuration so far: highest mean score among configs evaluated on at
# least `min_reps` trials (robust to a lucky single trial); falls back to the
# plain best if none reach min_reps yet.
incumbent_config <- function(agg, min_reps = 2) {
  robust <- agg |> dplyr::filter(is.finite(mean_score), n_ok >= min_reps)
  pick <- if (nrow(robust)) robust else dplyr::filter(agg, is.finite(mean_score))
  if (!nrow(pick)) return(NULL)
  best <- pick |> dplyr::arrange(dplyr::desc(mean_score)) |> head(1)
  list(config = config_from_json(best$config_json),
       mean_score = best$mean_score, n_ok = best$n_ok)
}

# This worker's position in the pool, for spreading a shared backlog across workers. Worker ids
# are free-form, so anything unparseable degrades to 1.
.worker_index <- function(settings) {
  w <- suppressWarnings(as.integer(settings$worker_id %||% "1"))
  if (!is.finite(w) || w < 1L) 1L else w
}

# Decide the trial for the next evaluation.
#
# A trial must reach settings$trial_replication DISTINCT configurations before the search
# spends steps on fresh trials. This is the PRECONDITION for trial_id as a surrogate feature,
# not a refinement of it -- LESSONS #19. Revisits are also cheaper: that trial's observations
# and dosage matrices are already cached.
#
# `cfg_hash` is the configuration this trial is for, when the caller has one; trials it has
# already been run on are then excluded.
choose_trial <- function(con, settings, conn = NULL, cfg_hash = NULL) {
  r <- suppressWarnings(as.integer(settings$trial_replication %||% 1L))
  if (!is.finite(r)) r <- 1L
  if (r <= 1L && is.null(cfg_hash)) return(sample_trial(settings, conn))

  td <- if (isTRUE(settings$simulate)) NULL else settings$target_domain
  evals <- read_evals(con) |>
             filter_evals_to_domain(td) |>
             filter_evals_to_scheme(settings$optimize_scheme) |>
             filter_evals_to_build(settings$build %||% OPTIMIZER_BUILD)
  seen <- if (is.null(cfg_hash)) character()
          else unique(as.character(evals$trial_id[evals$config_hash == cfg_hash]))
  if (r <= 1L || !nrow(evals)) return(.sample_unseen_trial(settings, conn, seen))

  # Trials that have been started but not yet replicated to `r` distinct configurations.
  backlog <- evals |>
    dplyr::group_by(trial_id) |>
    dplyr::summarise(n_cfg = dplyr::n_distinct(config_hash), .groups = "drop") |>
    dplyr::filter(n_cfg >= 1L, n_cfg < r, !(trial_id %in% seen)) |>
    dplyr::arrange(trial_id)
  if (!nrow(backlog)) return(.sample_unseen_trial(settings, conn, seen))

  # Offset by worker: without it every worker revisits the same trial and the backlog drains
  # one trial at a time. Wrapping is right here -- a one-entry backlog wants every worker on it.
  id <- backlog$trial_id[((.worker_index(settings) - 1L) %% nrow(backlog)) + 1L]

  if (isTRUE(settings$simulate)) return(.sim_trial(id))
  tryCatch(build_trial_descriptor(id, conn, settings),
           error = function(e) {
             message("  revisit of trial ", id, " failed (", conditionMessage(e),
                     ") -- sampling a fresh trial")
             .sample_unseen_trial(settings, conn, seen)
           })
}

# A fresh trial, preferring one not in `seen`: the real pipeline is deterministic in
# (config, trial, scheme), so repeating a pair recomputes a known score and double-counts it.
# After `tries` take what we have rather than fail a step.
.sample_unseen_trial <- function(settings, conn, seen = character(), tries = 5L) {
  tr <- sample_trial(settings, conn)
  for (i in seq_len(max(0L, tries - 1L))) {
    if (!(tr$id %in% seen)) break
    tr <- sample_trial(settings, conn)
  }
  tr
}

# Decide the next configuration to evaluate. Returns list(cfg, source[, ei]).
choose_config <- function(con, settings) {
  # The store is global; this optimization only learns from its own target domain.
  # Filtering here scopes EVERYTHING downstream -- seeds, the done/avoid set,
  # aggregation, the surrogate, and the incumbent -- to the in-domain slice, so a
  # config evaluated only in some OTHER domain counts as untried here and gets run
  # to build in-domain evidence.
  # In simulate mode there is no real target domain (synthetic trials carry no
  # program/location/year), so no domain filtering applies. The scheme filter always
  # applies: this run optimizes exactly one scheme, and the store may hold rows for
  # the other from a previous run.
  td          <- if (isTRUE(settings$simulate)) NULL else settings$target_domain
  evals       <- read_evals(con) |>
                   filter_evals_to_domain(td) |>
                   filter_evals_to_scheme(settings$optimize_scheme) |>
                   filter_evals_to_build(settings$build %||% OPTIMIZER_BUILD)
  done_hashes <- unique(evals$config_hash)

  # Phase 1: seed the five submissions.
  # Prediction5 submitted a different model per scheme, so the seeds are scheme-specific.
  #
  # Each worker takes a DIFFERENT un-done seed. This selection knows only what has FINISHED,
  # never what is in flight, so "the first un-done seed" is the same answer for every worker
  # that asks before the first one stores anything -- and an evaluation runs far longer than
  # run_workers.sh's 20 s launch stagger. Six workers then spent six evaluations on
  # Prediction1 before any other seed was touched. Offsetting by worker id gets all five
  # seeds running at once and still accumulates reps, one per seed per round.
  seeds  <- seed_configs(settings$optimize_scheme)
  hashes <- vapply(seeds, config_hash, character(1))
  undone <- seeds[!(hashes %in% done_hashes)]
  if (length(undone)) {
    i <- ((.worker_index(settings) - 1L) %% length(undone)) + 1L  # wraps past the seed count
    return(list(cfg = undone[[i]], source = paste0("seed:", names(undone)[i])))
  }

  # Phase 2: replication. A configuration gets settings$config_replication evaluations before
  # the search moves on -- LESSONS #19. The backlog counts EVALUATIONS, not distinct trials, so
  # it always drains; choose_trial() is what spreads them over distinct trials. Counting trials
  # would stall wherever only one exists, as under sim_fixed_trial.
  cr <- suppressWarnings(as.integer(settings$config_replication %||% 1L))
  if (is.finite(cr) && cr > 1L && nrow(evals)) {
    backlog <- evals |>
      dplyr::count(config_hash, name = "n_eval") |>
      dplyr::filter(n_eval < cr) |>
      dplyr::arrange(config_hash)
    # Workers past the end of the backlog explore instead of wrapping onto it: wrapping would
    # put all 22 of them on a one-entry backlog and hand that config 22 evaluations.
    w <- .worker_index(settings)
    if (w <= nrow(backlog)) {
      json <- evals$config_json[match(backlog$config_hash[w], evals$config_hash)]
      return(list(cfg = config_from_json(json), source = "replicate"))
    }
  }

  agg <- aggregate_scores(evals)
  n_scored <- sum(is.finite(agg$mean_score))

  # Phase 3: random exploration until the surrogate has enough to learn from.
  if (n_scored < settings$n_random_init) {
    return(list(cfg = fresh_random(done_hashes), source = "random_init"))
  }

  # Phase 4: surrogate-guided. Which fit is safe depends on the DESIGN of the store:
  #   BLOCKED  one row per (config, trial), trial_id a factor. Removes the trial variance from
  #            every config's estimate and keeps the replication weighting that collapsing to
  #            means throws away -- LESSONS #19.
  #   POOLED   one row per config mean, when the design cannot support trial_id
  #            (trial_feature_usable).
  scored  <- agg |> dplyr::filter(is.finite(mean_score))
  rows    <- dplyr::filter(evals, is.finite(score))
  # Drop rows whose trial has not yet been replicated, keeping trial_id in the model for the
  # rest. (Measured: unreplicated trials do NOT in fact damage the fit -- this is the
  # conservative reading, and trial_replication drains the backlog promptly, so few rows go.)
  if (nrow(rows)) {
    n_cfg <- tapply(rows$config_hash, rows$trial_id, function(x) length(unique(x)))
    rows  <- rows[n_cfg[as.character(rows$trial_id)] >= 2L, , drop = FALSE]
  }
  blocked <- nrow(rows) > 0 && dplyr::n_distinct(rows$trial_id) >= 2L
  if (blocked) {
    feats  <- configs_to_features(lapply(rows$config_json, config_from_json))
    feats$trial_id <- factor(rows$trial_id)
    surro  <- fit_surrogate(feats, rows$score,
                            ntree = settings$ntree, min_obs = settings$n_random_init)
    if (is.null(surro)) blocked <- FALSE          # too few rows after dropping -> pooled path
  }
  if (!blocked) {
    cfgs   <- lapply(scored$config_json, config_from_json)
    feats  <- configs_to_features(cfgs)
    surro  <- fit_surrogate(feats, scored$mean_score,
                            ntree = settings$ntree, min_obs = settings$n_random_init)
  }
  if (is.null(surro)) return(list(cfg = fresh_random(done_hashes), source = "random_fallback"))

  elites  <- get_elites(agg, settings$n_elites)
  cands   <- propose_candidates(agg, elites, done_hashes,
                                settings$n_cross, settings$n_mut, settings$n_rand)
  cfeats  <- configs_to_features(cands)
  # Score a candidate as its EXPECTED value over the observed trials -- the objective is
  # "best on a random trial", so trial_id is marginalised away rather than chosen.
  pr      <- if (blocked) predict_surrogate_marginal(surro, cfeats, levels(feats$trial_id))
             else         predict_surrogate(surro, cfeats)
  best    <- max(scored$mean_score)
  ei      <- expected_improvement(pr$mean, pr$sd, best, settings$ei_xi)

  pick <- which.max(ei)
  list(cfg = cands[[pick]], source = if (blocked) "acquisition_blocked" else "acquisition",
       ei = ei[pick], pred_mean = pr$mean[pick], pred_sd = pr$sd[pick])
}
