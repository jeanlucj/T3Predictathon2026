# optimizer.R
#
# The model-based + evolutionary search. One call to choose_config() decides
# which pipeline configuration to evaluate next, given everything tried so far:
#
#   1. While any of the five seeds is unevaluated, run that seed (baseline).
#   2. While there are too few scored configurations to fit a surrogate, run
#      fresh random configurations (exploration).
#   3. Otherwise: fit the random-forest surrogate on the per-configuration mean
#      scores, generate candidates by CROSSOVER and MUTATION of the current
#      elites plus fresh RANDOM draws, and pick the candidate with the highest
#      Expected Improvement.
#   4. With small probability, instead re-evaluate the incumbent on a new trial
#      to reduce the noise on the current leader.
#
# The driver (run_optimizer.R) supplies the trial; this module only chooses the
# configuration. Aggregating to a per-configuration mean before fitting the
# surrogate denoises the across-trial variability so the model sees each pipeline
# once with its average quality.

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
       affects = function(cfg) isTRUE(cfg$kernel.method %in% c("vanRaden_single", "rkhs_gaussian")))
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
aggregate_scores <- function(evals, min_n_test = 10L) {
  if (!nrow(evals)) {
    return(tibble::tibble(config_hash = character(), mean_score = numeric(),
                          unweighted = numeric(), n = integer(), n_ok = integer(),
                          config_json = character()))
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
  evals |>
    dplyr::group_by(config_hash) |>
    dplyr::summarise(
      mean_score  = if (any(.usable)) tanh(sum(.w[.usable] * .z[.usable]) / sum(.w[.usable]))
                    else NA_real_,
      unweighted  = mean(score[is.finite(score)]),
      n           = dplyr::n(),
      n_ok        = sum(.usable),
      config_json = dplyr::first(config_json),
      .groups = "drop")
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
  seeds <- seed_configs(settings$optimize_scheme)
  for (nm in names(seeds)) {
    if (!(config_hash(seeds[[nm]]) %in% done_hashes)) {
      return(list(cfg = seeds[[nm]], source = paste0("seed:", nm)))
    }
  }

  agg <- aggregate_scores(evals)
  n_scored <- sum(is.finite(agg$mean_score))

  # Phase 2: random exploration until the surrogate has enough to learn from.
  if (n_scored < settings$n_random_init) {
    return(list(cfg = fresh_random(done_hashes), source = "random_init"))
  }

  # Occasional re-evaluation of a top config on a new trial to denoise the leaders.
  # We re-run a random ELITE, not only the incumbent: each step scores ONE scheme
  # (one rep), so a config needs several selections to reach incumbent_min_reps. If
  # only the incumbent were ever re-run, no challenger could accumulate the reps to
  # overtake it and the incumbent would freeze on whichever config was denoised first.
  # Sampling across the elite pack lets challengers reach the reps bar and win.
  re_elites <- get_elites(agg, settings$n_elites)
  if (length(re_elites) && stats::runif(1) < settings$reeval_prob) {
    return(list(cfg = re_elites[[sample(length(re_elites), 1)]], source = "reeval_elite"))
  }

  # Phase 3: surrogate-guided. Fit on per-config mean scores.
  scored  <- agg |> dplyr::filter(is.finite(mean_score))
  cfgs    <- lapply(scored$config_json, config_from_json)
  feats   <- configs_to_features(cfgs)
  surro   <- fit_surrogate(feats, scored$mean_score,
                           ntree = settings$ntree, min_obs = settings$n_random_init)
  if (is.null(surro)) return(list(cfg = fresh_random(done_hashes), source = "random_fallback"))

  elites  <- get_elites(agg, settings$n_elites)
  cands   <- propose_candidates(agg, elites, done_hashes,
                                settings$n_cross, settings$n_mut, settings$n_rand)
  cfeats  <- configs_to_features(cands)
  pr      <- predict_surrogate(surro, cfeats)
  best    <- max(scored$mean_score)
  ei      <- expected_improvement(pr$mean, pr$sd, best, settings$ei_xi)

  pick <- which.max(ei)
  list(cfg = cands[[pick]], source = "acquisition",
       ei = ei[pick], pred_mean = pr$mean[pick], pred_sd = pr$sd[pick])
}
