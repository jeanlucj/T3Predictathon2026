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

# Restrict a set of evaluations to the run's trial universe. evals.sqlite is a shared archive
# of EVERY evaluation ever run (across target domains); a specific optimization only wants
# evidence from its own, on the premise that the best pipeline for a domain is the one to
# deploy on a new trial FROM that domain. The domain is DEFINED on the catalogue by
# .apply_target_domain() and resolved once per run into trial ids (pin_trial_universe), so
# membership here is a set test on an immutable id and nothing a trial's metadata does later
# can move a stored row in or out. A NULL universe (simulate mode) is no constraint.
filter_evals_to_universe <- function(evals, universe) {
  if (!length(universe) || !nrow(evals)) return(evals)
  dplyr::filter(evals, as.character(trial_id) %in% as.character(universe))
}

# Restrict evaluations to a single CV scheme. evals.sqlite is a shared archive that
# may hold rows for BOTH schemes (e.g. from an earlier run that optimized the other
# one); a given optimization targets exactly one, since CV0 and CV00 are distinct
# tasks that can have different optimal pipelines. Composes with filter_evals_to_universe.
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
                          pooled = numeric(), unweighted = numeric(), se = numeric(),
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
  out$se <- NA_real_                    # only the random-effects fit supplies one
  attr(out, "estimator") <- "pooled"
  # The median of the weights the fit actually uses. lmer reports sd_resid at unit weight, and
  # this is what converts it to a per-evaluation figure; computed here because only this
  # function knows which rows are usable and at what weight.
  attr(out, "median_weight") <-
    suppressWarnings(stats::median(evals$.w[evals$.usable], na.rm = TRUE))

  # Trial-adjusted estimate. `pooled` is confounded with WHICH trials a config drew, and
  # trials vary more than configs -- LESSONS #19 -- so a config that drew easy trials
  # outranks a better one that drew hard ones. A two-way random-effects fit removes that, and
  # shrinks the config BLUPs in proportion to replication, so a single lucky draw cannot post
  # an extreme value.
  #
  # Falls back to `pooled` whenever lme4 is unavailable or the design cannot support the fit.
  # Every failure path must return a value: this cannot take a run down. Whichever way it goes,
  # `estimator_note` carries the explanation, because a silent fallback also empties `se` --
  # and .contenders() drops every config without one, which disables contender replication.
  if (isTRUE(adjust_trial)) {
    # A WARNING IS NOT A FAILURE. lmer warns routinely here: sd_config is small against
    # sd_trial (LESSONS #19), which is exactly the regime that produces a boundary or
    # convergence warning. Discarding the fit for it threw away good estimates. Muffle and
    # record; only an error gives up on the fit.
    warn <- character()
    adj <- tryCatch(
      withCallingHandlers(
        .blup_scores(evals),
        warning = function(w) { warn <<- c(warn, conditionMessage(w))
                                invokeRestart("muffleWarning") }),
      error = function(e) { warn <<- c(warn, conditionMessage(e)); NULL })
    # One line: lme4's convergence text is multi-line and would break the report's bullet.
    note <- if (length(warn))
              trimws(gsub("\\s+", " ", paste(unique(warn), collapse = "; "))) else NULL
    if (!is.null(adj)) {
      i <- match(out$config_hash, adj$config_hash)
      ok <- !is.na(i) & is.finite(adj$blup[i])
      out$mean_score[ok] <- adj$blup[i][ok]
      out$se[ok] <- adj$se[i][ok]
      attr(out, "estimator")  <- "blup"
      attr(out, "var_comps")  <- attr(adj, "var_comps")
      attr(out, "estimator_note") <- note        # NULL when the fit ran clean
    } else {
      attr(out, "estimator_note") <-
        .blup_skip_reason(evals) %||% note %||% "the fit declined for an unrecorded reason"
    }
  }
  out
}

# Why .blup_scores() would decline, checked in the order it checks -- so the report can name
# the reason rather than print nothing. NULL means the fit can run.
.blup_skip_reason <- function(evals) {
  if (!requireNamespace("lme4", quietly = TRUE)) return("lme4 is not installed")
  d <- evals[evals$.usable, , drop = FALSE]
  if (!nrow(d)) return("no usable rows (a finite score and n_test >= min_n_test)")
  if (dplyr::n_distinct(d$trial_id) < 2L) return("only 1 distinct trial")
  # Identifiability: with one config per trial the trial effect and the residual are the same
  # term. lme4 errors on this ("number of levels ... must be < number of observations"), so
  # check rather than rely on the error.
  if (max(tapply(d$config_hash, d$trial_id, function(x) length(unique(x)))) < 2L)
    return("no trial has 2 configs -- the trial effect is not separable from the residual")
  if (dplyr::n_distinct(d$config_hash) < 2L) return("only 1 distinct config")
  NULL
}

# Config BLUPs from a two-way random-effects fit, plus the variance components. NULL when the
# design or the packages cannot support it -- the caller then keeps its pooled estimate.
.blup_scores <- function(evals) {
  if (!is.null(.blup_skip_reason(evals))) return(NULL)
  d <- evals[evals$.usable, , drop = FALSE]
  m <- suppressMessages(lme4::lmer(.z ~ 1 + (1 | trial_id) + (1 | config_hash),
                                   data = d, weights = .w))
  # condVar = TRUE attaches each level's posterior variance, which is what the contender rule
  # needs; it is on the z scale, as `blup` is before the tanh below.
  re <- lme4::ranef(m, condVar = TRUE)$config_hash
  se <- sqrt(as.numeric(attr(re, "postVar")))
  v  <- as.data.frame(lme4::VarCorr(m))
  sd_of <- function(g) { x <- v$sdcor[v$grp == g]; if (length(x)) x[1] else NA_real_ }
  structure(
    tibble::tibble(config_hash = rownames(re),
                   blup = tanh(as.numeric(lme4::fixef(m)[1]) + re[, 1]),
                   se   = se),
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

# Config hashes owed another evaluation, in a stable order, split into two tiers:
#   $base   short of settings$config_replication. Owed unconditionally -- this is the floor
#           every configuration gets, and throttling it would make the setting a suggestion.
#   $extra  a CONTENDER at or above the floor, earning one more each time this is called, so
#           the target rises for exactly the configs still in question. Rationed by the caller.
# Counts EVALUATIONS, not distinct trials, so it always drains. `universe` (eligible trial ids,
# NULL in simulate mode) removes a config already run on every eligible trial.
.replication_backlog <- function(evals, agg, settings, universe = NULL) {
  none <- list(base = character(), extra = character())
  cr <- suppressWarnings(as.integer(settings$config_replication %||% 1L))
  if (!is.finite(cr) || !nrow(evals)) return(none)
  cand <- .contenders(agg, settings$contender_z %||% 1, k = 8L)
  per <- evals |>
    dplyr::group_by(config_hash) |>
    dplyr::summarise(n_eval = dplyr::n(), n_trial = dplyr::n_distinct(trial_id), .groups = "drop")
  if (length(universe)) per <- dplyr::filter(per, n_trial < length(universe))
  per <- dplyr::arrange(per, config_hash)
  list(base  = per$config_hash[per$n_eval < cr],
       extra = per$config_hash[per$n_eval >= cr & per$config_hash %in% cand])
}

# Configs that cannot yet be ruled out as the best: their optimistic bound still reaches the
# leader's estimate. Returns config hashes, at most `k`, always including the leader -- so a
# result of length 1 means the field has been decided at this `z`. Empty when the estimator
# supplied no standard errors (the pooled fallback on a store too thin to fit).
.contenders <- function(agg, z = 1, k = 8L) {
  a <- agg[is.finite(agg$mean_score) & is.finite(agg$se %||% NA_real_), , drop = FALSE]
  if (!nrow(a)) return(character())
  ucb <- a$mean_score + z * a$se
  in_play <- ucb >= max(a$mean_score)          # >=: a bound that exactly reaches is not ruled out
  utils::head(a$config_hash[in_play][order(ucb[in_play], decreasing = TRUE)], k)
}

# Distinct configs a trial must reach before fresh trials are drawn. Rises with the store: the
# sqrt keeps the number of NEW trials growing as sqrt(n) instead of flattening to a constant.
.trial_target <- function(settings, n_scored) {
  base <- suppressWarnings(as.integer(settings$trial_replication %||% 1L))
  if (!is.finite(base)) base <- 1L
  if (base <= 1L) return(base)
  base + as.integer(floor(sqrt(max(0, n_scored)) / 5))
}

# This worker's position in the pool, for spreading a shared backlog across workers. Worker ids
# are free-form, so anything unparseable degrades to 1.
.worker_index <- function(settings) {
  w <- suppressWarnings(as.integer(settings$worker_id %||% "1"))
  if (!is.finite(w) || w < 1L) 1L else w
}

# Trials started but not yet replicated to `r` distinct configurations, in the order they are
# served. Restricted to the run's universe, so a trial that is no longer in it cannot hold the
# backlog open -- a non-empty backlog stops fresh trials being drawn at all.
.trial_backlog <- function(evidence, seen = character(), universe = NULL, r = 2L) {
  if (!nrow(evidence)) return(character())
  bl <- evidence |>
    dplyr::group_by(trial_id) |>
    dplyr::summarise(n_cfg = dplyr::n_distinct(config_hash), .groups = "drop") |>
    dplyr::filter(n_cfg >= 1L, n_cfg < r, !(as.character(trial_id) %in% as.character(seen)))
  if (length(universe))
    bl <- dplyr::filter(bl, as.character(trial_id) %in% as.character(universe))
  sort(as.character(bl$trial_id))
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
#
# `universe` is the run's pinned trial ids (pin_trial_universe); NULL is unconstrained.
choose_trial <- function(con, settings, conn = NULL, cfg_hash = NULL, universe = NULL) {
  base <- suppressWarnings(as.integer(settings$trial_replication %||% 1L))
  if (!is.finite(base)) base <- 1L
  if (base <= 1L && is.null(cfg_hash)) return(sample_trial(settings, conn))

  # Two views of the store, for two different questions. `computed` answers "has this pair
  # already been run", and is keyed exactly as claim_eval() keys it -- if the two disagree, a
  # pair can be proposed forever and claimed never. `evidence` answers "what does this run know
  # about the domain", and is the one the replication target is measured against.
  computed <- read_evals(con) |>
                filter_evals_to_scheme(settings$optimize_scheme) |>
                filter_evals_to_build(settings$build %||% OPTIMIZER_BUILD)
  evidence <- filter_evals_to_universe(computed, universe)
  # Trials this config has already been run on, PLUS the ones other workers are running it on
  # right now. `computed` holds only finished work, so without the claims a one-entry backlog
  # handed every worker the same trial and they all evaluated the identical pair. This is the
  # advisory half; claim_eval() in optimizer_step is the half that is actually atomic.
  seen <- if (is.null(cfg_hash)) character()
          else {
            done <- unique(as.character(computed$trial_id[computed$config_hash == cfg_hash]))
            live <- tryCatch(active_claims(con, settings$optimize_scheme),
                             error = function(e) NULL)
            union(done, if (is.null(live)) character()
                        else as.character(live$trial_id[live$config_hash == cfg_hash]))
          }
  if (base <= 1L || !nrow(evidence)) return(.sample_unseen_trial(settings, conn, seen))
  r <- .trial_target(settings, dplyr::n_distinct(evidence$config_hash[is.finite(evidence$score)]))

  backlog <- .trial_backlog(evidence, seen, universe, r)
  if (!length(backlog)) return(.sample_unseen_trial(settings, conn, seen))

  # Offset by worker: without it every worker revisits the same trial and the backlog drains
  # one trial at a time. Wrapping is right here -- a one-entry backlog wants every worker on it.
  id <- backlog[((.worker_index(settings) - 1L) %% length(backlog)) + 1L]

  if (isTRUE(settings$simulate)) return(.sim_trial(id))
  tryCatch(build_trial_descriptor(id, conn, settings),
           error = function(e) {
             message("  revisit of trial ", id, " failed (", conditionMessage(e),
                     ") -- sampling a fresh trial")
             .sample_unseen_trial(settings, conn, seen)
           })
}

# A trial not in `seen`: the real pipeline is deterministic in (config, trial, scheme), so
# repeating a pair recomputes a known score and double-counts it in that config's mean. Raises
# trials_exhausted() when every draw was already seen -- the caller then moves to another config.
.sample_unseen_trial <- function(settings, conn, seen = character(), tries = 5L) {
  for (i in seq_len(max(1L, tries))) {
    tr <- sample_trial(settings, conn)
    if (!(tr$id %in% seen)) return(tr)
  }
  trials_exhausted(length(seen))
}

# Decide the next configuration to evaluate. Returns list(cfg, source[, ei]).
#
# `universe`: eligible trial ids for the target domain (NULL in simulate mode), passed on to
# .replication_backlog(). `replicate = FALSE` skips the replication phase, for a caller whose
# replication pick turned out to have no trial left.
choose_config <- function(con, settings, universe = NULL, replicate = TRUE) {
  # The store is global; this optimization only learns from its own universe. Filtering here
  # scopes EVERYTHING downstream -- seeds, the done/avoid set, aggregation, the surrogate, and
  # the incumbent -- to that slice, so a config evaluated only in some OTHER domain counts as
  # untried here and gets run to build in-domain evidence. The scheme filter always applies:
  # this run optimizes exactly one scheme, and the store may hold rows for the other.
  evals       <- read_evals(con) |>
                   filter_evals_to_universe(universe) |>
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

  agg <- aggregate_scores(evals)
  n_scored <- sum(is.finite(agg$mean_score))

  # Phase 2: replication -- .replication_backlog() decides what is owed to whom. The base tier
  # is served whenever it is non-empty; only the contender tier is rationed, one worker in
  # `every`, keyed on the store's row count so a single worker alternates over time and N
  # workers split at any instant without lockstep.
  every <- max(1L, suppressWarnings(as.integer(settings$replicate_every %||% 1L)))
  w <- .worker_index(settings)
  if (isTRUE(replicate) && nrow(evals)) {
    bl <- .replication_backlog(evals, agg, settings, universe)
    pool <- if (length(bl$base)) bl$base
            else if ((nrow(evals) + w) %% every == 0L) bl$extra else character()
    if (length(pool)) {
      i <- ((w - 1L) %% length(pool)) + 1L
      json <- evals$config_json[match(pool[i], evals$config_hash)]
      return(list(cfg = config_from_json(json), source = "replicate"))
    }
  }

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
