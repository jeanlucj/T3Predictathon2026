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
       affects = function(cfg) identical(cfg$pheno_prep.ge_weighting, "env_gaussian")),
  list(build = "0.8.7",
       what  = "best_single_project requires focal coverage, as .best_panel does",
       affects = function(cfg) identical(cfg$geno_select.method, "best_single_project")),
  # Panel choice, so it reaches every config that picks ONE panel out of several: the
  # single-panel kernels, plus best_single_project, which applies the same rule to projects.
  list(build = "0.8.8",
       what  = paste("panel selection maximises the harmonic mean of focal and training",
                     "coverage, not their sum"),
       affects = function(cfg)
         isTRUE(cfg$kernel.method %in% c("vanRaden_single", "rkhs_gaussian")) ||
         identical(cfg$geno_select.method, "best_single_project"))
  # 0.9.0 deliberately has NO entry. It rewrote SELECTION -- the contender rule and its scale,
  # the parent pool, the reference-population sd_config and its pin, run identity, surrogate
  # case weights -- and none of that touches what a configuration computes on a trial. The one
  # 0.9.0 change that can: min_test_acc 5 -> 10, which is a floor in .best_panel_index(). An
  # entry predicated like 0.8.8's would retire 83% of the store to repair a case that bites only
  # when the winning panel covered 5-9 focal accessions AND another panel clears both floors --
  # and where it does not bite, the no-panel-clears fallback returns the same panel anyway. The
  # trade was judged not worth it; docs/SEARCH_ALGORITHM.md section 13 records the reasoning.
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

# The Fisher-z transform and the inverse-variance weights every score-level estimate uses,
# attached as .n_test / .usable / .z / .w.
#
# A score is a CORRELATION, so its precision depends on how many ACCESSIONS it was computed
# over -- n_test ranges from 10 to well over 1000. A plain mean gave those equal weight, so one
# lucky small trial could carry a configuration into get_elites. Weight on the Fisher-z scale
# instead, where the sampling variance is 1/(n-3) whatever the true correlation -- but CAPPED by
# settings$weight_cap, because sampling error is not the only error (see .eval_weight()).
#
# Factored out because write_report() has to describe the fit -- "the noise of one evaluation at
# a typical accession count" is sd_resid rescaled by these weights -- and restating the rule
# there would let the description drift from the fit it describes.
# settings$weight_cap, the ceiling on one evaluation's weight. No default here: a local one is
# how feasibility and usability came to disagree (LESSONS #34), and this number now governs two
# estimators, so it can only come from the setting.
.weight_cap <- function(settings) {
  x <- suppressWarnings(as.numeric(settings$weight_cap %||% NA_real_))
  if (!length(x) || is.na(x) || x <= 0) Inf else x
}

# The weight ONE evaluation carries, from its test-set size. Two independent errors corrupt a
# score and they behave differently as the trial grows, so the precisions add in series:
#
#   1/w = 1/(n_test - 3)  +  1/cap
#         \___________/      \___/
#         sampling error      config x trial interaction: does NOT shrink with the trial, so it
#         on the z scale      caps the weight. cap is the weight an unbounded trial would get.
#
# `cap = Inf` is the textbook Fisher weight n_test - 3, which over-weights large trials here:
# fitted unit-weight sd_resid should be 1 if n_test - 3 were the whole variance and it reaches
# 4.7 (dev/calibrate_weights.R, docs/LESSONS.md #36). NA n_test predates the column -- its score
# is real, we just cannot weight it -- so it takes the minimum weight rather than being dropped.
.eval_weight <- function(n_test, cap) {
  raw <- pmax(ifelse(is.na(n_test), 1, n_test - 3), 1)
  if (!is.finite(cap)) raw else 1 / (1 / raw + 1 / cap)
}

.fisher_weights <- function(evals, min_test, cap) {
  nt <- suppressWarnings(as.numeric(evals$n_test %||% NA_real_))
  dplyr::mutate(evals,
    .n_test = ifelse(is.na(nt), 0, nt),
    # Drop scores too small to inform anything. `min_test` is settings$min_test_acc, the same
    # floor score_predictions() refuses to produce a row below, so in a store written by this
    # build nothing is dropped here -- the test still matters for ARCHIVED rows, written when
    # feasibility stopped at 5 while this stopped at 10. NA n_test predates the column, so keep
    # it (its score is real, we just cannot weight it) with the minimum weight.
    .usable = is.finite(score) & (is.na(nt) | .n_test >= min_test),
    .z      = .to_z(score),
    .w      = .eval_weight(nt, cap))
}

# The z scale both estimators work on. Clamped because atanh(+-1) is Inf and a short test set
# can produce a perfect correlation; tanh() is the inverse, for reporting.
.to_z <- function(score) atanh(pmin(pmax(score, -0.999), 0.999))

# Mean score per distinct configuration (over the trials it was run on, within the
# already domain/scheme-filtered slice), with the parsed config attached. Failed
# evaluations (NA score) count toward n but not the mean; a config that only ever
# failed gets mean_score = NA.
# Minimum rows for a surrogate fit. Deliberately NOT settings$n_random_init: that is a
# search-policy budget ("explore at random until this many configs are scored"), this is a
# statistical floor on a regression fit, and tying them meant raising one raised the other.
.surrogate_min_obs <- function(settings) {
  x <- suppressWarnings(as.integer(settings$surrogate_min_obs %||% NA_integer_))
  if (!length(x) || is.na(x) || x < 1L) 25L else x
}

aggregate_scores <- function(evals, settings = NULL, min_test = NULL, weight_cap = NULL,
                             adjust_trial = TRUE) {
  # The floor comes from settings$min_test_acc, never from a default here: the same number gates
  # feasibility in score_predictions(), and a local default is how the two came to differ.
  # `min_test` is accepted directly so a caller probing the weighting (the tests) need not build
  # a settings object.
  if (is.null(min_test)) {
    if (is.null(settings))
      stop("aggregate_scores() needs `settings` (for min_test_acc), or an explicit min_test",
           call. = FALSE)
    min_test <- .min_test(settings)
  }
  if (is.null(weight_cap)) {
    if (is.null(settings))
      stop("aggregate_scores() needs `settings` (for weight_cap), or an explicit weight_cap",
           call. = FALSE)
    weight_cap <- .weight_cap(settings)
  }
  if (!nrow(evals)) {
    return(tibble::tibble(config_hash = character(), mean_score = numeric(),
                          pooled = numeric(), unweighted = numeric(), se = numeric(),
                          z_mean = numeric(), n = integer(), n_ok = integer(),
                          w_tot = numeric(), config_json = character()))
  }
  # `unweighted` is kept alongside `pooled` so the two can be compared on a real run.
  evals <- .fisher_weights(evals, min_test, weight_cap)
  out <- evals |>
    dplyr::group_by(config_hash) |>
    dplyr::summarise(
      pooled      = if (any(.usable)) tanh(sum(.w[.usable] * .z[.usable]) / sum(.w[.usable]))
                    else NA_real_,
      unweighted  = mean(score[is.finite(score)]),
      n           = dplyr::n(),
      n_ok        = sum(.usable),
      # Total weight behind the config, so a caller fitting on `mean_score` (the pooled
      # surrogate) can weight it the way the pooling did.
      w_tot       = sum(.w[.usable]),
      config_json = dplyr::first(config_json),
      .groups = "drop")
  out$mean_score <- out$pooled
  out$se <- NA_real_                    # only the random-effects fit supplies one
  # The z-scale companion to mean_score, supplied only by the random-effects fit. `.in_play()`
  # needs an estimate on the same scale as `se`; it falls back to atanh(mean_score) when this is
  # absent, which is exact for a blup row and the best available otherwise.
  out$z_mean <- NA_real_
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
    # sd_config from the reference population rather than from everything the search has
    # evaluated, when that can be estimated at all -- see .reference_sd_config(). NA keeps the
    # fit's own estimate, which is also what "fitted" asks for explicitly.
    pin <- if (identical(as.character((settings %||% list())$sd_config_source %||%
                                        "reference")[1], "reference"))
             .reference_sd_config(evals, settings) else NA_real_
    adj <- tryCatch(
      withCallingHandlers(
        .blup_scores(evals, sd_config_pin = if (is.finite(pin)) pin else NULL),
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
      out$z_mean[ok] <- adj$z_mean[i][ok]
      attr(out, "estimator")  <- "blup"
      attr(out, "var_comps")  <- attr(adj, "var_comps")
      attr(out, "within_config") <- attr(adj, "within_config")
      attr(out, "estimator_note") <- note        # NULL when the fit ran clean
      attr(out, "pin_route") <- attr(adj, "pin_route")
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
  if (!nrow(d)) return("no usable rows (a finite score and n_test >= min_test_acc)")
  if (dplyr::n_distinct(d$trial_id) < 2L) return("only 1 distinct trial")
  # Identifiability: with one config per trial the trial effect and the residual are the same
  # term. lme4 errors on this ("number of levels ... must be < number of observations"), so
  # check rather than rely on the error.
  if (max(tapply(d$config_hash, d$trial_id, function(x) length(unique(x)))) < 2L)
    return("no trial has 2 configs -- the trial effect is not separable from the residual")
  if (dplyr::n_distinct(d$config_hash) < 2L) return("only 1 distinct config")
  NULL
}

# ---- the reference population -------------------------------------------------------------
# sd_config is the SD of the configuration effects IN THE POPULATION CURRENTLY IN THE FIT, and
# the search narrows that population on purpose: late in a run most evaluated configurations are
# recombinations of a few parents, so the number measures the spread of the neighbourhood the
# search settled into rather than of the configuration space. It then does three jobs it is
# wrong for -- it is the shrinkage prior, it sets the unit the contender test is denominated in,
# and it feeds the homogenization that lowered it (docs/SEARCH_ALGORITHM.md section 11).
#
# So estimate it instead from configurations drawn WITHOUT reference to the surrogate.
# `evals.source` records which phase of choose_config() proposed each configuration, which makes
# those rows nameable: the Phase 3 random draws, the Phase 3b persistent draws and the Phase 5
# random fallback. Seeds are excluded -- five hand-picked submissions are not a sample of the
# space, and their spread would inflate the estimate. `replicate` rows are excluded as a SOURCE
# but not as EVIDENCE: a reference configuration keeps every row it ever earned, replication
# rows included.
#
# "random_init" and "random_reference" mean the same thing statistically and are kept apart only
# so the log shows which phase spent the step: the first is the one-time n_random_init budget,
# the second the persistent reference_every rate that keeps the population growing.
.REFERENCE_SOURCES <- c("random_init", "random_reference", "random_fallback")

# Hashes of the reference population, read from the FIRST evaluation of each configuration --
# later rows for the same configuration arrive through the replication phase and would otherwise
# be read as search-chosen.
#
# `tagged_only` decides what happens on rows written before `source` was stored. The two callers
# want different answers, and the difference is between a retrospective estimate and a forward
# spend:
#   FALSE  .reference_sd_config(): infer the population from arrival order. The alternative is
#          no estimate at all on every row written before 2026-11-03, and the inference is exact
#          for a store written by one run (Phase 3 completes before Phase 4 begins).
#   TRUE   .replication_backlog(): require a recorded tag. Spending an evaluation is an action,
#          and it should rest on a fact rather than an inference -- and the inference fails in
#          the direction that matters here, because on a small or untagged store "the first
#          n_random_init configurations" is most of the store, which would route nearly
#          everything into the reference-replication tier.
.reference_hashes <- function(evals, settings = NULL, tagged_only = FALSE) {
  if (!nrow(evals)) return(character())
  st <- settings %||% list()
  e  <- evals[order(as.character(evals$ts %||% rep("", nrow(evals))),
                    suppressWarnings(as.numeric(evals$id %||% rep(NA, nrow(evals))))), ,
              drop = FALSE]
  first <- e[!duplicated(e$config_hash), , drop = FALSE]
  src   <- as.character(first$source %||% rep(NA_character_, nrow(first)))
  keep  <- !is.na(src) & src %in% .REFERENCE_SOURCES
  if (any(keep) || isTRUE(tagged_only)) return(first$config_hash[keep])
  # `source` has only been STORED since 2026-11-03 (docs/LESSONS.md #38) although it was logged
  # from the first build, so for older rows the phase has to be inferred from arrival order:
  # Phase 3 runs to completion before Phase 4 starts, so the first n_random_init non-seed
  # configurations of a run ARE the random draws. Exact for a store written by one run; a store
  # holding several runs over one domain gets the earliest run's draws, which is a sample of the
  # same space.
  n <- suppressWarnings(as.integer(st$n_random_init %||% 25L))
  if (!length(n) || is.na(n) || n < 1L) n <- 25L
  seeds <- tryCatch(vapply(seed_configs(st$optimize_scheme %||% "CV0"),
                           config_hash, character(1)), error = function(e) character())
  utils::head(setdiff(first$config_hash, seeds), n)
}

# sd_config estimated on the reference population alone. NA whenever that cannot be done, and
# the caller then keeps the whole-population estimate. Four ways it declines:
#
#   too few reference configurations         settings$reference_min_configs
#   a design the fit cannot support          .blup_skip_reason(), e.g. one trial
#   a boundary estimate of 0                 carries no information, and pinning it would shrink
#                                            every BLUP exactly onto the grand mean -- strictly
#                                            worse than the fitted value it replaced
#   AN ESTIMATE THE POPULATION CANNOT        see below
#   SUPPORT
#
# The last one is not optional. This number becomes a FIXED prior, with no allowance anywhere
# for its own uncertainty, and it is fitted on as few as 20-odd configurations with one or two
# rows each -- so it can be large, confident and meaningless. Measured on the real runs by
# refitting on two thirds of the reference configurations: three runs are stable (`9a2664962c01`
# 0.0246, IQR 0.0216-0.0278; `f92e28ec6415` 0.0251, IQR 0.0240-0.0308; `34e15bd03838` 0.0560,
# IQR 0.0517-0.0630, none collapsing) while `e7517df121d5` returns 0.0517 from 37 usable rows
# over 24 configurations -- 1.5 rows each -- with an IQR of 0.0000-0.0924 and a collapse to the
# boundary in 9 of 25 subsamples. Pinning that resolved its field to 2 configurations out of
# 300, which the data cannot possibly support.
#
# So the estimate has to survive resampling. A row-count floor would not do: `80cfbd7c27e9` has
# 6.8 usable rows per reference configuration and is still borderline, because its 25
# configurations are spread over 24 trials and few of them share one. Resampling measures the
# thing that actually matters -- whether the answer depends on which configurations happened to
# be drawn -- rather than a proxy for it.
#
# Expects `evals` already carrying .usable/.z/.w from .fisher_weights(), as aggregate_scores()
# has them.
.reference_sd_config <- function(evals, settings) {
  st <- settings %||% list()
  m  <- suppressWarnings(as.integer(st$reference_min_configs %||% 20L))
  if (!length(m) || is.na(m) || m < 2L) m <- 20L
  ref <- .reference_hashes(evals, st)
  if (length(ref) < m) return(NA_real_)
  d <- evals[evals$config_hash %in% ref, , drop = FALSE]
  if (!nrow(d) || dplyr::n_distinct(d$config_hash[d$.usable]) < m) return(NA_real_)
  .fit_sd <- function(rows) {
    f <- tryCatch(suppressMessages(suppressWarnings(.blup_scores(rows))),
                  error = function(e) NULL)
    if (is.null(f)) return(NA_real_)
    v <- attr(f, "var_comps")
    if (!is.null(v) && "sd_config" %in% names(v)) as.numeric(v[["sd_config"]]) else NA_real_
  }
  sdc <- .fit_sd(d)
  if (!is.finite(sdc) || sdc <= 0) return(NA_real_)

  # The stability gate. Each draw refits on two thirds of the reference CONFIGURATIONS (not
  # rows: the question is whether the answer survives a different sample of the space). Rejected
  # when more than MAX_ZERO of the draws collapse to the boundary, which says the non-zero
  # estimate was an artefact of which configurations were in it. 0 draws disables the gate.
  #
  # Thresholds are statistical judgements rather than tuning knobs, so they live here. Measured
  # collapse rates on the seven real runs are 0-24%, so a quarter refuses only a reference
  # population that is genuinely too thin -- it is a guard, not a discriminator between those
  # runs (SEARCH_ALGORITHM.md section 11). The cost is a handful of extra fits on a few dozen
  # rows each, against an evaluation that takes tens of minutes.
  #
  # Two properties the gate needs and does not get for free:
  #   * K_FLOOR, not `m`. The subsample has to be able to drop enough configurations to probe
  #     anything. With 24 reference configurations and m = 20 it could drop 4, and the collapse
  #     rate it measured was not the one the decision was calibrated on.
  #   * A seed from the reference set itself, so the verdict is a function of the EVIDENCE. A
  #     re-rolled random decision flips the prior between iterations on unchanged data, which is
  #     the thrashing section 7 exists to stop.
  MAX_ZERO <- 0.25
  K_FLOOR  <- 12L
  B <- suppressWarnings(as.integer(st$reference_stability_draws %||% 25L))
  if (!length(B) || is.na(B) || B < 1L) return(sdc)
  cfgs <- sort(unique(d$config_hash[d$.usable]))
  k    <- max(K_FLOOR, ceiling(length(cfgs) * 2 / 3))
  if (k >= length(cfgs)) return(sdc)        # too few to resample from; the gate cannot run
  old_seed <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old_seed)) suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
          else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(strtoi(substr(rlang::hash(cfgs), 1, 7), 16L))
  draws <- vapply(seq_len(B), function(i)
    .fit_sd(d[d$config_hash %in% sample(cfgs, k), , drop = FALSE]), numeric(1))
  ok <- draws[is.finite(draws)]
  if (!length(ok) || mean(ok <= 0) > MAX_ZERO) NA_real_ else sdc
}

# Config BLUPs from a two-way random-effects fit, plus the variance components. NULL when the
# design or the packages cannot support it -- the caller then keeps its pooled estimate.
#
# `sd_config_pin` replaces the fit's own config variance with a value from elsewhere (the
# reference population above) and re-shrinks on it; NULL keeps lmer's estimate.
.blup_scores <- function(evals, sd_config_pin = NULL) {
  if (!is.null(.blup_skip_reason(evals))) return(NULL)
  d <- evals[evals$.usable, , drop = FALSE]
  m <- suppressMessages(lme4::lmer(.z ~ 1 + (1 | trial_id) + (1 | config_hash),
                                   data = d, weights = .w))
  # condVar = TRUE attaches each level's posterior variance, which is what the contender rule
  # needs; it is on the z scale, as `blup` is before the tanh below.
  re <- lme4::ranef(m, condVar = TRUE)$config_hash
  # `se` is the POSTERIOR SD of the configuration's random effect, not the standard error of a
  # mean. Its shape, APPROXIMATELY -- exact only for a balanced design with the crossed trial
  # term and the fixed effect ignored, so read it as the mechanism, not as a formula to check
  # digits against (it agrees closely when sd_config is small against sd_resid/sqrt(n), and runs
  # a few percent low otherwise):
  #
  #   se_blup(n) ~ sqrt( 1 / ( 1/sd_config^2 + n/sd_resid^2 ) )
  #
  # -- prior precision plus data precision. The fixed-effect quantity, which is what a reader
  # means by "how well do we know this configuration's accuracy", is the different
  #
  #   se_fixed(n) = sd_resid / sqrt(n)                (trial variance is already blocked out)
  #
  # They diverge sharply at this project's replication, because the prior dominates until n is
  # large against (sd_resid/sd_config)^2. With sd_config 0.015 and sd_resid 0.068 per
  # evaluation, the weight the BLUP puts on the DATA is 0.09 at n = 2 and only 0.54 at n = 24 --
  # so at the median configuration the reported se is ~90% prior. Shrinkage is right for
  # RANKING (it is what stops one lucky evaluation topping the list -- LESSONS #19) and
  # .contenders() should keep using it; write_report() shows both and labels which is which.
  #
  # Note also that the model has no config:trial term, so sd_resid confounds config x trial
  # interaction with measurement noise. The two are not separately identifiable without
  # replicating a (config, trial) cell, which the design avoids -- duplicate_cells counts those
  # as a defect. `within_config` below is the closest available window on that interaction.
  se <- sqrt(as.numeric(attr(re, "postVar")))
  v  <- as.data.frame(lme4::VarCorr(m))
  sd_of <- function(g) { x <- v$sdcor[v$grp == g]; if (length(x)) x[1] else NA_real_ }

  # Per-configuration spread across the trials it ran on, NOT pooled. `sd_adj` removes the
  # fitted trial effect: z_ij = mu + c_i + t_j + e_ij, so z - t_j leaves the configuration's own
  # scatter, which is what makes two configurations comparable when they cover different trial
  # subsets. Back-transformed so both columns are on the correlation scale.
  tre   <- lme4::ranef(m)$trial_id
  t_eff <- stats::setNames(tre[, 1], rownames(tre))
  tj    <- unname(t_eff[as.character(d$trial_id)]); tj[is.na(tj)] <- 0
  wc <- tibble::tibble(config_hash = d$config_hash, score = d$score,
                       adj = tanh(d$.z - tj)) |>
    dplyr::group_by(config_hash) |>
    dplyr::summarise(n_eval = dplyr::n(),
                     sd_raw = stats::sd(score),
                     sd_adj = stats::sd(adj), .groups = "drop")

  mu  <- as.numeric(lme4::fixef(m)[1])
  # z_mean is the configuration's estimate ON THE Z SCALE, which is the scale `se` is on. Both
  # are carried so that any comparison of an estimate with its own uncertainty can be made on
  # one scale and transformed afterwards, instead of adding a z-scale se to a tanh'd mean.
  # .in_play() is the consumer. docs/LESSONS.md #43.
  out <- tibble::tibble(config_hash = rownames(re),
                        z_mean = mu + re[, 1],
                        blup   = tanh(mu + re[, 1]),
                        se     = se)
  # sd_config here is always the FIT'S OWN estimate, so the number never changes meaning; the
  # extra entries appear only when a pin was applied, which keeps `var_comps` exactly the
  # three-element vector every existing reader expects in the unpinned case.
  vc  <- c(sd_trial = sd_of("trial_id"), sd_config = sd_of("config_hash"),
           sd_resid = sd_of("Residual"))
  s2r <- sd_of("Residual")^2
  sdf <- sd_of("config_hash")
  route <- NULL                         # which pin route ran, NULL when no pin was applied

  # Re-shrink at a PINNED config variance (the reference-population estimate). Both routes keep
  # the fit's own intercept, trial effects and sd_resid -- removing the trial effect is what the
  # crossed fit is for and that part is not in question -- and change only the config effects.
  #
  # Write prec_i(s) = 1/s^2 + W_i/sd_resid^2, the prior precision plus the data precision behind
  # configuration i, and lambda_i(s) = (W_i/sd_resid^2)/prec_i(s), the weight the BLUP puts on
  # the data. Since lambda_i(pin)/lambda_i(fit) = prec_i(fit)/prec_i(pin), ONE ratio
  #
  #   R_i = prec_i(fit) / prec_i(pin)
  #
  # governs both quantities: c_i scales by R_i and se_i by sqrt(R_i).
  #
  # ROUTE A (sd_config fitted > 0): rescale lme4's OWN output by R_i. This matters -- it is why
  # the switch is no longer confounded. The analytic one-way form ignores the uncertainty in the
  # trial effects and runs a few percent low against lme4's postVar, so computing it from
  # scratch changed the ESTIMATOR as well as the prior and the two could not be compared.
  # Rescaling inherits everything lme4 did and is exactly the identity at pin == fit.
  #
  # ROUTE B (the fitted sd_config is at or near the boundary): there is nothing to rescale.
  # Every c_i is 0 or a rounding artefact of 0, so rescaling multiplies no signal by a huge
  # factor -- and the test for this cannot be `sdf > 0`. A REML fit at the boundary returns
  # values like 6e-17 rather than an exact zero, which passes `> 0` while leaving the deviations
  # with no significant figure to carry: se recovers anyway (it is proportional either way) but
  # every estimate stays tied at the grand mean, so the rule still separates nothing and the
  # failure is silent. The condition is therefore on R itself. Past R_MAX the fitted prior has
  # shrunk the deviations by more than three orders of magnitude relative to the pin, so the
  # effects are recomputed from the trial-adjusted deviations instead -- exact, better
  # conditioned, and nearly continuous with route A at the threshold, since the two differ only
  # by the trial-effect-uncertainty correction route A inherits from lme4.
  #
  # This is the case the whole feature exists for: a boundary fit is exactly where the contender
  # rule has no `se` to read at all. The cost is that route B's se IS the analytic form, a few
  # percent low, because at the boundary lme4 offers nothing to calibrate against.
  if (!is.null(sd_config_pin) && is.finite(sd_config_pin) && sd_config_pin > 0 &&
      is.finite(s2r) && s2r > 0) {
    W <- tibble::tibble(config_hash = d$config_hash, w = d$.w, r = d$.z - mu - tj) |>
      dplyr::group_by(config_hash) |>
      dplyr::summarise(W = sum(w), dev = sum(w * r) / sum(w), .groups = "drop")
    R_MAX <- 1e3
    i    <- match(out$config_hash, W$config_hash)
    dprec <- W$W[i] / s2r                       # data precision behind each configuration
    ppin  <- 1 / sd_config_pin^2 + dprec
    R     <- if (is.finite(sdf) && sdf > 0) (1 / sdf^2 + dprec) / ppin else Inf
    # One route for every configuration in a fit, so `se` is on one estimator throughout: the
    # worst-conditioned configuration decides.
    route <- if (any(is.finite(R)) && max(R, na.rm = TRUE) <= R_MAX) "rescaled" else "recomputed"
    if (identical(route, "rescaled")) {
      ok <- !is.na(i) & is.finite(R)
      out$z_mean[ok] <- (mu + (out$z_mean - mu) * R)[ok]
      out$se[ok]     <- (out$se * sqrt(R))[ok]
    } else {
      ok <- !is.na(i) & is.finite(ppin) & ppin > 0
      out$z_mean[ok] <- (mu + (dprec / ppin) * W$dev[i])[ok]
      out$se[ok]     <- (1 / sqrt(ppin))[ok]
    }
    out$blup <- tanh(out$z_mean)
    vc["sd_config_ref"]  <- sd_config_pin
    vc["sd_config_used"] <- sd_config_pin
  }
  structure(out, var_comps = vc, within_config = wc, pin_route = route)
}

# Top-k configurations by mean score (only those with at least one success). A ranking of POINT
# estimates, asserting nothing about one configuration being better than another. It is no longer
# the parent pool -- get_parents() is (LESSONS #42) -- and is kept as the mean-ranked top-k that
# the churn and turnover measures are defined on (SEARCH_ALGORITHM.md section 9) and that the dev
# scripts compare the in-play set against. get_parents() reproduces this ordering inline for its
# floor rather than calling in, because it needs hashes and this returns parsed configurations.
get_elites <- function(agg, k = 8) {
  agg |>
    dplyr::filter(is.finite(mean_score)) |>
    dplyr::arrange(dplyr::desc(mean_score)) |>
    head(k) |>
    dplyr::pull(config_json) |>
    lapply(config_from_json)
}

# Parents for candidate generation: the configurations still in play -- the SAME set that earns
# contender replication. One definition of "potentially good" now drives both decisions, where
# parenthood used to be the top n_elites by point estimate and replication the in-play set, so a
# configuration could be judged good enough to deepen and not good enough to breed from, or the
# reverse, with no reason for the asymmetry.
#
# Capped at settings$n_parents and floored at settings$n_elites:
#   CAP    the in-play set is EVERYTHING when the estimator cannot separate configurations
#          (sd_config at the boundary), and unbounded parenthood is then recombination among
#          every configuration ever scored. Ordered by UCB, so the cap keeps the most promising.
#   FLOOR  the in-play set is EMPTY when the random-effects fit did not run, and an empty parent
#          pool turns propose_candidates() into pure random search without saying so. Topped up
#          from the mean ranking, which needs no `se`.
get_parents <- function(agg, settings) {
  pos <- function(x, default) {
    v <- suppressWarnings(as.integer(x %||% NA_integer_))
    if (!length(v) || is.na(v) || v < 1L) default else v
  }
  spec   <- .contender_spec(settings)
  n_min  <- pos((settings %||% list())$n_elites, 8L)
  # The cap cannot sit below the floor: the top-up would otherwise re-fill past it from the mean
  # ranking, which quietly turns n_parents into a no-op rather than a cap.
  hashes <- utils::head(.in_play(agg, spec$z, spec$delta, spec$rule),
                        max(n_min, pos((settings %||% list())$n_parents, 32L)))
  if (length(hashes) < n_min) {
    fin    <- agg[is.finite(agg$mean_score), , drop = FALSE]
    hashes <- utils::head(unique(c(hashes, fin$config_hash[order(fin$mean_score,
                                                                 decreasing = TRUE)])), n_min)
  }
  lapply(agg$config_json[match(hashes, agg$config_hash)], config_from_json)
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
#   $extra  earning one more each time this is called, so the target rises for exactly the
#           configs still in question. Rationed by the caller. Two kinds qualify:
#             - a CONTENDER at or above the floor (§7); and
#             - a REFERENCE configuration below settings$reference_replication. A variance
#               component needs replication within level, and at config_replication = 2 the
#               reference fit rests on two rows per configuration -- thin enough on one real run
#               to return a boundary 0 and be declined. Served here rather than in `base` on
#               purpose: `base` is unrationed and preempts everything, so a reference floor
#               there could starve the search.
# Counts EVALUATIONS, not distinct trials, so it always drains. `universe` (eligible trial ids,
# NULL in simulate mode) removes a config already run on every eligible trial.
.replication_backlog <- function(evals, agg, settings, universe = NULL) {
  none <- list(base = character(), extra = character())
  cr <- suppressWarnings(as.integer(settings$config_replication %||% 1L))
  if (!is.finite(cr) || !nrow(evals)) return(none)
  sp   <- .contender_spec(settings)
  cand <- .contenders(agg, sp$z, k = 8L, delta = sp$delta, rule = sp$rule)
  rr   <- suppressWarnings(as.integer(settings$reference_replication %||% 0L))
  if (!length(rr) || is.na(rr)) rr <- 0L
  ref  <- if (rr > cr) .reference_hashes(evals, settings, tagged_only = TRUE) else character()
  per <- evals |>
    dplyr::group_by(config_hash) |>
    dplyr::summarise(n_eval = dplyr::n(), n_trial = dplyr::n_distinct(trial_id), .groups = "drop")
  if (length(universe)) per <- dplyr::filter(per, n_trial < length(universe))
  per <- dplyr::arrange(per, config_hash)
  list(base  = per$config_hash[per$n_eval < cr],
       extra = per$config_hash[per$n_eval >= cr &
                               (per$config_hash %in% cand |
                                (per$n_eval < rr & per$config_hash %in% ref))])
}

# The contender test's three parameters, read once so that every call site agrees on them.
# All three are search POLICY: they change the trajectory, not the meaning of any stored row.
# Under the default settings$run_identity = "target" none of them is run-defining, so the whole
# selection rule can be retuned on a live run; under "experiment" contender_z is, and changing
# it rotates run_id. docs/SEARCH_ALGORITHM.md section 13.
.contender_spec <- function(settings) {
  num <- function(x, default) {
    v <- suppressWarnings(as.numeric(x %||% NA_real_))
    if (!length(v) || is.na(v) || v < 0) default else v
  }
  rule <- as.character((settings %||% list())$contender_rule %||% "racing")[1]
  if (!rule %in% c("racing", "point")) rule <- "racing"
  list(z     = num((settings %||% list())$contender_z, 0.25),
       delta = num((settings %||% list())$resolution_delta, 0),
       rule  = rule)
}

# Configurations that cannot yet be ruled out as the best, ALL of them, best first. Two parts:
#
#   THE CLAIM   ucb_i = mean_i + z*se_i, an optimistic view of configuration i.
#   THE BAR     "racing"  max_j(mean_j - z*se_j): the best value that can be DEFENDED about any
#                         configuration. Optimistic against pessimistic -- the comparison a
#                         racing algorithm makes, which eliminates a candidate only once its
#                         interval stops overlapping the leader's. Both sides are charged the
#                         same z, where the claim alone being inflated charges a thinly
#                         replicated candidate for its uncertainty and gives the leader a free
#                         pass on its own.
#               "point"   max_j(mean_j): optimistic against a point estimate. The bar is then
#                         the maximum of a set of noisy estimates -- the most upward-biased
#                         order statistic in the set -- so one lucky thinly replicated
#                         configuration sets it for everyone and moves it at every refit.
#                         Measured, the in-play set moves 2.4x to 4.6x further per checkpoint
#                         under this bar than under racing at the same z; kept for comparison
#                         (dev/sweep_contender_z.R).
#
# z does DOUBLE DUTY under racing: it inflates the claim AND lowers the bar, so the in-play set
# grows about twice as fast in z as under "point", and the two rules give comparable set sizes
# at z_racing ~ z_point / 2. At z = 0 the two are identical.
#
# `delta` is the RESOLUTION FLOOR: the difference in predictive accuracy actually worth
# resolving, chosen from OUTSIDE the run instead of measured in it. A configuration keeps its
# place when it misses the bar by less than delta, so membership flips only on evidence that
# moves further than the resolution anyone cares about -- a deadband. It is the only term in
# this test that does not shrink when the population homogenizes, which is the defect
# docs/SEARCH_ALGORITHM.md section 11 describes. Read it against sd_config: as delta approaches
# sd_config everything is in play, and that is the honest reading rather than a failure -- no
# configuration then differs from another by an amount worth having.
#
# ON SCALE. `mean_score` is on the correlation scale and `se` on the z scale, so the claim and
# the bar mix the two; tanh'(x) runs 0.99 to 0.79 over the range these runs occupy, so se is
# overstated by up to a quarter at the top of it. Left as it is deliberately: every z in
# dev/sweep_contender_z.R is calibrated against this comparison, so correcting the scale without
# re-sweeping would silently move the operating point. docs/LESSONS.md #40.
#
# Empty when the estimator supplied no standard errors (the pooled fallback, on a store too thin
# to fit) -- there are then no contenders at all, and get_parents() supplies the floor.
.in_play <- function(agg, z = 0.25, delta = 0, rule = "racing") {
  a <- agg[is.finite(agg$mean_score) & is.finite(agg$se %||% NA_real_), , drop = FALSE]
  if (!nrow(a)) return(character())
  # ONE SCALE. The estimate and its se are both on z, where the fit produced them, and the
  # interval endpoints are transformed to the correlation scale afterwards. tanh is monotone, so
  # this changes no in/out decision for z alone -- but `delta` is a difference in predictive
  # ACCURACY, which only means something on the correlation scale, and the comparison has to be
  # there for it to be added. `z_mean` comes from the random-effects fit; atanh(mean_score) is
  # the exact fallback for a caller that supplies only mean_score.
  zm <- if ("z_mean" %in% names(a)) as.numeric(a$z_mean) else rep(NA_real_, nrow(a))
  zm[!is.finite(zm)] <- .to_z(a$mean_score[!is.finite(zm)])
  ucb <- tanh(zm + z * a$se)
  bar <- if (identical(rule, "point")) max(tanh(zm)) else max(tanh(zm - z * a$se))
  # >=: a bound that exactly reaches the bar is not ruled out. The leader is in play under
  # either rule, since bar <= max(tanh(zm)) <= the ucb of whichever config holds that maximum.
  keep <- ucb + delta >= bar
  a$config_hash[keep][order(ucb[keep], decreasing = TRUE)]
}

# The first `k` of them, so a result of length 1 means the field has been decided at this spec.
.contenders <- function(agg, z = 0.25, k = 8L, delta = 0, rule = "racing") {
  utils::head(.in_play(agg, z, delta, rule), k)
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

  agg <- aggregate_scores(evals, settings)
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

  # Phase 3b: a PERSISTENT reference draw. n_random_init is a one-time gate on a cumulative
  # count, so once it is passed the reference population never grows again and the sd_config
  # estimate taken from it (.reference_sd_config) stops improving -- it stays frozen at whatever
  # the first n_random_init configurations could support. This rate keeps it growing: 1 step in
  # `reference_every` draws a fresh configuration however mature the run is.
  #
  # Placed after replication so it cannot preempt the config_replication floor, and before
  # acquisition so it is the acquisition budget it spends -- which is the honest accounting,
  # since its purpose is to buy an unbiased view of the space at the cost of exploitation.
  # Keyed on the store's row count plus the worker index, exactly as replicate_every is, so one
  # worker alternates over time and N workers split at any instant without lockstep.
  ref_every <- suppressWarnings(as.integer(settings$reference_every %||% 0L))
  if (length(ref_every) && !is.na(ref_every) && ref_every >= 2L &&
      (nrow(evals) + w) %% ref_every == 0L) {
    return(list(cfg = fresh_random(done_hashes), source = "random_reference"))
  }

  # Phase 4: surrogate-guided. Three ways to handle the trial effect, most structured first:
  #   MERF     one row per (config, trial), the trial effect a shrunken random effect fitted
  #            outside the forest. Costs one variance parameter rather than competing for tree
  #            structure, and predicts the trial-marginal value directly.
  #   BLOCKED  the same rows with trial_id as a forest feature, marginalised at prediction.
  #            Removes the trial variance from every config's estimate and keeps the
  #            replication weighting that collapsing to means throws away -- LESSONS #19.
  #   POOLED   one row per config mean, when the design cannot support a trial effect at all.
  scored  <- agg |> dplyr::filter(is.finite(mean_score))
  # The surrogate trains on the SAME rows the estimator ranks on, at the SAME weights: one
  # usability floor and one weight cap, both from settings, both applied by .fisher_weights().
  # It used to take every finite score at equal weight, which trained it on rows
  # aggregate_scores() discards and let a 900-accession correlation count no more than a
  # 12-accession one (LESSONS #36).
  rows    <- .fisher_weights(evals, .min_test(settings), .weight_cap(settings)) |>
               dplyr::filter(.usable)
  # Drop rows whose trial has not yet been replicated, keeping trial_id in the model for the
  # rest. (Measured: unreplicated trials do NOT in fact damage the fit -- this is the
  # conservative reading, and trial_replication drains the backlog promptly, so few rows go.)
  if (nrow(rows)) {
    n_cfg <- tapply(rows$config_hash, rows$trial_id, function(x) length(unique(x)))
    rows  <- rows[n_cfg[as.character(rows$trial_id)] >= 2L, , drop = FALSE]
  }
  # merf -> blocked -> pooled, each falling through when it cannot be fitted. settings$
  # surrogate_method names the most structured one to attempt; the trial effect needs at least
  # two trials either way.
  want    <- settings$surrogate_method %||% "merf"
  # How many rows a forest needs before it is worth fitting. This was settings$n_random_init
  # until 2026-11-03, which tied "how long to explore at random" to "how much data a regression
  # needs" -- two unrelated quantities that could not then be tuned apart.
  min_obs <- .surrogate_min_obs(settings)
  grouped <- nrow(rows) > 0 && dplyr::n_distinct(rows$trial_id) >= 2L
  mode    <- "pooled"
  surro   <- NULL
  # All three fit on the Z SCALE, as the estimator does. The trial effect that merf and blocked
  # remove is ADDITIVE, and additivity is incoherent on the correlation scale -- a +0.3 trial
  # shift applied to a configuration at 0.85 leaves the range. z is also where the weights are
  # valid, since 1/(n-3) is the sampling variance of z, not of r.
  if (grouped && identical(want, "merf")) {
    feats <- configs_to_features(lapply(rows$config_json, config_from_json))
    surro <- fit_surrogate_merf(feats, rows$.z, rows$trial_id, weights = rows$.w,
                                ntree = settings$ntree, min_obs = min_obs)
    if (!is.null(surro)) mode <- "merf"
  }
  if (is.null(surro) && grouped && want %in% c("merf", "blocked")) {
    feats <- configs_to_features(lapply(rows$config_json, config_from_json))
    feats$trial_id <- factor(rows$trial_id)
    surro <- fit_surrogate(feats, rows$.z, weights = rows$.w,
                           ntree = settings$ntree, min_obs = min_obs)
    if (!is.null(surro)) mode <- "blocked"
  }
  if (is.null(surro)) {
    feats <- configs_to_features(lapply(scored$config_json, config_from_json))
    surro <- fit_surrogate(feats, .to_z(scored$mean_score), weights = scored$w_tot,
                           ntree = settings$ntree, min_obs = min_obs)
    mode  <- "pooled"
  }
  if (is.null(surro)) return(list(cfg = fresh_random(done_hashes), source = "random_fallback"))

  elites  <- get_parents(agg, settings)
  cands   <- propose_candidates(agg, elites, done_hashes,
                                settings$n_cross, settings$n_mut, settings$n_rand)
  cfeats  <- configs_to_features(cands)
  # Score a candidate as its EXPECTED value over the observed trials -- the objective is
  # "best on a random trial", so trial_id is marginalised away rather than chosen.
  # merf already predicts the marginal -- its forest was fit to y - b, so the random effect is
  # zero at prediction. blocked has to average trial_id away explicitly.
  # One estimator, one scale, for both sides of the EI comparison: the candidates AND the
  # incumbent are the forest's own trial-marginal prediction.
  .marginal <- function(f) if (identical(mode, "blocked"))
                             predict_surrogate_marginal(surro, f, levels(feats$trial_id))
                           else predict_surrogate(surro, f)
  pr      <- .marginal(cfeats)
  # THE REFERENCE POINT IS A PREDICTION, NOT max(BLUP). A BLUP is shrunk toward the grand mean
  # in proportion to sd_config; a forest prediction is not shrunk at all. Comparing them put the
  # two sides of (mean_hat - best) on different estimation scales, and the gap grew as sd_config
  # fell: `best` sank, z = (mean_hat - best - xi)/sd rose, pnorm(z) -> 1, dnorm(z) -> 0, and
  # EI -> mean_hat - best - xi, a strictly increasing function of the predicted mean. The
  # sd*dnorm(z) term is the ONLY reason EI ever prefers a candidate the forest is unsure about,
  # so losing it makes acquisition plain greedy exploitation. The repair cannot run the other
  # way -- shrinking mean_hat needs the candidate's own replication, which is zero by
  # construction, so every candidate would collapse to the grand mean. LESSONS #37.
  #
  # Predicting at an evaluated configuration's own features is partly IN-SAMPLE, so `best` can
  # be optimistic; an out-of-bag reference would remove that. Measured against the unshrunk
  # pooled means, the forest smooths far harder than it memorises -- max(prediction) is 0.22 to
  # 0.52 where max(pooled z) is 0.70 to 0.93 -- so `best` lands well below the most extreme
  # thing observed, and an overstated `best` errs toward MORE weight on sd, which is the safe
  # direction for the failure being repaired. LESSONS #37 has the open question.
  ifeats  <- if (identical(mode, "pooled")) feats
             else configs_to_features(lapply(scored$config_json, config_from_json))
  best    <- max(.marginal(ifeats)$mean)
  # EI is NOT invariant under atanh -- it integrates the SIZE of the improvement, so on z a gain
  # at the top of the range counts for more than the same gain in the middle, which is the
  # intended reading. `ei_xi` is a z-scale margin.
  ei      <- expected_improvement(pr$mean, pr$sd, best, settings$ei_xi)

  pick <- which.max(ei)
  list(cfg = cands[[pick]],
       source = switch(mode, merf = "acquisition_merf", blocked = "acquisition_blocked",
                       "acquisition"),
       # pred_mean/pred_sd are on the z scale, as the fit is; pred_score back-transforms the
       # mean for anything reporting it beside a stored `score`.
       ei = ei[pick], pred_mean = pr$mean[pick], pred_sd = pr$sd[pick],
       pred_score = tanh(pr$mean[pick]))
}
