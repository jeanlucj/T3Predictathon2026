# evaluate.R
#
# Turns "run this configuration on a trial" into a number the optimizer can use:
# the Pearson correlation between predicted and observed per-accession grain-yield
# BLUEs on the held-out focal-trial accessions -- the Predictathon's own metric.
#
# Two modes, behind one interface (sample_trial / evaluate_config_on_trial):
#   * SIMULATE = TRUE  : a synthetic genomic-prediction world with a KNOWN
#     generative model in which some subtask choices are genuinely better than
#     others, and some are better only on certain kinds of trial (so crossover
#     and generalization matter). Runs offline in milliseconds; used by
#     tests/test_sim_loop.R to prove the optimizer actually learns before any
#     real compute is spent.
#   * SIMULATE = FALSE : the real pipeline (R/pipeline.R) on real T3 data
#     (R/data_access.R), scored against observed BLUEs.

library(tidyverse)

# ---------------------------------------------------------------------------
# Trial sampling. Returns a small descriptor list: id, plus covariates the
# scorer/pipeline use. In SIMULATE mode the covariates are drawn; in real mode
# they come from T3 metadata (see R/data_access.R::sample_real_trial).
# ---------------------------------------------------------------------------
sample_trial <- function(settings, conn = NULL) {
  if (settings$simulate) {
    # Fixed-trial mode: one constant descriptor, so a config's score stops depending on
    # which trial it happened to land on. This removes the LARGER of the two variance
    # sources (see settings$sim_fixed_trial). Mid-range values, so that no method branch
    # is degenerate: n_projects 3 keeps em_combine meaningful, env_var 0.5 keeps the G+E
    # and ge_weighting terms live.
    if (isTRUE(settings$sim_fixed_trial))
      return(list(id = "simtrial_fixed", n_projects = 3L, env_var = 0.5,
                  n_acc = 150L, heritability = 0.45))
    .sim_trial(paste0("simtrial_", sample.int(1e6, 1)))
  } else {
    sample_real_trial(settings, conn)               # defined in R/data_access.R
  }
}

# A simulated trial's attributes, DERIVED FROM ITS ID so the id actually identifies a trial.
# Needed because trial_replication revisits a trial by id: if the attributes were redrawn the
# "same" simulated trial would be a different trial each time, and neither the replication
# constraint nor the surrogate's trial_id feature could be exercised offline.
# The global RNG is saved and restored, so seeding here does not perturb the optimizer's own
# random stream (config sampling, candidate proposal).
.sim_trial <- function(id) {
  n <- suppressWarnings(as.integer(sub("^simtrial_", "", id)))
  if (!is.finite(n)) n <- 1L
  had <- exists(".Random.seed", envir = globalenv())
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(n)
  out <- list(
    id          = id,
    n_projects  = sample(1:4, 1),                 # how many genotyping projects cover it
    env_var     = stats::runif(1),                # environmental atypicality (0..1)
    n_acc       = sample(40:300, 1),              # number of accessions
    heritability = stats::runif(1, 0.2, 0.7)      # caps achievable accuracy
  )
  if (had) assign(".Random.seed", old, envir = globalenv())
  else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  out
}

# ---------------------------------------------------------------------------
# Evaluate one (config, trial, scheme). Returns list(score, n_test, status,
# reason, seconds). score is NA on failure; status/reason record why.
# ---------------------------------------------------------------------------
evaluate_config_on_trial <- function(cfg, trial, scheme, settings, conn = NULL) {
  t0 <- Sys.time()
  # Bracket this evaluation's memory peak (R/memory.R). Done before the pipeline runs and
  # read after it -- including after a failure, since an evaluation that dies of memory
  # pressure is exactly the one whose footprint we want on record.
  mem_reset()
  # Bound one evaluation's wall time (Inf = no bound; settings$max_eval_minutes explains why
  # that is the default). `transient = TRUE` restores the previous limit when this frame exits.
  # R checks the limit only at interpreter checkpoints, so a long call inside compiled code (a
  # BLAS crossproduct, a VCF parse) can overshoot it.
  cap_min <- settings$max_eval_minutes %||% Inf
  if (is.finite(cap_min) && cap_min > 0) {
    setTimeLimit(elapsed = cap_min * 60, transient = TRUE)
    on.exit(setTimeLimit(), add = TRUE)
  }
  # One error handler that branches on the condition class. (A single handler is
  # deliberate: re-raising a fatal with stop(e) from a multi-handler tryCatch gets
  # caught by that same tryCatch's sibling `error` handler, which would swallow
  # the halt; from a lone handler it propagates out as intended.)
  res <- tryCatch({
      if (settings$simulate) {
        .sim_evaluate(cfg, trial, scheme, settings)
      } else {
        preds_obs <- run_pipeline(cfg, trial, scheme, settings, conn)  # R/pipeline.R
        score_predictions(preds_obs$pred, preds_obs$obs)
      }
    },
    error = function(e) {
      # A run-level problem (bad settings, unimplemented method): re-raise so
      # run_optimizer.R halts instead of recording it as a per-config failure.
      if (inherits(e, "optimizer_fatal")) stop(e)
      # A demanding config this trial cannot satisfy: record it in the failure log
      # (status="infeasible", reason=code) and let the loop move on. EXPECTED.
      # A "suspect" infeasibility (funnel cliff -> data that should be visible is
      # not) is recorded under its own status so it can't hide among the real ones.
      if (inherits(e, "optimizer_infeasible"))
        return(list(score = NA_real_, n_test = 0L,
                    status = if (isTRUE(e$suspect)) "suspect" else "infeasible",
                    reason = e$code, detail = funnel_string(e$funnel)))
      # The wall-clock cap fired. Not a bug and not a data verdict -- its own status, so a
      # run that is merely too slow never hides among genuine crashes in the failure log.
      msg <- conditionMessage(e)
      if (grepl("reached elapsed time limit|reached CPU time limit", msg))
        return(list(score = NA_real_, n_test = 0L, status = "timeout",
                    reason = sprintf("exceeded max_eval_minutes = %g", cap_min)))
      # Anything unexpected: record it (so a stray bug does not kill a long run)
      # but tag it distinctly so it stands out in the log.
      list(score = NA_real_, n_test = 0L, status = "error", reason = msg)
    })
  res$detail <- res$detail %||% NA_character_
  res$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # peak_rss_mb is the honest cost (it includes BLAS/sommer allocations gc() cannot see);
  # peak_r_mb is kept because the ratio between them measures exactly that blind spot.
  res$peak_rss_mb <- mem_peak_rss_mb()
  res$peak_r_mb   <- mem_peak_mb()
  res$rss_mb      <- proc_rss_mb()
  res
}

# Score = correlation of predictions with observed BLUEs over shared accessions.
# Fails (NA) if fewer than 5 lines overlap or predictions/observations are
# constant -- matching the Predictathon's "task can fail" rule.
score_predictions <- function(pred, obs) {
  j <- dplyr::inner_join(
    tibble::tibble(germ = names(pred), pred = as.numeric(pred)),
    tibble::tibble(germ = names(obs),  obs  = as.numeric(obs)),
    by = "germ")
  m <- dplyr::filter(j, is.finite(pred), is.finite(obs))
  # Two unrelated failures used to share the `too_few_overlap` status, and telling them apart
  # after the fact was impossible: n_test is counted AFTER the finite filter, so both read as
  # "0 accessions". They mean opposite things --
  #   nothing JOINED          -> a coverage / name problem: predictions and observations are
  #                              about different accessions.
  #   joined but none FINITE  -> a modelling problem: the names lined up and the numbers are
  #                              NaN/Inf. Nothing about the data is missing.
  # Report the second separately, and say which side went non-finite.
  if (nrow(j) > 0 && nrow(m) == 0) {
    nap <- sum(!is.finite(j$pred)); nao <- sum(!is.finite(j$obs))
    return(list(score = NA_real_, n_test = 0L, status = "non_finite",
                reason = sprintf("%d joined, 0 finite (pred non-finite %d, obs non-finite %d)",
                                 nrow(j), nap, nao)))
  }
  if (nrow(m) < 5) return(list(score = NA_real_, n_test = nrow(m),
                               status = "too_few_overlap",
                               reason = sprintf("n<5 (%d joined, %d finite)", nrow(j), nrow(m))))
  if (stats::sd(m$pred) == 0 || stats::sd(m$obs) == 0) {
    return(list(score = NA_real_, n_test = nrow(m),
                status = "constant", reason = "no variance"))
  }
  list(score = stats::cor(m$pred, m$obs), n_test = nrow(m),
       status = "ok", reason = NA_character_)
}

# ===========================================================================
# SIMULATE world: a transparent stand-in for the real objective.
# .sim_true() is the deterministic expected predictive ability of a config on a
# trial; .sim_evaluate() adds single-trial sampling noise. The numbers encode
# plausible genomic-prediction wisdom AND deliberate interactions, so that (a)
# no single submitted seed is optimal, (b) recombining good blocks helps, and
# (c) some choices only pay off on certain trial types -- which is what makes
# "predict any random trial" a non-trivial target.
# ===========================================================================
.sim_true <- function(cfg, trial, scheme) {
  q <- 0

  # A. training-trial selection: targeted similarity beats both too-narrow and
  #    too-broad; k has an interior optimum around 20.
  q <- q + switch(cfg$train_select.method,
    top_k_similar     = 0.22 - 0.00035 * (cfg$train_select.k - 20)^2,
    # accession_overlap: the two-tier form (primary_only="no") beats primary-only.
    accession_overlap = if (identical(cfg$train_select.primary_only, "yes")) 0.10 else 0.17,
    same_program      = 0.14)

  # B. phenotype prep: shrinkage BLUP best; env weighting helps on atypical envs.
  q <- q + switch(cfg$pheno_prep.method,
    two_stage_blup = 0.14, blue_lm = 0.10, trial_center = 0.07, raw_mean = 0.03)
  if (identical(cfg$pheno_prep.ge_weighting, "env_gaussian")) {
    q <- q + 0.06 * trial$env_var - 0.02 * (1 - trial$env_var)  # pays off only on atypical envs
  }

  # C. genotyping data: one-hop with bridges is the sweet spot, and how strict the
  #    bridge requirement is has an interior optimum -- 1 admits panels hanging off a
  #    single shared line (a barely-identified cross-panel block), 5 turns away panels
  #    that carried real information.
  q <- q + switch(cfg$geno_select.method,
    focal_plus_onehop = 0.06 - 0.006 * (as.integer(cfg$geno_select.min_bridge) - 2)^2,
    best_single_project = 0.03, all_projects = 0.04)

  # D. kernel: EM-combine shines when many projects must be bridged; a single
  #    VanRaden GRM is best when one project covers the trial.
  q <- q + switch(cfg$kernel.method,
    em_combine      = 0.02 + 0.05 * (trial$n_projects - 1) / 3,
    vanRaden_single = 0.07 - 0.05 * (trial$n_projects - 1) / 3,
    rkhs_gaussian   = 0.04 - 0.01 * abs(log(cfg$kernel.rkhs_theta)))

  # E. model: multi-kernel G+E helps on atypical envs when E is on; RKHS only pays off
  #    with the RKHS kernel.
  q <- q + switch(cfg$model.method,
    gblup_sommer_GE = 0.04 + ifelse(identical(cfg$model.include_E, "yes"), 0.05 * trial$env_var, 0),
    gblup_rrblup    = 0.03,
    rkhs            = ifelse(cfg$kernel.method == "rkhs_gaussian", 0.06, 0.01))
  # ... and how lambda is chosen matters on its own: LOO is the robust all-rounder under
  # a misspecified model, REML is fine but trusts the model, and a fixed lambda is only
  # as good as the value (best near 1 here, worse by orders of magnitude either way).
  q <- q + switch(cfg$model.lambda_select,
    loo   = 0.03,
    reml  = 0.01,
    fixed = 0.02 - 0.012 * abs(log10(cfg$model.lambda_fixed)))

  # F. prediction post-processing: blending observed BLUE helps under CV0 (focal
  #    lines may recur in training) but HURTS under CV00 (they are masked out).
  w <- cfg$predict_post.blend_obs_w
  q <- q + if (scheme == "CV0") 0.10 * w - 0.12 * w^2 else -0.10 * w
  if (identical(cfg$predict_post.method, "cond_expectation")) q <- q + 0.02

  # Synergy bonus rewarding genuine recombination of strong blocks.
  if (cfg$train_select.method == "top_k_similar" &&
      cfg$kernel.method == "em_combine" &&
      cfg$model.method == "gblup_sommer_GE") q <- q + 0.04

  # Heritability caps achievable accuracy; squash to a plausible range.
  q <- q * (0.5 + trial$heritability)
  max(min(q, 0.85), -0.1)
}

.sim_evaluate <- function(cfg, trial, scheme, settings = NULL) {
  truth <- .sim_true(cfg, trial, scheme)
  # Observation noise on top of the deterministic truth. settings$sim_noise_sd = 0 makes an
  # evaluation exact, which (with sim_fixed_trial) turns the objective into a noiseless
  # black-box function -- the regime in which the search algorithm itself can be tested.
  sd    <- if (is.null(settings$sim_noise_sd)) 0.07 else settings$sim_noise_sd
  noise <- if (sd > 0) stats::rnorm(1, 0, sd) else 0
  list(score = max(min(truth + noise, 0.95), -0.3),
       n_test = trial$n_acc, status = "ok", reason = NA_character_)
}
