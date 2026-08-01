# surrogate.R
#
# The surrogate model at the heart of the SMAC-style optimizer: a random forest
# that predicts a configuration's score from its features, written transparently
# as a bag of regression trees (rpart) rather than pulled from a black-box
# package. Each tree sees a bootstrap sample of the evaluations and a random
# subset of the features (row + feature bagging = a random forest). The ensemble
# MEAN is the predicted score; the spread ACROSS trees is the model's
# uncertainty. Trees handle categorical method choices and NA (inapplicable)
# parameters natively, which is exactly why SMAC uses a forest instead of a
# Gaussian process for this kind of mixed search space.
#
# We then turn (mean, uncertainty) into an Expected Improvement acquisition score
# that the optimizer maximizes when choosing which configuration to run next.

library(tidyverse)

# Fit the surrogate. `features` is the data.frame from configs_to_features();
# `y` the observed scores (e.g. predictive ability). Returns a list of fitted
# trees plus the feature names, or NULL if there is not yet enough signal to fit.
fit_surrogate <- function(features, y,
                          ntree   = 300,
                          mtry    = max(2, floor(sqrt(ncol(features)))),
                          min_obs = 12) {
  ok <- is.finite(y)
  features <- features[ok, , drop = FALSE]
  y <- y[ok]
  if (length(y) < min_obs || stats::var(y) < .Machine$double.eps) return(NULL)

  feat_names <- names(features)
  ctrl <- rpart::rpart.control(minsplit = 6, minbucket = 3, cp = 0.001, xval = 0)

  trees <- lapply(seq_len(ntree), function(i) {
    rows <- sample.int(nrow(features), replace = TRUE)          # row bootstrap
    cols <- sample(feat_names, mtry)                            # feature bagging
    df <- features[rows, cols, drop = FALSE]
    df$.y <- y[rows]
    tryCatch(
      rpart::rpart(.y ~ ., data = df, method = "anova", control = ctrl),
      error = function(e) NULL)
  })
  trees <- trees[!vapply(trees, is.null, logical(1))]
  if (!length(trees)) return(NULL)

  structure(list(trees = trees, feat_names = feat_names, y_obs = y,
                 trial_levels = if ("trial_id" %in% feat_names) features$trial_id else NULL),
            class = "rf_surrogate")
}

# Is the eval slice replicated enough for trial_id to be a SAFE surrogate feature?
#
# rpart splits a many-level factor by ordering its levels on the response and cutting along
# that order -- so with one observation per trial the ordering IS the noise, and the split
# memorises it (in-sample R^2 0.98 against a pure-noise response, measured with this
# project's rpart.control; minbucket = 3 does not prevent it). With >= 2 configurations per
# trial the same fit recovers the real trial effect. Gate on that rather than trusting
# settings$trial_replication to have been in force for the whole store.
trial_feature_usable <- function(evals, min_cfg = 2L) {
  if (!nrow(evals) || !all(c("trial_id", "config_hash") %in% names(evals))) return(FALSE)
  n <- tapply(evals$config_hash, evals$trial_id, function(x) length(unique(x)))
  length(n) >= 2L && all(n >= min_cfg)
}

# Expected score of each candidate over the TRIAL POPULATION, when trial_id is a feature.
#
# The objective is "best on a random trial", so a candidate is scored by predicting it on
# every observed trial and averaging -- no unseen factor level is ever needed. Marginalise
# WITHIN each tree first, so the reported sd is the uncertainty of the marginal estimate
# across trees, not the trial-to-trial spread (which is a property of the trials, not of how
# well we know the configuration).
predict_surrogate_marginal <- function(model, newfeatures, trials) {
  trials <- factor(trials, levels = levels(model$trial_levels %||% factor(trials)))
  n_c <- nrow(newfeatures); n_t <- length(trials)
  grid <- newfeatures[rep(seq_len(n_c), times = n_t), , drop = FALSE]
  grid$trial_id <- rep(trials, each = n_c)
  per_tree <- vapply(model$trees, function(tr) {
    p <- as.numeric(stats::predict(tr, newdata = grid))
    rowMeans(matrix(p, nrow = n_c, ncol = n_t))     # mean over trials, inside this tree
  }, FUN.VALUE = numeric(n_c))
  if (is.null(dim(per_tree))) per_tree <- matrix(per_tree, nrow = n_c)
  list(mean = rowMeans(per_tree),
       sd   = apply(per_tree, 1, stats::sd))
}

# Predict mean and between-tree sd for new configurations (a feature data.frame).
predict_surrogate <- function(model, newfeatures) {
  # One column per tree, one row per candidate.
  preds <- vapply(model$trees, function(tr) {
    as.numeric(stats::predict(tr, newdata = newfeatures))
  }, FUN.VALUE = numeric(nrow(newfeatures)))
  if (is.null(dim(preds))) preds <- matrix(preds, nrow = nrow(newfeatures))

  mean_hat <- rowMeans(preds)
  sd_hat   <- apply(preds, 1, stats::sd)
  # Floor the uncertainty so a candidate the forest happens to agree on is not
  # treated as known with certainty (keeps exploration alive).
  sd_floor <- 0.25 * stats::sd(model$y_obs)
  sd_hat   <- pmax(sd_hat, sd_floor)
  list(mean = mean_hat, sd = sd_hat)
}

# Expected Improvement over the current best score (maximization). `xi` trades
# off exploration (larger -> more exploratory).
expected_improvement <- function(mean_hat, sd_hat, best, xi = 0.01) {
  z <- (mean_hat - best - xi) / sd_hat
  improvement <- (mean_hat - best - xi) * stats::pnorm(z) + sd_hat * stats::dnorm(z)
  improvement[sd_hat <= 0] <- 0
  pmax(improvement, 0)
}
