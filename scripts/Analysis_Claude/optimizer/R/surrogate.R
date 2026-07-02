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

  structure(list(trees = trees, feat_names = feat_names, y_obs = y),
            class = "rf_surrogate")
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
