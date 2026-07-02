# FW_regression_algorithms.R
#
# Finlay-Wilkinson-style analysis of prediction ALGORITHMS across prediction
# tasks. Each task = studyName:CV_type (focal trial x CV scheme); its "index" is
# the model-based EXPECTED accuracy for that task = sqrt(mean reliability) for
# Prediction3 (from prediction_expected_accuracy.csv, computed Friday). Each
# algorithm's accuracy y is modelled with a random slope on the index, asking
# whether algorithms differ in their sensitivity to task difficulty (the FW
# question: do some gain more on easy tasks / lose more on hard ones).
#
# Using an EXTERNAL index (rather than the task mean of y) avoids two problems of
# the classic FW index here: the circularity of regressing y on a mean of those
# same y's, and the collinearity of the index with the (1|task) intercept.
#
#   full:    y ~ 1 + (1 | task) + (1 + index | algorithm)
#   reduced: y ~ 1 + (1 | task) + (1 | algorithm)
# compared by likelihood-ratio test (ML fits).
#
# Accuracy source: core_germplasm_long_with_overlap_df.csv (column y);
# index source:    prediction_expected_accuracy.csv (expected_accuracy);
# task = studyName:CV_type, algorithm = submission_folder.

library(tidyverse)
library(lme4)

here::i_am("FW_regression_algorithms.R")

# Task index = Prediction3 expected accuracy (sqrt mean reliability), per task.
# Inputs live in the tracked, anonymized data/ folder (copied from the T3
# outputs; submission_folder already carries Prediction1..5).
task_index <- readr::read_csv(
  here::here("data", "prediction_expected_accuracy.csv"), show_col_types = FALSE) |>
  dplyr::select(studyName, CV_type, index = expected_accuracy)

dat <- readr::read_csv(
  here::here("data", "core_germplasm_long_with_overlap_df.csv"),
  show_col_types = FALSE) |>
  dplyr::filter(!is.na(y)) |>
  dplyr::left_join(task_index, by = c("studyName", "CV_type")) |>
  dplyr::transmute(
    task      = interaction(studyName, CV_type, sep = ":", drop = TRUE),
    algorithm = factor(submission_folder),
    y         = y,
    index     = index) |>
  dplyr::filter(!is.na(index))

cat(sprintf("rows: %d | tasks: %d | algorithms: %d\n",
            nrow(dat), dplyr::n_distinct(dat$task), dplyr::n_distinct(dat$algorithm)))

full    <- lme4::lmer(y ~ 1 + (1 | task) + (1 + index | algorithm),
                      data = dat, REML = FALSE)
reduced <- lme4::lmer(y ~ 1 + (1 | task) + (1 | algorithm),
                      data = dat, REML = FALSE)

cat("\n=== Reduced: y ~ 1 + (1|task) + (1|algorithm) ===\n")
print(summary(reduced)$varcor)
cat("\n=== Full: y ~ 1 + (1|task) + (1 + index|algorithm) ===\n")
print(summary(full)$varcor)
cat("isSingular(full):", lme4::isSingular(full), "\n")

cat("\n=== Likelihood-ratio test (reduced vs full) ===\n")
print(anova(reduced, full))

cat("\n=== Per-algorithm index sensitivity (random-slope BLUPs) ===\n")
print(round(lme4::ranef(full)$algorithm, 3))

# ---------------------------------------------------------------------------
# Main effect: does the index predict the TASK-MEAN accuracy? (the interaction
# test above sets this aside via (1|task), so check it directly: OLS over the
# 18 tasks of mean-across-algorithms accuracy on the index.)
# ---------------------------------------------------------------------------
task_means <- dat |>
  dplyr::group_by(task, index) |>
  dplyr::summarise(mean_acc = mean(y), n_alg = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(CV_type = sub(".*:", "", as.character(task)))

cat("\n=== Task-mean accuracy ~ index (OLS, all", nrow(task_means), "tasks) ===\n")
print(summary(lm(mean_acc ~ index, data = task_means)))
cat("Pearson r (all tasks):", round(cor(task_means$mean_acc, task_means$index), 3), "\n")

cat("\n=== within each CV scheme (index & accuracy both differ by scheme) ===\n")
for (cv in c("CV0", "CV00")) {
  d <- dplyr::filter(task_means, CV_type == cv)
  cat(sprintf("  %-4s  r = %+.3f  (n = %d)\n", cv, cor(d$mean_acc, d$index), nrow(d)))
}

# ---------------------------------------------------------------------------
# Classic Finlay-Wilkinson figure: per-(task, algorithm) accuracy vs the task
# mean (the FW environmental index = mean across algorithms), one regression per
# algorithm with its confidence band. Dashed line y = x is the "average
# algorithm" (FW slope 1). NB: this uses the task-mean index for the picture;
# the LRT above used the external expected-accuracy index.
# ---------------------------------------------------------------------------
# algorithm is already the anonymized Prediction1..5 label (submission_folder in
# the copied data/ CSV was anonymized on import).
dat_fw <- dat |>
  dplyr::group_by(task) |>
  dplyr::mutate(task_mean = mean(y)) |>
  dplyr::ungroup() |>
  dplyr::mutate(prediction = factor(as.character(algorithm),
                                    levels = sort(unique(as.character(algorithm)))))

fw_plot <- ggplot2::ggplot(
    dat_fw,
    ggplot2::aes(task_mean, y, colour = prediction, shape = prediction)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::scale_shape_manual(values = c(16, 17, 15, 18, 4)) +
  ggplot2::labs(x = "Task mean accuracy (across algorithms)",
                y = "Algorithm accuracy",
                colour = "Algorithm", shape = "Algorithm") +
  ggplot2::theme_bw(base_size = 13)

ggplot2::ggsave(here::here("output", "FW_algorithm_regression.png"),
                plot = fw_plot, width = 9, height = 6.5, dpi = 300)
cat("\nWrote output/FW_algorithm_regression.png\n")

# ---------------------------------------------------------------------------
# Fixed-effects FW heterogeneity-of-slopes test. Algorithm is a 5-level fixed
# factor, so there is NO variance-component / boundary / singular-fit issue;
# this is the classic Eberhart-Russell partition and the F-test on the
# index:algorithm interaction has a defined null.
#   common slope:        y ~ index + algorithm
#   heterogeneous slope: y ~ index * algorithm   (tests index:algorithm)
# Three indices: (a) task mean = classic FW index (but part-whole circular,
# each algorithm is 1/5 of its own index), (b) leave-one-out task mean (mean of
# the OTHER algorithms — removes the part-whole correlation), (c) exogenous
# expected accuracy (sqrt mean reliability; weak, R^2~0.15 vs task mean).
# ---------------------------------------------------------------------------
dat_fe <- dat_fw |>
  dplyr::group_by(task) |>
  dplyr::mutate(n_t = dplyr::n(), loo_mean = (sum(y) - y) / (n_t - 1)) |>
  dplyr::ungroup()

het_test <- function(d, idx) {
  d$xx <- d[[idx]]
  d <- d[!is.na(d$xx), ]
  red <- lm(y ~ xx + prediction, data = d)
  ful <- lm(y ~ xx * prediction, data = d)
  list(an = anova(red, ful), n = nrow(d))
}

for (idx in c("task_mean", "loo_mean", "index")) {
  lab <- c(task_mean = "task mean (classic FW, circular)",
           loo_mean  = "leave-one-out task mean (de-circularized)",
           index     = "exogenous expected accuracy (weak)")[idx]
  cat("\n=== Heterogeneity-of-slopes F-test  |  index =", lab, "===\n")
  res <- het_test(dat_fe, idx)
  print(res$an)
  cat("n =", res$n, "\n")
}

cat("\n=== Per-algorithm FW slopes (task-mean index, 95% CI) ===\n")
for (p in levels(dat_fw$prediction)) {
  d <- dplyr::filter(dat_fw, prediction == p)
  m <- lm(y ~ task_mean, data = d)
  ci <- confint(m)["task_mean", ]
  cat(sprintf("  %-12s slope = %+.2f  [%+.2f, %+.2f]  (n = %d)\n",
              p, coef(m)["task_mean"], ci[1], ci[2], nrow(d)))
}
