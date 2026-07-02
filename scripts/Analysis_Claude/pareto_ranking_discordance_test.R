# pareto_ranking_discordance_test.R
#
# Question: did observed accuracy rank the algorithms DIFFERENTLY from expected
# reliability, or is the apparent disagreement just sampling noise?
#
# "Different ranking" is recast as an ALGORITHM x METRIC INTERACTION, so that
# rejection supports the claim. For each trial we keep the algorithms scorable on
# BOTH metrics, rank them within the trial under reliability and under observed
# accuracy, and standardize each within-trial rank to [-1, 1] (best = +1, worst =
# -1) so trials with different numbers of algorithms are comparable. Working in
# within-trial ranks removes the incommensurable scales and the metric main
# effect, leaving the interaction cleanly defined.
#
# (1) OMNIBUS INTERACTION TEST. Per algorithm a, delta_{t,a} = (std reliability
#     rank) - (std accuracy rank) in trial t; S_a = sum_t delta_{t,a}. Statistic
#     T = sum_a S_a^2. Under H0 (no algorithm x metric interaction) the two
#     metric labels are exchangeable within each trial, which flips the sign of a
#     trial's whole delta vector -- an exact sign-flip test. With <= ~20 retained
#     trials we ENUMERATE all 2^n sign assignments for an exact p-value.
# (2) CONCORDANCE EFFECT SIZE. Per-trial Kendall's tau between the two rankings,
#     averaged over trials, with a trial-bootstrap 95% CI. tau near 0 / negative
#     is the honest quantitative form of "reliability did not predict the
#     accuracy ranking".
# (3) PREDICTION3 REVERSAL. Targeted test that Prediction3 ranked better on
#     accuracy than on reliability: exact sign-flip p on its own deltas (read from
#     the same enumeration) plus a paired Wilcoxon signed-rank cross-check.
#
# Trials where either metric is constant across the common set (no within-trial
# ranking information -- e.g. the degenerate 24Crk_AY2-3 panel, reliability all 0)
# are dropped. CV0 is the primary analysis; CV00 is reported but is too sparse
# (commonly-scorable sets of ~3, Prediction2 never scorable) to be informative.

library(tidyverse)

here::i_am("scripts/Analysis_Claude/pareto_ranking_discordance_test.R")

set.seed(20260624)
n_boot <- 20000L
algs   <- c("Pred1", "Pred2", "Pred3", "Pred4", "Pred5")

# Standardize within-trial ranks of a value vector to [-1, 1] (higher = better).
std_rank <- function(v) {
  k <- length(v)
  if (k < 2) return(rep(0, k))
  r <- rank(v, ties.method = "average")
  (r - (k + 1) / 2) / ((k - 1) / 2)
}

rel <- readr::read_csv(
  here::here("scripts", "Analysis_Claude", "output", "reliability",
             "reliability_pareto_flagged_genotyped.csv"),
  show_col_types = FALSE) |>
  transmute(cv, trial, alg, rel = sqrtrel_h2_0.5)

obs <- readr::read_csv(
  here::here("scripts", "Analysis_Claude", "output", "reliability",
             "observed_accuracy_vs_size.csv"),
  show_col_types = FALSE) |>
  transmute(cv, trial, alg,
            acc = if_else(task_status == "success", accuracy, NA_real_))

# Common-scorable cells, within-trial standardized ranks, and per-cell delta.
cells <- full_join(rel, obs, by = c("cv", "trial", "alg")) |>
  filter(!is.na(rel), !is.na(acc)) |>
  group_by(cv, trial) |>
  filter(n() >= 2, sd(rel) > 0, sd(acc) > 0) |>   # need ranking info in both
  mutate(rel_sr = std_rank(rel),
         acc_sr = std_rank(acc),
         delta  = rel_sr - acc_sr,
         tau_trial = cor(rel, acc, method = "kendall")) |>
  ungroup()

analyse_cv <- function(d, scheme) {
  d <- d |> filter(cv == scheme)
  trials <- sort(unique(d$trial))
  n_tr <- length(trials)

  # Trials x algorithms delta matrix (0 where an algorithm is not in the common set).
  D <- matrix(0, nrow = n_tr, ncol = length(algs),
              dimnames = list(trials, algs))
  for (i in seq_along(trials)) {
    sub <- d |> filter(trial == trials[i])
    D[i, match(sub$alg, algs)] <- sub$delta
  }

  # Exact sign-flip enumeration (all 2^n_tr metric-label swaps).
  signs  <- as.matrix(expand.grid(rep(list(c(-1, 1)), n_tr)))  # 2^n x n
  S_perm <- signs %*% D                                        # 2^n x 5
  S_obs  <- colSums(D)
  T_perm <- rowSums(S_perm^2)
  T_obs  <- sum(S_obs^2)
  p_omni <- mean(T_perm >= T_obs)

  p_alg <- vapply(seq_along(algs), function(j)
    mean(S_perm[, j]^2 >= S_obs[j]^2), numeric(1))

  # Concordance: mean per-trial Kendall tau, trial-bootstrap CI.
  taus <- d |> distinct(trial, tau_trial) |> pull(tau_trial)
  tau_mean <- mean(taus)
  boot_tau <- replicate(n_boot, mean(sample(taus, length(taus), replace = TRUE)))
  tau_ci   <- quantile(boot_tau, c(0.025, 0.975))

  # Per-algorithm rank means (for interpretation).
  means <- d |>
    group_by(alg) |>
    summarise(n_trials = n(),
              mean_rel_sr = mean(rel_sr),
              mean_acc_sr = mean(acc_sr),
              S = sum(delta), .groups = "drop")

  list(scheme = scheme, n_tr = n_tr, trials = trials,
       p_omni = p_omni, T_obs = T_obs,
       tau_mean = tau_mean, tau_ci = tau_ci, n_tau = length(taus),
       per_alg = means |>
         left_join(tibble(alg = algs, p_alg = p_alg), by = "alg"),
       D = D, d = d)
}

res_cv0  <- analyse_cv(cells, "CV0")
res_cv00 <- analyse_cv(cells, "CV00")

report <- function(r) {
  cat("\n=================== ", r$scheme,
      " : algorithm x metric rank interaction ===================\n", sep = "")
  cat("Retained trials (both metrics informative): ", r$n_tr, " (",
      paste(r$trials, collapse = ", "), ")\n", sep = "")
  cat(sprintf("Omnibus interaction: exact sign-flip p = %.4f  (T_obs = %.2f, %d sign assignments)\n",
              r$p_omni, r$T_obs, 2^r$n_tr))
  cat(sprintf("Concordance: mean per-trial Kendall tau = %.3f  95%% CI [%.3f, %.3f]  (n_tau = %d)\n",
              r$tau_mean, r$tau_ci[1], r$tau_ci[2], r$n_tau))
  cat("Per algorithm (std rank: +1 best, -1 worst; S = sum(rel - acc), >0 = better",
      "on reliability):\n")
  print(r$per_alg |>
          mutate(across(c(mean_rel_sr, mean_acc_sr, S), \(x) round(x, 3)),
                 p_alg = round(p_alg, 4)),
        n = Inf)
}

report(res_cv0)
report(res_cv00)

# ---- Targeted Prediction3 reversal (CV0) -------------------------------------
p3 <- cells |> filter(cv == "CV0", alg == "Pred3") |>
  transmute(trial, rel_sr, acc_sr, diff = acc_sr - rel_sr)   # >0 = better on accuracy
wt <- suppressWarnings(
  wilcox.test(p3$acc_sr, p3$rel_sr, paired = TRUE, exact = FALSE))
cat("\n=================== Prediction3 reversal (CV0) ===================\n")
cat("Per-trial standardized ranks (diff = accuracy - reliability):\n")
print(p3 |> mutate(across(c(rel_sr, acc_sr, diff), \(x) round(x, 3))), n = Inf)
cat(sprintf("Mean diff = %.3f (>0 => Prediction3 ranked higher on accuracy)\n",
            mean(p3$diff)))
cat(sprintf("Exact sign-flip p (from omnibus enumeration) = %.4f\n",
            res_cv0$per_alg$p_alg[res_cv0$per_alg$alg == "Pred3"]))
cat(sprintf("Paired Wilcoxon signed-rank p = %.4f\n", wt$p.value))

# ---- Save ---------------------------------------------------------------------
out <- bind_rows(
  res_cv0$per_alg  |> mutate(cv = "CV0",  omnibus_p = res_cv0$p_omni,
                             tau_mean = res_cv0$tau_mean,
                             tau_lo = res_cv0$tau_ci[1], tau_hi = res_cv0$tau_ci[2]),
  res_cv00$per_alg |> mutate(cv = "CV00", omnibus_p = res_cv00$p_omni,
                             tau_mean = res_cv00$tau_mean,
                             tau_lo = res_cv00$tau_ci[1], tau_hi = res_cv00$tau_ci[2])
) |>
  select(cv, alg, n_trials, mean_rel_sr, mean_acc_sr, S, p_alg,
         omnibus_p, tau_mean, tau_lo, tau_hi)

out_csv <- here::here("scripts", "Analysis_Claude", "output", "permutation",
                      "ranking_discordance_results.csv")
readr::write_csv(out, out_csv)
cat("\nWrote:", out_csv, "\n")
