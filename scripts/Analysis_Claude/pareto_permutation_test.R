# pareto_permutation_test.R
#
# Question: the five algorithms have different Pareto-frontier counts (how often
# each lay on the per-trial frontier for reliability and for observed accuracy,
# under CV0 and CV00). Are those differences a property of the algorithms, or
# just a snapshot of these particular focal trials?
#
# Two analyses, per panel (metric x CV scheme):
#   (1) OMNIBUS PERMUTATION TEST of the exchangeability null "all algorithms are
#       equally likely to be on the frontier". Frontier membership is competitive
#       (coupled within a trial) and blocked by trial, so we permute ALGORITHM
#       LABELS WITHIN EACH TRIAL, among the algorithms that were scorable there,
#       holding the per-trial number of frontier slots (the row sum) fixed. This
#       is the exact-permutation form of Cochran's Q and conditions away both the
#       competitive coupling (via the row sums) and unequal scorability (we only
#       permute among scorable competitors). Statistic = sum_j (C_j - E_j)^2,
#       where C_j is alg j's observed frontier count and E_j its permutation-mean
#       count. We also report a per-algorithm two-sided permutation p-value (is
#       this algorithm's count beyond what label-shuffling produces?).
#   (2) PER-ALGORITHM UNCERTAINTY over the choice of focal trials: a bootstrap
#       over trials (the resampling unit is the trial) giving a percentile CI on
#       each algorithm's frontier rate (count / scorable), plus a delete-one-trial
#       jackknife as a robustness check.
#
# Inputs (as instructed):
#   reliability:       output/reliability/reliability_pareto_flagged_genotyped.csv
#                      (y = sqrtrel_h2_0.5, x = n_train_geno)
#   observed accuracy: output/reliability/observed_accuracy_vs_size.csv
#                      (y = accuracy,       x = n_train_geno)
#
# Frontier is RECOMPUTED here (not read from the `pareto` flag) with a tie
# tolerance on y, so that algorithms whose points differ only by floating-point
# noise are both kept. This is what keeps Prediction1 and Prediction4 -- which are
# effectively identical for reliability CV0 -- on the frontier together rather
# than letting a 1e-16 difference put one on and the other off.

library(tidyverse)

here::i_am("scripts/Analysis_Claude/pareto_permutation_test.R")

set.seed(20260624)

n_perm <- 50000L   # omnibus / per-alg label permutations
n_boot <- 20000L   # trial bootstrap resamples
tol_y  <- 1e-6     # y differences <= tol_y are treated as ties (not domination)
algs   <- c("Pred1", "Pred2", "Pred3", "Pred4", "Pred5")

# ---- Pareto frontier (minimize x, maximize y) with a y-tie tolerance ----
# A dominates B iff A is no worse on both axes and strictly better on at least one,
# where "better on y" requires a gap larger than tol_y (x are integer counts and
# are compared exactly). Near-identical points therefore never dominate each other
# and are both retained on the frontier.
on_frontier <- function(x, y, tol = tol_y) {
  n <- length(x)
  keep <- logical(n)
  for (i in seq_len(n)) {
    dominated <- FALSE
    for (j in seq_len(n)) {
      if (j == i) next
      no_worse        <- (x[j] <= x[i]) && (y[j] >= y[i] - tol)
      strictly_better <- (x[j] <  x[i]) || (y[j] >  y[i] + tol)
      if (no_worse && strictly_better) { dominated <- TRUE; break }
    }
    keep[i] <- !dominated
  }
  keep
}

# ---- Build a tidy membership table: one row per (panel, trial, alg) ----
# scorable = TRUE/FALSE; frontier = TRUE/FALSE/NA (NA when not scorable).
build_membership <- function(df) {
  df |>
    group_by(metric, cv, trial) |>
    mutate(frontier = if (sum(scorable) > 0) {
      out <- rep(NA, n())
      out[scorable] <- on_frontier(x[scorable], y[scorable])
      out
    } else NA) |>
    ungroup()
}

rel <- readr::read_csv(
  here::here("scripts", "Analysis_Claude", "output", "reliability",
             "reliability_pareto_flagged_genotyped.csv"),
  show_col_types = FALSE) |>
  transmute(metric = "reliability", cv, trial, alg,
            x = n_train_geno, y = sqrtrel_h2_0.5,
            scorable = !is.na(y))

obs <- readr::read_csv(
  here::here("scripts", "Analysis_Claude", "output", "reliability",
             "observed_accuracy_vs_size.csv"),
  show_col_types = FALSE) |>
  transmute(metric = "observed_accuracy", cv, trial, alg,
            x = n_train_geno, y = accuracy,
            scorable = (task_status == "success") & !is.na(accuracy))

membership <- bind_rows(build_membership(rel), build_membership(obs))

# ---- Sanity check: Pred1 and Pred4 must share frontier membership for rel CV0 ----
p14 <- membership |>
  filter(metric == "reliability", cv == "CV0", alg %in% c("Pred1", "Pred4")) |>
  select(trial, alg, frontier) |>
  pivot_wider(names_from = alg, values_from = frontier)
stopifnot(all(p14$Pred1 == p14$Pred4))
message("Check passed: Pred1 and Pred4 have identical reliability-CV0 frontier ",
        "membership across all trials.")

# ---- Per-panel analysis -------------------------------------------------------
analyse_panel <- function(memb) {
  # memb: rows for one (metric, cv). Keep trials with >=1 scorable alg.
  trials <- memb |>
    filter(scorable) |>
    distinct(trial) |>
    pull(trial)

  # Per-trial: which alg indices are scorable, and how many frontier slots (R_i).
  per_trial <- lapply(trials, function(tr) {
    sub <- memb |> filter(trial == tr)
    idx <- match(sub$alg[sub$scorable], algs)
    won <- sub$alg[which(sub$frontier)]
    list(scorable_idx = idx, R = length(won))
  })

  # Observed counts and scorable denominators.
  C_obs <- integer(length(algs)); S <- integer(length(algs))
  for (pt in per_trial) {
    S[pt$scorable_idx] <- S[pt$scorable_idx] + 1L
  }
  for (alg_i in seq_along(algs)) {
    C_obs[alg_i] <- sum(memb$frontier[memb$alg == algs[alg_i]], na.rm = TRUE)
  }

  # Permutation-mean count E_j = sum over scorable trials of R_i / m_i.
  E <- numeric(length(algs))
  for (pt in per_trial) {
    m <- length(pt$scorable_idx)
    if (m > 0) E[pt$scorable_idx] <- E[pt$scorable_idx] + pt$R / m
  }

  # Permutation null: within each trial assign R_i slots at random among scorable.
  perm_counts <- matrix(0L, nrow = n_perm, ncol = length(algs))
  for (b in seq_len(n_perm)) {
    cnt <- integer(length(algs))
    for (pt in per_trial) {
      if (pt$R > 0) {
        winners <- if (length(pt$scorable_idx) == 1L) pt$scorable_idx
                   else sample(pt$scorable_idx, pt$R)
        cnt[winners] <- cnt[winners] + 1L
      }
    }
    perm_counts[b, ] <- cnt
  }

  # Omnibus statistic and p-value.
  T_obs  <- sum((C_obs - E)^2)
  T_perm <- rowSums(sweep(perm_counts, 2, E)^2)
  p_omni <- (1 + sum(T_perm >= T_obs)) / (1 + n_perm)

  # Per-algorithm two-sided permutation p (count beyond chance, given exposure).
  p_alg <- vapply(seq_along(algs), function(j) {
    if (S[j] == 0) return(NA_real_)
    pg <- (1 + sum(perm_counts[, j] >= C_obs[j])) / (1 + n_perm)
    pl <- (1 + sum(perm_counts[, j] <= C_obs[j])) / (1 + n_perm)
    min(1, 2 * min(pg, pl))
  }, numeric(1))

  # Trial bootstrap: resample trials, recount. Membership within a trial is fixed.
  M_front <- matrix(NA, nrow = length(trials), ncol = length(algs),
                    dimnames = list(trials, algs))
  M_score <- matrix(FALSE, nrow = length(trials), ncol = length(algs),
                    dimnames = list(trials, algs))
  for (ti in seq_along(trials)) {
    sub <- memb |> filter(trial == trials[ti])
    M_front[ti, match(sub$alg, algs)] <- sub$frontier
    M_score[ti, match(sub$alg, algs)] <- sub$scorable
  }
  boot_rate  <- matrix(NA_real_, nrow = n_boot, ncol = length(algs))
  boot_count <- matrix(0L,       nrow = n_boot, ncol = length(algs))
  nt <- length(trials)
  for (b in seq_len(n_boot)) {
    rs    <- sample.int(nt, nt, replace = TRUE)
    front <- M_front[rs, , drop = FALSE]
    score <- M_score[rs, , drop = FALSE]
    cnt   <- colSums(front, na.rm = TRUE)
    den   <- colSums(score)
    boot_count[b, ] <- cnt
    boot_rate[b, ]  <- ifelse(den > 0, cnt / den, NA_real_)
  }
  rate_lo <- apply(boot_rate, 2, quantile, 0.025, na.rm = TRUE)
  rate_hi <- apply(boot_rate, 2, quantile, 0.975, na.rm = TRUE)
  cnt_lo  <- apply(boot_count, 2, quantile, 0.025, na.rm = TRUE)
  cnt_hi  <- apply(boot_count, 2, quantile, 0.975, na.rm = TRUE)

  # Delete-one-trial jackknife on the rate (count / scorable).
  jack_min <- integer(length(algs)); jack_max <- integer(length(algs))
  jack_se  <- numeric(length(algs))
  rate_obs <- ifelse(S > 0, C_obs / S, NA_real_)
  for (alg_i in seq_along(algs)) {
    loo <- vapply(seq_along(trials), function(drop_i) {
      keep  <- setdiff(seq_along(trials), drop_i)
      cnt   <- sum(M_front[keep, alg_i], na.rm = TRUE)
      den   <- sum(M_score[keep, alg_i])
      if (den > 0) cnt / den else NA_real_
    }, numeric(1))
    jc <- vapply(seq_along(trials), function(drop_i) {
      sum(M_front[setdiff(seq_along(trials), drop_i), alg_i], na.rm = TRUE)
    }, numeric(1))
    jack_min[alg_i] <- min(jc); jack_max[alg_i] <- max(jc)
    loo <- loo[is.finite(loo)]
    g <- length(loo)
    jack_se[alg_i] <- if (g > 1)
      sqrt((g - 1) / g * sum((loo - mean(loo))^2)) else NA_real_
  }

  tibble(
    metric = memb$metric[1], cv = memb$cv[1], alg = algs,
    scored = S, count = C_obs, rate = round(rate_obs, 3),
    exp_count = round(E, 2),
    boot_rate_lo = round(rate_lo, 3), boot_rate_hi = round(rate_hi, 3),
    boot_count_lo = cnt_lo, boot_count_hi = cnt_hi,
    jack_count_min = jack_min, jack_count_max = jack_max,
    jack_rate_se = round(jack_se, 3),
    p_alg_2sided = round(p_alg, 4),
    omnibus_p = round(p_omni, 4)
  )
}

panels <- membership |>
  group_by(metric, cv) |>
  group_split()

results <- map_dfr(panels, analyse_panel)

# ---- Report ------------------------------------------------------------------
omni <- results |> distinct(metric, cv, omnibus_p)
cat("\n================ OMNIBUS PERMUTATION TEST (exchangeability null) ===========\n")
cat("H0: all algorithms equally likely to be on the frontier (labels exchangeable\n")
cat("    within trials, among scorable competitors).  ", n_perm, " permutations.\n\n", sep = "")
print(omni, n = Inf)

cat("\n================ PER-ALGORITHM SUMMARY =====================================\n")
cat("count = frontier count; scored = scorable trials; rate = count/scored.\n")
cat("exp_count = mean count under label permutation (exposure-adjusted expectation).\n")
cat("boot_* = 95% trial-bootstrap CI; jack_* = delete-one-trial range / SE.\n")
cat("p_alg_2sided = per-alg two-sided permutation p (not corrected for 20 tests).\n\n")
print(results |>
        select(metric, cv, alg, scored, count, rate, exp_count,
               boot_rate_lo, boot_rate_hi, jack_count_min, jack_count_max,
               jack_rate_se, p_alg_2sided),
      n = Inf)

out_csv <- here::here("scripts", "Analysis_Claude", "output", "permutation",
                      "pareto_permutation_results.csv")
readr::write_csv(results, out_csv)
cat("\nWrote:", out_csv, "\n")
