# reliability_summary.R
#
# Post-processes reliability_by_alg_trial.csv (written by reliability_pareto.R):
#   1. Recomputes the Pareto frontier (high sqrt-mean-reliability, small declared
#      training size) at h2 = 0.2 / 0.5 / 0.8 and tallies frontier membership per
#      algorithm -- a robustness check on the lambda (variance-ratio) assumption.
#   2. Draws the Pareto scatter (one panel per focal trial), frontier points ringed.
#
# LOCATION: Claude-generated; scripts/Analysis_Claude/. Reads/writes only under
# ./output/reliability/. Algorithms are labelled Prediction1..5 (no participant
# names) so outputs are PII-clean. Run after reliability_pareto.R.

library(tidyverse)
here::i_am("reliability_summary.R")
out_dir <- here::here("output", "reliability")

res <- readr::read_csv(file.path(out_dir, "reliability_by_alg_trial.csv"),
                       show_col_types = FALSE)

# Pareto frontier within a (trial, scheme): maximise rel, minimise size.
on_pareto <- function(rel, size) {
  ok <- !is.na(rel) & !is.na(size); out <- rep(NA, length(rel))
  out[ok] <- vapply(which(ok), function(i)
    !any(ok & rel >= rel[i] & size <= size[i] & (rel > rel[i] | size < size[i])),
    logical(1))
  out
}

# --- 1. h2-sensitivity of the frontier tally -------------------------------
h2cols <- c(h2_0.2 = "sqrtrel_h2_0.2", h2_0.5 = "sqrtrel_h2_0.5", h2_0.8 = "sqrtrel_h2_0.8")
tally_all <- imap_dfr(h2cols, function(col, h2) {
  res |>
    group_by(trial, cv) |>
    mutate(par = on_pareto(.data[[col]], size_declared)) |>
    ungroup() |>
    filter(!is.na(par)) |>
    group_by(alg, cv) |>
    summarise(n_pareto = sum(par), n_scored = n(), .groups = "drop") |>
    mutate(h2 = h2)
})
tally_wide <- tally_all |>
  mutate(cell = sprintf("%d/%d", n_pareto, n_scored)) |>
  select(alg, cv, h2, cell) |>
  pivot_wider(names_from = c(cv, h2), values_from = cell) |>
  arrange(alg)
readr::write_csv(tally_wide, file.path(out_dir, "pareto_tally_h2_sensitivity.csv"))
cat("=== Pareto frontier membership: n_on_frontier / n_scored_trials ===\n")
print(as.data.frame(tally_wide), row.names = FALSE)

# --- 2. Pareto scatter, one panel per trial (h2 = 0.5) ---------------------
lab <- c(Pred1="Prediction1", Pred2="Prediction2", Pred3="Prediction3",
         Pred4="Prediction4", Pred5="Prediction5")
plot_scheme <- function(scheme, file) {
  d <- res |>
    filter(cv == scheme, !is.na(sqrtrel_h2_0.5)) |>
    group_by(trial) |>
    mutate(frontier = on_pareto(sqrtrel_h2_0.5, size_declared),
           alglab = lab[alg]) |>
    ungroup()
  p <- ggplot(d, aes(size_declared, sqrtrel_h2_0.5)) +
    geom_point(data = ~filter(.x, frontier),
               shape = 21, size = 5, stroke = 1.1, colour = "grey20", fill = NA) +
    geom_point(aes(colour = alglab), size = 2.6) +
    facet_wrap(~ trial, scales = "free_x") +
    scale_x_log10() +
    labs(x = "training size", y = "√ mean reliability", colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file, p, width = 9, height = 7, dpi = 150)
  cat("wrote", file, "\n")
}
plot_scheme("CV0",  file.path(out_dir, "pareto_scatter_CV0.png"))
plot_scheme("CV00", file.path(out_dir, "pareto_scatter_CV00.png"))
