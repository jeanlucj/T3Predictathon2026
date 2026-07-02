# grm_pareto_genotyped.R
#
# Re-does the GRM reliability/Pareto frontier with the x-axis = number of
# GENOTYPED accessions in the training population (n_train_geno) instead of the
# DECLARED training size (size_declared) used by reliability_summary.R. This makes
# the GRM Pareto share an x-axis with the observed-accuracy variant
# (observed_accuracy_pareto.R), so the two are directly comparable.
#
# y-axis is unchanged: sqrt-mean expected reliability (h2=0.5 primary; 0.2/0.8
# sensitivity). Axes per (trial, scheme): HIGH reliability (good), SMALL genotyped
# training size (good). Reads reliability_by_alg_trial.csv (already written by the
# full-coverage reliability_pareto.R run). Writes ONLY *_genotyped.* files, so the
# declared-size outputs are preserved.
#
# LOCATION: Claude-generated; scripts/Analysis_Claude/. Read-only inputs; writes
# under ./output/reliability/.

library(tidyverse)
here::i_am("grm_pareto_genotyped.R")
out_dir <- here::here("output", "reliability")

res <- readr::read_csv(file.path(out_dir, "reliability_by_alg_trial.csv"),
                       show_col_types = FALSE)

# Pareto within (trial, scheme): maximise reliability, minimise GENOTYPED size.
on_pareto <- function(rel, size) {
  ok <- !is.na(rel) & !is.na(size); out <- rep(NA, length(rel))
  out[ok] <- vapply(which(ok), function(i)
    !any(ok & rel >= rel[i] & size <= size[i] & (rel > rel[i] | size < size[i])),
    logical(1))
  out
}

# --- flagged table + tally (h2 = 0.5) on genotyped size ---------------------
flagged <- res |>
  group_by(trial, cv) |>
  mutate(pareto = on_pareto(sqrtrel_h2_0.5, n_train_geno)) |>
  ungroup()
readr::write_csv(flagged, file.path(out_dir, "reliability_pareto_flagged_genotyped.csv"))

tally <- flagged |>
  filter(!is.na(pareto)) |>
  group_by(alg, cv) |>
  summarise(n_pareto = sum(pareto), n_scored = n(), .groups = "drop") |>
  pivot_wider(names_from = cv, values_from = c(n_pareto, n_scored)) |>
  arrange(alg)
readr::write_csv(tally, file.path(out_dir, "pareto_tally_genotyped.csv"))

# --- h2-sensitivity of the genotyped-size tally -----------------------------
h2cols <- c(h2_0.2 = "sqrtrel_h2_0.2", h2_0.5 = "sqrtrel_h2_0.5", h2_0.8 = "sqrtrel_h2_0.8")
tally_h2 <- imap_dfr(h2cols, function(col, h2) {
  res |>
    group_by(trial, cv) |>
    mutate(par = on_pareto(.data[[col]], n_train_geno)) |>
    ungroup() |>
    filter(!is.na(par)) |>
    group_by(alg, cv) |>
    summarise(n_pareto = sum(par), n_scored = n(), .groups = "drop") |>
    mutate(h2 = h2)
}) |>
  mutate(cell = sprintf("%d/%d", n_pareto, n_scored)) |>
  select(alg, cv, h2, cell) |>
  pivot_wider(names_from = c(cv, h2), values_from = cell) |> arrange(alg)
readr::write_csv(tally_h2, file.path(out_dir, "pareto_tally_h2_sensitivity_genotyped.csv"))

cat("=== GRM Pareto on GENOTYPED size: frontier membership (h2 sensitivity) ===\n")
print(as.data.frame(tally_h2), row.names = FALSE)

# --- figures (separate filenames) -------------------------------------------
lab <- c(Pred1="Prediction1",Pred2="Prediction2",Pred3="Prediction3",
         Pred4="Prediction4",Pred5="Prediction5")
alg_pal <- setNames(scales::hue_pal()(5), unname(lab))  # fixed colour per algorithm

plot_scheme <- function(scheme, file) {
  d <- res |>
    filter(cv == scheme, !is.na(sqrtrel_h2_0.5), !is.na(n_train_geno), n_train_geno > 0) |>
    group_by(trial) |>
    mutate(frontier = on_pareto(sqrtrel_h2_0.5, n_train_geno), alglab = lab[alg]) |>
    ungroup()
  p <- ggplot(d, aes(n_train_geno, sqrtrel_h2_0.5)) +
    geom_point(data = ~filter(.x, frontier), shape = 21, size = 5, stroke = 1.1,
               colour = "grey20", fill = NA) +
    geom_point(aes(colour = alglab), size = 2.6)
  # Prediction1 and Prediction4 select the SAME accessions and build the SAME
  # VanRaden GRM, so their points coincide EXACTLY only under CV0; mark that shared
  # point with a split half-disc: left (U+25D0)=Prediction1, right (U+25D1)=
  # Prediction4 (size tuned to the single-point diameter). Under CV00 their CV-
  # masked training sets differ, so they render as two separate points (no glyph).
  if (scheme == "CV0") {
    pair <- d |>
      filter(alg %in% c("Pred1", "Pred4")) |>
      group_by(trial) |>
      filter(dplyr::n_distinct(alg) == 2) |>
      summarise(n_train_geno = exp(mean(log(n_train_geno))),
                sqrtrel_h2_0.5 = mean(sqrtrel_h2_0.5), .groups = "drop")
    p <- p +
      geom_text(data = pair, label = "◐", colour = alg_pal["Prediction1"], size = 3.07) +
      geom_text(data = pair, label = "◑", colour = alg_pal["Prediction4"], size = 3.07)
  }
  p <- p +
    facet_wrap(~ trial, scales = "free_x") +
    scale_x_log10() +
    scale_colour_manual(values = alg_pal) +
    labs(x = "genotyped training size", y = "√ mean reliability", colour = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file, p, width = 9, height = 7, dpi = 150)
  cat("wrote", file, "\n")
}
plot_scheme("CV0",  file.path(out_dir, "pareto_scatter_genotyped_CV0.png"))
plot_scheme("CV00", file.path(out_dir, "pareto_scatter_genotyped_CV00.png"))
