# observed_accuracy_pareto.R
#
# A simpler variant of the GRM reliability/Pareto analysis (reliability_pareto.R).
# Instead of the model-derived sqrt-mean reliability, the y-axis is the OBSERVED
# predictive accuracy each submission actually achieved on the Predictathon,
# taken directly from the anonymized copy at
#   data/core_germplasm_long_with_overlap_df.csv
# (column `y` = per-task accuracy on the core germplasm; per studyName x
# submission_folder x CV_type, submission_folder = Prediction1..5).
#
# x-axis = number of GENOTYPED accessions going into the training population
# (n_train_geno), reused from reliability_by_alg_trial.csv so it matches the GRM
# analysis's notion of "training size that actually carries genomic information".
# NOTE: the GRM Pareto plotted DECLARED training size (size_declared); this variant
# deliberately uses the genotyped count instead, per request.
#
# Pareto axes per (trial, CV scheme): HIGH observed accuracy (good), SMALL
# genotyped training size (good). Tally frontier membership per algorithm.
#
# LOCATION: Claude-generated; scripts/Analysis_Claude/. Reads the observed-accuracy
# table and reliability_by_alg_trial.csv (read-only); writes only under
# ./output/reliability/. No BrAPI, no heavy compute.

library(tidyverse)
here::i_am("observed_accuracy_pareto.R")
out_dir <- here::here("output", "reliability")

# observed accuracy; submission_folder in the copied data/ CSV is already the
# anonymized Prediction1..5 -- map to the Pred1..5 short labels used elsewhere.
sub2alg <- c(Prediction1 = "Pred1", Prediction2 = "Pred2", Prediction3 = "Pred3",
             Prediction4 = "Pred4", Prediction5 = "Pred5")
obs <- readr::read_csv(here::here("data", "core_germplasm_long_with_overlap_df.csv"),
                       show_col_types = FALSE) |>
  dplyr::transmute(trial = studyName, alg = sub2alg[submission_folder],
                   cv = CV_type, accuracy = y, task_status)

# genotyped training size (n_train_geno) from the full-coverage GRM run
size <- readr::read_csv(file.path(out_dir, "reliability_by_alg_trial.csv"),
                        show_col_types = FALSE) |>
  dplyr::select(alg, trial, cv, n_train_geno)

d <- dplyr::left_join(obs, size, by = c("alg", "trial", "cv"))

# Pareto frontier within (trial, scheme): maximise accuracy, minimise n_train_geno.
# Cells with no observed accuracy (failed tasks) or no size are not scorable.
on_pareto <- function(acc, size) {
  ok <- !is.na(acc) & !is.na(size); out <- rep(NA, length(acc))
  out[ok] <- vapply(which(ok), function(i)
    !any(ok & acc >= acc[i] & size <= size[i] & (acc > acc[i] | size < size[i])),
    logical(1))
  out
}

d <- d |>
  group_by(trial, cv) |>
  mutate(pareto = on_pareto(accuracy, n_train_geno)) |>
  ungroup() |>
  arrange(trial, cv, dplyr::desc(accuracy))
readr::write_csv(d, file.path(out_dir, "observed_accuracy_vs_size.csv"))

tally <- d |>
  filter(!is.na(pareto)) |>
  group_by(alg, cv) |>
  summarise(n_pareto = sum(pareto), n_scored = n(), .groups = "drop") |>
  pivot_wider(names_from = cv, values_from = c(n_pareto, n_scored)) |>
  arrange(alg)
readr::write_csv(tally, file.path(out_dir, "observed_accuracy_pareto_tally.csv"))

cat("=== Observed-accuracy Pareto tally: n_on_frontier / n_scored ===\n")
tally |>
  mutate(CV0 = sprintf("%d/%d", n_pareto_CV0, n_scored_CV0),
         CV00 = sprintf("%d/%d", n_pareto_CV00, n_scored_CV00)) |>
  select(alg, CV0, CV00) |> as.data.frame() |> print(row.names = FALSE)

# Pareto scatter, one panel per trial (failed/NA-accuracy points dropped; any
# n_train_geno == 0 cannot sit on a log axis and is dropped from the plot only).
lab <- c(Pred1="Prediction1",Pred2="Prediction2",Pred3="Prediction3",
         Pred4="Prediction4",Pred5="Prediction5")
# Fixed colour per algorithm (= ggplot's default 5-hue palette) so each panel/
# figure uses the SAME colour for a given algorithm even when some are absent
# (e.g. Prediction2 has no CV00 score). Named values keep colours stable; with no
# `limits`, the legend shows only algorithms actually present in the figure.
alg_pal <- setNames(scales::hue_pal()(5), unname(lab))
plot_scheme <- function(scheme, file) {
  dd <- d |> filter(cv == scheme, !is.na(accuracy), !is.na(n_train_geno), n_train_geno > 0) |>
    mutate(alglab = lab[alg])
  p <- ggplot(dd, aes(n_train_geno, accuracy)) +
    geom_point(data = ~filter(.x, pareto), shape = 21, size = 5, stroke = 1.1,
               colour = "grey20", fill = NA) +
    geom_point(aes(colour = alglab), size = 2.6) +
    facet_wrap(~ trial, scales = "free_x") +
    scale_x_log10() +
    scale_colour_manual(values = alg_pal) +
    labs(x = "genotyped training size", y = "observed accuracy", colour = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file, p, width = 9, height = 7, dpi = 150)
  cat("wrote", file, "\n")
}
plot_scheme("CV0",  file.path(out_dir, "observed_accuracy_pareto_CV0.png"))
plot_scheme("CV00", file.path(out_dir, "observed_accuracy_pareto_CV00.png"))
