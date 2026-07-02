# predict_cv0_cv00.R
#
# Re-runs Prediction3's submitted model - the multikernel GBLUP (sommer::mmer
# with a genomic kernel G and an environmental kernel E) - for CV0 and CV00,
# using the BUG-FIXED per-trial GRMs from build_per_trial_grms.R.
#
# To isolate the effect of the GRM dosage-transposition fix, this REUSES the
# participant's per-scheme training-trial and accession SELECTION from
# ../Prediction3/Predictions/<trial>/CV{0,00}_{Trials,Accessions}.csv (CV0 and
# CV00 select different training-trial sets - see analysis notes). Only the GRM
# changes. The model itself mirrors multi_GBLUP() in GBLUP_predictions.Rmd:
#   y ~ 1 + Z_g u_g + Z_e u_e + e,  u_g ~ N(0, sigma2_g G),  u_e ~ N(0, sigma2_e E)
# test accessions get NA phenotypes and are predicted via the G (and E) kernels.
#
# LOCATION: Claude-generated; in scripts/Analysis_Claude/. Reads participant
# inputs from ../Prediction3/ and the per-trial GRMs from ./output/; writes
# predictions to ./output/predictions/<trial>/CV{0,00}_Predictions.csv (it does
# NOT touch the participant's Predictions/). Run from this directory.

library(tidyverse)
library(sommer)

here::i_am("predict_cv0_cv00.R")

pred_dir  <- here::here("..", "Prediction3", "Predictions")
data_dir  <- here::here("..", "Prediction3", "data")
analysis_dir <- here::here("..", "Prediction3", "Analysis")
grm_dir   <- here::here("output", "per_trial_grm")
out_root  <- here::here("output", "predictions")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

YIELD <- "grain_yield_kg_ha_co_321_0001218"

# Phenotypes (per trial x accession BLUEs) and the trial-wise environmental
# similarity matrix (the E kernel), indexed by study_name.
phenos <- readr::read_csv(file.path(data_dir, "genotyped_trial_phenotypes.csv"),
                          show_col_types = FALSE) |>
  dplyr::select(study_name, germplasm_name, dplyr::all_of(YIELD))
W <- as.matrix(readr::read_csv(file.path(analysis_dir, "Trial_W_mat.csv"),
                               show_col_types = FALSE))
rownames(W) <- colnames(W)

rd <- function(trial, f, col) {
  readr::read_csv(file.path(pred_dir, trial, f), show_col_types = FALSE)[[col]]
}

# Multikernel GBLUP for one focal trial + CV scheme -> test-accession predictions.
predict_one <- function(trial, cv) {           # cv in c("CV0", "CV00")
  test_acc     <- rd(trial, paste0(cv, "_Predictions.csv"), "germplasmName")
  train_trials <- rd(trial, paste0(cv, "_Trials.csv"), "studyName")
  train_acc    <- rd(trial, paste0(cv, "_Accessions.csv"), "germplasmName")

  G <- readRDS(file.path(grm_dir, paste0(trial, "_grm.rds")))
  g_acc <- rownames(G)

  # Training rows: selected training accessions in their selected training
  # trials, restricted to genotyped lines with an observed yield.
  subset_pheno <- phenos |>
    dplyr::filter(study_name %in% train_trials,
                  germplasm_name %in% train_acc,
                  germplasm_name %in% g_acc,
                  !is.na(.data[[YIELD]]))

  # Test rows: genotyped focal-trial accessions, NA phenotype.
  test_in <- intersect(test_acc, g_acc)
  test_pheno <- tibble::tibble(study_name = trial, germplasm_name = test_in)
  test_pheno[[YIELD]] <- NA_real_

  combined <- dplyr::bind_rows(subset_pheno, test_pheno)
  envs <- intersect(unique(combined$study_name), rownames(W))
  combined <- combined |> dplyr::filter(study_name %in% envs)
  gens <- unique(combined$germplasm_name)

  Gm <- G[gens, gens, drop = FALSE]
  Em <- W[envs, envs, drop = FALSE]

  combined$germplasm_name <- factor(combined$germplasm_name, levels = rownames(Gm))
  combined$study_name     <- factor(combined$study_name, levels = rownames(Em))

  # Fit on observed (training) rows only; test rows (focal env) are held out and
  # predicted through the kernels.
  train_data <- combined |> dplyr::filter(study_name != trial) |> as.data.frame()

  fit <- sommer::mmer(
    stats::as.formula(paste(YIELD, "~ 1")),
    random  = ~ sommer::vs(germplasm_name, Gu = Gm) +
                sommer::vs(study_name, Gu = Em),
    data    = train_data,
    verbose = FALSE)

  mu     <- fit$Beta$Estimate[1]
  g_blup <- fit$U$`u:germplasm_name`[[YIELD]]
  e_blup <- fit$U$`u:study_name`[[YIELD]]

  preds <- tibble::tibble(germplasmName = test_in) |>
    dplyr::mutate(
      g = g_blup[germplasmName],
      e = e_blup[trial],
      g = dplyr::if_else(is.na(g), 0, g),
      e = dplyr::if_else(is.na(e), 0, e),
      prediction = mu + g + e) |>
    dplyr::select(germplasmName, prediction)

  out_dir <- file.path(out_root, trial)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(preds, file.path(out_dir, paste0(cv, "_Predictions.csv")))
  message(sprintf("  %-28s %-4s  train rows %d, gens %d, envs %d, predicted %d",
                  trial, cv, nrow(train_data), length(gens), length(envs), nrow(preds)))
  invisible(preds)
}

# Run all trials by default, or only the trials named on the command line (used
# to fan the 9 trials across several Rscript processes). e.g.
#   Rscript predict_cv0_cv00.R 24Crk_AY2-3 STP1_2025_MCG
all_trials <- list.dirs(pred_dir, full.names = FALSE, recursive = FALSE)
args <- commandArgs(trailingOnly = TRUE)
trials <- if (length(args)) intersect(args, all_trials) else all_trials

for (trial in trials) {
  for (cv in c("CV0", "CV00")) {
    tryCatch(predict_one(trial, cv),
             error = function(e) warning(trial, " ", cv, ": ", conditionMessage(e)))
  }
}
message("\nDone. Predictions written under ", out_root)
