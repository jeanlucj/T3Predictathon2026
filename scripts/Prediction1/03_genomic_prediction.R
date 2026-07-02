# ============================================================
# 03_genomic_prediction.R
# GBLUP prediction: CV0 + CV00
# Training set defined by Training_Trial_Info.rds (primary + secondary)
# ============================================================

if (!requireNamespace("rrBLUP",    quietly = TRUE)) install.packages("rrBLUP")
if (!requireNamespace("AGHmatrix", quietly = TRUE)) install.packages("AGHmatrix")

library(rrBLUP)
library(AGHmatrix)
library(tidyverse)

here::i_am("scripts/Prediction1/03_genomic_prediction.R")

SNP_DIR         <- here::here("scripts", "Prediction1", "data", "SNP_matrices")
PHENO_FILE      <- here::here("scripts", "Prediction1", "data",
                              "Phenotype_processed", "pheno_blues.rds")
TRIAL_INFO_FILE <- here::here("scripts", "Prediction1", "data",
                              "Phenotype_processed", "Training_Trial_Info.rds")
OUTPUT_DIR      <- here::here("scripts", "Prediction1", "data", "Predictions")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

PRED_VCF_MAP <- list(
  "2025_AYT_Aurora"            = "2025-12-16_15_59_59_fileP7Mi",
  "24Crk_AY2-3"                = "2025-12-11_18_36_07_fileY6KQ",
  "25_Big6_SVREC_SVREC"        = "2025-12-16_15_37_13_filegz0x",
  "CornellMaster_2025_McGowan" = "2025-12-15_16_43_21_fileY4ZA",
  "YT_Urb_25"                  = "2026-01-12_15_02_24_file7kBP",
  "AWY1_DVPWA_2024"            = "2026-02-25_15_45_12_fileGxN4",
  "OHRWW_2025_SPO"             = "2025-12-16_14_44_02_fileoUdn",
  "TCAP_2025_MANKS"            = "2025-12-12_19_42_54_file8tAQ",
  "STP1_2025_MCG"              = "2025-12-16_19_25_20_filek_yE"
)

PRED_TRIALS <- names(PRED_VCF_MAP)

TRIAL_TO_STUDY <- c(
  "2025_AYT_Aurora"            = "10673",
  "24Crk_AY2-3"                = "10674",
  "25_Big6_SVREC_SVREC"        = "10675",
  "CornellMaster_2025_McGowan" = "10676",
  "YT_Urb_25"                  = "10677",
  "AWY1_DVPWA_2024"            = "10678",
  "OHRWW_2025_SPO"             = "10679",
  "TCAP_2025_MANKS"            = "10680",
  "STP1_2025_MCG"              = "10681"
)

# ── Load phenotype BLUEs ──────────────────────────────────
# pheno_blues columns:
#   focal_study_id | trial_name | study_db_id (training study) | source_type | germplasmName | BLUE
cat("Loading phenotype BLUE data...\n")
pheno_blues <- readRDS(PHENO_FILE)
cat(sprintf("  %d records, %d training studies, %d accessions\n",
            nrow(pheno_blues),
            n_distinct(pheno_blues$study_db_id),
            n_distinct(pheno_blues$germplasmName)))

# ── Load Training_Trial_Info ──────────────────────────────
cat("\nLoading Training_Trial_Info.rds...\n")
if (!file.exists(TRIAL_INFO_FILE))
  stop("Training_Trial_Info.rds not found. Please run 02_prepare_phenotypes.R first.")

training_info <- readRDS(TRIAL_INFO_FILE)
training_map  <- setNames(training_info,
                           sapply(training_info, `[[`, "study_id"))
cat(sprintf("  Loaded info for %d focal trials\n", length(training_map)))

# ── Load all SNP matrices ─────────────────────────────────
cat("\nLoading SNP matrices...\n")
all_snp <- list()
for (f in list.files(SNP_DIR, pattern = "_snp_matrix\\.rds$", full.names = TRUE)) {
  label <- gsub("_snp_matrix\\.rds$", "", basename(f))
  all_snp[[label]] <- readRDS(f)
  cat(sprintf("  %s: %d x %d\n", label,
              nrow(all_snp[[label]]), ncol(all_snp[[label]])))
}

# ── Build GRM ─────────────────────────────────────────────
# Convert -1/0/1 encoding to 0/1/2 for Gmatrix(), add ridge term
build_grm <- function(snp_mat) {
  mat <- pmin(pmax(round(snp_mat + 1), 0L), 2L)
  G   <- Gmatrix(SNPmatrix = mat, method = "VanRaden",
                 missingValue = NA, maf = 0.05)
  G + diag(1e-4, nrow(G))
}

# ── GBLUP model ───────────────────────────────────────────
run_gblup <- function(y_train, G, test_lines) {
  train_lines <- intersect(names(y_train), rownames(G))
  test_lines  <- intersect(test_lines, rownames(G))

  if (length(train_lines) < 10) {
    warning("Fewer than 10 training accessions in GRM; skipping GBLUP")
    return(NULL)
  }

  # Standardize y for numerical stability
  y_mu <- mean(y_train[train_lines], na.rm = TRUE)
  y_sd <- sd(y_train[train_lines],   na.rm = TRUE)
  if (is.na(y_sd) || y_sd == 0) y_sd <- 1
  y_s  <- (y_train[train_lines] - y_mu) / y_sd

  fit <- tryCatch(
    mixed.solve(y = y_s, K = G[train_lines, train_lines], SE = FALSE),
    error = function(e) { cat(sprintf("    %s\n", e$message)); NULL }
  )
  if (is.null(fit)) return(NULL)

  u <- fit$u; names(u) <- train_lines

  # Predict test accessions: u_test = G21 * G11^-1 * u_train
  Ginv <- tryCatch(
    solve(G[train_lines, train_lines]),
    error = function(e) {
      warning("G matrix singular; using Moore-Penrose pseudoinverse")
      MASS::ginv(G[train_lines, train_lines])
    }
  )
  u_test <- as.vector(G[test_lines, train_lines] %*% Ginv %*% u)
  names(u_test) <- test_lines

  # Back-transform to original scale
  u_test * y_sd + y_mu
}

# ── Select training set ───────────────────────────────────
# CV0:  use all primary + secondary training studies for this focal trial
# CV00: same as CV0, additionally exclude genotyped focal accessions
get_train_df <- function(focal_trial, cv, test_lines) {
  study_id <- TRIAL_TO_STUDY[focal_trial]
  info     <- training_map[[study_id]]

  if (is.null(info)) {
    warning(sprintf("No training info for study %s; using all other focal trial data", study_id))
    train_df <- pheno_blues %>% filter(focal_study_id != study_id)
  } else {
    pri_ids <- info$primary_study_db_ids[!is.na(info$primary_study_db_ids)]
    sec_ids <- info$secondary_study_db_ids[!is.na(info$secondary_study_db_ids)]
    all_ids <- unique(c(pri_ids, sec_ids))

    train_df <- pheno_blues %>%
      filter(focal_study_id == study_id,
             study_db_id    %in% all_ids)

    if (nrow(train_df) == 0) {
      warning(sprintf(
        "No phenotype data for training studies of study %s.\nEnsure 00_download_training_phenotypes.R has been run.",
        study_id
      ))
      train_df <- pheno_blues %>% filter(focal_study_id != study_id)
    } else {
      cat(sprintf("    [%s] primary+secondary -> %d records, %d studies, %d accessions\n",
                  cv, nrow(train_df),
                  n_distinct(train_df$study_db_id),
                  n_distinct(train_df$germplasmName)))
    }
  }

  if (cv == "CV00")
    train_df <- train_df %>% filter(!germplasmName %in% test_lines)

  train_df
}

# ── Write submission files ────────────────────────────────
write_submission <- function(preds, train_df, trial_out, cv) {
  write_csv(tibble(germplasmName = names(preds),
                   prediction    = as.numeric(preds)),
            file.path(trial_out, paste0(cv, "_Predictions.csv")))
  write_csv(tibble(studyName = as.character(unique(na.omit(train_df$study_db_id)))),
            file.path(trial_out, paste0(cv, "_Trials.csv")))
  write_csv(tibble(germplasmName = unique(train_df$germplasmName)),
            file.path(trial_out, paste0(cv, "_Accessions.csv")))
}

# ════════════════════════════════════════════════════════════
cat("\nStarting predictions...\n"); cat(strrep("=", 55), "\n")

for (focal_trial in PRED_TRIALS) {
  cat(sprintf("\n> %s\n", focal_trial))

  trial_out <- file.path(OUTPUT_DIR, focal_trial)
  dir.create(trial_out, recursive = TRUE, showWarnings = FALSE)

  pred_label <- PRED_VCF_MAP[[focal_trial]]
  if (!pred_label %in% names(all_snp)) {
    cat(sprintf("  SNP matrix not found: %s\n", pred_label)); next
  }

  pred_snp   <- all_snp[[pred_label]]
  test_lines <- rownames(pred_snp)
  cat(sprintf("  Test accessions: %d\n", length(test_lines)))

  # ── Build combined SNP matrix for GRM ──────────────────
  # Find SNP matrices that contain training accessions (>=5 matches)
  train_lines_all <- pheno_blues %>%
    filter(focal_study_id == TRIAL_TO_STUDY[focal_trial]) %>%
    pull(germplasmName) %>% unique()

  useful_snp <- names(all_snp)[sapply(names(all_snp), function(nm) {
    sum(train_lines_all %in% rownames(all_snp[[nm]])) >= 5
  })]
  cat(sprintf("  SNP matrices containing training accessions: %d\n", length(useful_snp)))

  # Merge by shared SNP positions (threshold: >= 500 shared SNPs)
  shared   <- colnames(pred_snp)
  combined <- pred_snp
  for (nm in useful_snp) {
    if (nm == pred_label) next
    common <- intersect(colnames(all_snp[[nm]]), shared)
    if (length(common) >= 500) {
      combined <- rbind(combined[, common, drop = FALSE],
                        all_snp[[nm]][, common, drop = FALSE])
      shared   <- common
    }
  }

  # If too few shared SNPs, lower threshold to 100
  if (ncol(combined) < 500) {
    shared   <- colnames(pred_snp)
    combined <- pred_snp
    for (nm in useful_snp) {
      if (nm == pred_label) next
      common <- intersect(colnames(all_snp[[nm]]), shared)
      if (length(common) >= 100) {
        combined <- rbind(combined[, common, drop = FALSE],
                          all_snp[[nm]][, common, drop = FALSE])
        shared   <- common
      }
    }
  }

  combined <- combined[!duplicated(rownames(combined)), ]
  cat(sprintf("  GRM: %d accessions x %d SNPs\n", nrow(combined), ncol(combined)))
  G <- build_grm(combined)

  # ── CV0 / CV00 predictions ─────────────────────────────
  for (cv in c("CV0", "CV00")) {
    cat(sprintf("  [%s] ", cv))

    train_df <- get_train_df(focal_trial, cv, test_lines)

    y_vec <- train_df %>%
      group_by(germplasmName) %>%
      summarise(y = mean(BLUE, na.rm = TRUE), .groups = "drop") %>%
      { setNames(.$y, .$germplasmName) }

    preds <- run_gblup(y_vec, G, test_lines)

    if (!is.null(preds)) {
      write_submission(preds, train_df, trial_out, cv)
      cat(sprintf("Done: %d accessions\n", length(preds)))
    } else {
      # Fallback: global mean imputation
      global_mean    <- mean(y_vec, na.rm = TRUE)
      preds_fallback <- setNames(rep(global_mean, length(test_lines)), test_lines)
      write_submission(preds_fallback, train_df, trial_out, cv)
      cat(sprintf("Fallback (mean imputation): %.1f\n", global_mean))
    }
  }
}

cat(sprintf("\n%s\nDone. Results saved to: %s\n", strrep("=", 55), OUTPUT_DIR))
