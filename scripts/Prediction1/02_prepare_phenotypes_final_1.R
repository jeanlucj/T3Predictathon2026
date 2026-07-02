# ============================================================
# 02_prepare_phenotypes.R
# 1. Auto-download latest Training_Trial_Info.rds from GitHub
# 2. Read primary + secondary phenotype CSVs and compute BLUEs
# ============================================================

library(tidyverse)
library(data.table)
library(dplyr)

here::i_am("scripts/Prediction1/02_prepare_phenotypes_final_1.R")

PHENO_DIR  <- here::here("scripts", "Prediction1", "data", "Phenotype")
OUTPUT_DIR <- here::here("scripts", "Prediction1", "data", "Phenotype_processed")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Download latest Training_Trial_Info.rds from GitHub ──
GITHUB_RDS_URL <- paste0(
  "https://raw.githubusercontent.com/jeanlucj/",
  "T3_predictathon_find_training_trials/main/",
  "output/Training_Trial_Info.rds"
)
LOCAL_TRIAL_INFO <- file.path(OUTPUT_DIR, "Training_Trial_Info.rds")

cat("Downloading latest Training_Trial_Info.rds from GitHub...\n")
tryCatch({
  download.file(GITHUB_RDS_URL, destfile = LOCAL_TRIAL_INFO,
                mode = "wb", quiet = TRUE, timeout = 60)
  cat("  Download successful\n")
}, error = function(e) {
  if (file.exists(LOCAL_TRIAL_INFO)) {
    cat(sprintf("  Network error, using local cache: %s\n", e$message))
  } else {
    stop("Download failed and no local cache found. Check network or download manually:\n  ",
         GITHUB_RDS_URL)
  }
})

training_info <- readRDS(LOCAL_TRIAL_INFO)
training_map  <- setNames(training_info,
                           sapply(training_info, `[[`, "study_id"))
cat(sprintf("  Loaded training info for %d focal trials\n\n", length(training_map)))

for (nm in names(training_map)) {
  info  <- training_map[[nm]]
  n_pri <- length(info$primary_study_db_ids[!is.na(info$primary_study_db_ids)])
  n_sec <- length(info$secondary_study_db_ids[!is.na(info$secondary_study_db_ids)])
  cat(sprintf("  study%s: primary=%d trials, secondary=%d trials\n",
              nm, n_pri, n_sec))
}
cat("\n")

# ── Read all phenotype CSVs ───────────────────────────────
# Folder structure (created by script 00):
#   Phenotype/study10673/study10673/primary_phenotypes_6513.csv
#   Phenotype/study10673/study10673/secondary_phenotypes_5376.csv
#
# CSV columns (janitor::clean_names() applied):
#   observation_variable_name, germplasm_name, value, study_db_id, rep, block

STUDY_TRIAL_MAP <- c(
  "10673" = "2025_AYT_Aurora",
  "10674" = "24Crk_AY2-3",
  "10675" = "25_Big6_SVREC_SVREC",
  "10676" = "CornellMaster_2025_McGowan",
  "10677" = "YT_Urb_25",
  "10678" = "AWY1_DVPWA_2024",
  "10679" = "OHRWW_2025_SPO",
  "10680" = "TCAP_2025_MANKS",
  "10681" = "STP1_2025_MCG"
)

cat("Reading phenotype CSVs...\n")
all_pheno <- list()

study_dirs <- list.dirs(PHENO_DIR, recursive = FALSE)
study_dirs <- study_dirs[
  gsub("study", "", basename(study_dirs)) %in% names(STUDY_TRIAL_MAP)
]

for (d in study_dirs) {
  study_id <- gsub("study", "", basename(d))

  # recursive=TRUE to handle nested subfolder (study10673/study10673/...)
  csvs  <- list.files(d, pattern = "^(primary|secondary)_phenotypes_.*\\.csv$",
                      full.names = TRUE, recursive = TRUE)
  n_pri <- length(list.files(d, pattern = "^primary_phenotypes_",   recursive = TRUE))
  n_sec <- length(list.files(d, pattern = "^secondary_phenotypes_", recursive = TRUE))
  cat(sprintf("  study%s: primary=%d, secondary=%d CSVs\n", study_id, n_pri, n_sec))

  if (length(csvs) == 0) {
    cat("    No CSV files found. Please run 00_download_training_phenotypes.R first.\n")
    next
  }

  df <- bind_rows(lapply(csvs, function(f) {
    tryCatch({
      tmp <- fread(f, data.table = FALSE)
      tmp$source_type <- ifelse(grepl("^primary", basename(f)), "primary", "secondary")
      tmp
    }, error = function(e) NULL)
  }))
  if (nrow(df) == 0) next

  df <- df %>%
    filter(grepl("grain_yield|Grain yield|CO_321:0001218",
                 observation_variable_name, ignore.case = TRUE)) %>%
    mutate(yield = suppressWarnings(as.numeric(value))) %>%
    filter(!is.na(yield), is.finite(yield)) %>%
    mutate(
      focal_study_id = study_id,
      trial_name     = STUDY_TRIAL_MAP[study_id],
      study_db_id    = as.integer(study_db_id)
    ) %>%
    select(focal_study_id, trial_name, study_db_id, source_type,
           germplasm_name, rep, block, yield)

  all_pheno[[study_id]] <- df
  cat(sprintf("    -> %d records, %d accessions, from %d training studies\n",
              nrow(df), n_distinct(df$germplasm_name), n_distinct(df$study_db_id)))
}

all_pheno_df <- bind_rows(all_pheno)
cat(sprintf("\nTotal: %d records, %d training studies, %d accessions\n",
            nrow(all_pheno_df),
            n_distinct(all_pheno_df$study_db_id),
            n_distinct(all_pheno_df$germplasm_name)))

# ── Compute BLUEs per training study ─────────────────────
cat("\nComputing BLUEs...\n")

compute_blues <- function(df) {
  # Outlier removal: |z| > 3
  df <- df %>%
    mutate(z = abs(yield - mean(yield, na.rm = TRUE)) /
                    sd(yield, na.rm = TRUE)) %>%
    filter(z < 3 | is.na(z)) %>%
    select(-z)

  tryCatch({
    has_rep   <- n_distinct(df$rep,   na.rm = TRUE) > 1
    has_block <- n_distinct(df$block, na.rm = TRUE) > 1

    formula <- if (has_rep && has_block) yield ~ germplasm_name + rep + block
               else if (has_rep)         yield ~ germplasm_name + rep
               else                      yield ~ germplasm_name

    fit   <- lm(formula, data = df)
    coefs <- coef(fit)

    intercept  <- coefs["(Intercept)"]
    germ_coefs <- coefs[grepl("^germplasm_name", names(coefs))]
    names(germ_coefs) <- gsub("^germplasm_name", "", names(germ_coefs))

    base_germ <- setdiff(unique(df$germplasm_name), names(germ_coefs))[1]
    blues <- c(setNames(intercept, base_germ), intercept + germ_coefs)

    tibble(germplasmName = names(blues), BLUE = as.numeric(blues))

  }, error = function(e) {
    # Fallback: per-accession mean if lm fails
    df %>%
      group_by(germplasm_name) %>%
      summarise(BLUE = mean(yield, na.rm = TRUE), .groups = "drop") %>%
      rename(germplasmName = germplasm_name)
  })
}

pheno_blues <- all_pheno_df %>%
  group_by(focal_study_id, trial_name, study_db_id, source_type) %>%
  group_modify(~ compute_blues(.x)) %>%
  ungroup()

cat(sprintf("Done: %d BLUE records\n", nrow(pheno_blues)))

# ── Save ──────────────────────────────────────────────────
saveRDS(all_pheno_df,  file.path(OUTPUT_DIR, "pheno_raw.rds"))
saveRDS(pheno_blues,   file.path(OUTPUT_DIR, "pheno_blues.rds"))
write_csv(pheno_blues, file.path(OUTPUT_DIR, "pheno_blues.csv"))

cat(sprintf("\nSaved to: %s\n", OUTPUT_DIR))
cat("  pheno_raw.rds\n")
cat("  pheno_blues.rds / pheno_blues.csv\n")
cat("  Training_Trial_Info.rds  (from GitHub)\n")
cat("\nNext step: run 03_genomic_prediction.R\n")
