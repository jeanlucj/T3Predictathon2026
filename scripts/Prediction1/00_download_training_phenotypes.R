# ============================================================
# 00_download_training_phenotypes.R
# Download all historical training phenotype data from T3/Wheat via BrAPI
#
# Required packages (not on CRAN, install manually):
#   remotes::install_github("IntegratedBreedingPlatform/BrAPI")
#   remotes::install_github("jeanlucj/T3BrapiHelpers")
#
# Prerequisite:
#   Training_Trial_Info.rds must have been downloaded by script 02
#
# Estimated runtime: 30-60 min per focal trial, several hours total
# To run a subset, set RUN_ONLY_STUDY_IDS below
# ============================================================

library(tidyverse)
library(BrAPI)
library(T3BrapiHelpers)
library(purrr)

here::i_am("scripts/Prediction1/00_download_training_phenotypes.R")

PHENO_DIR       <- here::here("scripts", "Prediction1", "data", "Phenotype")
TRIAL_INFO_FILE <- here::here("scripts", "Prediction1", "data",
                              "Phenotype_processed", "Training_Trial_Info.rds")
dir.create(PHENO_DIR, recursive = TRUE, showWarnings = FALSE)

# To run specific focal trials only, provide study IDs here
# e.g. RUN_ONLY_STUDY_IDS <- c("10679", "10680", "10681")
RUN_ONLY_STUDY_IDS <- c()

# ── Check dependencies ────────────────────────────────────
for (pkg in c("BrAPI", "T3BrapiHelpers", "janitor")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "Missing package %s. Install with:\n  remotes::install_github('%s')\n",
      pkg,
      switch(pkg,
             "BrAPI"          = "IntegratedBreedingPlatform/BrAPI",
             "T3BrapiHelpers" = "jeanlucj/T3BrapiHelpers",
             pkg)
    ))
  }
}

# ── Load Training_Trial_Info ──────────────────────────────
if (!file.exists(TRIAL_INFO_FILE))
  stop("Training_Trial_Info.rds not found. Please run 02_prepare_phenotypes.R first.")

training_info <- readRDS(TRIAL_INFO_FILE)
cat(sprintf("Loaded training info for %d focal trials\n\n", length(training_info)))

if (length(RUN_ONLY_STUDY_IDS) > 0) {
  training_info <- training_info[
    sapply(training_info, `[[`, "study_id") %in% RUN_ONLY_STUDY_IDS
  ]
  cat(sprintf("Running %d specified focal trials only\n\n", length(training_info)))
}

# ── Connect to T3/Wheat ───────────────────────────────────
cat("Connecting to T3/Wheat (enter credentials in the login window)...\n")
wheat_conn <- BrAPI::getBrAPIConnection("T3/Wheat")
cat("Connected successfully\n\n")

# ── Download yield phenotypes for a single study ──────────
# observationVariableDbIds=84527 is the T3 variable ID for Grain yield
save_yield_phenotypes <- function(study_id, prefix) {
  file_name_yield <- file.path(PHENO_DIR,
                               paste0(prefix, "_phenotypes_", study_id, ".csv"))

  if (file.exists(file_name_yield)) {
    cat(sprintf("    Skipping (already exists): %s\n", basename(file_name_yield)))
    return(file_name_yield)
  }

  extract_obs_unit_coord <- function(search_res) {
    plot_column <- search_res$observationUnitPosition$positionCoordinateX %||% NA_integer_
    plot_row    <- search_res$observationUnitPosition$positionCoordinateY %||% NA_integer_
    obs_tags <- search_res$observationUnitPosition$observationLevelRelationships |>
      bind_rows() |>
      dplyr::select(-levelOrder) |>
      tidyr::pivot_wider(names_from = levelName, values_from = levelCode) |>
      dplyr::mutate(row = plot_row, column = plot_column,
                    observationUnitDbId = search_res$observationUnitDbId)
    return(obs_tags)
  }

  yield_phenotypes <- tryCatch({
    wheat_conn$search(
      "observations",
      body = list(studyDbIds             = c(study_id),
                  observationVariableDbIds = c(84527))
    )$combined_data |>
      purrr::map(function(sl) { sl$season <- sl$season[[1]]; sl }) |>
      dplyr::bind_rows()
  }, error = function(e) {
    cat(sprintf("    study %s request failed: %s\n", study_id, e$message))
    return(data.frame())
  })

  if (nrow(yield_phenotypes) == 0) {
    cat(sprintf("    study %s: no yield data\n", study_id))
    return(NULL)
  }

  yield_phenotypes <- yield_phenotypes |>
    dplyr::select(studyDbId, season,
                  germplasmDbId, germplasmName,
                  observationUnitDbId, observationUnitName,
                  observationVariableName, value)

  resp <- wheat_conn$search(
    "observationunits",
    body = list(observationUnitDbIds = yield_phenotypes$observationUnitDbId)
  )
  plot_coords <- lapply(resp$combined_data, extract_obs_unit_coord) |>
    dplyr::bind_rows()

  yield_phenotypes <- merge(yield_phenotypes, plot_coords) |>
    dplyr::select(studyDbId, season,
                  germplasmDbId, germplasmName,
                  observationUnitDbId, observationUnitName,
                  rep, block, plot, row, column,
                  observationVariableName, value) |>
    dplyr::arrange(plot) |>
    janitor::clean_names()

  readr::write_csv(yield_phenotypes, file_name_yield)
  cat(sprintf("    study %s: %d records saved\n", study_id, nrow(yield_phenotypes)))
  return(file_name_yield)
}

# ── Download all training studies for one focal trial ─────
download_training_trials <- function(info) {
  study_id <- info$study_id
  subdir   <- file.path(PHENO_DIR, paste0("study", study_id))
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

  primary_ids <- info$primary_study_db_ids[!is.na(info$primary_study_db_ids)]
  cat(sprintf("  [primary] downloading %d studies...\n", length(primary_ids)))
  purrr::map(primary_ids, purrr::safely(save_yield_phenotypes),
             prefix = paste0("study", study_id, "/primary"), .progress = TRUE)

  secondary_ids <- info$secondary_study_db_ids[!is.na(info$secondary_study_db_ids)]
  if (length(secondary_ids) > 0) {
    cat(sprintf("  [secondary] downloading %d studies...\n", length(secondary_ids)))
    purrr::map(secondary_ids, purrr::safely(save_yield_phenotypes),
               prefix = paste0("study", study_id, "/secondary"), .progress = TRUE)
  }

  invisible(NULL)
}

# ── Main loop ─────────────────────────────────────────────
cat(strrep("=", 55), "\n")
cat("Starting download (~30-60 min per focal trial)\n")
cat(strrep("=", 55), "\n")

for (i in seq_along(training_info)) {
  info     <- training_info[[i]]
  study_id <- info$study_id
  n_pri    <- length(info$primary_study_db_ids[!is.na(info$primary_study_db_ids)])
  n_sec    <- length(info$secondary_study_db_ids[!is.na(info$secondary_study_db_ids)])
  cat(sprintf("\n[%d/%d] focal study%s  primary=%d, secondary=%d\n",
              i, length(training_info), study_id, n_pri, n_sec))
  cat(strrep("-", 40), "\n")
  tryCatch(download_training_trials(info),
           error = function(e) cat(sprintf("  Failed: %s\n", e$message)))
}

# ── Summary ───────────────────────────────────────────────
cat(sprintf("\n%s\nDownload complete. File counts:\n", strrep("=", 55)))
for (info in training_info) {
  sid    <- info$study_id
  subdir <- file.path(PHENO_DIR, paste0("study", sid))
  n_p    <- length(list.files(subdir, pattern = "^primary_phenotypes_.*\\.csv$",   recursive = TRUE))
  n_s    <- length(list.files(subdir, pattern = "^secondary_phenotypes_.*\\.csv$", recursive = TRUE))
  cat(sprintf("  study%s: primary=%d, secondary=%d files\n", sid, n_p, n_s))
}
cat("\nNext step: run 02_prepare_phenotypes.R\n")
