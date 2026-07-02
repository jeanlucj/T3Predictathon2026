# ============================================================
# 01_vcf_to_snp_matrix.R
# Convert all VCF files to SNP matrices
# - Auto-deduplication: prefer .vcf.gz over .vcf for same base name
# - Supports both .vcf and .vcf.gz
# - Already-processed files are skipped automatically
# ============================================================

library(vcfR)
library(tidyverse)

here::i_am("scripts/Prediction1/01_vcf_to_snp_matrix_final.R")

VCF_DIR    <- here::here("scripts", "Prediction1", "data", "Genotype")
OUTPUT_DIR <- here::here("scripts", "Prediction1", "data", "SNP_matrices")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

MAF_THRESHOLD  <- 0.05
MISS_THRESHOLD <- 0.20

# ── Deduplication: prefer .vcf.gz when both versions exist ─
all_files  <- list.files(VCF_DIR, pattern = "\\.vcf(\\.gz)?$",
                         full.names = TRUE, ignore.case = TRUE)
base_names <- tools::file_path_sans_ext(
               tools::file_path_sans_ext(basename(all_files)))

vcf_df <- tibble(path = all_files, base = base_names) %>%
  mutate(is_gz = grepl("\\.gz$", path, ignore.case = TRUE)) %>%
  group_by(base) %>%
  arrange(desc(is_gz)) %>%   # .gz first
  slice(1) %>%               # keep only one per base name
  ungroup()

cat(sprintf("Found %d VCF files; processing %d after deduplication\n\n",
            length(all_files), nrow(vcf_df)))

# ── Core function: VCF -> SNP matrix (additive -1/0/1 encoding) ──
process_vcf <- function(vcf_path, label) {
  out_path <- file.path(OUTPUT_DIR, paste0(label, "_snp_matrix.rds"))

  if (file.exists(out_path)) {
    cat(sprintf("  Skipping (already done): %s\n", label))
    return(invisible(out_path))
  }

  cat(sprintf("  Processing: %s\n", basename(vcf_path)))

  vcf <- tryCatch(
    read.vcfR(vcf_path, verbose = FALSE),
    error = function(e) {
      cat(sprintf("     Failed to read: %s\n", e$message))
      return(NULL)
    }
  )
  if (is.null(vcf)) return(NULL)

  # Deduplicate SNPs by CHROM:POS key
  # (Using the ID column is wrong: most VCFs have "." for all IDs,
  #  which would cause all but the first SNP to be dropped)
  snp_key <- paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], sep = ":")
  vcf     <- vcf[!duplicated(snp_key), ]

  cat(sprintf("     Raw: %d SNPs x %d samples\n",
              nrow(vcf@gt), ncol(vcf@gt) - 1))

  # Vectorized genotype encoding (replaces slow nested for+sapply loop)
  # extract.gt(as.numeric=TRUE) returns 0/1/2 (number of ALT alleles)
  # subtract 1 to get additive encoding: -1/0/1
  gt_num <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
  gt_num <- gt_num - 1L

  # Transpose to samples x SNPs
  geno <- t(gt_num)

  # Filter samples with missingness > threshold
  geno <- geno[rowMeans(is.na(geno)) <= MISS_THRESHOLD, , drop = FALSE]

  # Filter SNPs by missingness and MAF
  snp_miss <- colMeans(is.na(geno))
  freq     <- (colMeans(geno, na.rm = TRUE) + 1) / 2
  maf      <- pmin(freq, 1 - freq)
  keep     <- snp_miss <= MISS_THRESHOLD & !is.na(maf) & maf >= MAF_THRESHOLD
  geno     <- geno[, keep, drop = FALSE]

  # Impute missing values with rounded per-SNP mean
  # (rounding preserves integer encoding and avoids bias in downstream
  #  round(snp + 1) conversion in build_grm)
  col_means_rounded <- round(colMeans(geno, na.rm = TRUE))
  for (j in seq_len(ncol(geno))) {
    na_idx <- is.na(geno[, j])
    if (any(na_idx)) geno[na_idx, j] <- col_means_rounded[j]
  }

  cat(sprintf("     After filtering: %d samples x %d SNPs\n", nrow(geno), ncol(geno)))
  saveRDS(geno, out_path)
  cat(sprintf("     Saved.\n"))
  return(invisible(out_path))
}

# ── Batch processing ──────────────────────────────────────
cat("Starting...\n")
cat(rep("=", 50), "\n", sep = "")

success <- 0
fail    <- 0

for (i in seq_len(nrow(vcf_df))) {
  cat(sprintf("\n[%d/%d] ", i, nrow(vcf_df)))
  label  <- gsub("[^A-Za-z0-9_-]", "_", vcf_df$base[i])
  result <- process_vcf(vcf_df$path[i], label)
  if (!is.null(result)) success <- success + 1
  else                   fail    <- fail + 1
}

# ── Summary ───────────────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 50)))
done <- list.files(OUTPUT_DIR, pattern = "_snp_matrix\\.rds$")
cat(sprintf("Completed: %d SNP matrices\n", length(done)))
cat(sprintf("Failed:    %d\n", fail))
cat(sprintf("Output:    %s\n", OUTPUT_DIR))
cat("\nNext step: run 02_prepare_phenotypes.R\n")
