# reliability_pareto.R
#
# For each algorithm (Prediction1..5) and each of the 9 focal trials, reconstruct
# the relationship matrix (GRM/kernel) THAT ALGORITHM would build given its own
# parameters, then score how well its training-trial / accession selection plus
# its relationship estimation support prediction of the focal trial.
#
# Metric (per focal accession j, "non-normalized reliability", the GBLUP /
# selection-index quantity the user worked with for Prediction3):
#     r2_j = c' P^-1 c ,   P = Var(x) = K_tt + lambda I ,   c = Cov(x, G_j) = K_{t,j}
# with sigma2_g = 1 and lambda = sigma2_e/sigma2_g a COMMON variance ratio across
# all algorithms (so the score isolates SELECTION + RELATIONSHIP ESTIMATION, not
# per-algorithm variance-component fitting). lambda is set from a nominal
# heritability h2 via lambda = (1 - h2)/h2; we report h2 = 0.5 (primary) and
# 0.2 / 0.8 for sensitivity. Per-trial summary is the square-root mean reliability
#     sqrt( mean_j r2_j )   over a COMMON focal set (genotyped focal accessions),
# so every algorithm is scored on the same j's within a trial.
#
# The two Pareto axes per (trial, scheme): HIGH sqrt-mean-reliability (good),
# SMALL declared training-population size (good). An algorithm is on the Pareto
# frontier if no other algorithm has reliability >= and size <= (one strict). We
# tally frontier membership across the 9 trials per algorithm.
#
# GRM per algorithm == its own relationship-estimation method, built OFFLINE from
# the cached per-project GRMs (no BrAPI):
#   * vanRaden_single (P1, P2, P4): a single merged VanRaden GRM. The cached
#     projects use largely DISJOINT marker sets (distinct genotyping platforms),
#     so a real union-of-VCFs merge estimates within-platform relationships and
#     drives cross-platform pairs to the prior (0) under mean-imputation -- i.e. a
#     BLOCK-DIAGONAL assembly of the per-project GRMs. No cross-platform bridging.
#   * em_combine     (P3, P5): per-project GRM -> Wishart-EM combine, which DOES
#     bridge platforms through shared accessions (recovers cross-platform
#     relationships). This single difference -- bridge vs no-bridge -- is the
#     relationship-estimation quality the analysis isolates, alongside selection.
# Faithful-but-approximate (same philosophy as optimizer/R/seeds.R): the dominant
# levers -- which accessions/projects get connected, and single-merge vs EM-combine
# -- are reproduced. We do NOT re-derive each team's exact VCF, and the per-project
# GRMs use rrBLUP::A.mat defaults, so per-algorithm MAF / max-missing thresholds
# (a minor, second-order lever) are not differentiated.
#
# PROJECT UNIVERSE (expanded 2026-06-16). A project is included for an algorithm
# iff it genotyped one of that algorithm's needed accessions. The available
# universe is the FULL BrAPI per-accession genotyping-project map (76 projects),
# not just Prediction3's 31-project selection: expand_project_universe.R maps the
# 19,624-accession universe to its real genotyping projects and fetch_new_projects.R
# downloads + caches the additional ones, so each algorithm's projects are no
# longer bounded by P3's choice. A second pass (reextract_original31_full.R,
# 2026-06-17) re-extracted the original 31 projects at FULL sample coverage from
# the on-disk raw VCFs -- build_grm_for_cv00.R had sample-restricted the large
# ones to P3's accessions, hiding ~3,200 genotyped training accessions. Residual
# caveats: 12 covering projects had no usable archived VCF; two projects (8159,
# 10670) store samples under preliminary line names (GTC_*/F6*) that were since
# captured as germplasm SYNONYMS in T3, so the synonym-canonicalization step
# (build_synonym_map.R) now recovers them; and Prediction1/4 still sit at ~52%
# UNCOVERED because most of
# their "all-with-yield" lines are genuinely ungenotyped in any T3 project (a T3
# diagnostic: 500 uncovered P1 accessions link to only 4 genotyping projects).
# Effect on results: the expansion alone left the Pareto tally unchanged, but the
# FULL re-extraction raised Prediction1/4's YT_Urb reliability (0.821 -> 0.849,
# above Prediction5), changing the CV0 tally to Pred1 4, Pred2 8, Pred3 2,
# Pred4 4, Pred5 5 (was Pred1 3 / Pred4 2). Prior-state results are kept as
# output/reliability/*_31proj.csv (original) and *_subset31.csv (expanded, pre-
# full-extraction). See writing/notes/Algorithm_genotyping_project_matching.md.
#
# LOCATION: Claude-generated; lives in scripts/Analysis_Claude/. READS the five
# submissions' per-trial accession/trial lists from T3_predictathon_scripts/ and
# the dosage/GRM caches from ./output/; WRITES only under ./output/reliability/.
# The dosage/GRM caches are populated OFFLINE before this runs: build_grm_for_cv00.R
# (Prediction3's 31 projects) plus expand_project_universe.R + fetch_new_projects.R
# (the 33 additional fetched projects). This script itself needs no BrAPI.
# It does NOT touch any Prediction1..5 participant folder. Resumable: per-algorithm
# GRMs are cached; rerun picks up where it stopped. Run from this directory.

library(tidyverse)

here::i_am("reliability_pareto.R")
source(here::here("em_covariance_combiner.R"))

dosage_dir <- here::here("output", "dosage_cache")
grm_dir    <- here::here("output", "grm_cache")
# Each algorithm's submission lives in scripts/PredictionN/submission/ (the
# anonymized copies); resolve the per-algorithm submission dir.
pred_submission <- function(pred) here::here("..", pred, "submission")
# Genotype VCFs are large local data, not version-controlled; point vcf_dir at
# wherever you keep them locally (this script needs them to run).
vcf_dir    <- here::here("..", "..", "T3_predictathon_scripts", "data",
                         "recommended_vcfs")
out_dir    <- here::here("output", "reliability")
grm_cache_out <- file.path(out_dir, "grm")
dir.create(grm_cache_out, showWarnings = FALSE, recursive = TRUE)

# Synonym-aware matching: the recommended focal VCF read in focal_common() may
# carry an accession under a now-demoted preliminary/synonym name. Canonicalize
# its sample IDs to the primary name so they intersect covered_acc (whose
# rownames are already primary once the dosage caches are rebuilt). Fail-soft to
# identity if the map is absent; build it with build_synonym_map.R.
map_path     <- file.path(out_dir, "synonym_map.rds")
alias_lookup <- if (file.exists(map_path)) readRDS(map_path) else character()
canon <- function(s) { hit <- alias_lookup[s]; ifelse(is.na(hit), s, unname(hit)) }

# --- algorithm <-> submission folder, and each algorithm's GRM parameters -----
# (kernel.method / MAF / max_missing taken from optimizer/R/seeds.R)
ALGOS <- tribble(
  ~alg,    ~submission,     ~kernel,           ~maf,  ~max_missing,
  "Pred1", "Prediction1",   "vanRaden_single", 0.05,  0.20,
  "Pred2", "Prediction2",   "vanRaden_single", 0.01,  0.50,
  "Pred3", "Prediction3",   "em_combine",      0.05,  0.20,
  "Pred4", "Prediction4",   "vanRaden_single", 0.01,  0.50,
  "Pred5", "Prediction5",   "em_combine",      0.05,  0.20
)

TRIALS <- list.dirs(pred_submission("Prediction1"),
                    full.names = FALSE, recursive = FALSE)
H2 <- c(h2_0.2 = 0.2, h2_0.5 = 0.5, h2_0.8 = 0.8)  # lambda = (1-h2)/h2

# ---------------------------------------------------------------------------
# Project caches: per-project dosage matrices and per-project A.mat GRMs, keyed
# by project id (digits after dosage_/grm_). Loaded lazily; rownames give the
# accessions each project genotyped.
# ---------------------------------------------------------------------------
dosage_files <- list.files(dosage_dir, pattern = "^dosage_.*[.]rds$", full.names = TRUE)
proj_ids     <- stringr::str_extract(basename(dosage_files), "(?<=dosage_)\\d+")
names(dosage_files) <- proj_ids

.proj_acc_env <- new.env(parent = emptyenv())   # id -> accession names (cheap)
proj_accessions <- function(id) {
  if (is.null(.proj_acc_env[[id]]))
    .proj_acc_env[[id]] <- rownames(readRDS(dosage_files[[id]]))
  .proj_acc_env[[id]]
}
proj_dosage <- function(id) readRDS(dosage_files[[id]])
proj_grm <- function(id) {                       # cached full-project A.mat
  f <- file.path(grm_dir, paste0("grm_", id, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  g <- rrBLUP::A.mat(proj_dosage(id)); saveRDS(g, f); g
}
# accession -> sample size of its largest covering project (EM degrees of freedom)
PROJ_N <- vapply(proj_ids, function(id) length(proj_accessions(id)), integer(1))

# Read all project accession lists once (small) so we can find covering projects.
PROJ_ACC <- lapply(proj_ids, proj_accessions); names(PROJ_ACC) <- proj_ids
covered_acc <- unique(unlist(PROJ_ACC))

# ---------------------------------------------------------------------------
# Submission readers
# ---------------------------------------------------------------------------
rd_col <- function(alg, trial, file, col) {
  sub <- ALGOS$submission[ALGOS$alg == alg]
  p <- file.path(pred_submission(sub), trial, file)
  if (!file.exists(p)) return(character())
  v <- readr::read_csv(p, show_col_types = FALSE)[[col]]
  as.character(v)
}

# Common focal set for a trial = genotyped focal accessions from the recommended
# VCF (same j's for every algorithm). Cache the VCF sample names.
.vcf_env <- new.env(parent = emptyenv())
focal_common <- function(trial) {
  if (is.null(.vcf_env[[trial]])) {
    con <- file(file.path(vcf_dir, paste0(trial, ".vcf")), "r"); on.exit(close(con))
    samp <- character()
    repeat {
      l <- readLines(con, 1L); if (!length(l)) break
      if (startsWith(l, "#CHROM")) { samp <- canon(strsplit(l, "\t")[[1]][-(1:9)]); break }
    }
    .vcf_env[[trial]] <- intersect(samp, covered_acc)
  }
  .vcf_env[[trial]]
}

# ---------------------------------------------------------------------------
# Build the algorithm's GRM over genotyped(need): its own kernel method + QC.
#   need = union(genotyped training accessions, common focal set)
# Returns a relationship matrix K (rows/cols = accessions), scaled to mean diag~1.
# ---------------------------------------------------------------------------
std_grm <- function(g) g / mean(diag(g))

qc_markers <- function(X, maf, max_missing) {
  X <- X[, colMeans(is.na(X)) <= max_missing, drop = FALSE]
  af  <- colMeans(X, na.rm = TRUE) / 2
  X <- X[, pmin(af, 1 - af) >= maf, drop = FALSE]
  if (ncol(X) < 50) return(NULL)
  for (j in which(colMeans(is.na(X)) > 0)) {
    mu <- mean(X[, j], na.rm = TRUE); X[is.na(X[, j]), j] <- round(mu)
  }
  X
}

build_algo_grm <- function(alg, trial, cv) {
  cf <- file.path(grm_cache_out, paste0(alg, "_", trial, "_", cv, ".rds"))
  if (file.exists(cf)) return(readRDS(cf))

  par <- ALGOS[ALGOS$alg == alg, ]
  train_acc <- rd_col(alg, trial, paste0(cv, "_Accessions.csv"), "germplasmName")
  need <- union(intersect(train_acc, covered_acc), focal_common(trial))
  relevant <- proj_ids[vapply(PROJ_ACC, function(a) any(need %in% a), logical(1))]
  if (!length(relevant)) { saveRDS(NULL, cf); return(NULL) }

  if (par$kernel == "vanRaden_single") {
    # Faithful single-merge ("union of VCFs" -> one VanRaden GRM). The cached
    # projects use largely DISJOINT marker sets (different genotyping platforms;
    # many pairwise marker overlaps are 0 -- see the marker-overlap diagnostic),
    # so a real union merge estimates WITHIN-platform relationships while
    # CROSS-platform pairs are driven to the prior (0) by mean-imputation of the
    # other platform's sites. That is exactly a BLOCK-DIAGONAL assembly of the
    # per-project (standardized) GRMs. Unlike em_combine it does NOT bridge
    # platforms through shared accessions -- the core relationship-estimation
    # difference this analysis measures. Every needed accession (incl. all focal
    # lines) is placed, so a focal line with no same-platform training simply gets
    # reliability 0 (the algorithm's mean fallback), not a dropped/NA score.
    proj_of <- new.env(parent = emptyenv())          # accession -> best project id
    for (id in relevant) for (a in intersect(PROJ_ACC[[id]], need)) {
      if (is.null(proj_of[[a]]) || PROJ_N[[id]] > PROJ_N[[proj_of[[a]]]]) proj_of[[a]] <- id
    }
    accs <- ls(proj_of)
    if (!length(accs)) { saveRDS(NULL, cf); return(NULL) }
    K <- matrix(0, length(accs), length(accs), dimnames = list(accs, accs))
    by_proj <- split(accs, vapply(accs, function(a) proj_of[[a]], character(1)))
    for (id in names(by_proj)) {
      g <- std_grm(proj_grm(id))
      a <- intersect(by_proj[[id]], rownames(g))
      K[a, a] <- g[a, a]
    }

  } else {  # em_combine
    parts <- list(); dfs <- integer(0)
    for (id in relevant) {
      g <- proj_grm(id); keep <- rownames(g) %in% need
      if (sum(keep) > 0) {
        parts[[length(parts) + 1L]] <- std_grm(g[keep, keep, drop = FALSE])
        dfs <- c(dfs, PROJ_N[[id]])
      }
    }
    if (!length(parts)) { saveRDS(NULL, cf); return(NULL) }
    if (length(parts) == 1L) {
      K <- parts[[1]]
    } else {
      combined_names <- unique(unlist(lapply(parts, colnames)))
      idx <- lapply(parts, function(g) match(colnames(g), combined_names))
      res <- EMCovarianceCombiner(partial_covs = parts, var_indices = idx,
                                  degrees_freedom = dfs)
      K <- res$psi[-1, -1]
      rownames(K) <- colnames(K) <- combined_names
    }
  }
  K <- std_grm(K)
  saveRDS(K, cf)
  K
}

# ---------------------------------------------------------------------------
# Reliability r2_j = K_{j,t} (K_tt + lambda I)^-1 K_{t,j} for every common focal
# j, vectorized: solve once, then column sums. Returns sqrt(mean_j r2_j) for each
# heritability, plus bookkeeping (declared size, genotyped/connected counts).
# ---------------------------------------------------------------------------
reliability <- function(alg, trial, cv) {
  train_acc_all <- rd_col(alg, trial, paste0(cv, "_Accessions.csv"), "germplasmName")
  size_declared <- length(unique(train_acc_all))
  fc <- focal_common(trial)

  K <- tryCatch(build_algo_grm(alg, trial, cv), error = function(e) {
    message("  GRM fail ", alg, " ", trial, " ", cv, ": ", conditionMessage(e)); NULL })
  if (is.null(K)) {
    return(tibble(alg = alg, trial = trial, cv = cv, size_declared = size_declared,
                  n_focal = length(fc), n_train_geno = NA_integer_,
                  !!!setNames(as.list(rep(NA_real_, length(H2))),
                              paste0("sqrtrel_", names(H2)))))
  }

  knames <- rownames(K)
  t_acc <- intersect(intersect(train_acc_all, covered_acc), knames)  # genotyped training in K
  j_acc <- intersect(fc, knames)
  # Under CV0 a focal line may also sit in training via another trial; that is the
  # legitimate CV0 advantage, so we leave such j in t (the index handles it).
  if (length(t_acc) < 2 || !length(j_acc)) {
    return(tibble(alg = alg, trial = trial, cv = cv, size_declared = size_declared,
                  n_focal = length(j_acc), n_train_geno = length(t_acc),
                  !!!setNames(as.list(rep(NA_real_, length(H2))),
                              paste0("sqrtrel_", names(H2)))))
  }

  Ktt <- K[t_acc, t_acc, drop = FALSE]
  Ktj <- K[t_acc, j_acc, drop = FALSE]
  out <- list()
  for (nm in names(H2)) {
    lambda <- (1 - H2[[nm]]) / H2[[nm]]
    A <- Ktt; diag(A) <- diag(A) + lambda
    R <- tryCatch(solve(A, Ktj), error = function(e) MASS::ginv(A) %*% Ktj)
    r2 <- colSums(Ktj * R)               # r2_j for each focal j
    out[[paste0("sqrtrel_", nm)]] <- sqrt(mean(pmax(r2, 0)))
  }
  tibble(alg = alg, trial = trial, cv = cv, size_declared = size_declared,
         n_focal = length(j_acc), n_train_geno = length(t_acc), !!!out)
}

# ---------------------------------------------------------------------------
# Run: all algorithms x trials x schemes (or a subset of trials from the CLI).
# ---------------------------------------------------------------------------
args   <- commandArgs(trailingOnly = TRUE)
trials <- if (length(args)) intersect(args, TRIALS) else TRIALS

grid <- expand_grid(alg = ALGOS$alg, trial = trials, cv = c("CV0", "CV00"))
res <- pmap_dfr(grid, function(alg, trial, cv) {
  r <- reliability(alg, trial, cv)
  message(sprintf("  %-6s %-26s %-4s size=%d  nTrainGeno=%s  nFocal=%d  sqrtRel(h2=.5)=%s",
                  alg, substr(trial, 1, 26), cv, r$size_declared,
                  ifelse(is.na(r$n_train_geno), "NA", r$n_train_geno),
                  r$n_focal,
                  ifelse(is.na(r$sqrtrel_h2_0.5), "NA", formatC(r$sqrtrel_h2_0.5, digits = 3, format = "f"))))
  r
})

readr::write_csv(res, file.path(out_dir, "reliability_by_alg_trial.csv"))

# ---------------------------------------------------------------------------
# Pareto frontier per (trial, scheme): maximize sqrt-mean-reliability (h2=0.5),
# minimize declared size. Tally frontier membership per algorithm.
# ---------------------------------------------------------------------------
on_pareto <- function(rel, size) {
  ok <- !is.na(rel) & !is.na(size)
  res <- rep(NA, length(rel))
  res[ok] <- vapply(which(ok), function(i) {
    dominated <- any(ok & rel >= rel[i] & size <= size[i] &
                       (rel > rel[i] | size < size[i]))
    !dominated
  }, logical(1))
  res
}

pareto <- res |>
  group_by(trial, cv) |>
  mutate(pareto = on_pareto(sqrtrel_h2_0.5, size_declared)) |>
  ungroup()
readr::write_csv(pareto, file.path(out_dir, "reliability_pareto_flagged.csv"))

tally <- pareto |>
  filter(!is.na(pareto)) |>
  group_by(alg, cv) |>
  summarise(n_pareto = sum(pareto), n_trials = n(), .groups = "drop") |>
  pivot_wider(names_from = cv, values_from = c(n_pareto, n_trials))
readr::write_csv(tally, file.path(out_dir, "pareto_tally.csv"))

cat("\n================ Pareto frontier tally (out of", length(trials), "trials) ================\n")
print(as.data.frame(tally))
cat("\nWrote:\n  ", file.path(out_dir, "reliability_by_alg_trial.csv"),
    "\n  ", file.path(out_dir, "reliability_pareto_flagged.csv"),
    "\n  ", file.path(out_dir, "pareto_tally.csv"), "\n")
