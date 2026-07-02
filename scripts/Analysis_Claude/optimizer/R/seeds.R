# seeds.R
#
# The five Predictathon submissions, expressed in the six-subtask genome. These
# seed the optimizer so the submitted algorithms are an explicit baseline the
# search must beat, and so crossover has good building blocks from iteration one.
# Each seed sets the methods and the parameters that apply to them, plus the
# always-on parameters; repair_config() then NA-fills anything inapplicable and
# enforces canonical key order. These are faithful-but-approximate encodings: the
# point is to start the search near the submissions, not to reproduce them
# bit-for-bit.

library(tidyverse)

# Build a seed from a sparse list of overrides on top of safe always-on defaults.
.make_seed <- function(overrides) {
  base <- list(
    pheno_prep.ge_weighting = "none",
    pheno_prep.ge_bandwidth = 1,
    pheno_prep.standardize  = "none",
    geno_select.marker_thin = "1",
    kernel.maf              = 0.05,
    kernel.max_missing      = 0.20,
    kernel.ridge            = 1e-4,
    kernel.impute           = "mean_round",
    predict_post.blend_obs_w = 0,
    predict_post.min_overlap = 10
  )
  cfg <- modifyList(base, overrides)
  repair_config(cfg)
}

seed_configs <- function() {
  list(
    # P1 (GBLUP, rrBLUP): overlap-tier training trials; per-trial BLUE then
    # average; VanRaden GRM on merged markers; rrBLUP GBLUP; conditional
    # expectation; standardized BLUEs; mean fallback below 10 overlap.
    Prediction1 = .make_seed(list(
      train_select.method       = "accession_overlap",
      train_select.primary_min  = 4,
      train_select.secondary_min = 12,
      train_select.primary_only = "no",
      pheno_prep.method         = "blue_lm",
      pheno_prep.z_thr          = 3,
      pheno_prep.standardize    = "per_trial_z",
      geno_select.method        = "focal_plus_onehop",
      geno_select.min_bridge    = 1,
      kernel.method             = "vanRaden_single",
      model.method              = "gblup_rrblup",
      model.lambda_select       = "fixed",
      predict_post.method       = "cond_expectation"
    )),

    # P2 (end-to-end GBLUP, global union GRM): broad training set; BLUE; all
    # projects -> one global VanRaden GRM; rrBLUP GBLUP; direct BLUP.
    Prediction2 = .make_seed(list(
      train_select.method       = "accession_overlap",
      train_select.primary_only = "yes",
      train_select.primary_min  = 2,
      train_select.secondary_min = 12,
      pheno_prep.method         = "blue_lm",
      pheno_prep.z_thr          = 3,
      geno_select.method        = "all_projects",
      kernel.method             = "vanRaden_single",
      kernel.maf                = 0.01,
      kernel.max_missing        = 0.50,
      model.method              = "gblup_rrblup",
      model.lambda_select       = "fixed",
      predict_post.method       = "direct_blup"
    )),

    # P3 (multi-kernel GBLUP, EM-combined GRM): top-k genomically+environmentally
    # similar trials; BLUE; per-project GRMs combined by Wishart-EM; sommer G+E
    # multi-kernel; direct BLUP.
    Prediction3 = .make_seed(list(
      train_select.method       = "top_k_similar",
      train_select.k            = 15,
      train_select.similarity   = "both",
      pheno_prep.method         = "blue_lm",
      pheno_prep.z_thr          = 3,
      geno_select.method        = "focal_plus_onehop",
      geno_select.min_bridge    = 1,
      kernel.method             = "em_combine",
      model.method              = "gblup_sommer_GE",
      model.include_E           = "yes",
      predict_post.method       = "direct_blup"
    )),

    # P4 (environmentally-weighted GBLUP with blending): trial-centered,
    # G×E-weighted phenotypes; VanRaden GRM; LOO-ridge GBLUP; blend observed BLUE
    # with the genomic prediction.
    Prediction4 = .make_seed(list(
      train_select.method       = "accession_overlap",
      train_select.primary_only = "yes",
      train_select.primary_min  = 2,
      train_select.secondary_min = 12,
      pheno_prep.method         = "trial_center",
      pheno_prep.z_thr          = 3,
      pheno_prep.ge_weighting   = "env_gaussian",
      pheno_prep.ge_bandwidth   = 1,
      geno_select.method        = "all_projects",
      kernel.method             = "vanRaden_single",
      kernel.maf                = 0.01,
      kernel.max_missing        = 0.50,
      model.method              = "gblup_loo_ridge",
      model.lambda_select       = "loo",
      predict_post.method       = "direct_blup",
      predict_post.blend_obs_w  = 0.70
    )),

    # P5 (GBLUP/RKHS with environmental covariates): same-program training trials;
    # BLUE; EM-combined GRM across genotyping projects; RKHS with environmental
    # effects; direct BLUP.
    Prediction5 = .make_seed(list(
      train_select.method       = "same_program",
      train_select.prog_cap     = 20,
      pheno_prep.method         = "blue_lm",
      pheno_prep.z_thr          = 3,
      geno_select.method        = "focal_plus_onehop",
      geno_select.min_bridge    = 2,
      kernel.method             = "em_combine",
      model.method              = "rkhs",
      model.include_E           = "yes",
      predict_post.method       = "direct_blup"
    ))
  )
}
