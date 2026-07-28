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
    model.lambda_select     = "reml",
    model.lambda_fixed      = 1,
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

# The five submissions as configurations, reconciled against each submission's source code and
# writing/notes/algorithm_characterization_by_subtask.md. `scheme` matters because Prediction5
# submitted a DIFFERENT model per cross-validation scenario; the other four ran one algorithm
# under both.
#
# Residual gaps -- what the search space cannot express, so these seeds only approximate:
#   * P2 and P4 applied no trial-level exclusion at all. The closest expressible encoding is
#     accession_overlap with both tiers at their loosest thresholds. (For P4 this is in fact
#     faithful: it consumed the challenge-provided training packages, which are the same
#     germplasm-overlap selection P1's curated list was derived from.)
#   * P3 took the UNION of the top 15 genomically-similar and top 15 environmentally-similar
#     trials; `similarity = "both"` here ranks on the SUM of the two scores and takes the top k,
#     which is a different set. No union option exists.
#   * P3 chose genotyping projects by greedy set-cover over training-accession coverage; the
#     closest available method is focal_plus_onehop.
#   * P2 solved its ridge at lambda = 1e-5; `lambda_fixed` bottoms out at 1e-2, so the seed
#     sits at the floor and is more regularized than the submission.
#   * P5 fitted location and year as FIXED effects under CV0; no subtask-E method does that.
seed_configs <- function(scheme = c("CV0", "CV00")) {
  scheme <- match.arg(scheme)
  list(
    # P1 -- GBLUP (rrBLUP). Primary+secondary overlap tiers from the curated training-trial
    # list; per-study BLUEs after |z|>3 outlier removal, averaged; VanRaden GRM (AGHmatrix,
    # ridge 1e-4); mixed.solve on standardized BLUEs; conditional expectation for test lines.
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
      model.lambda_select       = "reml",     # mixed.solve estimates the ratio by REML
      predict_post.method       = "cond_expectation"
    )),

    # P2 -- end-to-end GBLUP on a global union GRM. NO trial-level selection (a broad
    # historical phenotype database plus the challenge packages), so the loosest tiers here.
    # Phenotypes are z-scored WITHIN each trial and used directly as the target -- no BLUE
    # step -- hence raw_mean + per_trial_z. Closed-form ridge at a fixed lambda. Predictions
    # are mu + G_um alpha, i.e. the conditional expectation, over all genotyped lines.
    Prediction2 = .make_seed(list(
      train_select.method       = "accession_overlap",
      train_select.primary_only = "no",
      train_select.primary_min  = 2,
      train_select.secondary_min = 4,
      pheno_prep.method         = "raw_mean",
      pheno_prep.standardize    = "per_trial_z",
      geno_select.method        = "all_projects",
      kernel.method             = "vanRaden_single",
      kernel.maf                = 0.01,
      kernel.max_missing        = 0.50,
      model.method              = "gblup_rrblup",
      model.lambda_select       = "fixed",
      model.lambda_fixed        = 0.01,       # submission used 1e-5; this is the range floor
      predict_post.method       = "cond_expectation"
    )),

    # P3 -- multi-kernel GBLUP on an EM-combined GRM. Top-15 genomic and top-15 environmental
    # similar trials; per-trial BLUEs; per-project GRMs combined by Wishart EM; sommer G+W;
    # breeding values read straight off the fitted object.
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

    # P4 -- environmentally weighted GBLUP with blending. Consumed the challenge training
    # packages whole (the same overlap-based selection behind P1's list, hence both tiers);
    # trial means removed, then environment-weighted per-accession means; focal VCF plus
    # supplementary training VCFs; LOO-tuned ridge; predictions for every line in the focal
    # VCF (conditional expectation), then blended with the observed BLUE where one exists.
    Prediction4 = .make_seed(list(
      train_select.method       = "accession_overlap",
      train_select.primary_only = "no",
      train_select.primary_min  = 4,
      train_select.secondary_min = 12,
      pheno_prep.method         = "trial_center",
      pheno_prep.z_thr          = 3,
      pheno_prep.ge_weighting   = "env_gaussian",
      pheno_prep.ge_bandwidth   = 1,
      geno_select.method        = "focal_plus_onehop",
      geno_select.min_bridge    = 1,
      kernel.method             = "vanRaden_single",
      kernel.maf                = 0.01,
      kernel.max_missing        = 0.50,
      model.method              = "gblup_rrblup",
      model.lambda_select       = "loo",      # log-grid leave-one-out search
      predict_post.method       = "cond_expectation",
      predict_post.blend_obs_w  = 0.70
    )),

    # P5 -- SCHEME-DEPENDENT, and the reason this function takes `scheme`. Same-program trial
    # selection; per-environment BLUEs with within-environment standardization after dropping
    # suspect environments; the single best-covering genotyping project per trial; per-project
    # GRMs combined by covariance_combiner.
    #   CV0  : GBLUP with location and year as fixed effects.
    #   CV00 : the combined GRM turned into a Gaussian kernel, fitted as RKHS with BRR
    #          environmental effects.
    Prediction5 = .make_seed(c(list(
      train_select.method       = "same_program",
      train_select.prog_cap     = 20,
      pheno_prep.method         = "blue_lm",
      pheno_prep.z_thr          = 3,
      pheno_prep.standardize    = "per_trial_z",
      geno_select.method        = "best_single_project",
      predict_post.method       = "direct_blup"),
      if (scheme == "CV00")
        list(kernel.method = "rkhs_gaussian", kernel.rkhs_theta = 1,
             model.method = "rkhs", model.include_E = "yes")
      else
        list(kernel.method = "em_combine",
             model.method = "gblup_rrblup", model.lambda_select = "reml")
    ))
  )
}
