# config_space.R
#
# The configuration space of a genomic-prediction pipeline: the six subtasks, the methods each
# can use (drawn from the five Predictathon submissions and the literature), and each method's
# tunable parameters. The single source of truth for the search space -- sampling, feature
# encoding, crossover and mutation are all defined here in terms of SUBTASKS.
#
# To add a method: add its name to a subtask's `methods`, add any new parameters to that
# subtask's `params` (tagging which methods use them), and add a matching branch in the
# dispatcher in R/pipeline.R. Nothing else needs to change.
#
# A configuration is a FLAT named list: one "<subtask>.method" key per subtask plus one key per
# parameter, e.g.
#   list(train_select.method = "accession_overlap",
#        train_select.primary_min = 4, train_select.secondary_min = 12, ...)
# Parameters inapplicable to the chosen method are NA.

library(tidyverse)

# ---------------------------------------------------------------------------
# The search space. For each subtask: the candidate methods, and the parameters
# (each tagged with the methods it applies to; `methods = NULL` means "all").
# Parameter types: "int" (integer range), "real" (continuous range; scale "log"
# samples log-uniformly), "cat" (categorical values).
# ---------------------------------------------------------------------------
SUBTASKS <- list(

  # A. Select training trials from T3 for a focal trial.
  train_select = list(
    methods = c("accession_overlap", "top_k_similar", "same_program"),
    params = list(
      primary_min   = list(type = "int", range = c(2, 8),   methods = "accession_overlap"),
      secondary_min = list(type = "int", range = c(4, 20),  methods = "accession_overlap"),
      # primary_only="yes": use only the primary tier (direct germplasm-overlap neighbours
      # of the focal trial; secondary_min ignored). "no": also add the secondary,
      # friends-of-friends tier at secondary_min.
      primary_only  = list(type = "cat", values = c("yes", "no"), methods = "accession_overlap"),
      k             = list(type = "int", range = c(5, 40),  methods = "top_k_similar"),
      similarity    = list(type = "cat", values = c("genomic", "environmental", "both"),
                           methods = "top_k_similar"),
      prog_cap      = list(type = "int", range = c(10, 40), methods = "same_program")
    )
  ),

  # B. Preprocess training phenotypes into one target value per accession.
  pheno_prep = list(
    methods = c("blue_lm", "trial_center", "two_stage_blup", "raw_mean"),
    params = list(
      z_thr        = list(type = "real", range = c(2, 5),  methods = c("blue_lm", "trial_center")),
      ge_weighting = list(type = "cat",  values = c("none", "env_gaussian"), methods = NULL),
      ge_bandwidth = list(type = "real", range = c(0.3, 3), scale = "log", methods = NULL),
      standardize  = list(type = "cat",  values = c("none", "per_trial_z"), methods = NULL)
    )
  ),

  # C. Select which genotyping data (projects / VCFs) to use.
  #
  # Marker density is NOT a search parameter: it is set by settings$dosage_budget_bytes at
  # parse time, and a per-config request could only subset that cache (LESSONS.md #16).
  geno_select = list(
    methods = c("focal_plus_onehop", "best_single_project", "all_projects"),
    params = list(
      # Admit a protocol group only if it shares at least `min_bridge` accessions with the
      # group(s) genotyping the FOCAL trial. Those shared lines are the rows em_combine stitches
      # panels through: 1 admits a panel hanging off a single line, 5 demands a real anchor.
      min_bridge  = list(type = "int", range = c(1, 5),  methods = "focal_plus_onehop")
    )
  ),

  # D. Turn markers into a relationship matrix / kernel.
  kernel = list(
    methods = c("vanRaden_single", "em_combine", "rkhs_gaussian"),
    params = list(
      maf         = list(type = "real", range = c(0.01, 0.10), methods = NULL),
      max_missing = list(type = "real", range = c(0.10, 0.50), methods = NULL),
      ridge       = list(type = "real", range = c(1e-5, 1e-2), scale = "log", methods = NULL),
      impute      = list(type = "cat",  values = c("mean_round", "mean"), methods = NULL),
      rkhs_theta  = list(type = "real", range = c(0.2, 5), scale = "log", methods = "rkhs_gaussian")
    )
  ),

  # E. Train the prediction model.
  #
  # `lambda_select` is how the ridge (variance ratio lambda = sigma2_e / sigma2_u) is CHOSEN,
  # which is a separate question from which model is fit -- hence a knob, not a method name.
  # It applies to every method: all of them reach the same GBLUP backbone when the G+E path
  # is off or unavailable.
  model = list(
    methods = c("gblup_rrblup", "gblup_sommer_GE", "rkhs"),
    params = list(
      include_E     = list(type = "cat", values = c("yes", "no"),
                           methods = c("gblup_sommer_GE", "rkhs")),
      lambda_select = list(type = "cat", values = c("reml", "fixed", "loo"), methods = NULL),
      # Used only when lambda_select = "fixed"; sampled but ignored otherwise.
      lambda_fixed  = list(type = "real", range = c(1e-2, 1e2), scale = "log", methods = NULL)
    )
  ),

  # F. Predict on the focal trial and post-process.
  predict_post = list(
    methods = c("direct_blup", "cond_expectation"),
    params = list(
      blend_obs_w = list(type = "real", range = c(0, 0.95), methods = NULL),
      min_overlap = list(type = "int",  range = c(3, 15),   methods = NULL)
    )
  )
)

# ---------------------------------------------------------------------------
# Helpers over the spec.
# ---------------------------------------------------------------------------

# Does parameter `p` apply to method `m`? (NULL `methods` => applies to all.)
.param_applies <- function(p_spec, method) {
  is.null(p_spec$methods) || method %in% p_spec$methods
}

# Sample one value for a parameter spec.
.sample_param <- function(p_spec) {
  switch(p_spec$type,
    int  = sample(seq.int(p_spec$range[1], p_spec$range[2]), 1),
    cat  = sample(p_spec$values, 1),
    real = {
      if (!is.null(p_spec$scale) && p_spec$scale == "log") {
        exp(stats::runif(1, log(p_spec$range[1]), log(p_spec$range[2])))
      } else {
        stats::runif(1, p_spec$range[1], p_spec$range[2])
      }
    },
    fatal(paste("unknown param type", p_spec$type), "bug_unknown_param_type")
  )
}

# Key name for a subtask's chosen method, and for a parameter.
.method_key <- function(subtask) paste0(subtask, ".method")
.param_key  <- function(subtask, param) paste0(subtask, ".", param)

# Sample method + applicable parameters for ONE subtask -> a named list (the
# block for that subtask). Inapplicable parameters are set to NA.
sample_block <- function(subtask) {
  spec   <- SUBTASKS[[subtask]]
  method <- sample(spec$methods, 1)
  block  <- stats::setNames(list(method), .method_key(subtask))
  for (pname in names(spec$params)) {
    p_spec <- spec$params[[pname]]
    block[[.param_key(subtask, pname)]] <-
      if (.param_applies(p_spec, method)) .sample_param(p_spec) else NA
  }
  block
}

# Sample a complete random configuration across all six subtasks.
sample_config <- function() {
  cfg <- list()
  for (subtask in names(SUBTASKS)) cfg <- c(cfg, sample_block(subtask))
  cfg
}

# For each subtask, take the WHOLE block (method + its params) from parent a or parent b: the
# "recombine submitted algorithms" operator, e.g. P1's trial selection with P5's kernel.
crossover <- function(a, b) {
  child <- list()
  for (subtask in names(SUBTASKS)) {
    keys   <- c(.method_key(subtask),
                paste0(subtask, ".", names(SUBTASKS[[subtask]]$params)))
    parent <- if (stats::runif(1) < 0.5) a else b
    child  <- c(child, parent[keys])
  }
  child
}

# Mutation: resample `n_subtasks` randomly chosen subtask blocks of a config.
mutate_config <- function(cfg, n_subtasks = 1) {
  subtasks <- sample(names(SUBTASKS), n_subtasks)
  for (subtask in subtasks) {
    keys <- c(.method_key(subtask),
              paste0(subtask, ".", names(SUBTASKS[[subtask]]$params)))
    cfg[keys] <- NULL
    cfg <- c(cfg, sample_block(subtask))
  }
  # Preserve canonical key order so hashing/encoding are stable.
  cfg[canonical_keys()]
}

# Canonical ordering of all config keys (method then params, per subtask).
canonical_keys <- function() {
  unlist(lapply(names(SUBTASKS), function(subtask) {
    c(.method_key(subtask),
      paste0(subtask, ".", names(SUBTASKS[[subtask]]$params)))
  }), use.names = FALSE)
}

# Stable hash of a configuration, for de-duplication in the store. Hashing the canonical
# full-precision JSON rather than the R list makes it immune to int-vs-double typing and
# key-order differences, so a config and its JSON round-trip hash identically.
config_hash <- function(cfg) {
  rlang::hash(config_to_json(cfg))
}

# Coerce a config so inapplicable params are NA and key order is canonical. Run after crossover,
# where a param may be irrelevant to the method inherited from the other parent.
repair_config <- function(cfg) {
  for (subtask in names(SUBTASKS)) {
    method <- cfg[[.method_key(subtask)]]
    for (pname in names(SUBTASKS[[subtask]]$params)) {
      key <- .param_key(subtask, pname)
      if (!.param_applies(SUBTASKS[[subtask]]$params[[pname]], method)) {
        cfg[[key]] <- NA
      } else if (length(cfg[[key]]) == 0 || isTRUE(is.na(cfg[[key]]))) {
        # Applies but absent (a seed config omitted it) or NA (crossover brought in a method
        # needing a param the donor block lacked): fill with a fresh draw.
        cfg[[key]] <- .sample_param(SUBTASKS[[subtask]]$params[[pname]])
      }
    }
  }
  cfg[canonical_keys()]
}

# ---------------------------------------------------------------------------
# Feature encoding for the surrogate. rpart consumes a data.frame of factors (methods,
# categorical params) and numerics, with NA left in place for inapplicable parameters -- rpart
# handles those via surrogate splits. The fixed schema makes every config encode to the same
# columns in the same order.
# ---------------------------------------------------------------------------

# Each column's type and, for factors, its levels: numeric params -> "numeric"; methods and
# categorical params -> factor over the spec's allowed values.
feature_schema <- function() {
  schema <- list()
  for (subtask in names(SUBTASKS)) {
    spec <- SUBTASKS[[subtask]]
    schema[[.method_key(subtask)]] <- list(kind = "factor", levels = spec$methods)
    for (pname in names(spec$params)) {
      p_spec <- spec$params[[pname]]
      key <- .param_key(subtask, pname)
      if (p_spec$type == "cat") {
        schema[[key]] <- list(kind = "factor", levels = p_spec$values)
      } else {
        schema[[key]] <- list(kind = "numeric")
      }
    }
  }
  schema
}

# Encode a LIST of configurations into a single data.frame matching the schema.
configs_to_features <- function(cfg_list) {
  schema <- feature_schema()
  cols <- lapply(names(schema), function(key) {
    numeric_col <- schema[[key]]$kind == "numeric"
    vals <- vapply(cfg_list, function(cfg) {
      v <- cfg[[key]]
      if (is.null(v) || length(v) == 0) v <- NA
      if (numeric_col) as.numeric(v) else as.character(v)  # NA -> typed NA
    }, FUN.VALUE = if (numeric_col) numeric(1) else character(1))
    if (numeric_col) as.numeric(vals) else factor(vals, levels = schema[[key]]$levels)
  })
  names(cols) <- names(schema)
  as.data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)
}
