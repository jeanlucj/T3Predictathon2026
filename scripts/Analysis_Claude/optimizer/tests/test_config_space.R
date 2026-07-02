# test_config_space.R
#
# Invariants of the genome: sampling produces well-formed configs, encoding is
# schema-stable, JSON round-trips preserve NA, and crossover/mutation keep
# configs valid. Run: Rscript tests/test_config_space.R

library(tidyverse)
here::i_am("tests/test_config_space.R")
source(here::here("settings.R"))
invisible(lapply(list.files(here::here("R"), "[.]R$", full.names = TRUE), source))

ok <- 0L; fail <- 0L
check <- function(cond, msg) {
  if (isTRUE(cond)) { ok <<- ok + 1L }
  else { fail <<- fail + 1L; cat("  FAIL:", msg, "\n") }
}

set.seed(1)

# 1. A sampled config has exactly the canonical keys, methods are valid, and
#    inapplicable params are NA while applicable ones are not.
for (i in 1:200) {
  cfg <- sample_config()
  check(identical(names(cfg), canonical_keys()), "sample_config key order")
  for (st in names(SUBTASKS)) {
    m <- cfg[[paste0(st, ".method")]]
    check(m %in% SUBTASKS[[st]]$methods, paste("valid method", st))
    for (p in names(SUBTASKS[[st]]$params)) {
      v <- cfg[[paste0(st, ".", p)]]
      applies <- is.null(SUBTASKS[[st]]$params[[p]]$methods) ||
                 m %in% SUBTASKS[[st]]$params[[p]]$methods
      if (applies) check(!is.na(v), paste("applicable param set", st, p))
      else check(is.na(v), paste("inapplicable param NA", st, p))
    }
  }
}

# 2. Encoding is schema-stable: same columns regardless of which configs.
f1 <- configs_to_features(replicate(3, sample_config(), simplify = FALSE))
f2 <- configs_to_features(list(sample_config()))
check(identical(names(f1), names(f2)), "feature columns stable")
check(identical(names(f1), names(feature_schema())), "feature columns = schema")

# 3. JSON round-trip preserves the config (including NA).
for (i in 1:100) {
  cfg <- sample_config()
  rt  <- config_from_json(config_to_json(cfg))
  check(config_hash(cfg) == config_hash(rt), "json round-trip hash")
}

# 4. Crossover/mutation yield valid (repaired) configs.
for (i in 1:100) {
  a <- sample_config(); b <- sample_config()
  child <- repair_config(crossover(a, b))
  check(identical(names(child), canonical_keys()), "crossover canonical keys")
  for (st in names(SUBTASKS)) {
    m <- child[[paste0(st, ".method")]]
    for (p in names(SUBTASKS[[st]]$params)) {
      v <- child[[paste0(st, ".", p)]]
      applies <- is.null(SUBTASKS[[st]]$params[[p]]$methods) ||
                 m %in% SUBTASKS[[st]]$params[[p]]$methods
      check(applies == !is.na(v), paste("crossover param consistency", st, p))
    }
  }
  mut <- mutate_config(a, 2)
  check(identical(names(mut), canonical_keys()), "mutate canonical keys")
}

# 5. The five seeds are well-formed configs.
for (nm in names(seed_configs())) {
  s <- seed_configs()[[nm]]
  check(identical(names(s), canonical_keys()), paste("seed canonical", nm))
}

cat(sprintf("\nconfig_space tests: %d passed, %d failed\n", ok, fail))
if (fail > 0) quit(status = 1)
