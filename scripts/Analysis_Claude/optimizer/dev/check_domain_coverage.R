# dev/check_domain_coverage.R
#
# Why is a named trial missing from the resolved universe?
#
# pin_trial_universe() aborts a run with `domain_not_covered` when a name in
# settings$target_domain$trials is absent from .eligible_trials(). That message says WHICH
# names, not WHY -- and there are four different reasons, with four different fixes. This walks
# each missing name down the same ladder .eligible_trials() applies and names the step that
# dropped it.
#
# The trap it exists for: .apply_target_domain() filters in sequence -- programs, years,
# locations, THEN trials. So `trials` intersects with the others rather than overriding them,
# and a named trial outside your `years` is gone before the name filter runs. A domain that
# names trials AND constrains years is easy to write and self-inconsistent.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript dev/check_domain_coverage.R
#
# Read-only: one catalogue fetch (usually cached), no store, no genotypes, no locks. Safe
# beside anything.

# The optimizer ROOT, not this script's directory: `.Renviron` -- and so the T3
# credentials -- is read from the WORKING DIRECTORY only, with no parent walk. `settings.R`
# is the marker for that root. here::i_am() below also halts from the wrong place, but only
# when it cannot find the project at all, and a cwd inside the project is not one of those
# cases: here() still resolves while .Renviron silently does not.
if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript dev/check_domain_coverage.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("dev/check_domain_coverage.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

s    <- optimizer_settings()
td   <- s$target_domain
want <- as.character(td$trials %||% character())
if (!length(want))
  stop("settings$target_domain$trials is empty -- nothing to check.\n",
       "  This script answers `domain_not_covered`, which only arises when trials are named.")

conn <- t3_connect(s)
cat("=== catalogue ===\n")
cand0 <- trial_catalog(conn, s) |> dplyr::filter(!is.na(study_db_id))
cat(sprintf("  %d trials, refreshed %s\n", nrow(cand0), .trial_catalog_ts(s) %||% "unknown"))
cat(sprintf("  focal_trait_db_id: %s%s\n", s$focal_trait_db_id %||% "(unset)",
            if (!is.null(s$focal_trait_db_id))
              "  <- the catalogue ALREADY holds only trials with this trait" else ""))

# One filter at a time, in .apply_target_domain's order, so `dropped_by` is the FIRST step that
# removed a name -- which is the one to fix.
steps <- list(
  programs  = function(d) if (is.null(td$programs))  d else dplyr::filter(d, program_name %in% td$programs),
  years     = function(d) if (is.null(td$years))     d else dplyr::filter(d, suppressWarnings(as.integer(year)) %in% as.integer(td$years)),
  locations = function(d) if (is.null(td$locations)) d else dplyr::filter(d, location_name %in% td$locations))

cat("\n=== target_domain ===\n")
for (k in c("programs", "years", "locations"))
  cat(sprintf("  %-10s %s\n", k,
              if (is.null(td[[k]])) "(unset -- no constraint)" else paste(td[[k]], collapse = ", ")))
cat(sprintf("  %-10s %d name(s)\n", "trials", length(want)))

# Absent from the trait-filtered catalogue is TWO different problems with two different fixes,
# so separate them against the unfiltered wheat trial list.
#
# Do NOT get that by calling trial_catalog() with focal_trait_db_id nulled: it caches under the
# singleton key "trial_catalog" with no identifier, so that call either returns the cached
# FILTERED table (useless) or overwrites the run's catalogue with an unfiltered one (worse).
# get_all_trial_meta_data() is the same fetch trial_catalog() makes in its no-trait branch, and
# calling it here goes nowhere near the cache.
absent0 <- setdiff(want, as.character(cand0$study_name))
untrait <- NULL
if (length(absent0)) {
  cat("\n  ", length(absent0), " name(s) absent from the trait-filtered catalogue;",
      " fetching the unfiltered trial list to say why ...\n", sep = "")
  untrait <- tryCatch(
    tibble::as_tibble(T3BrapiHelpers::get_all_trial_meta_data(conn, s$crop_name)) |>
      janitor::clean_names(),
    error = function(e) { cat("  (unfiltered fetch failed: ", conditionMessage(e),
                              " -- absent names cannot be split)\n", sep = ""); NULL })
}

alive   <- cand0
verdict <- setNames(rep(NA_character_, length(want)), want)
if (is.null(untrait)) {
  verdict[absent0] <- "absent (could not split: no unfiltered list)"
} else {
  all_wheat <- as.character(untrait$study_name)
  verdict[intersect(absent0, all_wheat)]   <- "no focal-trait data"
  verdict[setdiff(absent0, all_wheat)]     <- "name not found in T3"
}
for (k in names(steps)) {
  if (is.null(td[[k]])) next
  before <- as.character(alive$study_name)
  alive  <- steps[[k]](alive)
  gone   <- setdiff(intersect(want, before), as.character(alive$study_name))
  verdict[gone[is.na(verdict[gone])]] <- paste0("dropped by target_domain$", k)
}
verdict[is.na(verdict)] <- "PRESENT"

cat("\n=== verdict ===\n")
tbl <- tibble::tibble(trial = want, why = unname(verdict)) |> dplyr::arrange(why != "PRESENT", why, trial)
print(tbl, n = nrow(tbl))
cat(sprintf("\n  %d of %d present; %d missing\n",
            sum(tbl$why == "PRESENT"), nrow(tbl), sum(tbl$why != "PRESENT")))

# A name absent from the catalogue is usually a typo or a variant, not a missing trial, so give
# the nearest candidates rather than just "absent".
absent <- tbl$trial[tbl$why == "name not found in T3"]
if (length(absent)) {
  cat("\n=== nearest catalogue names for the absent ones ===\n")
  norm <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(x)))
  allnm <- as.character(cand0$study_name); alln <- norm(allnm)
  for (a in absent) {
    exact <- allnm[alln == norm(a)]
    near  <- if (length(exact)) exact else
      allnm[order(utils::adist(norm(a), alln))][seq_len(min(3, length(allnm)))]
    cat(sprintf("  %-40s -> %s%s\n", a, paste(near, collapse = " | "),
                if (length(exact)) "   <- SAME after case/punctuation, i.e. a formatting difference" else ""))
  }
}

cat("\nFixes, by verdict:\n")
cat("  no focal-trait data         the trial exists in T3 but records no observation of\n")
cat("                              focal_trait_db_id (", s$focal_trait_db_id %||% "unset",
    "). It can never be scored, so remove it\n", sep = "")
cat("                              from target_domain$trials -- this is a property of the data,\n")
cat("                              not a misconfiguration.\n")
cat("  name not found in T3        a typo, a formatting difference, or a trial not visible to\n")
cat("                              this login -- see the near-matches above.\n")
cat("  dropped by target_domain$X  the named trial falls outside constraint X. Either widen X\n")
cat("                              or drop the trial from target_domain$trials: naming a trial\n")
cat("                              does NOT exempt it from the other filters.\n")
