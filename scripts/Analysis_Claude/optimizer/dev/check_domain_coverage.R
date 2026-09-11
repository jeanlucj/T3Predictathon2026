# dev/check_domain_coverage.R
#
# Is every trial in target_domain actually USABLE? (Not just: does its name resolve.)
#
# TWO questions, one answer. `domain_not_covered` says which names failed to RESOLVE, not why --
# and resolving is not the same as being usable. A trial can resolve, be pinned into the universe,
# be counted in the domain size, and still never be evaluated: sample_real_trial() rejects on
# accession count and on an empty observation set with a bare `next`, silently. Trial 10154
# (23Ros_AY1-3Yield) did exactly that for 2,471 evaluations. So this reports one verdict per
# trial covering both questions.
#
# The trap it exists for: .apply_target_domain() filters in sequence -- programs, years,
# locations, THEN trials. So `trials` intersects with the others rather than overriding them,
# and a named trial outside your `years` is gone before the name filter runs. A domain that
# names trials AND constrains years is easy to write and self-inconsistent.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript dev/check_domain_coverage.R                     # the configured domain
#   Rscript dev/check_domain_coverage.R --trials=a,b,c       # a CANDIDATE domain
#
# Read-only: catalogue, per-trial accession and observation caches (30 days each), and the local
# project index. No store, no locks. Safe beside anything.
#
# HOW LONG. Everything is a per-trial cache read, so a domain the optimizer has already run is
# seconds; a fresh candidate domain of 15-30 trials is a minute or two, once, and then cached.
#
# ONE PREREQUISITE. If the projects index is not prewarmed, projects_for_accessions() falls back
# to the BrAPI wizard and the `bridged` column becomes the slow one. Fill it first:
#   Rscript tools/prepare_indices.R --only=projects     # ~110 calls, about a minute
# A run whose log says "project discovery: LOCAL" is already on the fast path.

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
# --trials= screens a CANDIDATE domain without editing settings.local.R first, which is the
# point of a pre-flight: you want the verdict before you commit the domain, not after.
args   <- commandArgs(trailingOnly = TRUE)
o_tr   <- grep("^--trials=", args, value = TRUE)
if (length(o_tr)) td$trials <- trimws(strsplit(sub("^--trials=", "", o_tr[1]), ",")[[1]])
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

# --- usability, for the trials that resolved --------------------------------
# Resolving is not the same question as being usable, and reading the first as the second is how
# a domain gets launched with trials the optimizer can never sample. Every gate below is one the
# pipeline already applies SILENTLY: sample_real_trial() rejects on accessions and on an empty
# observation set with a bare `next`, and a trial whose yields are all identical scores
# `constant` after a full evaluation per configuration.
#
# Costs are per-trial caches (cache/acc, cache/obs, both 30 days) plus the local project index,
# so a domain that has been run is seconds and a fresh one is a minute or two.
present <- want[verdict == "PRESENT"]
usable  <- setNames(rep(NA_character_, length(want)), want)
det     <- setNames(rep("", length(want)), want)

if (length(present)) {
  pid <- as.character(cand0$study_db_id[match(present, as.character(cand0$study_name))])

  cat("\n=== usability ===\n")
  n_acc <- purrr::map_int(pid, function(i)
    length(tryCatch(get_trial_accessions(i, conn, s), error = function(e) character())),
    .progress = "Trial accessions")

  # ONE call for every trial: get_observations() already carries its own progress bar and its
  # result carries study_id, so grouping gives the per-trial numbers without a loop.
  obs <- tryCatch(get_observations(pid, conn, s), error = function(e) NULL)
  ostat <- if (is.null(obs) || !nrow(obs))
    tibble::tibble(study_id = character(), n_obs = integer(), obs_sd = numeric())
  else obs |> dplyr::group_by(study_id) |>
    dplyr::summarise(n_obs = dplyr::n(),
                     obs_sd = stats::sd(value), .groups = "drop")
  n_obs  <- ostat$n_obs[match(pid, ostat$study_id)];  n_obs[is.na(n_obs)] <- 0L
  obs_sd <- ostat$obs_sd[match(pid, ostat$study_id)]

  # Genotyping structure. `bridged` asks whether any project that covers enough of THIS trial's
  # accessions also covers enough accessions belonging to the OTHER domain trials -- because
  # under CV00 masking makes the training set disjoint from the focal accessions, so a project
  # holding only this trial's germplasm can never supply training rows. Other domain trials are
  # a proxy for the training set, which train_select actually chooses per configuration: FALSE
  # is decisive, TRUE only means "not ruled out".
  acc_of <- purrr::map(pid, function(i)
    tryCatch(get_trial_accessions(i, conn, s), error = function(e) character()))
  names(acc_of) <- pid

  # Coverage counts come from the INVERTED INDEX (accession -> projects), never from
  # get_project_dosage(): a dosage matrix is hundreds of MB and loading one per project per
  # trial would take hours. The index is a local RDS and answers the same question exactly.
  idx  <- tryCatch(.project_index(s), error = function(e) NULL)
  univ <- tryCatch(.all_project_ids(conn, s), error = function(e) character())
  have_idx <- .index_covers(idx, univ, s, "proj_acc")
  if (!have_idx)
    cat("  ** the projects index does not cover the crop, so `bridged` cannot be computed\n",
        "     locally. Run: Rscript tools/prepare_indices.R --only=projects  (~1 min) **\n", sep = "")
  # project -> how many of `a` it holds
  covered <- function(a) {
    hit <- unlist(idx[intersect(a, names(idx))], use.names = FALSE)
    if (!length(hit)) integer(0) else table(as.character(hit))
  }
  mt <- .min_test(s); mn <- .min_train(s)
  geo <- purrr::map_dfr(seq_along(pid), function(j) {
    if (!have_idx) return(tibble::tibble(n_proj = NA_integer_, max_cov = NA_integer_,
                                         bridged = NA))
    cf <- covered(acc_of[[j]])
    # setdiff is LOAD-BEARING. Under CV00 mask_cv() removes every focal-trial accession from the
    # training set, so another trial's accessions count as training ONLY where they are not also
    # this trial's. Trials 10714 (132 accessions) and 10725 (196, a superset) are the proof: with
    # 10714 focal, 10725's extra 64 lines remain as training in the same panel and it succeeds;
    # with 10725 focal, 10714 contributes nothing at all and it fails 100%. Without the setdiff
    # this column would call 10725 bridged and be wrong.
    co <- covered(setdiff(unique(unlist(acc_of[-j], use.names = FALSE)), acc_of[[j]]))
    q  <- names(cf)
    tibble::tibble(n_proj  = length(cf),
                   max_cov = if (length(cf)) max(cf) else 0L,
                   # A project that holds enough of THIS trial and enough of the others: under
                   # CV00 the training set is disjoint from the focal accessions, so a project
                   # holding only this trial's germplasm can never supply training rows.
                   bridged = any(cf[q] >= mt &
                                 ifelse(q %in% names(co), co[q], 0L) >= mn))
  })

  u <- ifelse(n_acc < s$min_trial_acc,
         sprintf("%d accessions < min_trial_acc %d", n_acc, s$min_trial_acc),
       ifelse(n_obs == 0L, "no focal-trait OBSERVATIONS (the catalogue filters by id, get_observations by name)",
       ifelse(!is.finite(obs_sd) | obs_sd == 0, "all observations identical -- every evaluation would score `constant`",
       ifelse(is.na(geo$bridged), "USABLE (bridged UNCHECKED -- projects index is cold)",
       ifelse(!geo$bridged, sprintf("no genotyping project holds >=%d of its accessions AND >=%d from other domain trials", mt, mn),
              "USABLE")))))
  usable[present] <- u
  det[present] <- sprintf("acc %4d  obs %5d  sd %6s  proj %3s  max_cov %4s  bridged %s",
                          n_acc, n_obs,
                          ifelse(is.finite(obs_sd), sprintf("%.2f", obs_sd), "--"),
                          ifelse(is.na(geo$n_proj), "?", geo$n_proj),
                          ifelse(is.na(geo$max_cov), "?", geo$max_cov),
                          ifelse(is.na(geo$bridged), "?", ifelse(geo$bridged, "yes", "NO")))
}

cat("\n=== verdict ===\n")
# One verdict column, not two: "resolved" and "usable" are different questions, and showing them
# separately is what let a PRESENT-but-unusable trial through.
final <- ifelse(verdict != "PRESENT", unname(verdict), unname(usable))
tbl <- tibble::tibble(trial = want, verdict = final, detail = unname(det)) |>
  dplyr::arrange(verdict != "USABLE", verdict, trial)
print(tbl, n = nrow(tbl), width = Inf)
cat(sprintf("\n  %d of %d USABLE; %d resolved but unusable; %d did not resolve\n",
            sum(final == "USABLE"), length(final),
            sum(final != "USABLE" & verdict == "PRESENT"),
            sum(verdict != "PRESENT")))
if (any(final == "USABLE"))
  cat("\n  target_domain$trials for the usable set:\n    ",
      paste(sprintf('"%s"', want[final == "USABLE"]), collapse = ", "), "\n", sep = "")

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
cat("  USABLE                      resolves, is samplable, has varying focal-trait data, and at\n")
cat("                              least one genotyping project can supply both focal and\n")
cat("                              training rows. `bridged` uses the other domain trials as a\n")
cat("                              proxy for the training set, so this is a necessary condition,\n")
cat("                              not a sufficient one.\n")
cat("  < min_trial_acc             sample_real_trial() skips it SILENTLY, for ever -- the trial\n")
cat("                              is pinned into the universe and never evaluated.\n")
cat("  no focal-trait OBSERVATIONS the catalogue search and get_observations() use DIFFERENT\n")
cat("                              trait filters (id vs name); this trial satisfies one, not the\n")
cat("                              other. Also a silent skip.\n")
cat("  all observations identical  scores `constant` after a full evaluation per configuration.\n")
cat("  no genotyping project ...   under CV00 the training set is disjoint from the focal\n")
cat("                              accessions, so a project holding only this trial's germplasm\n")
cat("                              can never supply training rows. Remove it from the domain.\n")
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
