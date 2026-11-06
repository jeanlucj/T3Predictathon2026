# tools/plan_multiple.R
#
# Read multiple_settings.R, validate it, and emit one row per setting for a shell to loop over.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript tools/plan_multiple.R --check     # human summary + the wall-clock arithmetic
#   Rscript tools/plan_multiple.R             # TSV on stdout: file <TAB> days <TAB> holdouts
#
# Options:
#   --file=<path>   the plan; default multiple_settings.R in the optimizer root
#   --check         validate and describe, emit no TSV
#   --days=<n>      override days_per_setting (a short shakeout of a real plan)
#
# WHY A SEPARATE SCRIPT. container/t3multi_ceres.sh is bash and the plan is R, so something has
# to read one and hand it to the other. Putting it here rather than in the job script means the
# plan is checkable BEFORE a three-week allocation starts, from a login node, in a second -- and
# the checks below are the ones worth running early, because every one of them describes work
# that would otherwise fail on day four or, worse, succeed and mean nothing.
#
# On stdout goes ONLY the TSV; everything else is a message, i.e. stderr. A caller doing
# `plan=$(Rscript tools/plan_multiple.R)` must not get a banner in its first row.
#
# Read-only: no store, no network, no credentials.

if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript tools/plan_multiple.R --check\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("tools/plan_multiple.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", n, "="), "", h[1]) }
o_file  <- opt("file", here::here("multiple_settings.R"))
o_check <- "--check" %in% args
o_days  <- suppressWarnings(as.numeric(opt("days", NA)))

# How long a holdout test is allowed before the driver stops starting new settings. Matches the
# --time in container/submit_holdout.sh, which is the only measurement there is of what one
# costs. Read by the driver, printed here so the two cannot disagree silently.
HOLDOUT_RESERVE_HOURS <- 24

# ---- read ------------------------------------------------------------------
# A plan is R source, so an error in it is an R error. Evaluate it in its own environment: the
# plan must not be able to redefine a setting or a function this script just sourced.
if (!file.exists(o_file))
  stop("no plan at ", o_file, "\n  copy multiple_settings.R.example and edit it.", call. = FALSE)
e <- new.env(parent = globalenv())
tryCatch(sys.source(o_file, envir = e),
         error = function(err)
           stop(basename(o_file), " does not parse: ", conditionMessage(err), call. = FALSE))

days_global <- get0("days_per_setting", envir = e, ifnotfound = NULL)
if (is.finite(o_days)) days_global <- o_days
plan_list   <- get0("setting_file_names_holdouts", envir = e, ifnotfound = NULL)

if (is.null(days_global) || length(days_global) != 1L ||
    !is.finite(suppressWarnings(as.numeric(days_global))) || as.numeric(days_global) <= 0)
  stop(basename(o_file), " must define `days_per_setting` as one positive number", call. = FALSE)
if (!is.list(plan_list) || !length(plan_list))
  stop(basename(o_file), " must define `setting_file_names_holdouts` as a non-empty list()",
       call. = FALSE)

# ---- validate --------------------------------------------------------------
# Each check names the failure it prevents, because "invalid plan" on day four is not actionable.
problems <- character()
note     <- function(...) problems <<- c(problems, paste0(...))

rows <- purrr::imap(plan_list, function(el, i) {
  where <- sprintf("setting %d", i)
  if (!is.list(el)) { note(where, ": not a list()"); return(NULL) }
  sf <- el$setting_file_name
  ho <- el$holdouts
  if (!is.character(sf) || length(sf) != 1L || !nzchar(sf)) {
    note(where, ": `setting_file_name` must be one non-empty string"); return(NULL)
  }
  if (!is.character(ho) || !length(ho) || any(!nzchar(ho))) {
    note(where, " (", sf, "): `holdouts` must be a non-empty character vector. A setting with ",
         "no holdouts cannot be tested, which is the whole point of the sequence."); return(NULL)
  }
  d <- el$days %||% days_global
  if (length(d) != 1L || !is.finite(suppressWarnings(as.numeric(d))) || as.numeric(d) <= 0) {
    note(where, " (", sf, "): `days` must be one positive number"); return(NULL)
  }

  path <- if (startsWith(sf, "/")) sf else here::here(sf)
  if (!file.exists(path)) {
    note(where, ": no settings file at ", path, "\n      ",
         "settings.local.*.R is gitignored -- copy it to this checkout.")
    return(NULL)
  }
  # The same reader run_optimizer.R and holdout_test.R use, so a file that passes here cannot
  # fail there: it enforces that a domain file sets only keys that can identify a run.
  ov <- tryCatch(read_domain_settings(path), error = function(err) {
    note(where, " (", sf, "): ", conditionMessage(err)); NULL })
  if (is.null(ov)) return(NULL)

  # THE INVARIANT. A holdout inside its own target domain means the optimizer trains on the
  # trial the test then "holds out", and the comparison with the submissions is worthless.
  # Checkable offline only when the domain names its trials; programs/years/locations resolve
  # against the live catalogue, so those are reported as unchecked rather than assumed clean.
  td   <- ov$target_domain
  bad  <- if (!is.null(td$trials)) intersect(ho, as.character(td$trials)) else character()
  if (length(bad))
    note(where, " (", sf, "): holdout(s) ", paste(bad, collapse = ", "),
         " are IN target_domain$trials. Remove them from the domain or drop them as holdouts; ",
         "as written the test would hold out a trial the optimizer trained on.")
  unchecked <- is.null(td$trials) &&
    any(!vapply(td[c("programs", "years", "locations")], is.null, logical(1)))

  tibble::tibble(settings_file = sf, days = as.numeric(d),
                 holdouts = paste(ho, collapse = ","), n_holdouts = length(ho),
                 scheme = as.character(ov$optimize_scheme %||% NA_character_),
                 n_domain_trials = length(td$trials %||% character()),
                 unchecked = unchecked)
}) |> purrr::compact() |> dplyr::bind_rows()

# One file per setting: the holdout store is named from it (holdout_XYZ.sqlite), so two settings
# sharing a file would write one store and the second would be read as the first.
if (nrow(rows)) {
  dup <- rows$settings_file[duplicated(rows$settings_file)]
  if (length(dup))
    note("duplicate settings file(s): ", paste(unique(dup), collapse = ", "),
         ". Each setting needs its own file -- the holdout store is named after it.")
}

if (length(problems))
  stop("plan ", basename(o_file), " has ", length(problems), " problem(s):\n  - ",
       paste(problems, collapse = "\n  - "), call. = FALSE)
if (!nrow(rows)) stop("plan ", basename(o_file), " yielded no settings", call. = FALSE)

# ---- emit ------------------------------------------------------------------
if (o_check) {
  message("plan      : ", o_file)
  message("settings  : ", nrow(rows), "   days_per_setting: ", days_global,
          if (is.finite(o_days)) "  (overridden by --days)" else "")
  print(as.data.frame(rows[c("settings_file", "scheme", "days", "n_domain_trials",
                             "n_holdouts", "holdouts")]), row.names = FALSE)
  if (any(rows$unchecked))
    message("\nNOTE: ", sum(rows$unchecked), " setting(s) define the domain by ",
            "program/year/location rather than by trial name, so the holdout-not-in-domain ",
            "check could not run offline for them. Confirm with ",
            "dev/check_domain_coverage.R against today's catalogue.")
  # Holdout names are resolved to ids by holdout_test.R against the LIVE catalogue, so a trial
  # T3 has since renamed fails there, hours in, not here. The offline check cannot see that.
  message("\nNOT checked here (needs T3): that every holdout name still resolves. ",
          "Run dev/check_trial_renames.R if any of these trials are older than a few months.")
  opt_d  <- sum(rows$days)
  hold_d <- nrow(rows) * HOLDOUT_RESERVE_HOURS / 24
  message(sprintf("\nwall clock: %.1f days optimizing + up to %.1f days of holdout tests = %.1f days",
                  opt_d, hold_d, opt_d + hold_d))
  message(sprintf("            the job's --time must be at least that; 21-00:00:00 is the ceres maximum%s",
                  if (opt_d + hold_d > 21) "  ** THIS PLAN EXCEEDS IT **" else ""))
  message("\nrun it with: container/submit_multiple.local.sh")
} else {
  # TSV only. No header: the driver reads it with `while IFS=$'\t' read -r f d h`.
  for (i in seq_len(nrow(rows)))
    cat(sprintf("%s\t%g\t%s\n", rows$settings_file[i], rows$days[i], rows$holdouts[i]))
}
