# tools/report_replication.R
#
# How much work is left before every configuration reaches config_replication, and how long will
# it take? Ask this BEFORE concluding that a run is stuck.
#
# WHY IT MATTERS. choose_config() serves the two replication tiers unequally:
#
#   pool <- if (length(bl$base)) bl$base                       # every config short of the floor
#           else if ((nrow(evals) + w) %% every == 0L) bl$extra  # the contender tier, rationed
#           else character()
#
# The contender tier is reached ONLY when the base backlog is empty. So a run whose contenders
# hold steady without gaining trials is not stalled -- it is paying down the base backlog, and
# that is the number to look at. A contender cannot accumulate evaluations until it clears.
#
# Read it as a TREND, not a level. New configurations are created exactly when the backlog
# empties and re-enter it immediately, so a momentary zero means little; two readings an hour
# apart mean a lot.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript tools/report_replication.R                    # this run's domain, 12 workers assumed
#   Rscript tools/report_replication.R --workers=22       # the wall-clock estimate uses this
#   Rscript tools/report_replication.R state/old.sqlite   # a specific store (or --store=)
#
# Scoped to the current target domain, like tools/report_memory.R: the store outlives a domain
# and replication owed in a retired one is not work this run will do.
#
# SAFETY: opens a COPY with its -wal/-shm sidecars, never the live store. Read-only, seconds,
# safe beside a running job.

# The optimizer ROOT, not this script's directory: `.Renviron` -- and so the T3
# credentials -- is read from the WORKING DIRECTORY only, with no parent walk. `settings.R`
# is the marker for that root. here::i_am() below also halts from the wrong place, but only
# when it cannot find the project at all, and a cwd inside the project is not one of those
# cases: here() still resolves while .Renviron silently does not.
if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript tools/report_replication.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("tools/report_replication.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
nw <- suppressWarnings(as.integer(sub("^--workers=", "", grep("^--workers=", args, value = TRUE))))
if (!length(nw) || !is.finite(nw) || nw < 1L) nw <- 12L
pos <- args[!grepl("^--", args)][1]

s   <- optimizer_settings()
con <- open_store(.copy_store_with_sidecars(resolve_read_store(pos, "tools/report_replication.R"),
                                            file.path(tempdir(), "replication.sqlite")))
all_e <- read_evals(con)
rid   <- run_id_for(s, s$build %||% OPTIMIZER_BUILD)
rr    <- tryCatch(read_run(con, rid), error = function(e) NULL)
close_store(con)
uni <- run_universe(rr)

e <- all_e |> filter_evals_to_scheme(s$optimize_scheme) |>
              filter_evals_to_build(s$build %||% OPTIMIZER_BUILD)
if (length(uni)) e <- filter_evals_to_universe(e, uni)
if (!nrow(e)) stop("no rows in this run's domain/scheme/build slice.")
if (!length(uni))
  cat("** no run row for ", rid, ": not restricted to the pinned universe, only to scheme\n",
      "   and build. Configurations from another domain may be counted below. **\n", sep = "")

cr <- as.integer(s$config_replication %||% 1L)
cat(sprintf("\nslice: %d eval(s), %d config(s); universe %d trial(s)\n",
            nrow(e), dplyr::n_distinct(e$config_hash), length(uni %||% character())))
cat(sprintf("config_replication %d   replicate_every %d   trial_replication %d\n",
            cr, s$replicate_every %||% 1L, s$trial_replication %||% 1L))

# Counted exactly as .replication_backlog() does, including its universe drop: a configuration
# that has already covered every trial is owed nothing and leaves the backlog.
per <- e |> dplyr::group_by(config_hash) |>
  dplyr::summarise(n_eval = dplyr::n(), n_trial = dplyr::n_distinct(trial_id), .groups = "drop")
covered <- if (length(uni)) sum(per$n_trial >= length(uni)) else 0L
if (length(uni)) per <- dplyr::filter(per, n_trial < length(uni))

cat("\nevaluations per configuration (backlog-eligible only):\n")
print(dplyr::count(per, n_eval, name = "configs"), n = Inf)
if (covered) cat(sprintf("  (%d configuration(s) have covered the whole domain and are owed nothing)\n", covered))

short <- dplyr::filter(per, n_eval < cr)
owed  <- sum(cr - short$n_eval)
secs  <- suppressWarnings(stats::median(e$seconds[is.finite(e$seconds)], na.rm = TRUE))
cat(sprintf("\nBASE backlog: %d configuration(s) short of %d, %d evaluation(s) owed\n",
            nrow(short), cr, owed))
if (is.finite(secs)) {
  cat(sprintf("median evaluation: %.0f min (mean %.0f)\n", secs / 60,
              mean(e$seconds[is.finite(e$seconds)]) / 60))
  cat(sprintf("-> %.0f h of compute; %.1f h wall at %d workers\n",
              owed * secs / 3600, owed * secs / 3600 / nw, nw))
} else cat("no timing recorded, so no estimate.\n")

cat("\nThe contender tier is served only while this is EMPTY, so a non-zero backlog is why\n")
cat("contenders hold steady without gaining trials. Take two readings an hour apart: what\n")
cat("matters is whether the backlog TRENDS down, not whether it is momentarily zero.\n")
