# dev/check_trial_renames.R
#
# Did a trial's study_name CHANGE while the store was being written? Domain membership is
# decided by study_name, but the replication backlog groups by trial_id, so an id carrying two
# names is present enough to block the backlog and invisible enough never to satisfy it.
#
# Read-only, no network, seconds. Works on a copy, never the live store.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   ./container/run_in_container.sh exec dev/check_trial_renames.R

# The optimizer ROOT, not this script's directory: `.Renviron` -- and so the T3
# credentials -- is read from the WORKING DIRECTORY only, with no parent walk. `settings.R`
# is the marker for that root. here::i_am() below also halts from the wrong place, but only
# when it cannot find the project at all, and a cwd inside the project is not one of those
# cases: here() still resolves while .Renviron silently does not.
if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript dev/check_trial_renames.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("dev/check_trial_renames.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

if (.on_login_node()) stop(.login_node_message("dev/check_trial_renames.R"), call. = FALSE)
s  <- optimizer_settings()
hr <- function(t) cat("\n", t, "\n", strrep("-", nchar(t)), "\n", sep = "")

src <- if (!is.null(s$db_backup_path) && file.exists(s$db_backup_path)) s$db_backup_path else s$db_path
if (!file.exists(src)) stop("no store and no backup to read: ", src, call. = FALSE)
tmp <- tempfile(fileext = ".sqlite"); file.copy(src, tmp); on.exit(unlink(tmp), add = TRUE)
con <- open_store(tmp); on.exit(close_store(con), add = TRUE)

ev <- read_evals(con)
cat(sprintf("source : %s\nrows   : %d, %s to %s\n", src, nrow(ev),
            min(ev$ts, na.rm = TRUE), max(ev$ts, na.rm = TRUE)))

# --- THE question: every name trial 10262 was ever stored under ------------
hr("trial 10262: rows per study_name")
ev |> filter(as.character(trial_id) == "10262") |>
  group_by(study_name) |>
  summarise(n = n(), n_cfg = n_distinct(config_hash),
            first = min(ts), last = max(ts), .groups = "drop") |>
  arrange(first) |> print(n = 20)

# --- and the same question asked of the whole store ------------------------
# If renaming is a general hazard rather than one trial's bad luck, it shows up here.
hr("every trial_id stored under MORE THAN ONE study_name")
multi <- ev |> group_by(trial_id) |>
  summarise(n_names = n_distinct(study_name),
            names = paste(sort(unique(study_name)), collapse = " | "),
            n = n(), .groups = "drop") |>
  filter(n_names > 1) |> arrange(desc(n))
if (!nrow(multi)) cat("  none -- every trial_id has exactly one name in the store\n")
print(multi, n = 50)

# --- the 8 orphan rows under the old name ----------------------------------
hr("rows stored under 25_Big6_MASON_MASMI")
ev |> filter(study_name == "25_Big6_MASON_MASMI") |>
  group_by(trial_id) |> summarise(n = n(), first = min(ts), last = max(ts), .groups = "drop") |>
  print(n = 20)

# --- what the catalogue CACHE says now, and when it last refreshed ----------
# Read straight off disk: trial_catalog() would go to the network if the cache had expired,
# and the refresh TIME is half the evidence.
hr("catalogue cache: current name for 10262, and cache age")
cat_file <- file.path(s$cache_dir, "trial_catalog", "trial_catalog.rds")
if (!file.exists(cat_file)) {
  cat("  no cached catalogue at ", cat_file, "\n", sep = "")
} else {
  cat(sprintf("  file    : %s\n  modified: %s\n", cat_file,
              format(file.info(cat_file)$mtime, "%Y-%m-%d %H:%M:%S")))
  tc <- readRDS(cat_file)
  if (is.list(tc) && !is.data.frame(tc)) tc <- tc[[which(map_lgl(tc, is.data.frame))[1]]]
  tc |> filter(as.character(study_db_id) == "10262") |>
    select(any_of(c("study_db_id", "study_name", "program_name", "location_name", "year"))) |>
    print()
  hit <- tc |> filter(study_name %in% s$target_domain$trials) |> nrow()
  cat(sprintf("  catalogue rows matching target_domain$trials: %d of %d (want 32)\n",
              hit, nrow(tc)))
}
