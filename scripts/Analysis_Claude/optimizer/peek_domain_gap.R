# peek_domain_gap.R
#
# Why are stored evaluations OUT of the current slice, and when did the last in-slice one
# land? peek_backup.R reports the slice as a single number ("8004 row(s) are out of slice");
# this splits that number by cause and by date, which is what distinguishes "the run stopped"
# from "the run kept going but into a domain this run cannot see".
#
# Read-only, no network, seconds. Works on a copy, never the live store.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   ./container/run_in_container.sh exec peek_domain_gap.R

if (!file.exists("peek_domain_gap.R"))
  stop("run this FROM the optimizer directory:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("peek_domain_gap.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

if (.on_login_node()) stop(.login_node_message("peek_domain_gap.R"), call. = FALSE)
s  <- optimizer_settings()
hr <- function(t) cat("\n", t, "\n", strrep("-", nchar(t)), "\n", sep = "")

# The backup is the durable record; the node-local store is whatever this node happens to
# hold. Prefer the backup, fall back to the store, and say which one answered.
src <- if (!is.null(s$db_backup_path) && file.exists(s$db_backup_path)) s$db_backup_path else s$db_path
if (!file.exists(src)) stop("no store and no backup to read: ", src, call. = FALSE)
tmp <- tempfile(fileext = ".sqlite"); file.copy(src, tmp); on.exit(unlink(tmp), add = TRUE)
con <- open_store(tmp); on.exit(close_store(con), add = TRUE)

ev <- read_evals(con) |> mutate(day = substr(ts, 1, 10))
cat(sprintf("source : %s\n", src))
cat(sprintf("rows   : %d, %s to %s\n", nrow(ev), min(ev$ts, na.rm = TRUE), max(ev$ts, na.rm = TRUE)))
cat(sprintf("build  : %s   scheme: %s\n", s$build %||% OPTIMIZER_BUILD, s$optimize_scheme))

# The real predicates, not a re-implementation, so this cannot drift from what the run does.
in_dom <- ev |> filter_evals_to_domain(s$target_domain) |> pull(id)
ev <- ev |> mutate(
  dom    = if_else(id %in% in_dom, "in-domain", "out-of-domain"),
  sch    = if_else(scheme == s$optimize_scheme, "scheme-ok", "scheme-no"),
  bld    = if_else(build == (s$build %||% OPTIMIZER_BUILD), "build-ok", "build-no"),
  slice  = if_else(dom == "in-domain" & sch == "scheme-ok" & bld == "build-ok", "IN SLICE", "out"))

# --- what excludes each row ------------------------------------------------
hr("why rows are out of slice")
ev |> count(dom, sch, bld, slice) |> arrange(desc(n)) |> print(n = 50)

# --- when did each group last produce anything -----------------------------
hr("first and last row per group")
ev |> group_by(slice, dom) |> summarise(n = n(), first = min(ts), last = max(ts), .groups = "drop") |>
  arrange(desc(n)) |> print(n = 50)

# --- the domain columns themselves -----------------------------------------
# An out-of-domain row is either a trial from a DIFFERENT domain or a row whose domain
# columns were never recorded. Those two have completely different fixes, and only the
# study_name column tells them apart.
hr("study_name on out-of-domain rows")
oo <- ev |> filter(dom == "out-of-domain")
cat(sprintf("  rows with study_name NA or empty : %d of %d\n",
            sum(is.na(oo$study_name) | !nzchar(oo$study_name)), nrow(oo)))
oo |> filter(!is.na(study_name), nzchar(study_name)) |>
  count(study_name, trial_id, sort = TRUE) |> print(n = 60)

hr("target_domain$trials, and whether the store holds rows for each")
tibble(study_name = s$target_domain$trials %||% character()) |>
  left_join(ev |> count(study_name, name = "rows_in_store"), by = "study_name") |>
  mutate(rows_in_store = replace_na(rows_in_store, 0L)) |> print(n = 60)

# --- the calendar ----------------------------------------------------------
# The question two weeks away asks: did work keep landing, and did it land where this run
# can see it?
hr("rows per day, last 30 days present")
ev |> count(day, slice) |> pivot_wider(names_from = slice, values_from = n, values_fill = 0L) |>
  arrange(desc(day)) |> head(30) |> print(n = 30)
