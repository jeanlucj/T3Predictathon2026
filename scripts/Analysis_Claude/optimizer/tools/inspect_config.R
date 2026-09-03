# tools/inspect_config.R
#
# Print the stored configuration of specific evals. Read-only, no network, seconds.
#
# Written to settle whether trial 11033's `non_finite` evals used the G x E weighting path
# that build 0.7.3 fixed: `.trial_similarity` returned NA for a candidate trial with an
# unknown year (28% of the T3 catalogue), which became an NA weight, an NA target, and a
# non-finite prediction. That bug fires ONLY when pheno_prep.ge_weighting == "env_gaussian",
# so that one field is the confirmation.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript tools/inspect_config.R --ids=4,23
#
# Options:
#   --ids=4,23      eval ids to show. Default: every eval whose status is not `ok`.
#   --keys=a,b      config fields to show. Default: the six <subtask>.method keys plus
#                   pheno_prep.ge_weighting. `--keys=all` prints the whole configuration.
#
# SAFETY: copies the database and its -wal/-shm sidecars to a temp directory and reads the
# copy, so the live store is never opened (see tools/inspect_failures.R).

# The optimizer ROOT, not this script's directory: `.Renviron` -- and so the T3
# credentials -- is read from the WORKING DIRECTORY only, with no parent walk. `settings.R`
# is the marker for that root. here::i_am() below also halts from the wrong place, but only
# when it cannot find the project at all, and a cwd inside the project is not one of those
# cases: here() still resolves while .Renviron silently does not.
if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript tools/inspect_config.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("tools/inspect_config.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
split_csv <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",")[[1]])

o_ids  <- split_csv(opt("ids"))
o_keys <- split_csv(opt("keys", paste0(c(paste0(names(SUBTASKS), ".method"),
                                         "pheno_prep.ge_weighting"), collapse = ",")))

store_path <- resolve_read_store(what = "tools/inspect_config.R")
tmp <- .copy_store_with_sidecars(store_path, file.path(tempdir(), "peekcfg.sqlite"))
con <- DBI::dbConnect(RSQLite::SQLite(), tmp)
e   <- tibble::as_tibble(DBI::dbReadTable(con, "evals"))
DBI::dbDisconnect(con)

sel <- if (is.null(o_ids)) {
  dplyr::filter(e, status != "ok")
} else {
  dplyr::filter(e, as.character(id) %in% o_ids)
}
cat("store: ", store_path, "\nshowing ", nrow(sel), " of ", nrow(e), " evals\n\n", sep = "")
if (!nrow(sel)) { cat("no matching evals\n"); quit(save = "no") }

for (i in seq_len(nrow(sel))) {
  cfg <- tryCatch(config_from_json(sel$config_json[i]), error = function(err) NULL)
  cat(sprintf("--- eval %s | trial %s | %s | %s | %s ---\n",
              sel$id[i], sel$trial_id[i], sel$scheme[i], sel$status[i],
              sel$reason[i] %||% "-"))
  if (is.null(cfg)) { cat("  (config_json unparseable)\n\n"); next }
  keys <- if (identical(o_keys, "all")) names(cfg) else intersect(o_keys, names(cfg))
  for (k in keys) cat(sprintf("  %-28s %s\n", k, as.character(cfg[[k]])))
  # The 0.7.3 verdict, stated rather than left to the reader.
  if ("pheno_prep.ge_weighting" %in% names(cfg)) {
    hit <- identical(cfg$pheno_prep.ge_weighting, "env_gaussian")
    verdict <- c("does NOT use env_gaussian: the 0.7.3 bug is NOT the cause here",
                 "USES env_gaussian: the 0.7.3 NA-similarity bug could produce this")
    cat("  -> ", verdict[hit + 1L], "\n", sep = "")
  }
  cat("\n")
}
