# tools/inspect_failures.R
#
# Read-only look at every non-`ok` eval in the store, with the stored funnel unpacked into
# columns and a verdict per row. Answers -- in seconds, with no network and no re-running --
# the question tools/diagnose_failures.R would spend hours on:
#
#   is `insufficient_geno_overlap` caused by .best_panel passing over a focal-covering
#   panel, or by a name/synonym mismatch hiding genotypes?
#
# The funnel is already recorded: evaluate.R stores `detail = funnel_string(e$funnel)`, i.e.
#   focal_acc; train_acc; geno_projects; geno_acc; train_in; test_in
#
# Verdicts:
#   (a) PANEL RULE      test_in = 0 while train_in > 0 -- the kernel covers training lines
#                       but no focal lines. .best_panel maximizes coverage of
#                       union(train, focal), so a training-heavy panel can win.
#   (b) NAME MISMATCH   train_in = 0 AND test_in = 0 while geno_acc is large -- the genotypes
#                       are there but under names nothing matches. Fixing .best_panel would
#                       not help; look at the synonym path instead.
#   neither cliff       both sides non-zero, just below the thresholds: a genuinely thin trial.
#
# SAFETY: this never opens the live store. The database keeps data in a `-wal` sidecar, so
# the main file can be 0 bytes on its own -- all three files are copied to a temp directory
# and the COPY is queried. (Copying evals.sqlite without its -wal/-shm loses everything.)
#
# Run from THIS directory (R reads .Renviron and here() from the working directory):
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript tools/inspect_failures.R

# The optimizer ROOT, not this script's directory: `.Renviron` -- and so the T3
# credentials -- is read from the WORKING DIRECTORY only, with no parent walk. `settings.R`
# is the marker for that root. here::i_am() below also halts from the wrong place, but only
# when it cannot find the project at all, and a cwd inside the project is not one of those
# cases: here() still resolves while .Renviron silently does not.
if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript tools/inspect_failures.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("tools/inspect_failures.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

store_path <- resolve_read_store(what = "tools/inspect_failures.R")

tmp <- .copy_store_with_sidecars(store_path, file.path(tempdir(), "peek.sqlite"))
con <- DBI::dbConnect(RSQLite::SQLite(), tmp)
e   <- tibble::as_tibble(DBI::dbReadTable(con, "evals"))
DBI::dbDisconnect(con)

cat("store: ", store_path, "\n", sep = "")
cat("bytes: ", file.info(store_path)$size,
    if (file.exists(paste0(store_path, "-wal")))
      sprintf(" (+ %d in -wal)", file.info(paste0(store_path, "-wal"))$size) else "",
    "\nrows:  ", nrow(e), "\n\n", sep = "")
print(dplyr::count(e, status, sort = TRUE), n = Inf)

bad <- dplyr::filter(e, status != "ok")
if (!nrow(bad)) {
  cat("\nno non-ok rows -- nothing to explain.\n")
} else {
  # Pull one key out of the "k=v; k=v" funnel string.
  fget <- function(d, k) unname(vapply(d, function(x) {
    if (is.na(x)) return(NA_real_)
    p <- strsplit(trimws(strsplit(x, ";", fixed = TRUE)[[1]]), "=", fixed = TRUE)
    v <- stats::setNames(sapply(p, `[`, 2), sapply(p, `[`, 1))
    suppressWarnings(as.numeric(v[k]))
  }, numeric(1)))
  kern <- unname(vapply(bad$config_json, function(j)
    tryCatch(as.character(jsonlite::fromJSON(j)[["kernel.method"]]),
             error = function(err) NA_character_), character(1)))
  ti <- fget(bad$detail, "train_in"); te <- fget(bad$detail, "test_in")

  cat("\n")
  print(data.frame(
    id = bad$id, trial = bad$trial_id, scheme = bad$scheme, status = bad$status,
    kernel = kern, n_test = bad$n_test,
    focal = fget(bad$detail, "focal_acc"), geno_acc = fget(bad$detail, "geno_acc"),
    train_in = ti, test_in = te,
    verdict = ifelse(is.na(te), "- (no funnel recorded)",
      ifelse(te == 0 & ti >  0, "(a) PANEL RULE -- training lines but no focal lines",
      ifelse(te == 0 & ti == 0, "(b) NAME MISMATCH -- covers neither side",
      ifelse(ti == 0 & te >  0, "(b?) training side empty",
             "neither cliff -- genuinely thin"))))),
    right = FALSE, row.names = FALSE)

  # Rows whose status comes from SCORING (non_finite, too_few_overlap, constant) or from an
  # uncaught error carry a message rather than a funnel. Printing only the funnel columns
  # renders them "- (no funnel recorded)", which reads like missing information when the
  # reason string is in fact the whole explanation. Show it.
  noft <- dplyr::filter(bad, is.na(detail) | !grepl("test_in", detail))
  if (nrow(noft)) {
    cat("\nrows carrying a reason rather than a funnel (the reason IS the explanation):\n")
    print(as.data.frame(noft[, c("id", "trial_id", "status", "reason")]),
          right = FALSE, row.names = FALSE)
  }
}
