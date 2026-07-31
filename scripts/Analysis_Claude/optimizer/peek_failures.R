# peek_failures.R
#
# Read-only look at every non-`ok` eval in the store, with the stored funnel unpacked into
# columns and a verdict per row. Answers -- in seconds, with no network and no re-running --
# the question diagnose_failures.R would spend hours on:
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
#   Rscript peek_failures.R

if (!file.exists("peek_failures.R"))
  stop("run this FROM the optimizer directory:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript peek_failures.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("peek_failures.R")
source(here::here("settings.R"))

s <- optimizer_settings()
if (!file.exists(s$db_path)) stop("no store at ", s$db_path)

# Copy the database AND its write-ahead-log sidecars, then read the copy.
tmp <- file.path(tempdir(), "peek.sqlite")
invisible(file.copy(s$db_path, tmp, overwrite = TRUE))
for (ext in c("-wal", "-shm"))
  if (file.exists(paste0(s$db_path, ext)))
    invisible(file.copy(paste0(s$db_path, ext), paste0(tmp, ext), overwrite = TRUE))

con <- DBI::dbConnect(RSQLite::SQLite(), tmp)
e   <- tibble::as_tibble(DBI::dbReadTable(con, "evals"))
DBI::dbDisconnect(con)

cat("store: ", s$db_path, "\n", sep = "")
cat("bytes: ", file.info(s$db_path)$size,
    if (file.exists(paste0(s$db_path, "-wal")))
      sprintf(" (+ %d in -wal)", file.info(paste0(s$db_path, "-wal"))$size) else "",
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
