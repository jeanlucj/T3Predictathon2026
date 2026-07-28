# report_timing.R
#
# Where does an evaluation's wall time actually go? Reads the results store and decomposes
# it, so a "the optimizer is slow" impression becomes a number attached to a subtask method.
#
# Every row of `evals` carries `seconds` -- the wall time INSIDE evaluate_config_on_trial(),
# i.e. the whole six-subtask pipeline -- and `ts`, written when the row is stored. So for
# consecutive rows,
#
#     gap_i     = ts_i - ts_(i-1)          (end of one evaluation to the end of the next)
#     inside_i  = seconds_i                (the pipeline)
#     outside_i = gap_i - seconds_i        (trial sampling + choose_config's surrogate fit
#                                           + the report/checkpoint write)
#
# A large `outside` says the cost is NOT in the pipeline, which points at trial sampling or
# the surrogate rather than at genotype handling -- worth knowing before optimizing the
# wrong thing.
#
# READ-ONLY: safe to run while the optimizer is going.
# Run on the server (R 4.6.1) from the optimizer directory:  Rscript report_timing.R

suppressMessages(library(tidyverse))
options(width = 200)   # keep the tables from wrapping

store_path <- local({
  p <- Sys.getenv("OPTIMIZER_PATH")
  cand <- if (nzchar(p)) file.path(p, "state", "evals.sqlite") else file.path("state", "evals.sqlite")
  if (!file.exists(cand) && file.exists(file.path("state", "evals.sqlite")))
    cand <- file.path("state", "evals.sqlite")
  cand
})
if (!file.exists(store_path))
  stop("no store at ", store_path,
       "\nRun from the optimizer directory, or set OPTIMIZER_PATH (check with ",
       "Sys.getenv(\"OPTIMIZER_PATH\")).")

con <- DBI::dbConnect(RSQLite::SQLite(), store_path)
e <- DBI::dbGetQuery(con, "SELECT id, config_json, trial_id, scheme, status, reason,
                                  score, seconds, ts FROM evals ORDER BY id")
DBI::dbDisconnect(con)
if (!nrow(e)) stop("the store is empty: ", store_path)

# Pull one field out of each stored configuration. config_to_json() writes with
# auto_unbox = TRUE, so these are plain scalars.
fld <- function(field) {
  vapply(e$config_json, function(j) {
    v <- tryCatch(jsonlite::fromJSON(j)[[field]], error = function(err) NA)
    if (is.null(v) || length(v) != 1L) NA_character_ else as.character(v)
  }, character(1), USE.NAMES = FALSE)
}

e$ts     <- as.POSIXct(e$ts, tz = "UTC")
e$gap    <- c(NA_real_, as.numeric(diff(e$ts), units = "secs"))
e$inside <- e$seconds / 60
# `outside` is only meaningful for consecutive evaluations of the SAME continuous run. A gap
# that spans a restart carries the downtime with it, and a gap shorter than the evaluation it
# follows is impossible -- both are dropped rather than allowed to distort the summary.
plausible <- is.finite(e$gap) & is.finite(e$seconds) & e$gap >= e$seconds &
             e$gap < 6 * 3600
e$outside <- ifelse(plausible, (e$gap - e$seconds) / 60, NA_real_)
n_dropped <- sum(is.finite(e$gap)) - sum(plausible)
e$kernel  <- fld("kernel.method")
e$geno    <- fld("geno_select.method")
e$model   <- fld("model.method")
e$train   <- fld("train_select.method")

hdr <- function(x) cat("\n", x, "\n", strrep("-", nchar(x)), "\n", sep = "")

cat(sprintf("store: %s\nrows: %d   span: %s -> %s   (%.1f h)\n",
            store_path, nrow(e), format(min(e$ts)), format(max(e$ts)),
            as.numeric(difftime(max(e$ts), min(e$ts), units = "hours"))))

hdr("MINUTES INSIDE the pipeline (median), kernel x geno_select")
print(round(tapply(e$inside, list(e$kernel, e$geno), median, na.rm = TRUE), 1))
hdr("...and how many evaluations sit in each cell")
print(table(e$kernel, e$geno))

hdr("MINUTES INSIDE the pipeline (median), by model and by train_select")
print(round(tapply(e$inside, e$model, median, na.rm = TRUE), 1))
print(round(tapply(e$inside, e$train, median, na.rm = TRUE), 1))

hdr("MINUTES OUTSIDE the pipeline  (gap - seconds: sampling + surrogate + checkpoint)")
ok <- is.finite(e$outside)
if (!any(ok)) {
  cat("  no usable consecutive pairs (every gap spans a restart or is implausible)\n")
} else {
  print(round(summary(e$outside[ok]), 1))
  # Both in minutes: outside is already /60, gap is still in seconds.
  share <- sum(e$outside[ok]) / (sum(e$gap[ok]) / 60)
  cat(sprintf("\n  median outside = %.1f min vs median inside = %.1f min\n",
              median(e$outside[ok]), median(e$inside, na.rm = TRUE)))
  cat(sprintf("  outside accounts for %.0f%% of wall time over those pairs\n", 100 * share))
  cat("  (a large share means the GRM path is NOT the bottleneck)\n")
}
if (n_dropped > 0)
  cat(sprintf("  [%d of %d gaps ignored: spanned a restart, exceeded 6 h, or ran shorter\n",
              n_dropped, sum(is.finite(e$gap))),
      "   than the evaluation they follow]\n", sep = "")

hdr("status / reason mix")
print(table(e$status, useNA = "ifany"))
if (any(!is.na(e$reason))) print(sort(table(e$reason[!is.na(e$reason)]), decreasing = TRUE))

hdr("slowest 10 evaluations")
print(e |>
        dplyr::arrange(dplyr::desc(inside)) |>
        dplyr::transmute(id, trial_id, kernel, geno, model, scheme, status,
                         inside_min = round(inside, 1)) |>
        head(10), row.names = FALSE)

hdr("last 15 evaluations")
print(e |>
        dplyr::transmute(id, kernel, geno, scheme, status, reason,
                         inside_min = round(inside, 1), outside_min = round(outside, 1)) |>
        tail(15), row.names = FALSE)
