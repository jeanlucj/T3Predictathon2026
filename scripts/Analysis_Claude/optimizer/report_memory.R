# report_memory.R
#
# What does one evaluation cost in MEMORY? The answer sets the two numbers that cannot be
# reasoned out from the code: how big `dosage_budget_bytes` can be on this machine, and how
# many workers fit in its RAM.
#
# Reads peak_r_mb (peak R heap for that evaluation, from R/memory.R) and rss_mb off the
# store. Both are NA for rows written before the instrumentation existed -- those are
# reported as a count, not silently dropped.
#
# The headline is the "workers that fit" table at the end. It is deliberately driven by the
# MAXIMUM observed peak rather than the median: a worker that dies of an out-of-memory kill
# takes its whole evaluation with it, so the tail is the number that matters, not the centre.
#
# READ-ONLY: safe to run while the optimizer is going (and in WAL mode a reader never blocks
# the writers). Run from the optimizer directory:
#   Rscript report_memory.R                 # the live store
#   Rscript report_memory.R state/old.sqlite   # an archived one
#
# For the machine as a whole in real time -- all workers plus everything else on the node --
# use monitor_memory.sh, which needs no restart to attach to a running job.

suppressMessages(library(tidyverse))
options(width = 200)

store_path <- local({
  arg <- commandArgs(trailingOnly = TRUE)
  if (length(arg) && nzchar(arg[1])) return(arg[1])
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
have <- DBI::dbListFields(con, "evals")
if (!("peak_r_mb" %in% have)) {
  DBI::dbDisconnect(con)
  stop("this store has no peak_r_mb column -- it predates the memory instrumentation.\n",
       "Run the optimizer once with the current code (open_store migrates the table in ",
       "place), then re-run this report.")
}
e <- DBI::dbGetQuery(con, "SELECT id, config_json, trial_id, scheme, status, score, seconds,
                                  peak_r_mb, rss_mb, worker, dosage_budget, ts
                           FROM evals ORDER BY id")
DBI::dbDisconnect(con)
if (!nrow(e)) stop("the store is empty: ", store_path)

# One field out of each stored configuration (config_to_json writes plain scalars).
fld <- function(field) {
  vapply(e$config_json, function(j) {
    v <- tryCatch(jsonlite::fromJSON(j)[[field]], error = function(err) NA)
    if (is.null(v) || length(v) != 1L) NA_character_ else as.character(v)
  }, character(1), USE.NAMES = FALSE)
}
e$kernel <- fld("kernel.method")
e$geno   <- fld("geno_select.method")
e$ts     <- as.POSIXct(e$ts, tz = "UTC")
e$peak_gb <- e$peak_r_mb / 1024

hdr <- function(x) cat("\n", x, "\n", strrep("-", nchar(x)), "\n", sep = "")

m <- dplyr::filter(e, is.finite(peak_r_mb))
n_missing <- nrow(e) - nrow(m)
cat(sprintf("store: %s\nrows: %d   with memory recorded: %d   span: %s -> %s\n",
            store_path, nrow(e), nrow(m), format(min(e$ts)), format(max(e$ts))))
if (n_missing > 0)
  cat(sprintf("  [%d row(s) predate the instrumentation and carry no memory figure]\n",
              n_missing))
if (!nrow(m))
  stop("no row has a memory figure yet -- let the optimizer complete an evaluation first.")

# Marker density is not a config parameter, so rows made under different dosage budgets are
# not comparable on memory (nor, strictly, on score). Say so before summarizing across them.
buds <- sort(unique(m$dosage_budget[is.finite(m$dosage_budget)]))
if (length(buds) > 1) {
  hdr("WARNING: several dosage budgets in this store")
  print(m |>
          dplyr::filter(is.finite(dosage_budget)) |>
          dplyr::group_by(dosage_budget_gb = round(dosage_budget / 1e9, 1)) |>
          dplyr::summarise(n = dplyr::n(),
                           median_peak_gb = round(median(peak_gb), 1),
                           max_peak_gb = round(max(peak_gb), 1), .groups = "drop") |>
          as.data.frame(), row.names = FALSE)
  cat("\n  Rows made at different marker densities are not comparable. Size the machine\n",
      "  from the rows matching the budget you intend to RUN at.\n", sep = "")
}

hdr("PEAK R HEAP PER EVALUATION (GB)")
print(round(quantile(m$peak_gb, c(0, .5, .75, .9, .95, .99, 1)), 2))

hdr("peak GB by kernel x geno_select (median / max)")
print(round(tapply(m$peak_gb, list(m$kernel, m$geno), median), 2))
cat("\nmax:\n")
print(round(tapply(m$peak_gb, list(m$kernel, m$geno), max), 2))

hdr("peak GB by outcome -- does a failure cost as much as a success?")
print(m |>
        dplyr::group_by(status) |>
        dplyr::summarise(n = dplyr::n(),
                         median_gb = round(median(peak_gb), 2),
                         max_gb = round(max(peak_gb), 2), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(max_gb)) |>
        as.data.frame(), row.names = FALSE)

if (any(is.finite(m$rss_mb))) {
  hdr("process RSS at end of evaluation (GB) -- heap plus what R has not returned to the OS")
  r <- m$rss_mb[is.finite(m$rss_mb)] / 1024
  print(round(quantile(r, c(0.5, 0.9, 1)), 2))
  cat(sprintf("\n  RSS/heap ratio at the median: %.2fx\n",
              median(r) / median(m$peak_gb)))
}

if (length(unique(m$worker[!is.na(m$worker)])) > 1) {
  hdr("by worker -- confirms the workers are sharing the store, and are alike")
  print(m |>
          dplyr::group_by(worker) |>
          dplyr::summarise(n = dplyr::n(),
                           median_gb = round(median(peak_gb), 2),
                           max_gb = round(max(peak_gb), 2),
                           median_min = round(median(seconds, na.rm = TRUE) / 60, 1),
                           .groups = "drop") |>
          as.data.frame(), row.names = FALSE)
}

hdr("hungriest 10 evaluations")
print(m |>
        dplyr::arrange(dplyr::desc(peak_gb)) |>
        dplyr::transmute(id, trial_id, kernel, geno, scheme, status,
                         peak_gb = round(peak_gb, 2),
                         min = round(seconds / 60, 1)) |>
        head(10), row.names = FALSE)

# ---------------------------------------------------------------------------
# The point of the report: how many workers fit.
#
# Sizing on the max, not the median, is deliberate -- workers do not take turns, so the
# machine has to survive several of them being in a heavy evaluation at once. The 30%
# headroom covers the OS, the page cache (large here: VCF and RDS traffic), and the
# fact that R returns memory to the OS lazily.
# ---------------------------------------------------------------------------
hdr("HOW MANY WORKERS FIT")
# Simulate-mode evaluations touch no genotypes at all, so their peak is just the baseline R
# heap. Sizing a machine from them would suggest hundreds of workers fit. Say so rather than
# print a number that is off by two orders of magnitude.
n_sim <- sum(grepl("^simtrial_", m$trial_id))
if (n_sim > 0)
  cat(sprintf(paste0("NOTE: %d of %d rows are SIMULATE-mode (simtrial_*), which allocate no\n",
                     "      genotype data. The figures below are meaningless unless they come\n",
                     "      from real evaluations.\n\n"), n_sim, nrow(m)))

p95  <- quantile(m$peak_gb, 0.95)
pmax_ <- max(m$peak_gb)
cat(sprintf("observed peak per evaluation: p95 %.1f GB, max %.1f GB (n = %d)\n\n",
            p95, pmax_, nrow(m)))
fits <- tibble::tibble(ram_gb = c(64, 128, 256, 512, 1024)) |>
  dplyr::mutate(usable_gb   = 0.7 * ram_gb,
                workers_p95 = pmax(1, floor(usable_gb / p95)),
                workers_max = pmax(1, floor(usable_gb / pmax_)))
print(as.data.frame(fits |> dplyr::mutate(usable_gb = round(usable_gb))), row.names = FALSE)
cat("\n  workers_max is the safe number; workers_p95 assumes ~1 evaluation in 20 is heavy\n",
    "  and that they rarely coincide. Prefer workers_max unless you are watching the run.\n",
    "  Raising dosage_budget_bytes raises both peaks -- re-run this after changing it.\n",
    "  dosage_total_budget_bytes is what actually BOUNDS the peak; if max above is well\n",
    "  under it, the cap is not binding and there is room to raise the density budget.\n",
    sep = "")
