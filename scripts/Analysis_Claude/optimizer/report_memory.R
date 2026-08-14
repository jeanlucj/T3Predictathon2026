# report_memory.R
#
# What does one evaluation cost in MEMORY? The answer sets the two numbers that cannot be
# reasoned out from the code: how big `dosage_budget_bytes` can be on this machine, and how
# many workers fit in its RAM.
#
# Sizes from peak_rss_mb -- the evaluation's true peak RSS (R/memory.R). peak_r_mb, R's own
# heap peak, is ALSO reported but must not be sized from: gc() cannot see BLAS, sommer or dist
# allocations, and on one real run it understated the requirement by 2.7x. Rows written before
# either column existed are reported as a count, not silently dropped.
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

# Which store to read is resolve_read_store()'s decision -- ASK IT, do not rebuild the path
# from OPTIMIZER_HOME. settings.local.R commonly moves db_path (onto /workdir, or onto a SLURM
# node's $TMPDIR), and guessing would silently report on a stale copy no worker is writing to.
# On a cluster it falls back to the durable backup, saying so, because the live store is
# node-local to the job that wrote it.
source(here::here("settings.R"))
source(here::here("R", "store.R"))
store_path <- resolve_read_store(commandArgs(trailingOnly = TRUE)[1], "report_memory.R")

con <- DBI::dbConnect(RSQLite::SQLite(), store_path)
have <- DBI::dbListFields(con, "evals")
if (!("peak_r_mb" %in% have)) {
  DBI::dbDisconnect(con)
  stop("this store has no peak_r_mb column -- it predates the memory instrumentation.\n",
       "Run the optimizer once with the current code (open_store migrates the table in ",
       "place), then re-run this report.")
}
# peak_rss_mb postdates peak_r_mb, so a store may have one and not the other.
has_rss_peak <- "peak_rss_mb" %in% have
e <- DBI::dbGetQuery(con, sprintf("SELECT id, config_json, trial_id, scheme, status, score,
                                          seconds, peak_r_mb, rss_mb, worker, dosage_budget, ts%s
                                   FROM evals ORDER BY id",
                                  if (has_rss_peak) ", peak_rss_mb" else ""))
if (!has_rss_peak) e$peak_rss_mb <- NA_real_
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

hdr <- function(x) cat("\n", x, "\n", strrep("-", nchar(x)), "\n", sep = "")

# SIZE FROM RSS. peak_r_mb is R's heap peak and cannot see BLAS/sommer/dist allocations, so it
# UNDER-states the requirement -- measured at 2.7x on one real run. Use it only when the honest
# figure is absent, and say so loudly, because this is the number a machine gets sized from.
use_rss  <- any(is.finite(e$peak_rss_mb))
e$peak_gb <- if (use_rss) e$peak_rss_mb / 1024 else e$peak_r_mb / 1024
peak_src <- if (use_rss) "peak_rss_mb (true process peak)" else "peak_r_mb (R heap -- AN UNDER-ESTIMATE)"

m <- dplyr::filter(e, is.finite(peak_gb))
n_missing <- nrow(e) - nrow(m)
cat(sprintf("store: %s\nrows: %d   with memory recorded: %d   span: %s -> %s\nsizing from: %s\n",
            store_path, nrow(e), nrow(m), format(min(e$ts)), format(max(e$ts)), peak_src))
if (!use_rss)
  cat("\n!! No peak_rss_mb in this store, so every figure below comes from R's heap peak and\n",
      "   UNDERSTATES real memory use -- by 2.7x on the one run where both were measured.\n",
      "   Do NOT size a worker count from this. Re-run the optimizer with current code on\n",
      "   Linux to record peak_rss_mb, and meanwhile use rss_max_mb from monitor_memory.sh.\n",
      sep = "")
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

hdr(if (use_rss) "PEAK RSS PER EVALUATION (GB) -- the sizing figure" else "PEAK R HEAP PER EVALUATION (GB) -- an under-estimate")
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

# The blind spot, quantified: how much of an evaluation's real footprint gc() cannot account
# for. A large ratio means the cost is in compiled code (GRM/BLAS/sommer), which
# dosage_total_budget_bytes does NOT bound -- that setting caps dosage bytes only.
both <- dplyr::filter(m, is.finite(peak_rss_mb), is.finite(peak_r_mb), peak_r_mb > 0)
if (nrow(both)) {
  hdr("what gc() CANNOT see: peak RSS / peak R heap")
  rr <- both$peak_rss_mb / both$peak_r_mb
  print(round(quantile(rr, c(0, .5, .9, 1)), 2))
  cat(sprintf("\n  median %.1fx -- of a %.1f GB evaluation, ~%.1f GB is invisible to peak_r_mb\n",
              median(rr), max(both$peak_rss_mb) / 1024,
              (max(both$peak_rss_mb) - max(both$peak_r_mb)) / 1024))
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
