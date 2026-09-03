# memory.R
#
# Peak memory for ONE evaluation, which is what decides dosage_budget_bytes and how many
# workers fit in RAM. The PEAK is the measure: dosage matrices are allocated, merged, QC'd and
# discarded within one evaluation, so a sample taken between evaluations sees almost nothing.
#
# SIZE A MACHINE FROM mem_peak_rss_mb(), NOT mem_peak_mb(). On one 8-worker run the two
# disagreed by 2.7x (33.6 GB heap vs 89.4 GB RSS): gc() sees only what R allocated itself, and
# misses the BLAS workspace behind tcrossprod, sommer's Armadillo objects, stats::dist, pages
# R has not returned, and the interpreter itself.
#
#   mem_peak_mb()       peak R heap since reset. A LOWER BOUND; kept because its ratio to the
#                       RSS peak measures what gc() cannot see.
#   mem_peak_rss_mb()   peak RSS since reset -- the honest figure, and what the OOM killer
#                       acts on. Linux only; NA elsewhere.
#   proc_rss_mb()       resident size now, not a peak. A cheap sanity check.
#   proc_peak_rss_mb()  the raw VmHWM read, monotonic over the process life unless reset.
#
# The per-evaluation RSS peak is possible because Linux lets a process reset its own high-water
# mark: writing "5" to /proc/self/clear_refs (CLEAR_REFS_MM_HIWATER_RSS, kernel >= 4.0) sets
# VmHWM back down to the current VmRSS. mem_reset() does that, so VmHWM read at the end of an
# evaluation is that evaluation's peak rather than the whole worker's.
#
# For the machine as a whole while a run is going, use tools/watch_memory.sh instead.

library(tidyverse)

# --- per-evaluation RSS peak (Linux) ---------------------------------------
# Ask the kernel to drop this process's high-water RSS to its current RSS. Returns whether the
# write went through; whether the kernel HONOURED it is a separate question (see below).
.rss_peak_reset <- function() {
  if (!file.exists("/proc/self/clear_refs")) return(FALSE)
  isTRUE(tryCatch({ cat("5\n", file = "/proc/self/clear_refs"); TRUE },
                  error = function(e) FALSE, warning = function(w) FALSE))
}

# Probed ONCE per session and cached (same in-RAM pattern as .vcf_download_fails in
# data_access.R). The probe is necessary rather than paranoid: a kernel that ignores the write
# leaves VmHWM at its old value, and VmHWM only ever climbs -- so every row would silently
# record the whole worker's high-water mark instead of the evaluation's, which looks like a
# plausible number and is not one. Verified by resetting and checking VmHWM actually fell to
# meet VmRSS.
.mem_env <- new.env(parent = emptyenv())

.rss_peak_resettable <- function() {
  if (!is.null(.mem_env$resettable)) return(.mem_env$resettable)
  ok <- FALSE
  if (.rss_peak_reset()) {
    hwm <- proc_peak_rss_mb(); rss <- proc_rss_mb()
    # Allow a little slack: RSS can move between the two reads.
    ok <- is.finite(hwm) && is.finite(rss) && hwm <= rss * 1.05 + 16
  }
  if (!ok)
    message("memory: this platform will not reset the peak-RSS counter -- ",
            "peak_rss_mb will be NA and peak_r_mb (an UNDER-estimate) is all that is recorded")
  .mem_env$resettable <- ok
  ok
}

# Start a new measurement window. `full = TRUE` forces a complete collection so the peak that
# follows is attributable to this evaluation rather than to garbage left by the previous one;
# the clear_refs write does the same for the kernel's RSS high-water mark.
mem_reset <- function() {
  if (.rss_peak_resettable()) .rss_peak_reset()
  invisible(gc(reset = TRUE, full = TRUE))
}

# Peak RSS in MB since the last mem_reset(). NA where the counter cannot be reset -- a
# monotonic run-level high-water mark would be worse than nothing, because it is indistinguishable
# from a per-evaluation peak once it is in the store.
mem_peak_rss_mb <- function() {
  if (!.rss_peak_resettable()) return(NA_real_)
  proc_peak_rss_mb()
}

# Peak R heap in MB since the last mem_reset(), summed over the two cell kinds (Ncells =
# language objects, Vcells = the vectors and matrices, which is where the dosages live).
#
# The column is found BY NAME: gc()'s matrix has a variable number of columns -- a
# "limit (Mb)" column appears only where a memory limit is set (it is there on macOS, absent
# on a stock Linux build), so a positional index silently reads the wrong column on one
# platform or the other. "max used" is always followed by its "(Mb)" companion.
mem_peak_mb <- function() {
  g <- tryCatch(gc(verbose = FALSE), error = function(e) NULL)
  if (is.null(g)) return(NA_real_)
  j <- which(colnames(g) == "max used")
  if (!length(j) || j[1] + 1L > ncol(g)) return(NA_real_)
  sum(g[, j[1] + 1L], na.rm = TRUE)
}

# Resident set size of THIS process, in MB. /proc reports kB directly, so there is no page-size
# assumption to get wrong. NA on a platform with neither /proc nor ps (nothing downstream
# requires it -- the heap peak is the number the budget is derived from).
.proc_status_kb <- function(field) {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  ln <- tryCatch(grep(paste0("^", field, ":"), readLines("/proc/self/status", warn = FALSE),
                      value = TRUE), error = function(e) character())
  if (!length(ln)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", ln[1]))
}

proc_rss_mb <- function() {
  kb <- .proc_status_kb("VmRSS")                                  # Linux
  if (is.na(kb) && .Platform$OS.type == "unix") {                 # macOS and friends
    kb <- tryCatch(as.numeric(system2("ps", c("-o", "rss=", "-p", Sys.getpid()),
                                      stdout = TRUE, stderr = FALSE)),
                   error = function(e) NA_real_, warning = function(w) NA_real_)
  }
  if (length(kb) != 1 || is.na(kb)) NA_real_ else kb / 1024
}

# Process high-water RSS in MB (Linux only; NA elsewhere). Monotonic -- a run-level figure.
proc_peak_rss_mb <- function() {
  kb <- .proc_status_kb("VmHWM")
  if (is.na(kb)) NA_real_ else kb / 1024
}

# Compact "1.4 GB" / "820 MB" for log lines and messages.
fmt_mb <- function(mb) {
  if (!length(mb) || !is.finite(mb)) return("NA")
  if (mb >= 1024) sprintf("%.1fGB", mb / 1024) else sprintf("%.0fMB", mb)
}
