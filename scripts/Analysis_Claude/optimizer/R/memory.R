# memory.R
#
# How much memory does ONE evaluation take? The answer decides two things that cannot be
# settled by reading the code: what `dosage_budget_bytes` this machine can afford, and how
# many workers fit in its RAM (see README "Running several workers").
#
# The measure that matters here is the PEAK, not the footprint at any instant: the dosage
# matrices are allocated, merged, QC'd and discarded inside a single evaluation, so a sample
# taken between evaluations sees almost nothing. gc(reset = TRUE) at the start of an
# evaluation and gc() at the end brackets exactly that peak.
#
# R-heap peak vs process RSS -- both are recorded because they answer different questions:
#   mem_peak_mb()      peak R heap since the last reset. PER-EVALUATION, attributable to a
#                      configuration, and the right number for sizing a budget (essentially
#                      all of the memory here is R-allocated matrices).
#   proc_rss_mb()      the process's resident size right now, including what R has allocated
#                      but not returned to the OS, plus BLAS/loaded libraries. Always >= the
#                      heap number, and it is what the OOM killer looks at.
#   proc_peak_rss_mb() the process's high-water RSS. MONOTONIC over the life of the process,
#                      so it is a RUN-level number, not a per-evaluation one -- do not
#                      difference it between evaluations.
# For the machine as a whole while a run is going, use monitor_memory.sh instead.

library(tidyverse)

# Start a new measurement window. `full = TRUE` forces a complete collection so the peak that
# follows is attributable to this evaluation rather than to garbage left by the previous one.
mem_reset <- function() invisible(gc(reset = TRUE, full = TRUE))

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
