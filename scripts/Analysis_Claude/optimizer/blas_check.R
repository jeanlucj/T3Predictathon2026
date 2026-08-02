# blas_check.R
#
# Is R's linear algebra actually running on multiple cores? Prints what R is linked against,
# what the thread environment says, and then MEASURES it -- because the first two can look
# right while the third is still single-threaded.
#
#   Rscript blas_check.R                      # as your shell has it configured
#   OMP_NUM_THREADS=8 Rscript blas_check.R    # to see whether the setting takes effect
#
# Reading the result -- "cores busy" is user CPU time over elapsed time on a compute-bound
# matrix multiply:
#   ~1.0            single-threaded. Either reference BLAS, or threads capped to 1.
#   ~N              N threads in use.
#
# Three things can hold it at 1, in the order worth checking:
#
#   1. THE LAUNCHER CAPS IT. run_workers.sh exports OMP/OPENBLAS/MKL/VECLIB _NUM_THREADS from
#      its second argument, which DEFAULTS TO 2 -- on purpose, so N workers do not each grab
#      every core and thrash. `./run_workers.sh 8` therefore uses 8 x 2 = 16 cores, and spare
#      cores are expected. Raise it deliberately: `./run_workers.sh 6 8`.
#   2. REFERENCE BLAS. A BLAS path containing "libRblas" is single-threaded by construction
#      and no environment variable changes that. On BioHPC, R builds from 4.4.3 onward are
#      compiled with OpenBLAS; 4.0.5-4.4.2 are not.
#   3. NOTHING TO THREAD. This pipeline is not BLAS-bound -- see RUNBOOK.md section 1.

cat("R version   :", R.version.string, "\n")
cat("BLAS        :", extSoftVersion()[["BLAS"]], "\n")
cat("LAPACK      :", La_library(), "\n")
if (requireNamespace("RhpcBLASctl", quietly = TRUE))
  cat("BLAS procs  :", RhpcBLASctl::blas_get_num_procs(), "\n")
cat("cores on box:", parallel::detectCores(), "\n")
for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS"))
  cat(sprintf("%-23s %s\n", v, Sys.getenv(v, "<unset>")))

# A COMPUTE-bound square multiply. The tall-thin products this pipeline actually does are
# memory-bandwidth-bound and would understate threading, so they are the wrong probe here.
n <- 3000
A <- matrix(rnorm(n * n), n); B <- matrix(rnorm(n * n), n)
invisible(A %*% B)                                   # warm up
tm <- system.time(for (i in 1:3) invisible(A %*% B))
busy <- tm[["user.self"]] / max(tm[["elapsed"]], 1e-9)
cat(sprintf("\n%dx%d matmul x3: elapsed %.2fs  user %.2fs  -> ~%.1f cores busy\n",
            n, n, tm[["elapsed"]], tm[["user.self"]], busy))
cat(if (busy < 1.5) {
  "  SINGLE-THREADED. Check the three causes in this file's header, in that order.\n"
} else {
  sprintf("  Threaded, roughly %.0f-way.\n", busy)
})
