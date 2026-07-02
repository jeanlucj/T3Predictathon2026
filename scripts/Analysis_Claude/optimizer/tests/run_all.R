# run_all.R  --  run every tests/test_*.R in a fresh process and aggregate.
#
#   Rscript tests/run_all.R          # fast test files only
#   Rscript tests/run_all.R --all    # also run the slow end-to-end sim loop
#
# Each test file exits non-zero on failure; this reports a single PASS/FAIL per
# file and exits non-zero if any failed.

suppressMessages(library(tidyverse))
here::i_am("tests/run_all.R")

run_slow <- "--all" %in% commandArgs(trailingOnly = TRUE)
slow_files <- c("test_sim_loop.R")     # minutes, not seconds

files <- list.files(here::here("tests"), pattern = "^test_.*[.]R$")
if (!run_slow) files <- setdiff(files, slow_files)
files <- sort(files)

cat("Running", length(files), "test file(s)",
    if (!run_slow) "(skipping slow; pass --all to include)" else "", "\n\n")

status <- vapply(files, function(f) {
  cat(strrep("=", 64), "\n", f, "\n", strrep("=", 64), "\n", sep = "")
  s <- system2("Rscript", here::here("tests", f), stdout = "", stderr = "")
  cat("\n[", if (s == 0) "PASS" else "FAIL", "]", f, "\n\n")
  s
}, FUN.VALUE = integer(1))

cat(strrep("=", 64), "\nSUMMARY\n", strrep("=", 64), "\n", sep = "")
for (f in names(status)) cat(sprintf("  %-26s %s\n", f, if (status[f] == 0) "PASS" else "FAIL"))
cat(sprintf("\n%d/%d test files passed\n", sum(status == 0), length(status)))
quit(status = if (any(status != 0)) 1L else 0L)
