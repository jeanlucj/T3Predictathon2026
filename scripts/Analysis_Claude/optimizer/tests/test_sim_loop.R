# test_sim_loop.R
#
# End-to-end proof that the optimization machinery actually LEARNS, run offline
# against the synthetic world in R/evaluate.R (no network, seconds to run). It
# runs the full loop, then checks that the discovered incumbent beats the best
# submitted seed and that scores improve from the random phase to the surrogate
# phase. Run: Rscript tests/test_sim_loop.R

library(tidyverse)
here::i_am("tests/test_sim_loop.R")
source(here::here("settings.R"))
invisible(lapply(list.files(here::here("R"), "[.]R$", full.names = TRUE), source))
source(here::here("run_optimizer.R"))   # defines run_optimizer(); does not auto-run when sourced

set.seed(42)

settings <- modifyList(optimizer_settings(), list(
  simulate    = TRUE,
  db_path     = tempfile(fileext = ".sqlite"),
  stop_file   = tempfile(),               # never created -> no early stop
  report_path = tempfile(fileext = ".md"),
  n_random_init = 25,
  max_iters   = 320,
  checkpoint_every = 50
))

run_optimizer(settings)

con   <- open_store(settings$db_path)
evals <- read_evals(con)
agg   <- aggregate_scores(evals)
inc   <- incumbent_config(agg, settings$incumbent_min_reps)

seed_hashes <- vapply(seed_configs(), config_hash, character(1))
best_seed <- agg |>
  dplyr::filter(config_hash %in% seed_hashes, is.finite(mean_score)) |>
  dplyr::summarise(m = max(mean_score)) |> dplyr::pull(m)

# Did the surrogate phase outscore the random phase? Compare mean score of the
# first n_random_init scored evals vs the last 50.
ok_ev <- evals |> dplyr::filter(is.finite(score)) |> dplyr::arrange(id)
early <- mean(head(ok_ev$score, settings$n_random_init))
late  <- mean(tail(ok_ev$score, 50))

cat("\n================ SIM LOOP RESULT ================\n")
cat(sprintf("evaluations stored:        %d (%d scored)\n", nrow(evals), nrow(ok_ev)))
cat(sprintf("best seed (submission):    %.3f\n", best_seed))
cat(sprintf("incumbent mean score:      %.3f over %d trials\n", inc$mean_score, inc$n_ok))
cat(sprintf("early random mean:         %.3f\n", early))
cat(sprintf("late surrogate mean:       %.3f\n", late))
cat("incumbent config:\n"); cat(format_config(inc$config), "\n")

fail <- 0L
if (!(inc$mean_score > best_seed)) {
  cat("\nFAIL: incumbent did not beat the best submitted seed\n"); fail <- 1L
}
if (!(late > early)) {
  cat("FAIL: surrogate phase did not improve on the random phase\n"); fail <- 1L
}
close_store(con)
if (fail == 0L) cat("\nPASS: optimizer beats submissions and improves over random search.\n")
quit(status = fail)
