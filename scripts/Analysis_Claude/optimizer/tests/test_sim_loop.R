# test_sim_loop.R
#
# End-to-end proof that the search ALGORITHM works, run offline against the synthetic world
# in R/evaluate.R (no network). Run: Rscript tests/test_sim_loop.R
#
# Two questions have to be kept apart here, and conflating them is what made an earlier
# version of this file unreliable:
#
#   1. Does the algorithm work?  A correctness property. Test it where a correct
#      implementation MUST win: a deterministic objective. If the surrogate cannot beat
#      random search when every evaluation is exact, the search is broken and no budget
#      rescues it. THIS FILE ASSERTS THAT.
#   2. How many evaluations does it need at a given noise level?  An empirical property of
#      the objective, not of the code. Measured below, asserted nowhere -- gating on it
#      makes the suite fail on the RNG instead of on a regression.
#
# .sim_true() is deterministic in (config, trial, scheme). Two things make an EVALUATION of
# it noisy, and both are turned off here (see settings$sim_noise_sd / $sim_fixed_trial):
# the observation noise .sim_evaluate() adds, and -- the larger source -- the fact that the
# objective is a mean over a heterogeneous trial population sampled one trial at a time.

library(tidyverse)
here::i_am("tests/test_sim_loop.R")
source(here::here("settings.R"))
invisible(lapply(list.files(here::here("R"), "[.]R$", full.names = TRUE), source))
source(here::here("run_optimizer.R"))   # defines run_optimizer(); does not auto-run when sourced

set.seed(42)

settings <- modifyList(optimizer_settings(), list(
  simulate    = TRUE,
  sim_noise_sd    = 0,      # exact observations
  sim_fixed_trial = TRUE,   # one trial, so the objective is a fixed function of the config
  config_replication = 1,   # with one trial and no noise a second eval is an exact duplicate
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

# The objective this run was optimizing, evaluated exactly.
trial   <- sample_trial(settings)                       # the fixed descriptor
quality <- function(cfg) .sim_true(cfg, trial, settings$optimize_scheme)

inc_q  <- quality(inc$config)
seed_q <- max(vapply(seed_configs(), quality, numeric(1)))

# What RANDOM SEARCH would deploy on the same budget. With exact evaluations its choice is
# simply the best config it drew, so the arm can be computed directly instead of run: the
# maximum true quality over ~300 random draws (the run evaluates 320, a few of them repeats
# and seeds). Averaged over 3 replicates at a fixed seed so the bar is reproducible.
# Verified against actually running the random arm through run_optimizer(): analytic
# 0.557 +- 0.019 vs measured 0.530 / 0.568 / 0.559 over three seeds.
set.seed(20260727)
rand_bar <- mean(replicate(3, max(vapply(seq_len(300),
                                         function(i) quality(sample_config()), numeric(1)))))

cat("\n================ SIM LOOP RESULT (deterministic regime) ================\n")
cat(sprintf("evaluations stored:        %d\n", nrow(evals)))
cat(sprintf("incumbent true quality:    %.3f over %d reps\n", inc_q, inc$n_ok))
cat(sprintf("random search would get:   %.3f   [ASSERTED: incumbent > this]\n", rand_bar))
cat(sprintf("best submitted seed:       %.3f   [ASSERTED: incumbent > this]\n", seed_q))
cat("incumbent config:\n"); cat(format_config(inc$config), "\n")

fail <- 0L
# Replication, at whatever level this regime asks for. Gating on incumbent_min_reps instead
# would fail here by construction, since config_replication is 1 in the deterministic regime.
if (!(nrow(evals) > 0 && !is.null(inc) && inc$n_ok >= settings$config_replication)) {
  cat("\nFAIL: no incumbent with the required replication after a full run\n"); fail <- 1L
}
# The correctness claim: with an exact objective the surrogate must beat equal-budget random
# search. Measured margin is ~+0.05 (0.606 vs 0.553 over seeds 42/7/2024), so a bare
# inequality has room without being brittle.
if (!(inc_q > rand_bar)) {
  cat(sprintf("\nFAIL: incumbent %.3f did not beat equal-budget random search %.3f\n",
              inc_q, rand_bar)); fail <- 1L
}
if (!(inc_q > seed_q)) {
  cat(sprintf("FAIL: incumbent %.3f did not beat the best submitted seed %.3f\n",
              inc_q, seed_q)); fail <- 1L
}

# ---------------------------------------------------------------------------------------
# Question 2, measured and NOT asserted. Same fair experiment (two arms, 320 evaluations
# each, both deploying the incumbent their own evaluations selected), three regimes, seeds
# 42/7/2024, quality scored on the objective each regime was optimizing:
#
#   regime                                 surrogate   random   diff    surrogate wins
#   deterministic (this file)                  0.606    0.553  +0.053            3/3
#   trial variation only (no obs. noise)       0.513    0.525  -0.012            1/3
#   full defaults (obs. noise + trial var)     0.527    0.528  -0.001            2/3
#
# So the algorithm is sound, and its advantage is spent entirely on variance at a
# 320-evaluation budget. Note the middle row: with observation noise switched OFF the search
# is still no better than random, so the binding constraint is TRIAL HETEROGENEITY -- the
# objective is a mean over a heterogeneous population estimated one draw at a time -- not
# measurement noise. Practically, replicating each config across more trials is what would
# help; making individual evaluations more precise would not. See BACKGROUND.md sec. 4.
# ---------------------------------------------------------------------------------------

close_store(con)

# ---------------------------------------------------------------------------------------
# Liveness: a run of iterations that store NOTHING must stop the worker.
#
# A skip is ordinary on its own -- that configuration had no trial left. An unbroken run of
# them means the search cannot place work anywhere, and it looks identical from the outside to
# a healthy run: same iteration rate, same backup age. Without this bound a worker spins until
# the scheduler kills it. optimizer_step is stubbed because the point is the loop's response,
# not any particular reason for the skip.
cat("\nliveness: consecutive skips halt the run\n")
real_step <- optimizer_step
optimizer_step <- function(con, settings, conn = NULL)
  list(source = "replicate", skipped = TRUE, trial = "simtrial_fixed")
msgs <- character()
withCallingHandlers(
  run_optimizer(modifyList(settings, list(
    db_path = tempfile(fileext = ".sqlite"), stop_file = tempfile(),
    report_path = tempfile(fileext = ".md"),
    max_iters = 500, max_consec_skip = 5, checkpoint_every = 1000))),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
optimizer_step <- real_step

iters <- suppressWarnings(as.integer(sub(".*iter ([0-9]+).*", "\\1",
                                         grep("iter [0-9]+", msgs, value = TRUE))))
if (!any(grepl("nothing evaluated in", msgs))) {
  cat("\nFAIL: an unbroken run of skips did not halt the worker\n"); fail <- 1L
}
if (!(length(iters) && max(iters, na.rm = TRUE) <= 6L)) {
  cat(sprintf("\nFAIL: ran %d iterations on a max_consec_skip of 5\n",
              if (length(iters)) max(iters, na.rm = TRUE) else 0L)); fail <- 1L
}

if (fail == 0L) cat("\nPASS: the surrogate search beats random search and the submissions on an exact objective.\n")
quit(status = fail)
