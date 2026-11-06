# settings.R -- every knob in one place.
#
# Machine-specific values go in settings.local.R (gitignored), never here.

# Stamped into the log, the report and every stored eval. Last digit for a minor change, middle
# for a major one. If a bump changes what a configuration COMPUTES, add a BUILD_CHANGES entry
# (R/optimizer.R) so the invalidated rows can be retired.
OPTIMIZER_BUILD <- "0.9.0"

# TRUE puts state, logs and the cache backup on durable storage instead of the work directory.
# OPTIMIZER_REMOTE decides if set; otherwise TRUE whenever OPTIMIZER_HOME is set.
.detect_remote_server <- function() {
  ovr <- tolower(trimws(Sys.getenv("OPTIMIZER_REMOTE")))
  if (nzchar(ovr)) return(ovr %in% c("1", "true", "yes", "on"))
  nzchar(Sys.getenv("OPTIMIZER_HOME"))
}
remote_server <- .detect_remote_server()

# Read settings.local.R, which defines `settings_override <- list(...)` layered over the
# defaults below.
.local_overrides <- function(path = here::here("settings.local.R")) {
  if (!file.exists(path)) return(list())
  e <- new.env(parent = globalenv())
  tryCatch(sys.source(path, envir = e),
           error = function(err) stop(basename(path), " failed to source: ", conditionMessage(err)))
  ov <- get0("settings_override", envir = e, ifnotfound = list())
  if (!is.list(ov)) stop(basename(path), " must define `settings_override` as a list()")
  ov
}

# Reject a path given the wrong KIND of thing. `*_dir` names a directory; `*_file` and `*_path`
# name a file, filename included.
.validate_paths <- function(s) {
  file_settings <- c("db_path", "db_backup_path", "report_path", "stop_file", "cache_ready_file")
  dir_settings  <- c("cache_dir", "log_dir", "cache_backup_dir")
  for (k in file_settings) {
    v <- s[[k]]
    if (is.null(v) || !nzchar(v)) next
    if (dir.exists(v))
      stop(sprintf(paste0("%s must name a FILE, not a directory -- got \"%s\".\n",
                          "  Append a filename, e.g. \"%s\"."),
                   k, v, file.path(v, basename(.default_basename(k)))), call. = FALSE)
    # basename() strips a trailing separator, so test the raw string.
    if (grepl("[/\\\\]$", v))
      stop(sprintf(paste0("%s must name a FILE, but \"%s\" ends in a path separator.\n",
                          "  Append a filename, e.g. \"%s%s\"."),
                   k, v, v, .default_basename(k)), call. = FALSE)
  }
  for (k in dir_settings) {
    v <- s[[k]]
    if (is.null(v) || !nzchar(v)) next
    if (file.exists(v) && !dir.exists(v))
      stop(sprintf("%s must name a DIRECTORY, not a file -- got \"%s\".", k, v), call. = FALSE)
  }
  invisible(TRUE)
}

# A SLURM head node: the scheduler's tools are on PATH but we hold no allocation. Nothing here
# belongs on one -- a run needs node-local disk for WAL, and even the read-only diagnostics
# load a multi-GB store and fit models.
#
# LIMITATION: inside the apptainer image `squeue` need not be on PATH, so this can return FALSE
# on a login node. container/run_in_container.sh makes the same test on the HOST, where squeue
# always is, and that is the check that covers the documented cluster path.
.on_login_node <- function() {
  nzchar(Sys.which("squeue")) && !nzchar(Sys.getenv("SLURM_JOB_ID"))
}

# The message every caller gives for it, so the recipe is written once.
.login_node_message <- function(what) {
  paste0(what, " must not run on a login node.\n",
         "  Get an allocation first:\n",
         "    salloc -N1 -n4 --mem=32G -t 1:00:00 -A <account>\n",
         "    module load apptainer\n",
         "    ./run_in_container.sh exec <script.R>\n",
         "  (find your account with: sacctmgr -Pns show user format=account,defaultaccount)")
}

# Store and cache paths for a run inside a SLURM job, for settings.local.R to splice in:
#
#   settings_override <- c(cluster_scratch_paths(), list(dosage_budget_bytes = 16e9))
#
# Must live in settings.R, not R/: .local_overrides() runs before R/ is sourced. Kept out of
# settings.local.R itself -- LESSONS #26.
cluster_scratch_paths <- function(subdir = paste0("t3opt_", Sys.info()[["user"]])) {
  tmp <- Sys.getenv("TMPDIR")
  # SLURM_JOB_ID, not the path, is what distinguishes a compute node's local disk from a login
  # node's shared storage.
  if (!nzchar(tmp))
    stop("TMPDIR is unset -- cannot place the store on node-local disk.")
  # Deliberately stricter than .on_login_node(): that one also requires squeue on PATH, and a
  # missing squeue must not be what decides a WAL store may go on shared storage.
  if (!nzchar(Sys.getenv("SLURM_JOB_ID")))
    stop("SLURM_JOB_ID is unset, so this is not running inside a job.\n",
         "  On a login node $TMPDIR is shared storage, and SQLite's WAL cannot work there.\n",
         "  Get an allocation first:\n",
         "    salloc -N1 -n4 --mem=16G -t 1:00:00 -A <account>\n",
         "  (find your account with: sacctmgr -Pns show user format=account,defaultaccount)")
  if (!nzchar(Sys.getenv("OPTIMIZER_HOME")))
    stop("OPTIMIZER_HOME is unset -- set it in .Renviron to your /project path.")

  # Per-user subdirectory: $TMPDIR may be shared with others on the node.
  base <- file.path(tmp, subdir)
  list(
    db_path   = file.path(base, "evals.sqlite"),
    cache_dir = file.path(base, "cache")
  )
}

# The filename each file setting defaults to, for the error messages above.
.default_basename <- function(k)
  switch(k, db_path = "evals.sqlite", db_backup_path = "evals_backup.sqlite",
         report_path = "report.md", stop_file = "STOP", cache_ready_file = ".cache_ready", "file")

# Layer overrides on the defaults; warn on an unknown key, usually a typo.
.apply_overrides <- function(defaults, ov) {
  if (!length(ov)) return(defaults)
  unknown <- setdiff(names(ov), names(defaults))
  if (length(unknown))
    warning("settings.local.R sets unknown setting(s) [typo?]: ", paste(unknown, collapse = ", "))
  modifyList(defaults, ov)
}

# local_overrides = FALSE returns the tracked defaults only; tests use it.
# `extra` (a list, e.g. from read_domain_settings()) is applied last, each key REPLACING the whole
# value: modifyList would drop NULL fields inside target_domain and change the run_id hash.
optimizer_settings <- function(local_overrides = TRUE, extra = NULL) {
  if (remote_server && !nzchar(Sys.getenv("OPTIMIZER_HOME")))
    stop("remote mode is on but OPTIMIZER_HOME is not set in the environment. Add it to ",
         ".Renviron and RESTART R (.Renviron is read only at startup), or set ",
         "OPTIMIZER_REMOTE=false. Check with Sys.getenv(\"OPTIMIZER_HOME\").")
  perm_dir <- if (remote_server) Sys.getenv("OPTIMIZER_HOME") else here::here()

  defaults <- list(
    build = OPTIMIZER_BUILD,

    # ---- mode -------------------------------------------------------------
    simulate = TRUE,             # TRUE = synthetic offline world; FALSE = real T3 pipeline

    # ---- simulate-mode noise (offline testing only) -----------------------
    sim_noise_sd    = 0.07,      # sd of the observation noise .sim_evaluate() adds
    sim_fixed_trial = FALSE,     # TRUE removes trial-to-trial heterogeneity

    # ---- CV schemes -------------------------------------------------------
    schemes = c("CV0", "CV00"),  # what the diagnostics sanity-check
    optimize_scheme = "CV00",     # the ONE scheme this run targets; must be in `schemes`

    # How the surrogate handles the trial effect. "merf" puts it outside the tree as a shrunken
    # random effect, "blocked" makes trial_id a forest feature marginalised at prediction,
    # "pooled" ignores trials and fits per-config means. Each falls back to the next when it
    # cannot be fitted. dev/surrogate_bakeoff.R scores all three on the current store.
    surrogate_method = "merf",

    # ---- search behaviour -------------------------------------------------
    # n_random_init is a BUDGET: draw configs at random until this many have a finite
    # mean_score, then let acquisition take over. It is a one-time gate on a cumulative count,
    # so raising it on a run that is already past it does nothing. `reference_every` below is
    # the persistent rate, and is the knob that keeps the reference population growing. This key
    # is run-defining only under run_identity = "experiment"; under "target" it is policy and
    # can be changed on a live run. See SEARCH_ALGORITHM.md §13.
    n_random_init       = 25,    # random configs before the surrogate takes over
    # surrogate_min_obs is a STATISTICAL FLOOR: rows a forest needs before it is worth fitting
    # at all. Until 2026-11-03 this was n_random_init itself, which meant "explore longer" and
    # "demand more data per tree" could not be set apart. Not run-defining.
    #
    # It counts whatever rows the chosen mode fits on, and the modes differ: merf and blocked
    # get one row per (config, trial), pooled gets one row per CONFIG. So a value above
    # n_random_init leaves pooled permanently unfittable -- phase 3 guarantees only
    # n_random_init scored configs before acquisition begins -- and the search falls through to
    # random_fallback instead. That is safe but wasteful, so keep this <= n_random_init unless
    # you intend to retire the pooled fallback. Defaulting both to 25 preserves the behaviour
    # the single coupled setting gave.
    surrogate_min_obs   = 25,    # rows required before a surrogate is fitted
    incumbent_min_reps  = 2,     # trials a config needs before it can be "incumbent"
    ei_xi               = 0.01,  # Expected-Improvement exploration constant
    # The parent pool for crossover/mutation is the CONTENDER set (get_parents()), so these two
    # are the cap and the floor on it, not the pool itself:
    #   n_parents  cap. The in-play set is every scored configuration when the estimator cannot
    #              separate them, and unbounded parenthood is then recombination among all of
    #              them. Ordered by UCB, so the cap keeps the most promising.
    #   n_elites   floor, filled from the mean ranking. Covers the case where the random-effects
    #              fit did not run at all: no `se`, so no contenders, so without this the
    #              candidate pool would be pure random draws and nothing would say so.
    n_elites            = 8,     # FLOOR on the parent pool (mean-ranked, needs no se)
    n_parents           = 32,    # CAP on the parent pool (UCB-ranked contenders)
    n_cross             = 60, n_mut = 60, n_rand = 60,  # candidate pool sizes
    ntree               = 300,   # trees in the random-forest surrogate

    # ---- budget -----------------------------------------------------------
    max_iters        = Inf,      # stop after this many evaluations
    checkpoint_every = 1,        # write a report snapshot every N iterations

    # ---- real-mode trial sampling (ignored when simulate = TRUE) ----------
    crop_name         = "Wheat",
    focal_trait       = "Grain yield - kg/ha|CO_321:0001218",  # as it appears in T3
    # observationVariableDbId for focal_trait; NULL samples the whole crop catalogue.
    focal_trait_db_id = "84527",

    # ---- genotype source combining (subtask C) ----------------------------
    # Fraction of the smaller marker panel two projects must share to count as the SAME
    # protocol and be combined before marker QC -- LESSONS #12.
    merge_containment     = 0.95,
    # Fraction of the smaller accession set that makes two projects in a group REDUNDANT.
    redundant_acc_overlap = 0.90,

    # em_combine degrees of freedom: a partial covariance's df is its relative WEIGHT, so only
    # ratios matter. Effective sample size sets the ordering, then values are re-centred on the
    # mean with the spread capped by the stdev. Keep the mean well above the stdev.
    em_df_mean  = 60,
    em_df_stdev = 15,

    # Memory budget for a SINGLE project's dense dosage (samples x markers x 4 B), and the sole
    # control on cached marker density: each project is parsed once, at the densest thinning
    # that fits. Changing it changes scores.
    dosage_budget_bytes = 2e9,
    # Cap on the COMBINED size of all projects loaded for one trial -- what actually bounds a
    # worker's peak. Over it, the largest are served at a coarser thin by column-subsetting the
    # cache. Inf disables. Sizing: RUNBOOK_INTERACTIVE.md §2.
    dosage_total_budget_bytes = 12e9,
    # Re-parse a cached dosage denser when this machine's budget affords it, at the cost of a
    # one-time re-download. FALSE keeps whatever thin is cached.
    dosage_redensify = TRUE,

    brapi_tries               = 4,   # attempts per BrAPI call; 1 disables retry
    # The startup connect + login, which every worker must clear before it can do anything.
    # Longer than brapi_tries because the trade-off is different: a worker that cannot connect
    # has nothing else to do, while a call inside the loop is holding a claim. 8 attempts at
    # 2 * 2^(n-1) s backoff rides out roughly 5 minutes of network trouble -- run_workers.sh
    # staggers starts by 20 s, so a shorter blip than that takes every worker at once.
    brapi_connect_tries       = 8,
    vcf_max_download_attempts = 3,   # per-session failed VCF downloads before skipping -- LESSONS #10

    min_trial_acc      = 30,     # skip trials with fewer genotyped accessions
    min_train_trials   = 3,      # skip focal trials with too few training TRIALS
    min_train_acc      = 20,     # genotyped overlap a GBLUP needs to fit, in ACCESSIONS
    # The single floor on the TEST overlap, in ACCESSIONS. It governs both halves of the
    # question, which is why there is one setting and not two:
    #   FEASIBILITY  below it, run_pipeline() and score_predictions() decline the task, so no
    #                evaluation is paid for and no row is stored.
    #   USABILITY    aggregate_scores() keeps a score only at or above it, because a Fisher-z
    #                weight of n_test - 3 cannot rescue a correlation this short.
    # Set them apart and the difference is pure waste: an evaluation that costs a full pipeline
    # run, succeeds, is stored, and is then invisible to every estimate -- 36% of one run's
    # successful evaluations, before this was unified (docs/LESSONS.md #34). The Predictathon's
    # own "task can fail" rule is 5, so 10 declines some tasks it would score; a correlation on
    # 5 accessions is not worth an evaluation either way.
    min_test_acc       = 10,
    # The CEILING on how much one evaluation can count, as a WEIGHT. Two independent errors
    # corrupt a score, and they behave differently as the trial grows:
    #   SAMPLING ERROR   variance 1/(n_test - 3) on the Fisher-z scale. Shrinks with the trial.
    #   CONFIG x TRIAL   the configuration really performing differently on THIS trial than on
    #   INTERACTION      an average one; variance 1/weight_cap. Does NOT shrink -- a bigger
    #                    trial measures that trial's idiosyncrasy more precisely, it does not
    #                    make the idiosyncrasy smaller.
    # So the weight is the two precisions in series, 1/w = 1/(n_test - 3) + 1/weight_cap, and
    # this number is the weight a trial of unbounded size would get. It therefore reads as a
    # trial size: 50 means "no trial counts for more than one with 53 usable accessions". Inf
    # restores the uncapped n_test - 3.
    # Without it the weight spans 7 to 1434 across a real store (205x) and the largest trials
    # decide everything; at 50 the span is 8x. Governs BOTH the aggregate estimator and the
    # surrogate, from one helper, because they estimate the same quantity from the same rows
    # (docs/LESSONS.md #36). dev/calibrate_weights.R re-derives it.
    weight_cap         = 50,
    # WHICH POPULATION sd_config is estimated on. It is the shrinkage prior and it sets the unit
    # the contender test is denominated in, but as a fitted quantity it is the spread of the
    # configurations the SEARCH CHOSE to evaluate -- which narrows as the search concentrates,
    # so the prior tightens because the search succeeded, and every BLUP collapses toward the
    # grand mean (SEARCH_ALGORITHM.md §11).
    #   "reference"  estimate it on the configurations drawn WITHOUT the surrogate -- the Phase
    #                3 random draws, read off evals.source -- and re-shrink on that. Seeds are
    #                excluded: five hand-picked submissions are not a sample of the space.
    #   "fitted"     the whole-population REML estimate, i.e. the behaviour before 2026-11-03.
    # Falls back to "fitted" silently whenever the reference estimate cannot be had, which the
    # report states; `source` has only been stored since 2026-11-03, so for older rows the
    # reference set is inferred from arrival order (.reference_hashes()). Not run-defining: it
    # changes an estimate derived from rows, not what any row means.
    sd_config_source      = "reference",
    reference_min_configs = 20,  # reference configs required before that estimate is trusted
    # The reference sd_config becomes a FIXED prior with no allowance for its own uncertainty,
    # and it is fitted on ~20-odd configurations with one or two rows each, so it can be large,
    # confident and meaningless. Each draw refits on two thirds of the reference
    # configurations; the estimate is rejected (and the fitted value kept) when too many draws
    # collapse to the boundary. Measured: this is what separates a trustworthy 0.0246 on
    # 9a2664962c01 from a 0.0517 on e7517df121d5 whose subsample IQR is 0.0000-0.0924. A few
    # extra fits on a few dozen rows, against an evaluation that takes tens of minutes. 0
    # disables the gate. The draws are seeded from the reference set itself, so the verdict is a
    # function of the evidence and does not flip between iterations on unchanged data.
    reference_stability_draws = 25,
    # A PERSISTENT exploration rate, which n_random_init is not: that is a one-time gate on a
    # cumulative count, so once it is passed the reference population never grows again and its
    # sd_config estimate stops improving. 1 step in this many draws a fresh random
    # configuration no matter how mature the run is, logged as source = "random_reference" so it
    # is distinguishable in the log while counting as a reference draw. 0 disables, restoring
    # the frozen reference population.
    reference_every       = 20,
    # Evaluations a REFERENCE configuration earns, against config_replication for everything
    # else. A variance component needs replication WITHIN level: at config_replication = 2 the
    # reference fit rests on about two rows per configuration, which was thin enough on one real
    # run (7d75f214584a) for it to return a boundary 0 and decline. Served through the RATIONED
    # tier of .replication_backlog(), never the unrationed floor, so a large value here cannot
    # starve the search. <= config_replication disables it.
    reference_replication = 4,
    # The two halves of LESSONS #19, both applying in simulate mode as well. 1 disables either.
    # Both are the BASE of a schedule that rises as the store grows -- see .trial_target() and
    # the replication backlog in R/optimizer.R.
    trial_replication  = 2,      # distinct configs a trial reaches before the search moves on
    config_replication = 2,      # evaluations every config gets, spread over distinct trials
    # WHAT COUNTS AS STILL IN PLAY. One definition, used for two things: which configurations
    # earn extra evaluations, and which are parents for candidate generation (get_parents()).
    #
    # contender_rule names what the optimistic bound mean_i + z*se_i is compared against:
    #   "racing"  max_j(mean_j - z*se_j), the best DEFENSIBLE value in the set. Optimistic
    #             against pessimistic, so both sides are charged the same z, and a candidate is
    #             eliminated only once its interval stops overlapping the leader's.
    #   "point"   max_j(mean_j). The bar is then the maximum of a set of noisy estimates, which
    #             one lucky thinly replicated configuration can set for everyone; measured, the
    #             in-play set moves 2.4x to 4.6x further per checkpoint than under racing at the
    #             same z (dev/sweep_contender_z.R).
    #
    # z does DOUBLE DUTY under racing -- it inflates the claim and lowers the bar -- so the two
    # rules give comparable set sizes at z_racing ~ z_point/2. 0.25 here is the racing
    # equivalent of 0.5 under "point", which is where the sweep found the most checkpoints with
    # a usable set size. Under run_identity = "experiment" contender_z is run-defining, so
    # changing it starts a new run_id; under the default "target" none of these three is, which
    # is the reading SEARCH_ALGORITHM.md §13 argues for.
    contender_rule     = "racing",
    contender_z        = 0.25,
    # The RESOLUTION FLOOR: the difference in predictive accuracy worth resolving, on the
    # correlation scale, chosen from outside the run instead of measured in it. A configuration
    # keeps its place when it misses the bar by less than this, so membership flips only on
    # evidence that moves further than the resolution anyone cares about. It is the only term in
    # the contender test that does NOT shrink as the evaluated population homogenizes, which is
    # the defect SEARCH_ALGORITHM.md §11 describes. Read it against sd_config, which the
    # reference population puts at 0.013 to 0.056 across this project's runs: as delta
    # approaches sd_config everything is in play, and that is the honest reading rather than a
    # failure -- no configuration then differs from another by an amount worth having. 0
    # switches it off.
    #
    # 0.002 is measured, not chosen (dev/sweep_contender_rule.R, racing z 0.25, reference
    # sd_config, pooled over the runs with a non-zero fitted sd_config):
    #
    #   delta   in_play (median)   % of scored   checkpoint-to-checkpoint jump
    #   0.000         7                2.4            0.66
    #   0.002        16                6.4            0.38
    #   0.005        33                12.8           0.31
    #   0.010        88                24.7           0.27
    #   0.020       117                63.7           0.24
    #
    # 0.002 cuts the thrashing by 42% against no floor while keeping the set comfortably above
    # the cap of 8 and inside a tenth of the population. Larger values keep buying stability --
    # a set of 117 can hardly move -- but by then the rule has stopped discriminating. Note that
    # SEARCH_ALGORITHM.md §11 suggested 0.02 before this was measured; at 0.02 two thirds of the
    # field is in play.
    resolution_delta   = 0.002,
    replicate_every    = 3,      # rations the CONTENDER tier only: 1 worker in this many takes
                                 # from it. The config_replication floor is never rationed.
    max_sample_fail    = 25,     # consecutive trial-sampling failures before halting
    max_consec_skip    = 50,     # consecutive iterations storing nothing before halting: the
                                 # liveness bound, not a tuning knob

    # Known-feasible trials used as a bug oracle by check_canaries(): failure here means the
    # code is hiding real data. NULL disables. The nine Predictathon focal trials, named
    # studyName -> studyDbId so diagnostics can find the submission folders.
    canary_trials    = c(`2025_AYT_Aurora`            = "10673",
                         `24Crk_AY2-3`                = "10674",
                         `25_Big6_SVREC_SVREC`        = "10675",
                         `CornellMaster_2025_McGowan` = "10676",
                         `YT_Urb_25`                  = "10677",
                         `AWY1_DVPWA_2024`            = "10678",
                         `OHRWW_2025_SPO`             = "10679",
                         `TCAP_2025_MANKS`            = "10680",
                         `STP1_2025_MCG`              = "10681"),
    canary_weak_trials = c("10674", "10678", "10681"),  # marginal: soft-warn, not ALARM
    oracle_trials      = c("10675", "10677"),           # data-rich, for sweep_rich_trials()
    brapi_host         = "wheat.triticeaetoolbox.org",

    # WHAT A RUN IS. The run signature decides when new evaluations join the existing evidence
    # and when they start a separate history (R/store.R .run_signature_keys):
    #   "target"      a run is a TARGET. Identified by optimize_scheme, target_domain,
    #                 focal_trait_db_id, min_trial_acc and simulate -- the keys that change what
    #                 a stored score MEANS. Search policy (trial_replication,
    #                 config_replication, replicate_every, contender_z, n_random_init) can then
    #                 be retuned without starting a new run, so evidence keeps accumulating on
    #                 one domain while the method improves.
    #   "experiment"  a run is an EXPERIMENT. All ten keys identify it, so every policy change
    #                 starts a fresh run. Use this to compare search policies against each other.
    # Switching between them rotates run_id ONCE, because the two hashes differ; after that the
    # chosen reading is stable. It does NOT change which rows the search reads -- that is the
    # universe, scheme and build filters (R/optimizer.R), none of which consult run_id -- so no
    # store migration is involved either way. SEARCH_ALGORITHM.md §13.
    run_identity     = "target",

    # ---- target domain: which trials the pipeline is optimized FOR --------
    # Focal trials are drawn from here. NULL places no constraint; fields combine with AND.
    # Training trials are unconstrained.
    target_domain = list(
      programs  = NULL,          # programName values, or NULL
      years     = NULL,          # years, or NULL
      locations = NULL,          # locationName values, or NULL
      trials    = c("25_MSU_AYT_MI_ALLMI", "25_MSU_AYT_MI_CENMI",
                    "25_MSU_AYT_MI_HURMI", "25_MSU_AYT_MI_MONMI",
                    "25_MSU_AYT_MI_SANMI", "25_MSU_PYT_MASMI", "25_MSU_PYT_SVREC",
                    "Big6_2023_FRA", "Big6_2023_FRE",
                    "Big6_2023_LAF", "Big6_2023_MAS",
                    "Big6_2023_PRN", "Big6_2023_SVREC", "Big6_2023_URB",
                    "Big6_2023_WOO", "Big6_2024_SVREC",
                    "Big6_2025_MAS", "Big6_2026_NEO",
                    "Big6_2026_URB", "YT_Addie_23",
                    "YT_Neo_23", "YT_Stp_23",
                    "YT_Urb_23", "YT_Urb_25")
    ),

    # ---- paths ------------------------------------------------------------
    db_path     = file.path(perm_dir, "state", "evals.sqlite"),
    stop_file   = file.path(perm_dir, "state", "STOP"),
    # Written by the leader once the cache is restored; other workers wait for it.
    cache_ready_file = file.path(perm_dir, "state", ".cache_ready"),
    report_path = file.path(perm_dir, "state", "report.md"),
    log_dir     = file.path(perm_dir, "logs"),
    # The regenerable cache stays on the work disk even in remote mode.
    cache_dir   = here::here("cache"),

    # ---- backups ----------------------------------------------------------
    # Additive rsync of the cache to durable storage; NULL disables. Restored to cache_dir at
    # startup when the work cache is empty.
    cache_backup_dir   = if (remote_server) file.path(perm_dir, "cache") else NULL,
    cache_sync_minutes = 30,     # 0 disables -- LESSONS #25
    # Where the live store is copied, by VACUUM INTO. NULL disables.
    db_backup_path     = if (remote_server) file.path(perm_dir,
                                                      "state", "evals_backup.sqlite") else NULL,

    # ---- several workers on one store -------------------------------------
    # N workers share one db_path and one cache_dir. db_path MUST BE ON LOCAL DISK: WAL
    # coordinates through an mmap'd -shm file that a network filesystem does not provide, and
    # on NFS the pragma fails silently. The leader (worker 1) restores the cache at startup and
    # signals cache_ready_file -- LESSONS #25.
    worker_id   = { w <- trimws(Sys.getenv("OPTIMIZER_WORKER")); if (nzchar(w)) w else "1" },
    is_leader   = !nzchar(trimws(Sys.getenv("OPTIMIZER_WORKER"))) ||
                  identical(trimws(Sys.getenv("OPTIMIZER_WORKER")), "1"),

    # ---- concurrent cache access ------------------------------------------
    # get_project_dosage holds a per-project lock across download -> parse -> cache-write; the
    # loser waits for the winner's result.
    lock_wait_minutes  = 60,
    lock_stale_minutes = 90      # break a lock older than this (its holder died)
  )
  s <- .apply_overrides(defaults, if (isTRUE(local_overrides)) .local_overrides() else list())
  if (length(extra)) s[names(extra)] <- extra
  # Validate the EFFECTIVE value, so a settings.local.R override is checked too.
  if (length(s$surrogate_method) != 1L ||
      !(s$surrogate_method %in% c("merf", "blocked", "pooled")))
    stop("surrogate_method must be one of merf, blocked, pooled -- got: ",
         paste(s$surrogate_method, collapse = ", "))
  if (length(s$surrogate_min_obs) != 1L ||
      !is.finite(suppressWarnings(as.integer(s$surrogate_min_obs))) ||
      as.integer(s$surrogate_min_obs) < 1L)
    stop("surrogate_min_obs must be a single integer >= 1 -- got: ",
         paste(s$surrogate_min_obs, collapse = ", "))
  if (length(s$run_identity) != 1L || !(s$run_identity %in% c("target", "experiment")))
    stop("run_identity must be one of target, experiment -- got: ",
         paste(s$run_identity, collapse = ", "))
  for (nm in c("reference_every", "reference_replication", "reference_stability_draws")) {
    v <- suppressWarnings(as.integer(s[[nm]]))
    if (length(v) != 1L || is.na(v) || v < 0L)
      stop(nm, " must be a single integer >= 0 -- got: ", paste(s[[nm]], collapse = ", "))
  }
  if (length(s$contender_rule) != 1L || !(s$contender_rule %in% c("racing", "point")))
    stop("contender_rule must be one of racing, point -- got: ",
         paste(s$contender_rule, collapse = ", "))
  if (length(s$sd_config_source) != 1L ||
      !(s$sd_config_source %in% c("reference", "fitted")))
    stop("sd_config_source must be one of reference, fitted -- got: ",
         paste(s$sd_config_source, collapse = ", "))
  for (nm in c("contender_z", "resolution_delta")) {
    v <- suppressWarnings(as.numeric(s[[nm]]))
    if (length(v) != 1L || !is.finite(v) || v < 0)
      stop(nm, " must be a single finite number >= 0 -- got: ",
           paste(s[[nm]], collapse = ", "))
  }
  for (nm in c("n_parents", "n_elites", "reference_min_configs")) {
    v <- suppressWarnings(as.integer(s[[nm]]))
    if (length(v) != 1L || is.na(v) || v < 1L)
      stop(nm, " must be a single integer >= 1 -- got: ", paste(s[[nm]], collapse = ", "))
  }
  if (length(s$optimize_scheme) != 1L || !is.character(s$optimize_scheme) ||
      !(s$optimize_scheme %in% s$schemes))
    stop("optimize_scheme must be a single scheme that is one of `schemes` (",
         paste(s$schemes, collapse = ", "), "); got: ",
         paste(s$optimize_scheme, collapse = ", "))
  .validate_paths(s)
  s
}
