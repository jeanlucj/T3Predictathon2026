# settings.R -- every knob in one place.
#
# Machine-specific values go in settings.local.R (gitignored), never here.

# Stamped into the log, the report and every stored eval. Last digit for a minor change, middle
# for a major one. If a bump changes what a configuration COMPUTES, add a BUILD_CHANGES entry
# (R/optimizer.R) so the invalidated rows can be retired.
OPTIMIZER_BUILD <- "0.8.6"

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
           error = function(err) stop("settings.local.R failed to source: ", conditionMessage(err)))
  ov <- get0("settings_override", envir = e, ifnotfound = list())
  if (!is.list(ov)) stop("settings.local.R must define `settings_override` as a list()")
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
optimizer_settings <- function(local_overrides = TRUE) {
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

    # ---- search behaviour -------------------------------------------------
    n_random_init       = 25,    # random configs before the surrogate takes over
    incumbent_min_reps  = 2,     # trials a config needs before it can be "incumbent"
    ei_xi               = 0.01,  # Expected-Improvement exploration constant
    n_elites            = 8,     # elites that seed crossover/mutation
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
    vcf_max_download_attempts = 3,   # per-session failed VCF downloads before skipping -- LESSONS #10

    min_trial_acc      = 30,     # skip trials with fewer genotyped accessions
    min_train_trials   = 3,      # skip focal trials with too few training TRIALS
    min_train_acc      = 20,     # genotyped overlap a GBLUP needs to fit, in ACCESSIONS
    min_test_acc       = 5,      # the Predictathon's own "task can fail" rule, in ACCESSIONS
    # The two halves of LESSONS #19, both applying in simulate mode as well. 1 disables either.
    # Both are the BASE of a schedule that rises as the store grows -- see .trial_target() and
    # the replication backlog in R/optimizer.R.
    trial_replication  = 2,      # distinct configs a trial reaches before the search moves on
    config_replication = 2,      # evaluations every config gets, spread over distinct trials
    contender_z        = 1,      # z in blup + z*se; a config stops earning reps once its
                                 # optimistic bound falls below the leader's estimate
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

    # ---- target domain: which trials the pipeline is optimized FOR --------
    # Focal trials are drawn from here. NULL places no constraint; fields combine with AND.
    # Training trials are unconstrained.
    target_domain = list(
      programs  = NULL,          # programName values, or NULL
      years     = NULL,          # years, or NULL
      locations = NULL,          # locationName values, or NULL
      trials    = c("23_Big6_MAS", "23_Big6_SVREC", "24_Big6_FHB_FHB",
                    "24_Big6_MASMI", "24_Big6_SVREC", "24_Jhon_FHB_CRD_FHB",
                    "25_Big6_Core_FHB_FHB", "25_Big6_MASON_MASMI",
                    "25_MSU_AYT_MI_ALLMI", "25_MSU_AYT_MI_CENMI", "25_MSU_AYT_MI_FHB",
                    "25_MSU_AYT_MI_HURMI", "25_MSU_AYT_MI_MONMI", "25_MSU_AYT_MI_SANMI",
                    "25_MSU_PYT_FHB", "25_MSU_PYT_MASMI", "25_MSU_PYT_SVREC",
                    "Big6_Fra_23", "Big6_Fre_23", "Big6_Laf_23", "Big6_Mas_23",
                    "Big6_Prn_23", "Big6_Scb_23", "Big6_Urb_23", "Big6_Woo_23",
                    "IL_Scb_23", "SI_Urb_22", "YT_Addie_23", "YT_Neo_23", "YT_Stp_23",
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
  # Validate the EFFECTIVE value, so a settings.local.R override is checked too.
  if (length(s$optimize_scheme) != 1L || !is.character(s$optimize_scheme) ||
      !(s$optimize_scheme %in% s$schemes))
    stop("optimize_scheme must be a single scheme that is one of `schemes` (",
         paste(s$schemes, collapse = ", "), "); got: ",
         paste(s$optimize_scheme, collapse = ", "))
  .validate_paths(s)
  s
}
