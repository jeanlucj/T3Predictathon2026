# settings.R
#
# Every knob in one place. The driver and the optimizer read from here, so you
# tune the run by editing this file (or overriding fields after sourcing it).

# remote_server: when TRUE, permanent files (state / logs / cache backup) go to
# durable HOME storage instead of the work directory (see the paths section and
# README "Running on a remote server"). It is DRIVEN BY THE ENVIRONMENT -- NOT
# hard-coded -- so the same tracked settings.R works unedited on your laptop and
# on the server, and `git pull` never conflicts on this file:
#   * OPTIMIZER_REMOTE, if set to a truthy value (1/true/yes/on) or falsy (0/false/no/off),
#   decides explicitly; otherwise
#   * it auto-detects: TRUE when OPTIMIZER_PATH is set in the environment, else FALSE.
# The server's .Renviron sets OPTIMIZER_PATH, so the server is remote
# automatically and a laptop (no OPTIMIZER_PATH) is local -- with no edit to
# this file on either machine.
.detect_remote_server <- function() {
  ovr <- tolower(trimws(Sys.getenv("OPTIMIZER_REMOTE")))
  if (nzchar(ovr)) return(ovr %in% c("1", "true", "yes", "on"))   # explicit override wins
  nzchar(Sys.getenv("OPTIMIZER_PATH"))                            # else auto-detect
}
remote_server <- .detect_remote_server()

# Machine-specific setting overrides live in an UNTRACKED settings.local.R (gitignored; see
# settings.local.R.example) so the tracked settings.R stays identical on every machine and
# `git pull` never conflicts on it. That file defines `settings_override <- list(...)`; those
# values are layered on top of the defaults below. Use it for anything you tune per machine --
# simulate, dosage_budget_bytes, max_hours, run_startup_canary, custom paths, ... (remote_server
# itself is env-driven, above -- set OPTIMIZER_PATH / OPTIMIZER_REMOTE for the local/remote
# split rather than overriding it here).
.local_overrides <- function(path = here::here("settings.local.R")) {
  if (!file.exists(path)) return(list())
  e <- new.env(parent = globalenv())
  tryCatch(sys.source(path, envir = e),
           error = function(err) stop("settings.local.R failed to source: ", conditionMessage(err)))
  ov <- get0("settings_override", envir = e, ifnotfound = list())
  if (!is.list(ov)) stop("settings.local.R must define `settings_override` as a list()")
  ov
}

# Layer overrides on top of the defaults; warn on an unknown key (usually a typo, since it
# would otherwise be silently added and never read).
.apply_overrides <- function(defaults, ov) {
  if (!length(ov)) return(defaults)
  unknown <- setdiff(names(ov), names(defaults))
  if (length(unknown))
    warning("settings.local.R sets unknown setting(s) [typo?]: ", paste(unknown, collapse = ", "))
  modifyList(defaults, ov)
}

optimizer_settings <- function() {
  # Permanent, resumable state (sqlite store, report, logs, STOP flag) goes under `perm_dir`:
  # durable HOME storage in remote mode, else the work dir. Fail loudly if remote mode is on but
  # OPTIMIZER_PATH is unset (only possible via an explicit OPTIMIZER_REMOTE=true) -- otherwise the
  # paths would silently resolve to the filesystem root ("/state/...").
  if (remote_server && !nzchar(Sys.getenv("OPTIMIZER_PATH")))
    stop("remote mode is on but OPTIMIZER_PATH is not set in the environment. Add it to ",
         ".Renviron and RESTART R (.Renviron is read only at startup), or set ",
         "OPTIMIZER_REMOTE=false. Check with Sys.getenv(\"OPTIMIZER_PATH\").")
  perm_dir <- if (remote_server) Sys.getenv("OPTIMIZER_PATH") else here::here()

  defaults <- list(
    # ---- mode -------------------------------------------------------------
    # TRUE  = synthetic offline world (fast; for verifying the machinery).
    # FALSE = real T3 pipeline on real data (network + heavy model fits).
    simulate = TRUE,

    # ---- simulate-mode variance (offline testing only) --------------------
    # .sim_true() -- the synthetic objective -- is DETERMINISTIC in (config, trial, scheme).
    # Two separate things make an *evaluation* of it noisy, and they are independently
    # controllable so that "does the search work at all?" can be tested apart from "how many
    # evaluations does a given noise level need?":
    #   sim_noise_sd     sd of the observation noise .sim_evaluate() adds (0 = none).
    #   sim_fixed_trial  TRUE returns the SAME trial descriptor every time, removing
    #                    trial-to-trial heterogeneity. This is the LARGER of the two sources
    #                    and the easier one to overlook: .sim_true ends with q*(0.5 + h) for
    #                    h ~ U(0.2, 0.7), which alone swings a score by about +-25%, and
    #                    env_var / n_projects additionally change WHICH methods are best.
    # The defaults reproduce the realistic, noisy regime; tests/test_sim_loop.R turns both
    # down to isolate the algorithm.
    sim_noise_sd    = 0.07,
    sim_fixed_trial = FALSE,

    # ---- CV schemes -------------------------------------------------------
    # `schemes` lists every CV scheme the DIAGNOSTICS sanity-check for code bugs
    # (check_canaries / sweep_rich_trials / diagnose_trial test all of them, since a
    # method can behave differently under masking).
    schemes = c("CV0", "CV00"),
    # The OPTIMIZER targets exactly ONE scheme per run: CV0 and CV00 are distinct
    # prediction tasks with (potentially) different optimal pipelines, so blending
    # them would only find a compromise. To optimize both, run the optimizer twice
    # -- once per scheme -- against the same store (evals.sqlite is a shared archive
    # scoped per scheme on read) and cache (the second run reuses cached genotypes).
    # Must be a single string and a member of `schemes` (validated below).
    optimize_scheme = "CV0",

    # ---- search behaviour -------------------------------------------------
    n_random_init       = 25,    # random configs before the surrogate takes over
    incumbent_min_reps  = 2,     # trials a config needs before it can be "incumbent"
    reeval_prob         = 0.15,  # chance an iteration re-evaluates a random elite (denoise leaders)
    ei_xi               = 0.01,  # Expected-Improvement exploration constant
    n_elites            = 8,     # elites that seed crossover/mutation
    n_cross             = 60, n_mut = 60, n_rand = 60,  # candidate pool sizes
    ntree               = 300,   # trees in the random-forest surrogate

    # ---- budget / background control -------------------------------------
    max_iters     = Inf,         # stop after this many evaluations
    max_hours     = Inf,         # ... or this much wall-clock, whichever first
    # Wall-clock cap on ONE evaluation, in minutes. Inf = no cap, and that is the default.
    #
    # PITFALL: a cap CENSORS NON-RANDOMLY. The slow configurations are the thorough ones
    # (em_combine over many panels, all_projects), which are also the ones that may predict
    # best. And it costs more than one score: a timeout stores NA, aggregate_scores() averages
    # only finite scores, so a config that always times out has mean_score = NA and can never
    # become incumbent however good it is. Set a cap too low and the search will "discover" that
    # cheap pipelines win -- an artifact of the cap.
    #
    # If you set one, derive it from the distribution of SUCCESSFUL run times (report_timing.R),
    # not from the median over all outcomes: successes run far longer than the median. Over-run
    # is recorded as status = "timeout" -- its own status, so the censoring stays visible. R can
    # only interrupt at interpreter checkpoints, so compiled calls may overshoot the cap.
    max_eval_minutes = Inf,
    checkpoint_every = 1,        # write report snapshot every N iterations

    # ---- real-mode trial sampling (ignored when simulate = TRUE) ----------
    crop_name        = "Wheat",
    # The trait to optimize prediction of, as it appears in T3 ("name|CO_id").
    # Not limited to yield -- set this to any numeric trait to optimize that.
    focal_trait      = "Grain yield - kg/ha|CO_321:0001218",
    # T3 observationVariableDbId for the focal trait. Used to fetch ONLY trials
    # that measured it (server-side filter on the studies search), instead of
    # the whole crop catalogue. Find it once with:
    #   conn$wizard("trials", list(traits = "<dbid>"))   # confirms it narrows trials
    # Grain yield - kg/ha is 84527 on wheat.triticeaetoolbox.org. Set NULL to
    # disable the filter and sample from all trials (the old, slower behaviour).
    focal_trait_db_id = "84527",
    # ---- genotype source combining (subtask C) ----------------------------
    # Two genotyping projects are treated as the SAME protocol -- and so are
    # combined before marker QC -- when they share at least this fraction of the
    # smaller project's marker panel. Protocol *ids* cannot be used for this: a
    # "V2"/"v2.1" protocol is the same protocol scored against a different
    # reference genome, so it carries a different id but (near-)identical markers.
    # On the T3 projects seen so far the two classes separate cleanly: same-panel
    # pairs share 100% of the smaller panel, everything else <= 75%.
    merge_containment = 0.95,
    # Within a protocol group, two projects are REDUNDANT when they share at least
    # this fraction of the smaller accession set (same lines re-called). Keep the
    # project with more accessions; if the accession sets are identical, keep the
    # one with more markers.
    redundant_acc_overlap = 0.90,

    # Memory budget for a SINGLE PROJECT's dense dosage matrix (samples x markers x 4 B).
    # This is the SOLE control on cached marker density: each project is downloaded and
    # parsed ONCE, at the densest thinning that fits this budget (thin 1 = full markers for
    # anything that fits; a 7.5M-marker GBS panel is thinned to fit). Marker density is NOT a
    # search parameter -- this setting is the whole story (LESSONS.md #16) -- so set it
    # deliberately: a bigger budget caches denser (better GRMs, more disk/RAM per project),
    # smaller thins large panels harder.
    #
    # PITFALL: this is NOT a cap on the process. choose_geno_sources loads EVERY project
    # covering a trial at once, so the peak is the SUM over those projects -- on the T3 wheat
    # archive, roughly 6x this number, and more on a trial covered by many panels. Use
    # dosage_total_budget_bytes below to bound the sum; use this one to set density.
    # Changing it changes marker density and therefore SCORES, and since density is not a
    # config parameter the store cannot tell the difference on its own -- which is why
    # every row records the budget it was produced under (the dosage_budget column).
    dosage_budget_bytes = 2e9,

    # Cap on the COMBINED dense size of all projects loaded for one trial -- the setting that
    # actually bounds a worker's peak, and hence how many workers fit in RAM. When the
    # covering projects exceed it, choose_geno_sources serves the largest ones at a coarser
    # marker thin (by column-subsetting the cache, so it costs no re-download) until the sum
    # fits, and records that it did so in the evaluation's detail.
    #
    # Peak per worker runs to roughly 1.5-2x this figure (marker QC and the merge each need a
    # copy), so size it at about a third of the RAM each worker may have:
    #
    #     workers x dosage_total_budget_bytes x 2  <  usable RAM
    #
    #   RAM     workers   dosage_budget_bytes   dosage_total_budget_bytes
    #   256 GB     8              4e9                    10e9
    #   512 GB     8              8e9                    20e9
    #   512 GB     4             16e9                    40e9
    #
    # These are starting points from the cached panel sizes; replace them with the measured
    # peak_r_mb distribution once a day of evaluations is in the store (report_memory.R).
    # Inf disables the cap (the pre-existing behaviour: no bound on the sum at all).
    dosage_total_budget_bytes = 12e9,

    # Re-densify a cached dosage when THIS machine's budget affords a denser parse than the
    # cache was thinned to (e.g. a laptop-thinned cache warmed onto a big-memory server): on
    # access, if dosage_budget_bytes now allows a smaller thin than the cached file, re-download
    # and re-parse denser (a one-time cost, only for projects that were thinned). FALSE keeps
    # whatever thin is cached (no re-downloads when you move to a bigger machine).
    dosage_redensify = TRUE,

    # How many times to attempt a BrAPI network call before giving up. A flaky T3
    # server (intermittent HTTP 500 / timeout) is ridden out by retrying with backoff, so
    # a transient error does not surface as a spurious infeasible/error. 1 disables retry.
    brapi_tries      = 4,

    # Per-SESSION cap on failed VCF-download attempts for one project before it is skipped
    # for the rest of the run (T3 is likely down for that archive). Effort also decreases on
    # each retry: the in-call retry count drops (brapi_tries, then -1 per prior failure) and
    # the backoff shortens. Resets on a new run, so a transiently-unavailable project is
    # retried fresh. Prevents a timing-out project from being re-stormed on every trial.
    vcf_max_download_attempts = 3,

    min_trial_acc    = 30,       # skip trials with fewer genotyped accessions
    min_train_trials = 3,        # skip focal trials we cannot assemble a training set for
    max_sample_fail  = 25,       # consecutive trial-sampling failures before halting
                                 # (guards against a too-restrictive domain or a
                                 #  down network spinning the loop forever)
    # Run the full check_canaries() bug oracle at run_optimizer() startup (real
    # mode only). It confirms the pipeline can still predict known-feasible
    # trials BEFORE spending hours, guarding the silent-failure mode where a
    # data-hiding bug looks like ordinary infeasibility. But it is the heaviest,
    # most network-bound phase and runs before the resumable loop, so set FALSE
    # to skip it once you have validated the pipeline (or when a flaky server
    # makes the canary downloads themselves the risk). You can always run
    # check_canaries(s, conn) by hand.
    run_startup_canary = FALSE,
    # Known-feasible trials used as a BUG ORACLE: check_canaries() runs the most
    # permissive config on these and expects success. A canary that comes back
    # infeasible can only mean the code is hiding real data. NULL = no canary
    # check. These are the nine Predictathon focal trials (studyDbId on
    # wheat.triticeaetoolbox.org), resolved by studyName:
    #   2025_AYT_Aurora 10673, 24Crk_AY2-3 10674, 25_Big6_SVREC_SVREC 10675,
    #   CornellMaster_2025_McGowan 10676, YT_Urb_25 10677, AWY1_DVPWA_2024 10678,
    #   OHRWW_2025_SPO 10679, TCAP_2025_MANKS 10680, STP1_2025_MCG 10681.
    # Named studyName -> studyDbId so diagnostics can map to the participant
    # submission folders (named by studyName) for the anchor.
    canary_trials    = c(`2025_AYT_Aurora`            = "10673",
                         `24Crk_AY2-3`                = "10674",
                         `25_Big6_SVREC_SVREC`        = "10675",
                         `CornellMaster_2025_McGowan` = "10676",
                         `YT_Urb_25`                  = "10677",
                         `AWY1_DVPWA_2024`            = "10678",
                         `OHRWW_2025_SPO`             = "10679",
                         `TCAP_2025_MANKS`            = "10680",
                         `STP1_2025_MCG`              = "10681"),
    # The three trials expected to be marginal for several methods (soft-warn,
    # not hard CANARY ALARM, when they fail).
    canary_weak_trials = c("10674", "10678", "10681"),
    # Data-rich trials for sweep_rich_trials(): the config-sweep CODE oracle runs every
    # method/branch here, where feasibility is a property of the code, not the data, so a
    # branch that cannot predict on EITHER is a bug. Big6 (10675) + YT_Urb (10677).
    oracle_trials = c("10675", "10677"),
    brapi_host       = "wheat.triticeaetoolbox.org",

    # ---- target domain: which trials the pipeline is being optimized FOR --
    # The optimizer maximizes mean accuracy over random focal trials drawn from
    # this domain. Leave a field NULL to place no constraint on it; supply a
    # vector to restrict random sampling so the resulting pipeline is tailored
    # to that subpopulation (e.g. one breeding program's recent trials at a few
    # locations). Only the FOCAL trial is constrained -- a pipeline may still
    # pull training trials from anywhere. `trials` pins the focal set to specific
    # studyName values directly (the most direct way to say "optimize for these
    # trials"); it combines (AND) with the other fields. Example:
    #   target_domain = list(programs  = c("Cornell", "OSU"),
    #                        years     = 2018:2025,
    #                        locations = c("Ithaca", "Wooster"),
    #                        trials    = c("2025_AYT_Aurora", "YT_Urb_25"))
    target_domain = list(
      programs  = NULL,          # character vector of programName values, or NULL
      years     = NULL,          # integer vector of years, or NULL
      locations = NULL,          # character vector of locationName values, or NULL
      # character vector of studyName (trial) values, or NULL
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
    # Durable, resumable state -> perm_dir (home on a remote server; see above).
    db_path     = file.path(perm_dir, "state", "evals.sqlite"),
    stop_file   = file.path(perm_dir, "state", "STOP"),
    report_path = file.path(perm_dir, "state", "report.md"),
    log_dir     = file.path(perm_dir, "logs"),
    # The large, regenerable cache stays on the WORK disk (fast, expendable) even in remote
    # mode -- it is backed up to durable storage periodically instead (cache_backup_dir).
    cache_dir   = here::here("cache"),

    # ---- cache backup -----------------------------------------------------
    # Periodically rsync the (regenerable) cache to durable storage so an abrupt kill on a
    # scratch/work disk does not force re-downloading. The sync is ADDITIVE (cache files are
    # write-once) and excludes the transient raw_project/ VCFs. On a remote server it defaults
    # to <OPTIMIZER_PATH>/cache; NULL disables it. Restored back to cache_dir at startup if the
    # work cache is empty (e.g. a fresh node after a scratch purge). An in-process sync cannot
    # survive a SIGKILL of R -- see the README for an external cron/systemd safety net.
    cache_backup_dir   = if (remote_server) file.path(Sys.getenv("OPTIMIZER_PATH"), "cache") else NULL,
    # Minimum minutes between in-process cache backups (0 disables). Throttled on wall time,
    # not iterations, so checkpoint_every = 1 doesn't rsync every step.
    cache_sync_minutes = 120,

    # ---- several workers on one store -------------------------------------
    # Multiple worker processes may run against ONE db_path and ONE cache_dir, each
    # evaluating a different configuration and all contributing to the same archive (see
    # README "Running several workers" and run_workers.sh). SQLite handles this fine here --
    # one small INSERT per multi-minute evaluation -- given the WAL + busy_timeout that
    # open_store sets. Two things follow:
    #
    #  1. db_path MUST BE ON LOCAL DISK. WAL coordinates through an mmap'd -shm file, which a
    #     network filesystem does not provide; on NFS the pragma silently fails (open_store
    #     warns) and concurrent writers risk corruption. So on a cluster point db_path at
    #     /workdir and set db_backup_path to durable storage.
    #  2. Only ONE worker does the shared-resource housekeeping (the report, the cache rsync,
    #     the store backup) -- eight processes rsyncing one tree and overwriting one
    #     report.md is pure waste. That is the "leader", worker 1.
    #
    # worker_id comes from the OPTIMIZER_WORKER environment variable that run_workers.sh
    # sets, so the same settings.R serves every worker; a lone run is worker "1" (the leader)
    # with no environment set up at all.
    worker_id   = { w <- trimws(Sys.getenv("OPTIMIZER_WORKER")); if (nzchar(w)) w else "1" },
    is_leader   = !nzchar(trimws(Sys.getenv("OPTIMIZER_WORKER"))) ||
                  identical(trimws(Sys.getenv("OPTIMIZER_WORKER")), "1"),
    # Where the leader copies the live store (VACUUM INTO; safe on a live database). This is
    # what makes a local-disk db_path safe to lose. NULL disables it. Defaults to the durable
    # location the store used to live at, so remote runs keep the same recovery story.
    db_backup_path    = if (remote_server) file.path(Sys.getenv("OPTIMIZER_PATH"),
                                                     "state", "evals_backup.sqlite") else NULL,
    db_backup_minutes = 30,

    # ---- concurrent cache access ------------------------------------------
    # Workers share cache_dir, so the expensive downloads are shared too -- but only if they
    # do not race. get_project_dosage holds a per-project lock across download -> parse ->
    # cache-write (the download deletes any stale copy first, so two workers on one project
    # would otherwise delete each other's in-flight file). A worker that finds the lock held
    # waits for the winner's result rather than duplicating a multi-GB download.
    lock_wait_minutes  = 60,   # how long to wait for another worker's download to finish
    lock_stale_minutes = 90    # break a lock older than this (its holder died)
  )
  s <- .apply_overrides(defaults, .local_overrides())
  # The optimizer targets exactly one scheme per run (validate the effective value,
  # so a settings.local.R override is checked too).
  if (length(s$optimize_scheme) != 1L || !is.character(s$optimize_scheme) ||
      !(s$optimize_scheme %in% s$schemes))
    stop("optimize_scheme must be a single scheme that is one of `schemes` (",
         paste(s$schemes, collapse = ", "), "); got: ",
         paste(s$optimize_scheme, collapse = ", "))
  s
}
