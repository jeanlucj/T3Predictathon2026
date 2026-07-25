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

  list(
    # ---- mode -------------------------------------------------------------
    # TRUE  = synthetic offline world (fast; for verifying the machinery).
    # FALSE = real T3 pipeline on real data (network + heavy model fits).
    simulate = TRUE,

    # ---- which CV schemes to score each evaluation under ------------------
    schemes = c("CV0", "CV00"),

    # ---- search behaviour -------------------------------------------------
    n_random_init       = 25,    # random configs before the surrogate takes over
    incumbent_min_reps  = 2,     # trials a config needs before it can be "incumbent"
    reeval_prob         = 0.15,  # chance an iteration re-evaluates the incumbent
    ei_xi               = 0.01,  # Expected-Improvement exploration constant
    n_elites            = 8,     # elites that seed crossover/mutation
    n_cross             = 60, n_mut = 60, n_rand = 60,  # candidate pool sizes
    ntree               = 300,   # trees in the random-forest surrogate

    # ---- budget / background control -------------------------------------
    max_iters     = Inf,         # stop after this many evaluations
    max_hours     = Inf,         # ... or this much wall-clock, whichever first
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

    # Memory budget for a single project's dense dosage matrix (samples x markers x 4 B).
    # This is the SOLE control on cached marker density: each project is downloaded and
    # parsed ONCE, at the densest thinning that fits this budget (thin 1 = full markers for
    # anything that fits; a 7.5M-marker GBS panel is thinned to fit). A config's requested
    # marker_thin is then derived by column-subsetting that one cache -- never re-parsed --
    # so set this deliberately: a bigger budget caches denser (better GRMs, more disk/RAM
    # per project), smaller thins large panels harder. Keep it well under total RAM, since
    # choose_geno_sources loads all covering projects and the GRM needs headroom too.
    dosage_budget_bytes = 2e9,

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
    cache_sync_minutes = 120
  )
}
