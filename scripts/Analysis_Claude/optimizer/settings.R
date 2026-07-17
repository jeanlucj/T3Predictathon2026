# settings.R
#
# Every knob in one place. The driver and the optimizer read from here, so you
# tune the run by editing this file (or overriding fields after sourcing it).

optimizer_settings <- function() {
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

    # Memory budget for a single project's dense dosage matrix. A project whose
    # samples x markers x 4 bytes would exceed this is auto-thinned (every k-th marker)
    # until it fits -- a 7.5M-marker GBS panel cannot be held densely and does not need
    # every marker for a GRM. The effective thinning is recorded in the dosage cache name.
    dosage_budget_bytes = 2e9,

    min_trial_acc    = 30,       # skip trials with fewer genotyped accessions
    min_train_trials = 3,        # skip focal trials we cannot assemble a training set for
    max_sample_fail  = 25,       # consecutive trial-sampling failures before halting
                                 # (guards against a too-restrictive domain or a
                                 #  down network spinning the loop forever)
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
      trials    = c("CornellMaster_2007_McGowan", "CornellMaster_2008_Helfer",
                    "CornellMaster_2008_Ketola", "CornellMaster_2008_Snyder",
                    "CornellMaster_2009_Helfer", "CornellMaster_2009_Ketola",
                    "CornellMaster_2009_McGowan", "CornellMaster_2010_Helfer",
                    "CornellMaster_2010_Ketola", "CornellMaster_2010_Snyder",
                    "CornellMaster_2011_Helfer", "CornellMaster_2011_Ketola",
                    "CornellMaster_2011_McGowan", "CornellMaster_2012_Helfer",
                    "CornellMaster_2012_Ketola", "CornellMaster_2012_Snyder",
                    "CornellMaster_2013_Helfer", "CornellMaster_2013_Ketola",
                    "CornellMaster_2013_McGowan", "CornellMaster_2014_Helfer",
                    "CornellMaster_2014_Ketola", "CornellMaster_2014_Snyder",
                    "CornellMaster_2015_Helfer", "CornellMaster_2015_Ketola",
                    "CornellMaster_2015_McGowan", "CornellMaster_2016_Helfer",
                    "CornellMaster_2016_Ketola", "CornellMaster_2016_Snyder",
                    "CornellMaster_2017_Ketola", "CornellMaster_2017_Snyder",
                    "CornellMaster_2018_Helfer", "CornellMaster_2018_Ketola",
                    "CornellMaster_2019_Ketola", "CornellMaster_2020_Helfer",
                    "CornellMaster_2021_Snyder", "CornellMaster_2022_Helfer",
                    "CornellMaster_2023_McGowan")
    ),

    # ---- paths ------------------------------------------------------------
    db_path     = here::here("state", "evals.sqlite"),
    stop_file   = here::here("state", "STOP"),
    report_path = here::here("state", "report.md"),
    log_dir     = here::here("logs"),
    cache_dir   = here::here("cache")
  )
}
