# A self-improving genomic-prediction pipeline optimizer

This subsystem takes the five Predictathon submissions, decomposes each into the
**six subtasks** every genomic-prediction algorithm must perform, and treats the
*choice of method and parameters for each subtask* as a configuration to be
optimized. It then searches that configuration space for a pipeline that
predicts grain yield in a **randomly chosen T3 trial** better, on average, than
any single submitted algorithm.

It is designed to run unattended in the background, accumulate evidence over
time, and get progressively better at predicting *any* trial -- not just the
nine trials used in the Predictathon.

---

## 1. The six subtasks (the configuration space)

Every submission, however it is written, does these six things. Each submission
does each one somewhat differently, and each has tunable knobs. We catalogue the
methods so they can be freely recombined (`R/config_space.R`).

| # | Subtask | Methods drawn from submissions / literature |
|---|---------|---------------------------------------------|
| A | **Select training trials** from T3 for a focal trial | accession_overlap with a `primary_only` toggle — direct germplasm-overlap neighbours of the focal trial, optionally plus a second "friends-of-friends" tier (P1/P2); top-k genomic/environmental similarity (P3); same breeding program (P5) |
| B | **Preprocess** training phenotypes into per-accession targets | per-trial BLUE `y~germ+rep+block` then average (P1), trial-centering `y - trial_mean + grand_mean` (P4), two-stage BLUP shrinkage, raw mean (baseline); optional G×E environmental weighting (P4) |
| C | **Select genotyping data** (which projects / VCFs) | focal+one-hop with bridges — the focal trial's panel plus any panel sharing `min_bridge` accessions with it (Analysis_Claude), best single project (P5), all projects / global union (P2) |
| D | **Build the relationship / kernel** from markers | VanRaden GRM on a merged matrix (P1/P2/P4), Wishart-EM combine of per-project GRMs (P3/P5), RKHS Gaussian kernel (P5 CV00) |
| E | **Train the prediction model** | rrBLUP GBLUP (P1), sommer multi-kernel G+E (P3), RKHS (P5); orthogonally, how the ridge is chosen — `lambda_select` ∈ REML (P1) / fixed (P2) / leave-one-out grid (P4) |
| F | **Predict** on the focal (test) trial | direct BLUP, conditional expectation `G21 G11^-1 u` (P1), blend with observed BLUE (P4), mean fallback when overlap is too small |

A **configuration** is one choice of method + parameters for each subtask -- a
point in a mixed categorical/continuous space. `R/config_space.R` is the single
source of truth: it lists the methods and parameter ranges, can sample a random
valid configuration, encode one as features for the surrogate, and recombine two
configurations (crossover) or perturb one (mutation).

To extend the system with a new idea from the literature, you add a method name
and its parameters to `config_space.R` and a branch in the corresponding
dispatcher in `R/pipeline.R`. Nothing else changes.

## 2. The objective: predict *any* trial

We do **not** optimize accuracy on the nine Predictathon trials. We optimize the
**expected predictive ability over randomly sampled T3 trials** drawn from a
configurable **target domain** (`settings$target_domain`). With the domain
unconstrained, the result is the best all-round pipeline; constraining it to a
list of breeding programs, years, locations, and/or specific trial names tailors the
pipeline to that subpopulation (e.g. the best pipeline for one program's recent trials). Only the
focal trial is constrained -- a pipeline may still pull training trials from
anywhere. The trait being predicted is also configurable
(`settings$focal_trait`); it defaults to grain yield but can be any numeric trait.

One *evaluation* is: sample a random T3 trial that has grain-yield data and
enough genotyped lines; run the configuration's full pipeline to predict that
trial under a cross-validation scheme (CV0 or CV00); score it by the Pearson
correlation between predicted and observed per-accession BLUEs on the held-out
focal-trial accessions. This is exactly the metric the Predictathon used
(`T3_predictathon_scripts/analysis/30_All_accuracies...Rmd`).

Because every evaluation uses a *fresh* random trial, a configuration that
scores well has to generalize. The optimization target is the mean score across
trials, so the search is driven toward configurations that are robust across the
whole database, not tuned to nine specific environments.

This makes the objective **noisy and expensive**: noisy because a single trial
is a high-variance sample of a configuration's true quality; expensive because
each evaluation pulls phenotypes and genotypes and fits mixed models. The
optimizer is built around those two facts.

## 3. The optimizer: random-forest-surrogate Bayesian optimization, with evolutionary recombination

For a **noisy, expensive, mixed categorical+continuous** search space, the
literature consensus (SMAC; Hutter et al. 2011) is a **random-forest-surrogate
Bayesian optimization**: it handles categorical choices and scales with the
number of evaluations better than Gaussian-process BO, and is far more
sample-efficient than a pure genetic algorithm. We use that as the backbone and
fold in evolutionary operators because the user's framing -- *recombine subtask
methods from different algorithms* -- is literally crossover.

The loop (`R/optimizer.R`), one iteration:

1. **Surrogate fit.** Train a surrogate on all evaluations so far: a bagged
   ensemble of regression trees (`rpart`) mapping a configuration's features to
   its score. This is SMAC's random forest, written transparently
   (`R/surrogate.R`): each tree sees a bootstrap sample; the ensemble mean is the
   predicted score and the spread across trees is the uncertainty. Trees handle
   categoricals and missing (irrelevant) parameters natively.

2. **Candidate generation.** Propose many candidate configurations three ways,
   so the search both exploits good structure and keeps exploring:
   - **Crossover** of two elite configurations -- swap whole subtask blocks
     (e.g. take P1's training-trial selection, P5's kernel, P4's blending). This
     is the "recombine the submitted algorithms" operation.
   - **Mutation** of an elite -- resample one subtask's method/parameters.
   - **Random** sampling of fresh configurations.

3. **Acquisition.** Score every candidate by **Expected Improvement** under the
   surrogate (mean and uncertainty), which balances "predicted to be good" with
   "uncertain, worth probing." Pick the top candidate(s) to actually run.

4. **Evaluate** the chosen configuration on a freshly sampled random trial,
   append `(config, trial, scheme, score)` to the store, and update the
   incumbent (best configuration by surrogate-predicted mean over trials).

Early iterations (before there is enough data to fit a surrogate) are pure random
+ the five submitted configurations as seeds, so the submissions are an explicit
baseline the search must beat.

Why this and not a plain GA or plain grid search: grid search wastes the
expensive evaluations; a plain GA needs a large population per generation under
noise; the surrogate lets us learn from *every* evaluation -- including ones on
different trials -- and spend real evaluations only where the model says it is
worth it.

## 4. Running in the background, getting better over time

- **Persistent store** (`R/store.R`, SQLite via `DBI`/`RSQLite`): every
  evaluation is a row `(config_hash, config_json, trial_id, scheme, score,
  n_test, status, reason, detail, seconds, ts)`. `status` ∈
  `ok | infeasible | suspect | error`; `reason` is a machine code and `detail`
  holds the failure funnel. The optimizer's entire state is reconstructable from
  this table, so the process can be killed and resumed at any time with no loss.
- **Caches** (`cache/`): the trial catalogue, accession lists, observations,
  genotyping-project lists, archived VCFs, and VCF→dosage matrices are memoized to
  disk and reused across configurations and iterations, so the cost of touching a
  given trial is paid once (see the outputs table in §7).
- **Driver** (`run_optimizer.R`): a loop with a wall-clock budget, an iteration
  budget, and a stop-file (`state/STOP`) so you can halt it cleanly. It
  checkpoints after every evaluation and writes a heartbeat to `logs/`.
  Launch with `nohup Rscript run_optimizer.R &` (see README).
- **Reporting** (`R/report.R`): the learning curve (best score vs. evaluations),
  the current incumbent configuration, and a marginal analysis of which subtask
  *methods* most improve the score -- so the run is interpretable, not a black
  box.

## 5. Verifying it, and guarding against data-hiding bugs

Two verification layers, because the failure modes differ:

- **The math (offline).** `R/evaluate.R` has a **SIMULATE mode**: a synthetic
  genomic-prediction world with a known generative model in which some subtask
  choices are genuinely better than others. `tests/test_sim_loop.R` runs the full
  optimizer against it offline in seconds and checks the incumbent rises above the
  random/seed baseline. The Tier-1 unit tests (`tests/test_subtasks.R`,
  `test_config_space.R`) cover the deterministic subtask and configuration-space logic.
- **The plumbing (real data).** Because the optimizer *skips and continues* past a
  trial it cannot evaluate, a bug that hides real data looks identical to a trial
  that is genuinely infeasible. Two mechanisms surface that: the **failure-log
  condition system** (`R/conditions.R` — `infeasible` / `fatal` / `suspect`, with a
  per-failure data "funnel" so a many-in/zero-out *cliff* is visible), and the
  **canary bug-oracle** (`R/diagnostics.R`), which validates the live pipeline
  against an independent ground truth (the five teams' submission files) via a
  calibrate-then-freeze step and runs at startup in real mode.

See `EVALUATION.md` for how to run and read both layers.

## 6. File map

```
optimizer/
  README.md            <- how to run it (user-facing)
  EVALUATION.md        <- evaluation & validation runbook (arm_evaluation, tests, data-flow, tooling)
  DESIGN.md            <- this file (what the optimizer is + where everything lives)
  BACKGROUND.md        <- challenges (statistical + data-management) and how each is met
  run_optimizer.R      <- background entry point (budget, stop-file, checkpoint, canary startup)
  settings.R           <- one place for all knobs (budgets, paths, simulate, canary trials)
  R/
    config_space.R     <- the six-subtask config space: sample / encode / crossover / mutate
    pipeline.R         <- parameterized pipeline: run a config's six subtasks on a trial
    data_access.R      <- BrAPI + caching: sample trials, pull phenotypes & dosages
    evaluate.R         <- run a config on a trial under CV0/CV00 -> score (+ SIMULATE)
    conditions.R       <- infeasible / fatal / sample_failed signals + funnel encoding
    surrogate.R        <- bagged-rpart RF surrogate: mean, uncertainty, Expected Improvement
    optimizer.R        <- candidate generation (crossover/mutate/random) + acquisition
    store.R            <- SQLite results store (resumable state)
    report.R           <- learning curve, incumbent, method importance, failure log
    diagnostics.R      <- canary oracle (anchor / calibrate / check / coverage) + diagnose_trial
    seeds.R            <- the five submissions encoded as starting configurations
  tests/
    run_all.R          <- run all fast test files (--all also runs the sim loop)
    test_subtasks.R    <- Tier 1 deterministic core (subtask oracles)
    test_config_space.R<- config-space invariants (sample/encode/crossover round-trips)
    test_sim_loop.R    <- end-to-end optimizer run in SIMULATE mode (offline)
  cache/   state/   logs/   <- generated at runtime (see §8)
```

## 7. Function inventory (where everything is)

One line per function, grouped by file. `.name` = file-private helper.

**`settings.R`** — `optimizer_settings()`: all knobs in one list.

**`run_optimizer.R`** — `run_optimizer()`: the resumable loop (budget, STOP,
checkpoint, startup canary); `optimizer_step()`: pick config → sample trial →
evaluate under each scheme → store.

**`R/config_space.R`** (the configuration space) — `SUBTASKS` (the search-space spec);
`sample_block()`/`sample_config()`: draw a random config; `crossover()`,
`mutate_config()`: evolutionary operators; `repair_config()`: enforce
param-applies-to-method invariant; `canonical_keys()`, `config_hash()`;
`feature_schema()`/`configs_to_features()`: encode configs for the surrogate;
`.param_applies`/`.sample_param`/`.method_key`/`.param_key`.

**`R/seeds.R`** — `seed_configs()`: the five submissions as starting configs;
`.make_seed()`: build one from overrides + safe defaults.

**`R/pipeline.R`** (the parameterized pipeline) — `run_pipeline()`: the six subtasks
end to end; dispatchers `select_training_trials()`, `build_targets()`,
`choose_geno_sources()`, `build_kernel()`, `train_model()`, `predict_test()`;
`mask_cv()`: CV0/CV00 masking; helpers `.find_related`, `.trial_similarity`,
`.blue_per_trial`, `.per_acc_blue`, `.group_by_panel` (group projects into protocols
by marker overlap — protocol *ids* cannot do this, since a `V2`/`v2.1` protocol is the
same protocol against another reference genome), `.prune_redundant` (drop re-called
duplicate projects within a group), `.merge_markers`, `.qc_markers`, `.best_panel`,
`.onehop_filter` (subtask C's hop: keep the focal trial's panel plus any panel sharing
`min_bridge` accessions with it — one hop, not the transitive closure, which is what makes
`focal_plus_onehop` narrower than `all_projects`),
`.ridge_blup`/`.loo_lambda` (GBLUP at a given ridge; the PRESS-identity leave-one-out grid
search ported from Prediction4),
`.bridge_accessions` (accessions genotyped in >1 protocol group — the shared rows
`em_combine` stitches its per-panel GRMs through), `.vanraden` (VanRaden GRM over the
needed accessions, centred and scaled on the panel population's allele frequencies —
**not** `rrBLUP::A.mat`, which assumes `{-1,0,1}` coding and silently drops every
alt-major marker when fed `{0,1,2}` dosages), `.fit_sommer_GE`.

Subtask C returns one **full-population** dosage matrix per protocol group; subtask D
runs marker QC and estimates allele frequencies on that population and only then forms
the GRM over the accessions the trial needs. Because `K_ij` depends on other accessions
only through the allele frequency, this is exactly the `[need, need]` submatrix of the
population's GRM, at `O(n_need²·m)` rather than `O(n_pop²·m)`. Under `em_combine`, the
per-panel GRMs additionally retain the **bridge** accessions (genotyped in >1 protocol
group) so `covariance_combiner()` can stitch the panels through their shared rows; the
bridges are dropped from the combined GRM, which is subset back to `need`.

**`R/data_access.R`** (BrAPI + caching) — `cached()`: the on-disk memoizer;
`trial_catalog()`: focal-trait trials + metadata + coords; `sample_real_trial()`,
`build_trial_descriptor()`, `.trial_descriptor()`: focal trial descriptors;
`get_observations()`/`.obs_tibble()`: phenotypes; `get_trial_accessions()` (also the
substrate `.find_related`/`.trial_similarity` compute germplasm overlap from);
`projects_for_accessions()`: genotyping **project** ids; `get_project_dosage()` +
`.ensure_project_vcf()`/`.vcf_complete()`/`.vcf_header()`/`.vcf_stat()`/`.vcf_to_dosage()`:
VCF → dosage — a base-R **streaming, chunked** reader (bounded memory; skips malformed
variant lines; rejects transposed/non-VCF archives; gzip/BGZF transparent). Each project is
parsed **once**, at the densest thin its dense size fits `settings$dosage_budget_bytes`
(`.cache_thin()`), and that budget is the **only** control on marker density — the pipeline
no longer asks for a thin of its own (`marker_thin` was removed from the search space: a
request could only be served by subsetting this cache, so it was a silent no-op whenever the
budget had already thinned harder). `get_project_dosage()` keeps a `marker_thin` argument for
manual probes, served by column-subsetting (`.find_densest_dosage()` + read-time subsample)
rather than by re-downloading. A `stat_<id>` / `unparseable_<id>` cache lets a coarse-sizing decision
and an unparseable-archive verdict persist without the (deleted) VCF. A transient VCF-download
failure is tracked **per session** (`.vcf_download_plan`, in RAM): effort decreases per prior
failure and the project is skipped after `settings$vcf_max_download_attempts`, so a timing-out
archive is not re-stormed on every covering trial (resets next run).
`.focal_trait_parts`/`.matches_trait`/`.apply_target_domain`.

**`R/evaluate.R`** — `sample_trial()`: real or synthetic trial; `score_predictions()`:
the Predictathon metric; `evaluate_config_on_trial()`: run + score, branching on the
condition classes; `.sim_true`/`.sim_evaluate`: the SIMULATE objective.

**`R/conditions.R`** — `infeasible()`, `fatal()`, `sample_failed()`: the three
"didn't work" signals; `funnel_string()`: encode the stage-count funnel.

**`R/surrogate.R`** — `fit_surrogate()`/`predict_surrogate()`: bagged-rpart RF
(mean + uncertainty); `expected_improvement()`: the acquisition function.

**`R/optimizer.R`** — `aggregate_scores()`: per-config mean (denoise);
`get_elites()`, `fresh_random()`, `propose_candidates()`: candidate generation;
`incumbent_config()`: robust best-so-far; `choose_config()`: the phase logic
(seeds → replicate → random → acquisition); `choose_trial()`: the trial to pair
with it, excluding those that configuration has already been run on.

**`R/store.R`** — `open_store()`/`close_store()`; `store_eval()`/`read_evals()`;
`n_evals()`, `tried_hashes()`; `config_to_json()`/`config_from_json()`
(NA-preserving).

**`R/report.R`** — `write_report()`: the Markdown snapshot; `format_config()`,
`method_importance()`, `failure_summary()`.

**`R/diagnostics.R`** (validation) — `canary_anchor()`: independent per-trial counts
from submission files; `calibrate_canary_trials()`/`print_calibration()`: our counts
vs the anchor; `canary_configs()`: the frozen coverage configs; `canary_coverage()`:
coverage check; `check_canaries()`: the startup oracle; `diagnose_trial()`: per-trial
funnel + raw re-derivation; `canary_config()`/`.canary_filler()`: permissive configs.

## 8. Where outputs go

| Directory | Written by | Contents |
|---|---|---|
| `cache/` | `cached()` in `data_access.R` | `trial_catalog.rds`; `acc_<id>.rds`; `obs_<sid>.rds`; `proj_<hash>.rds`; `raw_project_<id>.vcf`; `dosage_<id>[_thin<k>]_sz<size>.rds`; `stat_<id>.rds`; `unparseable_<id>.rds`. Regenerable; the only expensive thing here is the download/extraction, paid once. |
| `state/` | `store.R`, `report.R` | `evals.sqlite` (the single source of truth), `report.md` (rewritten each checkpoint), `STOP` (touch to halt). |
| `logs/` | `run_optimizer.R` (via the launch redirect) | `run.out` — per-iteration heartbeat, the startup canary line, `CANARY ALARM` / `FATAL` / `step error` messages. |

Paths are all set in `settings.R` (`cache_dir`, `db_path`, `report_path`,
`stop_file`, `log_dir`) — see README.md for splitting them across disks on BioHPC.

## References

- Hutter, Hoos & Leyton-Brown 2011, *Sequential Model-Based Optimization for
  General Algorithm Configuration* (SMAC) -- RF-surrogate model-based search.
- Jarquin et al. 2017, *The Plant Genome* -- CV0/CV00 scenarios.
- VanRaden 2008 -- GRM. Endelman 2011 -- rrBLUP. Akdemir CovCombR -- Wishart-EM
  combine of partial relationship matrices.
