# Developer & validation guide

This is the document for **working on, evaluating, and validating** the optimizer
-- for me now and in six months. `README.md` is for running it; `DESIGN.md` is the
static architecture (every function, where outputs go); `BACKGROUND.md` is why it
is built this way and the challenges it must handle. This file connects those: it
maps the test suite, traces how one evaluation flows function-to-function, shows
how to step through any single function, documents how the code touches the
cache / state / logs, and explains the validation/debugging tooling.

The system is genuinely complex -- the task and the messy T3 data layer make that
unavoidable -- so the goal here is that **reviewing this document + the tests is
faster than re-reading all the code.**

---

## 1. Orient yourself in 60 seconds

- Code is in `R/`; one entry point `run_optimizer.R`; all knobs in `settings.R`.
- **Offline** everything runs in `simulate = TRUE` (a synthetic objective, no
  network). **Real mode** (`simulate = FALSE`) pulls from T3/Wheat over BrAPI.
- Two ways to gain confidence without a week of compute: the **test suite** (§2,
  offline, deterministic) and the **canary oracle** (§6, real data, checks the
  live pipeline against an independent ground truth).
- Durable state is two dirs: `state/` (SQLite store + report) and `cache/`
  (downloaded data). Everything is recomputed from `state/evals.sqlite` on restart.

To load the whole subsystem into an interactive R session:

```r
setwd("scripts/Analysis_Claude/optimizer")
library(tidyverse)
source("settings.R")
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
s <- optimizer_settings()
```

---

## 2. The test suite

### Principle: oracle tests, not snapshot tests

A **snapshot test** asserts "the function returns whatever it currently returns."
It re-encodes the implementation, passes even when the code is wrong, and takes as
long to audit as the code. We avoid those. Every test has an **independent oracle**
-- the expected answer is knowable without reading the implementation. Four kinds:

- **Closed-form** -- correlation of a vector with itself is 1; a GRM among
  genetically identical lines has identical rows; centering by trial mean removes a
  constant added to one trial.
- **Property / invariant** -- Expected Improvement is ≥ 0; a relationship matrix is
  symmetric; two-stage BLUP shrinks effects toward the mean; QC drops a monomorphic
  marker.
- **Metamorphic / differential** -- double the phenotypes → BLUEs double; permute
  accessions → predictions permute identically; **CV00 excludes exactly the
  focal-trial phenotypes CV0 keeps**; raising `blend_obs_w` changes CV0 predictions
  but is a no-op under CV00.
- **Tiny hand-computable** -- a 3-marker × 2-sample VCF whose dosage matrix you can
  write down by hand.

Because the tests *run*, a test that encodes a wrong expectation usually **fails
against correct code**, surfacing the disagreement to adjudicate rather than
silently lying. Each test file comments its oracle inline.

Unit tests cover the **deterministic, offline** logic -- most of the optimizer.
They do **not** verify that the live BrAPI responses have the shape
`R/data_access.R` assumes; that is what the canary oracle (§6) is for. *Unit tests
prove the math; the canary proves the plumbing.*

### Running, and expected outcomes

```bash
cd scripts/Analysis_Claude/optimizer
Rscript tests/run_all.R          # fast files; one PASS/FAIL summary
Rscript tests/run_all.R --all    # also the slow end-to-end sim loop
Rscript tests/test_subtasks.R    # a single file
```

Each file is self-contained: it sources the subsystem, runs hand-rolled `check()`
assertions (no `testthat`, to match project style), prints a pass/fail count, and
exits non-zero on failure. `run_all.R` runs each in a fresh process and aggregates.

| Command | Expected |
|---|---|
| `tests/test_config_space.R` | ~5 genome invariants checked across ~400 sampled/recombined configs → prints `config_space tests: 8007 passed, 0 failed` (8007 = individual assertions, not distinct cases) |
| `tests/test_subtasks.R` | `Tier 1 subtask tests: 45 passed, 0 failed` |
| `tests/run_all.R` | `2/2 test files passed` |
| `tests/test_sim_loop.R` (or `run_all.R --all`) | `PASS: optimizer beats submissions and improves over random search`, exit 0 |

If a count drifts after a change, that is the regression signal -- reconcile it
before moving on.

### Tier 1 — deterministic core *(implemented: `tests/test_subtasks.R`)*

| Component (function) | Oracle |
|---|---|
| `score_predictions` (`evaluate.R`) | cor(x, a·x+b)=1 for a>0, =−1 for a<0; only shared accession names join; <5 overlap → NA `too_few_overlap`; constant → NA `constant`. |
| `.vcf_to_dosage` (`data_access.R`) — guards the transposition bug | hand 3-marker × 2-sample VCF → **accessions × markers**, `0/0→0,0/1→1,1/1→2`, correct names; sample subset → one row; `thin=2` keeps every 2nd marker. |
| `.qc_markers` (`pipeline.R`) | monomorphic dropped (MAF=0); high-missing dropped; remaining NA imputed (mean vs mean_round); errors below 50 markers. |
| `.merge_markers` (`pipeline.R`) | keeps marker intersection; de-dups repeated accessions; falls back to largest project when <50 shared; single project as-is. |
| `build_kernel` (`pipeline.R`) | VanRaden K symmetric, clones give identical rows; `ridge` raises the diagonal by exactly `ridge`; RKHS K symmetric, diag 1, off-diag (0,1], clones→1; `em_combine` with one project equals the plain VanRaden GRM. |
| `build_targets` + `.blue_per_trial` (`pipeline.R`) | `raw_mean` = mean of per-trial means; `trial_center` invariant to a constant added to one trial; `blue_lm` recovers group means; outlier removal at `z_thr`; `two_stage_blup` shrinks spread; `per_trial_z` → mean 0 sd 1. |
| `predict_test` (`pipeline.R`) | conditional expectation with a clone test line → `mu + u[clone-source]`; `blend_obs_w` blends under CV0, no-op under CV00; fallback below `min_overlap` → all predictions `mean(targets)`. |
| `mask_cv` (`pipeline.R`) | CV0 keeps every row; CV00 drops exactly the focal-accession rows (CV0 set minus CV00 set = the focal rows). |

### Genome invariants *(implemented: `tests/test_config_space.R`)*

Sampling, encoding, JSON round-trips, and crossover/mutation all produce
well-formed configurations. Property-based: ~5 invariants (canonical key order,
valid method, param-applies⇔not-NA, JSON round-trip hash, crossover/mutation
validity) are asserted across ~400 randomly sampled and recombined configs, which
is why the counter prints `8007 passed` — those are individual assertions, not
8007 distinct test cases.

### End-to-end *(implemented: `tests/test_sim_loop.R`, slow)*

Runs the whole optimizer against the synthetic world (`.sim_true`/`.sim_evaluate`
in `evaluate.R`) and asserts the discovered incumbent beats the submitted seeds and
that the surrogate phase improves on random search. This is the regression test for
any change to the search space or optimizer; run it after such changes.

### Tiers 2–4 — planned

Specified for later; not yet implemented. Oracles already worked out:

- **Tier 2 — optimizer mechanics.** `expected_improvement` (≥0; ↑ with mean; ↑ with
  sd when mean<best; →0 as sd→0); `fit_surrogate`/`predict_surrogate` (NULL below
  `min_obs`; monotone signal → high rank-correlation; sd floor); `aggregate_scores`/
  `get_elites`/`incumbent_config` (mean excludes failures, `n` counts them; top-k;
  `min_reps` honored); `crossover`/`mutate_config`/`propose_candidates` (each child
  block matches one parent; mutation changes ≥1 block; candidates repaired);
  `choose_config` phases (empty store → seed; < `n_random_init` → `random_init`;
  enough → `acquisition`).
- **Tier 3 — store / data helpers / reporting.** `store_eval`↔`read_evals`
  round-trip incl. NA score and `detail`; `.matches_trait`/`.focal_trait_parts`/
  `.apply_target_domain`; `.obs_tibble`+`get_observations` on a fake response;
  `cached` write/read; `report.R` (`format_config`, `method_importance`,
  `failure_summary`, "BEATS submissions" only when incumbent > best seed).
- **Tier 4 — network subtasks via a mock `conn`.** Inject a fake `conn` (closures
  for `$wizard`/`$search`/`$vcf_archived`) + a fake `trial_catalog` and test the
  *logic* of `select_training_trials`, `.trial_similarity`, `choose_geno_sources`
  (most-similar ranks first; `same_program` respects `prog_cap`; `top_k_similar`
  returns `k`; `best_single_project` picks the most-covering project). Needs a
  `tests/fixtures.R`.

---

## 3. Data-flow trace: one evaluation, function by function

The unit of work is **one configuration on one trial under one CV scheme → a
score**. Following the output→input handoffs:

```
run_optimizer()                         [run_optimizer.R]  the loop; budget, STOP, checkpoint
└─ optimizer_step(con, settings, conn)  [run_optimizer.R]
   ├─ choice <- choose_config(con, settings)        [optimizer.R]   -> {cfg, source}
   │     read_evals(con) -> aggregate_scores()      [store.R/optimizer.R]
   │     phase: seeds (seed_configs) | random (fresh_random) |
   │            reeval (incumbent_config) | acquisition:
   │            fit_surrogate(configs_to_features(scored), mean_score)  [surrogate.R/config_space.R]
   │            propose_candidates(crossover/mutate_config/fresh_random) [optimizer.R/config_space.R]
   │            predict_surrogate() -> expected_improvement() -> pick    [surrogate.R]
   ├─ trial  <- sample_trial(settings, conn)         [evaluate.R]   -> trial descriptor
   │     real: sample_real_trial()                   [data_access.R]
   │        trial_catalog() (focal-trait + target_domain filter)
   │        get_trial_accessions(); get_observations() feasibility check
   │        .trial_descriptor()  -> {id, accessions, program, location, year, lat/long/elev}
   └─ for scheme in settings$schemes:
        ev <- evaluate_config_on_trial(cfg, trial, scheme, settings, conn)  [evaluate.R]
              real: preds_obs <- run_pipeline(cfg, trial, scheme, settings, conn)  [pipeline.R]
                    score_predictions(preds_obs$pred, preds_obs$obs)
              sim:  .sim_evaluate(cfg, trial, scheme)
        store_eval(con, cfg, trial$id, scheme, ev$score, ev$n_test,
                   ev$status, ev$reason, ev$detail, ev$seconds)              [store.R]
```

Inside `run_pipeline()` (`pipeline.R`), the six subtasks chain output→input:

```
get_observations(focal) -> .per_acc_blue()                 -> obs (named numeric, the truth)
select_training_trials(cfg, trial)                          -> train_ids
   ( .find_related() / .trial_similarity() / trial_catalog filter )
get_observations(train_ids) -> mask_cv(scheme) -> build_targets(cfg)  -> targets (named numeric)
   ( .blue_per_trial() ; optional G×E weighting via .trial_similarity )
choose_geno_sources(cfg, names(targets), focal_acc)         -> dosage_list (accessions×markers)
   ( projects_for_accessions() -> get_project_dosage() -> .ensure_project_vcf()/.vcf_to_dosage() )
build_kernel(cfg, dosage_list)                              -> K (relationship/kernel matrix)
   ( .merge_markers() , .qc_markers() ; em_combine via covariance_combiner )
train_model(cfg, targets[train_in], K, ...)                -> fit   ( .fit_sommer_GE() )
predict_test(cfg, fit, K, train_in, test_in, targets, scheme) -> pred (named numeric)
run_pipeline returns list(pred = pred, obs = obs)          -> score_predictions()
```

Any subtask that cannot complete for this (trial, config) raises an **`infeasible`**
condition (§6); `run_pipeline` records the **funnel** of stage counts so a cliff is
visible. A wiring/settings/code error raises **`fatal`** and halts the run.

---

## 4. Stepping through one function at a time

Load the subsystem (the snippet in §1). Then exercise pieces directly.

**Offline (no network), using SIMULATE:**

```r
s   <- modifyList(optimizer_settings(), list(simulate = TRUE))
cfg <- seed_configs()$Prediction3            # or sample_config()
tr  <- sample_trial(s)                        # synthetic trial descriptor
evaluate_config_on_trial(cfg, tr, "CV0", s)   # -> score/status via .sim_evaluate
# optimizer internals:
con <- open_store(tempfile(fileext = ".sqlite"))
choose_config(con, s)                         # which config the loop would pick now
```

**Real path, one cached trial (network only on first touch; then cached):**

```r
s    <- modifyList(optimizer_settings(), list(simulate = FALSE))
conn <- BrAPI::createBrAPIConnection(s$brapi_host, is_breedbase = TRUE)
trial <- build_trial_descriptor("10676", conn, s)   # a specific study by id
# walk the six subtasks one at a time:
tr_ids  <- select_training_trials(seed_configs()$Prediction1, trial, conn, s)
tr_obs  <- mask_cv(get_observations(tr_ids, conn, s), trial$accessions, "CV0")
targets <- build_targets(seed_configs()$Prediction1, tr_obs, trial, conn, s)
dl      <- choose_geno_sources(seed_configs()$Prediction1, names(targets), trial$accessions, conn, s)
K       <- build_kernel(seed_configs()$Prediction1, dl)
# or the whole thing:
po <- run_pipeline(seed_configs()$Prediction1, trial, "CV0", s, conn)
score_predictions(po$pred, po$obs)
```

Tips: run any real-mode probe with stdin redirected from `/dev/null` (in a script:
`Rscript … </dev/null`) so `conn$vcf_archived()` can never block on its interactive
file-pick prompt. To see where data is lost in a single trial, use
`diagnose_trial()` (§6) -- it prints the funnel beside an independent re-derivation
of the raw counts.

---

## 5. How the code touches cache / state / logs

**`cache/`** -- the on-disk memoizer is `cached(settings, name, expr, max_age_days)`
in `data_access.R`; every entry is `cache/<name>.rds` (or a raw `.vcf`). Touch a
trial once, reuse forever.

| Cache file | Written by | Holds | Max age |
|---|---|---|---|
| `trial_catalog.rds` | `trial_catalog()` | focal-trait trials + metadata + lat/long/elev | 7 d |
| `acc_<studyid>.rds` | `get_trial_accessions()` | accession names of a trial | 30 d |
| `obs_<studyid>.rds` | `get_observations()` | all numeric observations of a study | 30 d |
| `proj_<hash>.rds` | `projects_for_accessions()` | genotyping **project** ids covering an accession set | 30 d |
| `raw_project_<id>.vcf` | `.ensure_project_vcf()` | the archived VCF (validated complete) | ∞ (re-downloaded if truncated) |
| `dosage_<id>[_thin<N>]_sz<size>_<keephash>.rds` | `get_project_dosage()` | accessions×markers dosage; key includes VCF **size** + a hash of the requested samples so a re-downloaded VCF or a different sample set never returns a stale subset | ∞ |
| `calibration_lightweight.rds` | (diagnostics, when you save it) | the last calibration table | — |

**`state/`** -- the single source of truth.

- `evals.sqlite` (`store.R`): one table `evals` with columns
  `id, config_hash, config_json, trial_id, scheme, score, n_test, status, reason,
  detail, seconds, ts`. Every evaluation -- success **and** failure -- is a row;
  the incumbent, surrogate training data, and "what to try next" are all recomputed
  from this table, so the process is killable/resumable. `status` ∈
  `ok | infeasible | suspect | error`; `reason` is the machine code; `detail` is the
  failure funnel.
- `report.md` (`report.R::write_report`): rewritten every `checkpoint_every`
  iterations -- incumbent, subtask-method importance, the failure log
  (`failure_summary`), the "⚠ Suspected bugs" section, and the running-best curve.
- `STOP`: create it (`touch state/STOP`) to halt cleanly after the current eval.

**`logs/`** -- `run_optimizer.R` writes a per-iteration heartbeat to stderr; the
launch command redirects it to `logs/run.out`. The startup canary result and any
`CANARY ALARM` / `step error` / `FATAL` messages appear here.

---

## 6. Validation & debugging tooling

### Failure classification (`R/conditions.R`)

Three deliberately distinct "didn't work" signals, so the loop knows whether to
continue or stop:

- **`infeasible(code, detail, funnel, suspect)`** -- this *one* (trial, config,
  scheme) cannot be evaluated (too little genotyped overlap, no training trials,
  too few markers, ...). Expected. `evaluate_config_on_trial` records it
  (`status = "infeasible"`, `reason = code`, `detail = funnel`) and the loop samples
  a fresh trial + config next iteration. If the funnel shows a **cliff** (lots of
  data in, ~zero out), it is tagged `suspect` instead -- a likely bug, broken out in
  the report.
- **`fatal(message, code)`** -- the whole run cannot proceed (bad settings, an
  unimplemented method). Re-raised past `evaluate.R` so `run_optimizer()` halts
  instead of spinning. (Implementation note: `evaluate_config_on_trial` uses a
  *single* `error` handler that branches on class -- re-raising a fatal from a
  multi-handler `tryCatch` would be swallowed by the sibling `error` handler.)
- **`sample_failed(detail)`** -- no usable trial could be sampled at all. Not a
  config's fault; the loop tolerates a few in a row then stops
  (`settings$max_sample_fail`).

`failure_summary()` (`report.R`) turns the stored failures into: counts by status,
dominant reasons, and the **failure rate of each subtask method** -- i.e. which
configurations break most often.

### The canary bug-oracle (`R/diagnostics.R`)

The "skip and continue" design means a **bug that hides data looks identical to a
trial that is genuinely infeasible.** The oracle catches that, via a
**calibrate-then-freeze** bootstrap against ground truth that does **not** touch our
code -- the five Predictathon teams' submission files.

```
canary_anchor(settings)            -> per-trial true accession count (median of the
                                      five algorithms' CV0_Predictions.csv, read from
                                      scripts/PredictionN/submission/). No optimizer code.
calibrate_canary_trials(s, conn,    -> our pipeline's per-trial counts joined to the anchor,
                        deep=FALSE)    with a `divergent` flag. Cheap path = genotyping_projects
                                      wizard membership; deep=TRUE also extracts VCF dosage and
                                      flags a wizard-vs-dosage gap (the synonym-name signal).
print_calibration(cal)             -> the side-by-side table; investigate any `divergent` row.
canary_configs()                   -> the FROZEN coverage configs: one per focal trial, assigned
                                      so that across the canaries every method + branch-level runs.
canary_coverage()                  -> which methods/levels the configs cover (should be all).
check_canaries(s, conn)            -> run each trial under its config, every scheme. Hard
                                      CANARY ALARM on a strong trial; soft warning on a weak one
                                      (settings$canary_weak_trials). Auto-run at startup in real mode.
```

**Why this breaks the catch-22:** the oracle shares the suspect code, so it cannot
*first-time-validate* it alone. The independent anchor does that: calibration runs
the suspect functions **and compares each count to the anchor**; a bug rarely
corrupts two independent derivations identically, so agreement is positive evidence
and a divergence localizes the bug. Freeze configs only once calibration agrees
(human-reviewed). After that, the frozen oracle's job is *regression* detection,
which it does fine despite sharing code.

### `diagnose_trial(id, settings, conn)` (`R/diagnostics.R`)

Replays one trial and prints the data funnel stage by stage **next to an
independent re-derivation of the raw counts from T3** -- including, per genotyping
project, whether the dosage matrix rownames actually intersect the trial's
accession names. A line like `overlap with accessions = 0` on a non-empty VCF is
the decisive synonym/name-mismatch signal.

### Lessons baked in (don't re-learn the hard way)

- **Cache poisoning.** A partial/raced VCF download once cached a tiny dosage
  matrix *forever* (`get_project_dosage` used `max_age_days = Inf`). Now
  `.ensure_project_vcf()` validates completeness (`.vcf_complete`) and re-downloads
  truncated files, and the dosage cache key includes the VCF byte size + a hash of
  the requested samples. If genotype counts ever look low, clear `cache/dosage_*`
  (regenerable) and any 0-sample `cache/raw_project_*.vcf`.
- **Interactive prompt.** `conn$vcf_archived()` can prompt to pick a file and hang a
  non-interactive run -- always launch real-mode work with `</dev/null`.

---

## 7. Extending the search space

To add a literature method:

1. Add the method name (and any new parameters, tagged with the methods they apply
   to) to the relevant subtask in `R/config_space.R::SUBTASKS`.
2. Add a matching dispatcher branch in `R/pipeline.R` (the `select_training_trials`
   / `build_targets` / `choose_geno_sources` / `build_kernel` / `train_model` /
   `predict_test` switch for that subtask). An unhandled method must `fatal()`.
3. Re-run `Rscript tests/run_all.R` (genome invariants + subtasks) and, after a
   search-space change, `tests/test_sim_loop.R`.
4. Re-run `canary_coverage()` so the new method is exercised by a canary (extend a
   `canary_configs()` entry if needed).

Nothing else in the optimizer needs to change -- the surrogate encoding, crossover,
mutation, store, and report all derive from `SUBTASKS`.
