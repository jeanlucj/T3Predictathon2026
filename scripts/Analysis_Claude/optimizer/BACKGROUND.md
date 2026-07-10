# Background: sources and design rationale

`README.md` says *how* to run it, `DESIGN.md` says *what* it is and where every
function lives, and `EVALUATION.md` is the evaluation/validation runbook. This
document records *why* it is built the way it is — the sources consulted, the design
decisions (including alternatives considered and rejected), and the statistical and
data-management challenges it must handle and how each is addressed. It is meant as
a permanent reference so the reasoning is not lost.

---

## 1. What constrained the design

Three properties of this problem drove almost every decision:

- **The objective is expensive.** One evaluation pulls phenotypes and genotypes
  from T3, builds a relationship matrix, and fits a mixed model. Evaluations cost
  seconds-to-minutes, not microseconds, so the search must be *sample-efficient*:
  it cannot afford grid search or a large evolutionary population per generation.
- **The objective is noisy.** A configuration's quality is its expected
  predictive ability across trials, but any one trial is a high-variance sample of
  that expectation. The optimizer must tolerate noise and not chase a single lucky
  trial.
- **The search space is mixed and structured.** Each subtask is a categorical
  method choice with its own continuous/integer parameters, and many parameters
  are only meaningful for one method. The optimizer must handle categoricals and
  conditional parameters natively.

The user's framing — *recombine subtask methods from different algorithms* — adds
a fourth: the representation should make "take this algorithm's training-trial
selection and that algorithm's kernel" a first-class operation.

---

## 2. Sources consulted

### 2a. Optimization / AutoML

- **Hutter, Hoos & Leyton-Brown (2011), "Sequential Model-Based Optimization for
  General Algorithm Configuration" (SMAC), LION-5.** The foundational paper for
  this design. It introduced using a **random forest** as the surrogate model in
  Bayesian optimization (instead of a Gaussian process), precisely because random
  forests handle categorical and conditional hyperparameters and scale to many
  evaluations. Our surrogate (`R/surrogate.R`) is a deliberately transparent
  re-implementation of this idea.

- **SMAC3 — the maintained reference implementation.**
  <https://github.com/automl/SMAC3>. Confirms the modern, widely-used form of the
  method: RF surrogate + Expected Improvement, with mechanisms for racing
  configurations under noise. We borrowed the architecture, not the package.

- **AutoML.org HPO overview.**
  <https://www.automl.org/hpo-overview/hpo-tools/smac/>. A concise survey placing
  SMAC against the other two dominant Bayesian-optimization approaches —
  Gaussian-process + Expected Improvement, and the Tree-structured Parzen
  Estimator (TPE). Useful for the comparison in §3.

- **"An improved hyperparameter optimization framework for AutoML systems using
  evolutionary algorithms," Scientific Reports (2023).**
  <https://www.nature.com/articles/s41598-023-32027-3>. Representative of the
  evolutionary-algorithm branch of AutoML, and of hybrids (e.g. DEHB combines
  Differential Evolution with Hyperband; GA-PARSIMONY combines a GA with Bayesian
  optimization). It is why we fold evolutionary crossover/mutation into a
  model-based loop rather than choosing one camp.

### 2b. Genomic prediction (the methods being recombined and the metric)

These are the established references behind the subtasks themselves; several are
already cited in the repo's own documentation.

- **Jarquín et al. (2017), *The Plant Genome* 10(2).** Defines the CV0 / CV00
  (and CV1/CV2) cross-validation scenarios for predicting genotypes in untested
  environments. Our evaluation scores exactly these scenarios, matching the
  Predictathon specification.
- **VanRaden (2008), *J. Dairy Sci.* 91:4414–4423.** The genomic relationship
  matrix `G = XX'/p` used in the `vanRaden_single` kernel.
- **Endelman (2011), *The Plant Genome* 4:250–255.** The `rrBLUP` package /
  mixed-model GBLUP that is the reliable backbone of the model subtask.
- **CovCombR (Akdemir).** The Wishart-EM algorithm for combining partial
  relationship matrices estimated from different marker sets — the basis of the
  `em_combine` kernel, used here through `T3BrapiHelpers::covariance_combiner`.
- **Predictathon scoring** (`T3_predictathon_scripts/analysis/30_...Rmd`). The
  metric is the Pearson correlation between predicted and observed per-accession
  BLUEs on a "core" accession set, with a task counted as failed if too few
  accessions are predicted or predictions are constant. `score_predictions()`
  reproduces this.

---

## 3. Design decisions and rationale

### The objective: expected accuracy over *random* trials, not the nine focal trials
Optimizing on the nine Predictathon trials would overfit to nine environments and
would not satisfy the goal of predicting *any* trial. Defining the objective as
the mean score over freshly sampled trials makes generalization the thing being
optimized: a configuration only looks good if it works on trials it has never
seen. The cost is higher variance per evaluation, which the rest of the design
absorbs.

### The target domain (tailoring what "any trial" means)
The set of trials the objective averages over is itself a design choice. By
default it is "all T3 trials measuring the focal trait," but `settings$target_domain`
lets the user restrict random focal-trial sampling to a list of breeding programs,
years, locations, and/or specific trial names (studyName). This is what lets the same machinery produce a pipeline
tuned for one program rather than a global all-rounder, since the best pipeline
genuinely differs across populations. Term: we call this the **target domain**;
in the algorithm-configuration literature the analogous concept is the **instance
set** (the set of problem instances over which average performance is optimized).
Only the *focal* trial is constrained -- the pipeline may still select training
trials from anywhere, because that selection is part of what is being optimized.
The focal *trait* is likewise configurable (`settings$focal_trait`), so the system
is not specific to yield.

### Optimizer family: RF-surrogate Bayesian optimization (SMAC), not the alternatives
| Candidate | Why not (here) |
|-----------|----------------|
| **Grid search** | Wastes the expensive evaluations; infeasible over a mixed conditional space. |
| **Pure random search** | A reasonable baseline (and we use it for warm-up), but it never learns from past evaluations. |
| **Gaussian-process Bayesian optimization** | The textbook sample-efficient choice, but GPs need a kernel over the inputs; categorical + conditional parameters are awkward, and GPs scale poorly (≈O(n³)) as evaluations accumulate over a long background run. |
| **TPE** | Handles conditionals well, but models densities of good/bad configs rather than giving the calibrated mean+uncertainty we want for Expected Improvement, and is less transparent to re-implement. |
| **Pure genetic algorithm** | Matches the "recombine" framing, but under a noisy, expensive objective it needs a large population per generation and discards the information in every individual evaluation. |
| **RF-surrogate BO (chosen)** | Learns from every evaluation; handles categoricals/conditionals/NA natively; scales linearly in evaluations; robust to noise. This is SMAC. |

We then **add evolutionary operators on top** rather than choosing one camp,
because (a) crossover *is* the user's "recombine algorithms" operation, and (b)
crossover/mutation of elites is an effective candidate generator for the
acquisition step in a combinatorial space, supplying structured proposals the
surrogate then ranks. This hybrid is in the spirit of DEHB / GA-PARSIMONY.

### A self-contained surrogate instead of a package
The R environment had **no** Bayesian-optimization, GA, or random-forest package
installed (`mlrMBO`, `GA`, `ParBayesianOptimization`, `ranger`, `randomForest`
were all absent), and the explicit requirement was code the user can read and
understand. Rather than add a heavy dependency that must compile and that a
background process would silently depend on, the surrogate is a bag of `rpart`
regression trees (`rpart` ships with base R): each tree sees a bootstrap sample of
the rows and a random subset of the features (row + feature bagging = a random
forest); the ensemble mean is the prediction and the spread across trees is the
uncertainty. This is ~40 transparent lines and is exactly the surrogate SMAC
prescribes. Trees also handle the conditional parameters for free — an
inapplicable parameter is `NA`, and `rpart` routes `NA` through surrogate splits
rather than requiring imputation.

### Genome representation: a flat list of per-subtask blocks
A configuration is a flat named list where each subtask contributes a
`<subtask>.method` key plus its parameter keys. This makes the three operations
the search needs trivial and readable: **crossover** swaps whole subtask blocks
between two parents (so a child literally inherits P1's training-trial selection
and P5's kernel); **mutation** resamples one subtask's block; **encoding** maps
each key to one surrogate feature column. `repair_config()` enforces the
invariant that a parameter is present exactly when it applies to the chosen
method — necessary because crossover can pair a method from one parent with a
parameter block from the other.

### Seeding with the five submissions
The submissions are encoded as starting configurations (`R/seeds.R`) and run
first. This gives (a) an explicit, like-for-like baseline the search must beat —
reported as "best seed" — and (b) good building blocks for crossover from
iteration one, so the search starts near known-decent pipelines instead of from
noise.

### Denoising: aggregate to a per-configuration mean before fitting the surrogate
Because each trial is a noisy sample, the surrogate is trained on the **mean score
per distinct configuration** across all the trials and schemes it has been run on,
not on raw `(config, trial)` points. This denoises the signal the model learns
from. A more elaborate alternative — modelling `(config, trial)` jointly with
trial descriptors as extra features, so the surrogate learns *which configs suit
which trial types* — is deliberately deferred (see §4); the per-config mean is
simpler, and the optimum of the mean is precisely the "best pipeline for a random
trial" we want.

### Acquisition and keeping exploration alive
Candidates are ranked by **Expected Improvement** over the current best mean
score, which balances "predicted to be good" against "uncertain, worth probing."
Two safeguards prevent premature convergence: the between-tree uncertainty is
**floored** (`0.25 × sd(observed scores)`) so a config the trees happen to agree
on is never treated as known with certainty, and a fraction of the candidate pool
is always **fresh random** draws, not just elite recombinations.

### Re-evaluating the incumbent (the noise tax)
With a noisy objective, the apparent best config may just have had lucky trials.
Two mechanisms guard against crowning a fluke: an evaluation is occasionally spent
**re-running the current incumbent on a new trial** (`reeval_prob`) to tighten its
mean, and the reported **incumbent must have at least `incumbent_min_reps`
successful trials** before it can hold the title. This trades a little speed for a
trustworthy leader.

### Persistence and resumability: one SQLite table as the single source of truth
A background job that runs for days must survive being killed. Every evaluation —
including failures, with their reason — is a row in `state/evals.sqlite`, and the
optimizer's entire state (incumbent, surrogate training data, what to try next) is
recomputed from that table on startup. There is no separate checkpoint that could
drift out of sync with the results.

### Caching: pay for a trial once
Phenotype pulls, accession lists, genotyping-project lookups, and VCF→dosage
extractions are memoized to disk (`cache/`, plus the existing per-project dosage
and GRM caches). Touching a given trial is therefore expensive once and cheap
forever after — essential when the search revisits trials and reuses the same
trial across CV0 and CV00 within one iteration.

### Verify the machinery before spending real compute (SIMULATE mode)
A bug in a multi-day background job is costly. `R/evaluate.R` provides a synthetic
objective with a *known* generative structure (some subtask choices genuinely
better, some better only on certain trial types, plus realistic per-trial noise).
`tests/test_sim_loop.R` runs the entire loop against it offline in seconds and
asserts the discovered incumbent beats the seeds and that the surrogate phase
improves on random search. This is what lets us trust the optimizer before
`simulate = FALSE` is ever set. It also doubles as a regression test for any
future change to the search space or optimizer.

### Scoring and failure handling
`score_predictions()` reproduces the Predictathon metric (Pearson correlation of
predicted vs observed BLUEs) and its failure rules (fewer than 5 overlapping
accessions, or constant predictions → `NA` with a recorded reason). Failures are
*stored*, not discarded, so the surrogate can learn to avoid configurations — or
config/trial-type combinations — that reliably break, rather than blindly
re-proposing them. Failures are classified by a small condition system
(`R/conditions.R`): an **`infeasible`** signal means this one (trial, config,
scheme) cannot be evaluated (expected — recorded and the loop continues); a
**`fatal`** signal means the run cannot proceed (bad settings, an unimplemented
method — the loop halts); a **`sample_failed`** signal means no usable trial could
be drawn. The stored `status`/`reason`/`detail` columns feed `failure_summary()`,
which reports the failure rate of each subtask method — the data needed to see
which (often more demanding) configurations are fragile. See "Telling genuine
infeasibility from a data-hiding bug" below for why this classification matters.

### Reusing proven code from this repo
The genotype path (VCF download, dosage extraction, EM combine) reuses the
patterns already validated in `scripts/Analysis_Claude/build_grm_for_cv00.R`,
including its **dosage-transposition bugfix** (building the dosage matrix in the
VCF's native marker × sample orientation and transposing afterward). This avoids
re-introducing a known, subtle correctness bug.

---

## 3b. The data-management challenge: the T3 layer

The statistics above assume clean inputs. Getting clean inputs out of T3/Wheat over
BrAPI is itself a substantial part of the work — the responses are inconsistent and
several plausible-looking calls return wrong-but-non-erroring results. Each issue
below was hit during bring-up; the fix is in `R/data_access.R` unless noted.

- **Trait → trial filtering needs the numeric variable id, not the name.** Passing
  the human trait string (`"Grain yield - kg/ha|CO_321:0001218"`) to the breeder
  wizard's `traits` filter *silently returns every trial* with a "no matching trait"
  warning — the wizard's trait dimension only lists derived `COMP:` traits. We
  filter the studies search by `observationVariableDbIds` instead
  (`settings$focal_trait_db_id`, grain yield = `84527`), which narrows ~9030 wheat
  trials to ~7561 grain-yield ones. `focal_trait` (the name) is still used for
  matching within an observation set.
- **Trial metadata has no coordinates / no year.** `get_all_trial_meta_data`
  returns snake_case columns (janitor) but no `latitude`/`longitude`/`year`. We
  derive `year` from `start_date` and join `latitude`/`longitude`/`elevation` per
  location via `T3BrapiHelpers::get_lat_long_elev_from_location_vec`, cached with the
  catalogue — which is what enables the environmental-similarity and G×E paths.
- **janitor snake_case columns.** All `T3BrapiHelpers` outputs are
  `janitor::clean_names()`-ed (`study_db_id`, `program_name`, …), while raw
  `conn$search`/`conn$wizard` responses keep camelCase. The two are not
  interchangeable; using the wrong case throws "object not found" inside
  `dplyr::filter`.
- **Genotyping protocol id ≠ project id.** `conn$vcf_archived()` needs a downloadable
  genotyping **project** id; `get_geno_protocol_from_germ_vec` returns **protocol**
  ids. We get project ids from the `genotyping_projects` wizard
  (`projects_for_accessions()`), batching accessions.
- **`/observations` shape.** The records live in `resp$combined_data` (the
  auto-paginated list), not `resp$data` (search metadata); each record is a nested
  list, and rep/block are not returned (they would require an `observationUnits`
  query), so BLUE estimation degrades gracefully to per-germplasm means.
- **VCF sample synonyms.** Some archived VCFs name samples by old/preliminary
  synonyms, so a focal accession can be genotyped yet not match by name — a small
  residual tail after the fixes above (e.g. ~22 of 139 on one trial). This is the
  one item only partially handled; the canary calibration measures its size per
  trial.
- **Cache poisoning from partial downloads.** `get_project_dosage` caches dosage
  with `max_age_days = Inf`. A truncated/raced VCF download once produced a tiny
  dosage that was then cached *forever*, silently shrinking every GRM built from it
  (observed: one project's focal overlap collapsed 115 → 6). Fix: `.ensure_project_vcf()`
  validates the VCF is complete (`.vcf_complete`) and re-downloads a truncated one;
  the dosage cache key now includes the VCF byte size **and** a hash of the
  requested samples, so a re-downloaded VCF or a different sample set can never
  return a stale subset.
- **`vcf_archived` can prompt.** When a project has multiple archived files,
  `conn$vcf_archived()` prompts interactively and hangs a non-interactive run; always
  redirect stdin from `/dev/null`.

## 3c. Telling genuine infeasibility from a data-hiding bug

The optimizer must *skip and continue* past trials a configuration genuinely cannot
predict — but that creates a dangerous failure mode: a **bug that hides real data
produces the same "infeasible" verdict** as a genuinely unpredictable trial, so a
broken build can run for hours looking like it is working. Three defences:

- **Funnel-cliff detection.** Each `infeasible` records the data cardinalities at
  every pipeline stage. Genuine infeasibility shows a plausible, gradual funnel; a
  bug shows a *cliff* — large inputs collapsing to ~zero at a join/match step (the
  signature of a name/key mismatch). A cliff is tagged `suspect`, not `infeasible`,
  and broken out in the report.
- **The canary bug-oracle.** Nine known-feasible trials (the Predictathon focal
  trials) are each run under a frozen config; together the configs exercise every
  subtask method and branch-level (`canary_coverage()`), so a bug in *any* branch
  trips the oracle. It runs at startup; a `CANARY ALARM` means the build is hiding
  data. The catch — *an oracle built from the suspect code cannot first-time-validate
  that code* — is broken with an **independent anchor**: `canary_anchor()` derives
  each trial's true accession count from the five teams' submission files (which
  never touched our code), and calibration runs our functions **and compares to the
  anchor**. A bug rarely corrupts two independent derivations identically, so
  agreement is positive evidence and a divergence localizes the bug. Configs are
  frozen only once calibration agrees; thereafter the oracle does *regression*
  detection, which it can do despite sharing code.
- **`diagnose_trial()`.** Drills into one trial, printing the funnel beside an
  independent re-derivation of the raw counts (per project, whether VCF sample names
  intersect the accession names) — turning a flagged trial into a localized cause in
  seconds.

---

## 4. Known limitations and future directions

- **Trial-aware surrogate.** Adding trial descriptors (heritability proxy, number
  of genotyping projects, environmental atypicality, training-set size) as
  surrogate features would let the model learn which pipeline suits which trial
  type, and could support *adaptive* pipelines that branch on trial properties.
  The current design optimizes a single fixed pipeline.
- **Smarter budget allocation.** Successive-halving / Hyperband (or DEHB) could
  spend fewer full evaluations on clearly-bad configs by testing them on cheaper
  proxies (fewer markers, fewer training trials) first. The current loop gives
  every evaluation the full budget.
- **Parallel evaluation.** `future`/`furrr` are available; several trials/configs
  could be evaluated concurrently, with the surrogate refit between batches. The
  current loop is sequential for simplicity and reproducibility.
- **Multi-objective.** CV0 and CV00 are scored separately but optimized as one
  mean. If the two scenarios trade off, a Pareto or scalarized multi-objective
  treatment may be preferable.
- **Real run partially validated.** The real-data path has now been exercised
  against the live BrAPI server: the project-membership layer agrees with the
  independent submission-file anchor on all nine canary trials, the dosage
  cache-poisoning bug is fixed, and `check_canaries()` is green (all nine trials
  predict end-to-end under their coverage configs). Still open: honest cross-trial
  *accuracy* has not been measured (the canary scores test feasibility, not tuned
  accuracy); the VCF synonym tail is only partially handled; and some methods'
  feasibility on the broader (non-focal) trial population is unconfirmed.

---

## 5. References

1. Hutter F., Hoos H.H., Leyton-Brown K. (2011). *Sequential Model-Based
   Optimization for General Algorithm Configuration.* LION-5. (SMAC.)
2. SMAC3 reference implementation. <https://github.com/automl/SMAC3>
3. AutoML.org, *HPO overview — SMAC.*
   <https://www.automl.org/hpo-overview/hpo-tools/smac/>
4. *An improved hyperparameter optimization framework for AutoML systems using
   evolutionary algorithms.* Scientific Reports (2023).
   <https://www.nature.com/articles/s41598-023-32027-3>
5. Jarquín D. et al. (2017). *Increasing genomic-enabled prediction accuracy by
   modeling genotype × environment interactions in Kansas wheat.* The Plant
   Genome 10(2). doi:10.3835/plantgenome2016.12.0130
6. VanRaden P.M. (2008). *Efficient methods to compute genomic predictions.*
   J. Dairy Sci. 91:4414–4423.
7. Endelman J.B. (2011). *Ridge regression and other kernels for genomic
   selection with R package rrBLUP.* The Plant Genome 4:250–255.
8. Akdemir D. *CovCombR* — Wishart-EM combination of partial covariance/relationship
   matrices (used via `T3BrapiHelpers::covariance_combiner`).
