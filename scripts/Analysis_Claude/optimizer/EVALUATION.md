# Evaluation guide — stepping through the optimizer

This is the document for **evaluating and validating** the optimizer from the
RStudio console -- for me now and in six months. `README.md` is for running it;
`DESIGN.md` is the static architecture (every function, where outputs go);
`BACKGROUND.md` is why it is built this way and the challenges it must handle.
This file is the **runbook**: it lists the modules, gives a fast→slow order to
check them in, and walks each one step by step -- one console line at a time --
saying what each step should return, what to eyeball, and the subtle failure
signatures to hunt for.

The system is genuinely complex -- the task and the messy T3 data layer make that
unavoidable -- so the goal is that **working through this document is faster and
more thorough than re-reading all the code**, and that it maximizes the chance of
catching the bugs that make the optimizer *look* like it works without delivering
on its objective.

The objective, to keep in view the whole time: discover one pipeline
configuration whose **mean Pearson correlation between predicted and observed
per-accession BLUEs, over randomly sampled T3 grain-yield trials**, beats every
one of the five submitted algorithms. A bug that silently hides data or leaks the
answer will still produce a plausible-looking number -- that is what we are
guarding against.

---

## 1. The modules, in one screen

Eleven `R/` modules, one entry point, one settings file. The **group** column is
the name you pass to `arm_evaluation()` (§2) to `debug()` that module.

| Module (`R/…`) | Owns | Off/online | Group |
|---|---|---|---|
| `config_space.R` | the six-subtask genome: sample / encode / crossover / mutate | offline | `genome` |
| `seeds.R` | the five submissions as starting configs | offline | `genome` |
| `optimizer.R` | candidate generation + acquisition + phase logic + incumbent | offline | `engine` |
| `surrogate.R` | bagged-rpart RF surrogate: mean, uncertainty, Expected Improvement | offline | `engine` |
| `store.R` | SQLite results store (resumable state) | offline | `store` |
| `report.R` | learning curve, incumbent, method importance, failure log | offline | `store` |
| `evaluate.R` | run a config on a trial → score; the SIMULATE objective | offline\* | `scoring` |
| `conditions.R` | `infeasible` / `fatal` / `sample_failed` signals + funnel | offline | — |
| `data_access.R` | BrAPI + caching: trials, phenotypes, accessions, projects, dosages | **online** | `data` |
| `pipeline.R` | the parameterized six-subtask pipeline + dispatchers + `mask_cv` | **online** | `subtaskA…F`, `flow` |
| `diagnostics.R` | canary oracle (anchor / calibrate / check / coverage) + `diagnose_trial` | **online** | `diagnostics` |

\* `evaluate.R` is offline in `simulate = TRUE`; in real mode it calls
`run_pipeline()`, which is online.

The six subtasks inside `pipeline.R` (the recombinable "genome"): **A** select
training trials, **B** preprocess phenotypes → per-accession targets, **C** select
genotyping data, **D** build relationship/kernel, **E** train model, **F** predict.
Each is its own group (`subtaskA` … `subtaskF`) so you can arm exactly one.

---

## 2. How to evaluate: the plan

Three principles, in priority order:

1. **Offline before online.** Prove the math with `simulate = TRUE` (and the test
   suite, §7) before spending a single download. Levels L1–L4 are offline.
2. **Fast before slow.** Follow `EVAL_ORDER` (printed by `eval_groups()`); it runs
   the instant offline modules first and the network/heavy ones last, so a cheap
   bug is caught before an expensive one. Within the live layer, pick the
   *smallest* work first (§5).
3. **One group armed at a time.** `arm_evaluation("subtaskB")`, run the pipeline,
   step through *only* subtask B, `disarm_evaluation()`. Arming everything at once
   makes the debugger unusable.

### The three tools

| Tool (in `R/evaluation.R`) | Use it to | When |
|---|---|---|
| `arm_evaluation(group)` / `disarm_evaluation()` | drop into `debug()` for one module and watch the Global Environment change line by line | stepping through logic |
| `peek(x)` | print a one-line health summary of the object flowing between subtasks (shape, NA counts, degeneracy, rowname overlap) and pass it through | inspecting an intermediate |
| `diagnose_trial(id, …)` / `canary_anchor()` / `calibrate_canary_trials()` | compare the pipeline's own counts to an **independent** re-derivation / ground truth | confirming a real result is real |

`eval_groups()` prints the whole menu. The debugger keys, once: **`n`** next line,
**`s`** step into a call, **`c`** finish this function (run to its return),
**`Q`** quit the debugger. **Always `disarm_evaluation()` before a background run**
— otherwise every pipeline call blocks on the debugger.

`peek()` recognizes the five shapes that carry the pipeline's state and flags the
signatures that otherwise pass silently:
- **named numeric** (`targets`, `obs`, `pred`) → length, distinct count, NA count,
  range; warns if all finite values are identical.
- **matrix** (a dosage matrix, `K`) → dims, square?, symmetric?, diagonal range;
  pass `accessions=` to get **rowname overlap** (the synonym-mismatch signal).
- **list of matrices** (`dosage_list`) → per-project dims and overlap.
- **tibble** (`train_obs`) → rows×cols and **NA count per column** (flags `ALL NA`).

---

## 3. Bootstrap (paste once per session)

```r
setwd("scripts/Analysis_Claude/optimizer")   # if not already the working dir
library(tidyverse)
source("settings.R")
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
source("run_optimizer.R")                     # defines run_optimizer()/optimizer_step()
                                              # (its auto-run is guarded, so sourcing is safe)
s <- optimizer_settings()
eval_groups()                                 # confirm the tooling loaded
```

For the **online** levels (L5 on) also open the live connection once and reuse it:

```r
s    <- modifyList(optimizer_settings(), list(simulate = FALSE))
conn <- BrAPI::createBrAPIConnection(s$brapi_host, is_breedbase = TRUE)
```

> **Caveat.** `conn$vcf_archived()` can prompt interactively to pick a file and
> will hang a non-interactive run. When you script a real-mode probe with
> `Rscript`, redirect stdin: `Rscript probe.R </dev/null`. Interactively in
> RStudio you will simply see the prompt.

Durable state lives in two dirs and is regenerable: `state/` (the SQLite store +
`report.md`; `touch state/STOP` halts a run) and `cache/` (downloaded data). See §8.

---

## 4. Level-by-level walkthrough

Each level: **arm** the group, run the numbered console lines (each assigns a
Global-Environment variable), and check the three annotations —
**→ returns** (shape you should see), **🔍 eyeball** (the `peek()`/inspection and
what healthy looks like), **🚩 red flag** (the subtle failure to hunt), plus the
italic *Assumption:* the code is making at that step (test **it**, not just the
output). `disarm_evaluation()` when you finish a level.

### L1 — `genome` (offline, instant)

```r
arm_evaluation("genome")
cfg  <- sample_config()                        # step with n/c to watch each block fill
cfgs <- seed_configs()                         # the 5 submissions as configs
kid  <- crossover(cfgs$Prediction1, cfgs$Prediction5)
mut  <- mutate_config(cfg, n_subtasks = 1)
disarm_evaluation()
```

- **→ returns** `cfg` is a named list with one `*.method` per subtask plus its
  parameters; `cfgs` has 5 named entries; `kid`/`mut` are configs.
- **🔍 eyeball** `config_hash(cfg)` is stable across repeated calls;
  `str(cfg)` shows every subtask present; each block in `kid` matches *one* parent
  (`kid$kernel.method %in% c(cfgs$Prediction1$kernel.method, cfgs$Prediction5$kernel.method)`).
- **🚩 red flag** a parameter that does not apply to its method is non-`NA`
  (should be `NA`), or a required applicable parameter is `NA`; `mut` identical to
  `cfg` (mutation must change ≥1 block).
- *Assumption:* `config_space.R::SUBTASKS` is the single source of truth — every
  method named here must have a dispatcher branch in `pipeline.R` (§10).

### L2 — `engine` (offline)

Build a tiny store of synthetic evaluations, then watch the decision engine.

```r
s   <- modifyList(optimizer_settings(), list(simulate = TRUE))
con <- open_store(tempfile(fileext = ".sqlite"))
for (i in 1:30) { tr <- sample_trial(s); cfg <- sample_config()
  for (sc in s$schemes) { ev <- evaluate_config_on_trial(cfg, tr, sc, s)
    store_eval(con, cfg, tr$id, sc, ev$score, ev$n_test, ev$status, ev$reason, ev$detail %||% NA, ev$seconds) } }
arm_evaluation("engine")
ch <- choose_config(con, s)                    # step: which phase? seed/random/acquisition
disarm_evaluation()
```

- **→ returns** `ch` is `list(cfg, source, ei)`; `source` is one of
  `seed` / `random_init` / `reeval` / `acquisition`.
- **🔍 eyeball** with few evals `source` is `seed` then `random_init`; only past
  `n_random_init` does it become `acquisition`. Inside acquisition,
  `expected_improvement()` is always ≥ 0 and rises with predicted mean.
- **🚩 red flag** the surrogate predicts a flat score for every candidate (no
  signal learned); acquisition picks the same config forever (EI collapsed);
  `incumbent_config()` returns a config with fewer than `incumbent_min_reps` reps.
- *Assumption:* the surrogate trains on the **per-config mean** score
  (`aggregate_scores`), i.e. the best config for a *random* trial is the one with
  the best average — not tuned to any single trial.

### L3 — `store` (offline)

```r
arm_evaluation("store")
ev  <- read_evals(con)                          # the whole table as a tibble
write_report(con, s)                            # writes state/report.md
disarm_evaluation()
```

- **→ returns** `ev` has columns `config_hash, config_json, trial_id, scheme,
  score, n_test, status, reason, detail, seconds, ts`; failures are rows too
  (with `NA` score).
- **🔍 eyeball** round-trip a config: `config_from_json(config_to_json(cfg))`
  equals `cfg` **including `NA` params**; `report.md` renders incumbent + method
  importance + failure log.
- **🚩 red flag** `NA` scores dropped instead of stored (the failure log is how we
  see data-hiding bugs); JSON round-trip changes a parameter's type.
- *Assumption:* the entire optimizer state is reconstructable from this one table,
  so a killed run resumes with no loss.

### L4 — `scoring` (offline; the CV0/CV00 crux)

```r
arm_evaluation("scoring")
x  <- setNames(rnorm(20), paste0("g", 1:20))
score_predictions(x, 2 * x + 1)                 # must be +1 (cor of x with a*x+b, a>0)
tr <- sample_trial(s); cfg <- sample_config()
evaluate_config_on_trial(cfg, tr, "CV00", s)    # step through the sim objective
disarm_evaluation()
```

- **→ returns** `score_predictions` returns a scalar in [-1, 1] or `NA` (with a
  reason) when <5 accessions overlap or predictions are constant.
- **🔍 eyeball** `cor(x, x) == 1`; only shared accession *names* join; `mask_cv`
  under `"CV00"` drops **exactly** the focal-trial accessions that `"CV0"` keeps
  (`setdiff(names(cv0_targets), names(cv00_targets))` = the focal rows).
- **🚩 red flag** a suspiciously **high** score (leakage — the focal answer reached
  the training side); **CV0 and CV00 scoring identically** (the CV00 mask is a
  no-op — the whole challenge collapses).
- *Assumption:* this is the exact Predictathon metric
  (`T3_predictathon_scripts/analysis/30_All_accuracies…Rmd`).

> **Everything below is online.** Switch to the real connection (§3) and follow §5:
> smallest / cheapest work first.

### L5 — `data` (live BrAPI; start on the CACHED trial)

Trial `10676` is already in `cache/`, so this touches **no** network on the
happy path — the cheapest possible first live check.

```r
arm_evaluation("data")
trial <- build_trial_descriptor("10676", conn, s)   # step: watch the descriptor fill
obs   <- get_observations("10676", conn, s)
acc   <- get_trial_accessions("10676", conn, s)
disarm_evaluation()
peek(obs); peek(trial$accessions)
```

- **→ returns** `trial` is a list with `id, accessions, n_acc, program, location,
  year, lat/long/elev`; `obs` is a tidy tibble
  (`study_id, germplasm_name, value, unit_id, rep, block, col, row`).
- **🔍 eyeball** `peek(obs)` — `value` finite, and **`rep`/`block` NOT all `NA`**
  (they can be legitimately absent, but all-`NA` means BLUE silently degrades to
  germplasm means downstream — the canonical subtle bug); `trial$lat` finite if
  the location has coordinates.
- **🚩 red flag** `n_acc` far below the T3 web UI's count for the trial; a
  camelCase-vs-snake_case column surprise (raw `conn$search` responses are
  camelCase, `T3BrapiHelpers` outputs are `janitor::clean_names()`).
- *Assumption:* the `obs_<sid>` cache is trait-independent — the focal trait is
  filtered at *read* time, so changing `focal_trait` reuses it.

### L6 — `subtaskA` select training trials (start with FEW trials)

```r
arm_evaluation("subtaskA")
# HIGH overlap threshold first -> few neighbours -> few downloads later:
cfgA <- modifyList(seed_configs()$Prediction1,
                   list(train_select.method = "accession_overlap",
                        train_select.primary_only = "yes",
                        train_select.primary_min = 20))
tr_ids <- select_training_trials(cfgA, trial, conn, s)
disarm_evaluation()
```

- **→ returns** `tr_ids` a character vector of study ids (the focal id removed
  downstream in `run_pipeline`).
- **🔍 eyeball** a *high* `primary_min` returns a short list; lowering it (or
  turning off `primary_only`) lengthens it monotonically. `top_k_similar` returns
  exactly `k`; `same_program` respects `prog_cap`.
- **🚩 red flag** an empty list where the trial obviously has neighbours (a
  `.find_related` / id-space bug), or a list that ignores the threshold.
- *Assumption:* `.find_related` reads `T3BrapiHelpers`' tabyl columns
  (`other_study_db_id`, `n`); a shape change there silently empties the result.

### L7 — `subtaskB` preprocess phenotypes → targets

```r
arm_evaluation("subtaskB")
train_obs <- mask_cv(get_observations(tr_ids, conn, s), trial$accessions, "CV0")
targets   <- build_targets(cfgA, train_obs, trial, conn, s)
disarm_evaluation()
peek(train_obs); peek(targets)
```

- **→ returns** `targets` a named numeric, one value per training accession.
- **🔍 eyeball** `peek(train_obs)` — no column is `ALL NA`; `peek(targets)` —
  `distinct` ≫ 1 and a plausible yield range.
- **🚩 red flag** **targets all identical / degenerate** (`build_targets`
  collapsed — e.g. every trial centered to the same mean); `rep`/`block` all `NA`
  so `.blue_per_trial` fell back to plain means for every trial without your
  noticing.
- *Assumption:* per-trial BLUE `y ~ germ (+rep +block)` only adds `rep`/`block` to
  the model **when they vary** — so their absence is tolerated, not detected.

### L8 — `subtaskC` select genotyping data (the synonym check)

```r
arm_evaluation("subtaskC")
dl <- choose_geno_sources(cfgA, names(targets), trial$accessions, conn, s)
disarm_evaluation()
peek(dl, accessions = c(names(targets), trial$accessions))
```

- **→ returns** `dl` a named list of accessions×markers dosage matrices, with an
  `n_projects` attribute.
- **🔍 eyeball** `peek(dl, accessions=…)` — each matrix's **rowname overlap** with
  the accession set is > 0.
- **🚩 red flag** a **non-empty dosage matrix with overlap = 0** — the decisive
  synonym / name-mismatch signature (VCF samples under old preliminary names);
  `n_projects > 0` but `length(dl) == 0` (found projects, extracted nothing → a
  download/parse bug, flagged `suspect` by `run_pipeline`). Cross-check with
  `diagnose_trial("10676", s, conn)` which prints "overlap with accessions" from an
  independent re-derivation.
- *Assumption:* a genotyping **project** id (not protocol id) is what
  `vcf_archived()` accepts; the dosage cache is keyed by project + thin + VCF byte
  size, so a re-downloaded VCF of a different size invalidates a stale dosage.

### L9 — `subtaskD` build the kernel

```r
arm_evaluation("subtaskD")
K <- build_kernel(cfgA, dl)
disarm_evaluation()
peek(K)
```

- **→ returns** `K` a square relationship/kernel matrix over the genotyped
  accessions.
- **🔍 eyeball** `peek(K)` — **square, symmetric**, diagonal in a sane range
  (~1 for VanRaden after standardization; exactly 1 on an RKHS diagonal). Clones
  (identical dosage rows) give identical `K` rows.
- **🚩 red flag** `symmetric = FALSE`; a diagonal far from expectation; `em_combine`
  returning a matrix one row/col too big (the phantom-variance row must be
  dropped); fewer than 50 markers surviving QC (`.qc_markers` should have raised
  `too_few_markers`).
- *Assumption:* `ridge` is added to the diagonal *after* the base kernel; markers
  are the intersection across projects (with a largest-project fallback below 50
  shared markers).

### L10 — `subtaskE` + `subtaskF` train and predict

```r
arm_evaluation(c("subtaskE", "subtaskF"))
geno <- rownames(K); tin <- intersect(names(targets), geno); ten <- intersect(trial$accessions, geno)
fit  <- train_model(cfgA, targets[tin], K, tin, ten, trial, train_obs)
pred <- predict_test(cfgA, fit, K, tin, ten, targets, "CV0", s)
disarm_evaluation()
peek(pred)
```

- **→ returns** `fit` a model object (`kind = "gblup"`, a `u` random effect, a
  `mu`); `pred` a named numeric over the test accessions.
- **🔍 eyeball** `peek(pred)` — **not constant**, plausible range; under the
  fallback (too little overlap) it is deliberately `mean(targets)` for everyone.
- **🚩 red flag** predictions constant when overlap was adequate (the model
  collapsed, not the fallback); blending changing predictions under **CV00** (it
  must be a no-op there — those phenotypes were masked).
- *Assumption:* `predict_test`'s `cond_expectation` computes `mu + K21 K11⁻¹ u`;
  the observed-BLUE blend applies **only** under CV0.

### L11 — `flow` end to end (watch the funnel)

```r
arm_evaluation("flow")               # run_pipeline + optimizer_step only (no double-break)
po <- run_pipeline(cfgA, trial, "CV0", s, conn)   # c through it, watch stage counts
disarm_evaluation()
score_predictions(po$pred, po$obs)
```

- **→ returns** `po = list(pred, obs)`; the score is a plausible correlation
  (roughly in the range the Predictathon reported, not ~1 and not `NA`).
- **🔍 eyeball** the stage counts stay sizeable end to end — no **cliff** (many
  accessions in, ~zero out) between the phenotyped and genotyped sides.
- **🚩 red flag** a cliff `run_pipeline` tags `suspect`; a score too high (leakage)
  or always `NA` (silent infeasibility).
- *Assumption:* an unevaluable (trial, config) raises `infeasible` and is *recorded*
  (not a crash); only a wiring/settings error raises `fatal` and halts.

### L12 — `diagnostics` (the independent oracle)

```r
arm_evaluation("diagnostics")
anchor <- canary_anchor(s)                       # counts from the 5 teams' submission CSVs (no optimizer code)
cal    <- calibrate_canary_trials(s, conn, anchor)   # our counts vs the anchor, per trial
print_calibration(cal)
check_canaries(s, conn)                          # the startup bug oracle
disarm_evaluation()
```

- **→ returns** `cal` a per-trial table with a `divergent` flag; `check_canaries`
  prints `CANARY ALARM` (hard, strong trial) or a soft warning (weak trial).
- **🔍 eyeball** every non-`divergent` row means our pipeline's count **agrees**
  with the independent anchor — positive evidence the plumbing is right.
- **🚩 red flag** any `divergent` row (our count outside 0.5×–2× the anchor) — a
  bug localized to that trial; a `CANARY ALARM` on a strong trial means the
  "infeasible" verdicts elsewhere cannot be trusted.
- *Assumption:* the anchor shares **none** of the optimizer's code, so a bug rarely
  corrupts it and our derivation identically — agreement is real evidence. Freeze
  `canary_configs()` only once calibration agrees (human-reviewed).

---

## 5. Fast → slow within the live layer

Live work is dominated by VCF downloads and mixed-model fits. Order it so the
cheap paths clear first:

- **Cached trial first.** L5–L11 above use `10676` (already in `cache/`) so the
  first full end-to-end pass downloads nothing. Only then move to a fresh trial.
- **Few training trials before many.** In subtask A, start with a *high*
  `accession_overlap` `primary_min` (fewest neighbours → fewest downloads); lower
  it only once the high-threshold path is clean. The prompt's rule of thumb —
  overlap 20 before overlap 3.
- **Smallest canary before largest.** `calibrate_canary_trials()` /
  `canary_anchor()` print per-trial accession counts; check the smallest trial
  before the largest so a bug surfaces on cheap data.
- **Thin markers on the first touch of a big project.** Pass `marker_thin > 1`
  (e.g. via a config's `geno_select.marker_thin`) so the first parse of a large
  VCF is cheap; repeat un-thinned once the path is trusted.
- **Watch the progress bars.** The slow loops report progress: `"Load project
  dosages"` (per project in `choose_geno_sources`), `"Similarity: candidate
  trials"`, `"Observations from study ids"`, `"Genotyping projects for
  accessions"`, `"Canary trials"`, `"Calibrate canaries"`. A bar that stalls tells
  you *which* item is slow.

---

## 6. Subtle-bug catalogue

The failure modes that make a broken optimizer look like a working one, each with
the check that exposes it. This is the list to internalize — the whole point is
noticing something "a little funny" that the code's own checks did not.

| Signature | What it means | Expose it with |
|---|---|---|
| `rep`/`block` (or any column) **all `NA`** | BLUE silently degraded to germplasm means | `peek(train_obs)` NA-per-column |
| **targets all identical / degenerate** | `build_targets` collapsed | `peek(targets)` distinct count |
| **dosage rowname overlap = 0** on a non-empty matrix | synonym / name mismatch | `peek(dl, accessions=…)`; `diagnose_trial` |
| **K not symmetric / bad diagonal / clones differ** | kernel construction bug | `peek(K)` |
| **predictions constant → `NA`**, or **score ≈ 1** | model collapse, or leakage | `peek(pred)`; sanity-check the score |
| **CV0 and CV00 scores identical** | `mask_cv` no-op — the CV00 mask isn't excluding focal accessions | the L4 differential check |
| **funnel cliff** (many in, ~0 out), tagged `suspect` | data-hiding bug, not genuine infeasibility | the `flow` funnel; `report.md` ⚠ section |
| **incumbent never beats the seeds** after many iters | search-space / surrogate bug | `report.md` learning curve |

The meta-rule behind all of them: **trust agreement between two independent
derivations; distrust a single number.** Our per-trial count vs `canary_anchor`;
the pipeline funnel vs `diagnose_trial`'s raw re-derivation. One bug rarely
corrupts both the same way — so a match is evidence, and a mismatch localizes.

---

## 7. The test suite (offline, deterministic)

### Principle: oracle tests, not snapshot tests

A **snapshot test** asserts "the function returns whatever it currently returns."
It re-encodes the implementation, passes even when the code is wrong, and takes as
long to audit as the code. We avoid those. Every test has an **independent oracle**
— the expected answer is knowable without reading the implementation. Four kinds:

- **Closed-form** — correlation of a vector with itself is 1; a GRM among
  genetically identical lines has identical rows; centering by trial mean removes a
  constant added to one trial.
- **Property / invariant** — Expected Improvement is ≥ 0; a relationship matrix is
  symmetric; two-stage BLUP shrinks effects toward the mean; QC drops a monomorphic
  marker.
- **Metamorphic / differential** — double the phenotypes → BLUEs double; permute
  accessions → predictions permute identically; **CV00 excludes exactly the
  focal-trial phenotypes CV0 keeps**; raising `blend_obs_w` changes CV0 predictions
  but is a no-op under CV00.
- **Tiny hand-computable** — a 3-marker × 2-sample VCF whose dosage matrix you can
  write down by hand.

Because the tests *run*, a test that encodes a wrong expectation usually **fails
against correct code**, surfacing the disagreement to adjudicate rather than
silently lying. Each test file comments its oracle inline.

Unit tests cover the **deterministic, offline** logic — most of the optimizer.
They do **not** verify that the live BrAPI responses have the shape
`R/data_access.R` assumes; that is what the canary oracle (§9) is for. *Unit tests
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
| `tests/test_config_space.R` | ~5 genome invariants across ~400 sampled/recombined configs → `config_space tests: 8007 passed, 0 failed` (8007 = individual assertions) |
| `tests/test_subtasks.R` | `Tier 1 subtask tests: 45 passed, 0 failed` |
| `tests/run_all.R` | `2/2 test files passed` |
| `tests/test_sim_loop.R` (or `run_all.R --all`) | `PASS: optimizer beats submissions and improves over random search`, exit 0 |

If a count drifts after a change, that is the regression signal — reconcile it
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
validity) asserted across ~400 randomly sampled and recombined configs — hence the
`8007 passed` counter (individual assertions, not distinct cases).

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

## 8. How the code touches cache / state / logs

**`cache/`** — the on-disk memoizer is `cached(settings, name, expr, max_age_days)`
in `data_access.R`; every entry is `cache/<name>.rds` (or a raw `.vcf`). Touch a
trial once, reuse forever.

| Cache file | Written by | Holds | Max age |
|---|---|---|---|
| `trial_catalog.rds` | `trial_catalog()` | focal-trait trials + metadata + lat/long/elev | 7 d |
| `acc_<studyid>.rds` | `get_trial_accessions()` | accession names of a trial | 30 d |
| `obs_<studyid>.rds` | `get_observations()` | all numeric observations of a study | 30 d |
| `proj_<hash>.rds` | `projects_for_accessions()` | genotyping **project** ids covering an accession set | 30 d |
| `raw_project_<id>.vcf` | `.ensure_project_vcf()` | the archived VCF (validated complete). **Transient**: deleted once its full dosage is cached | until dosage cached |
| `dosage_<id>[_thin<N>]_sz<size>.rds` | `get_project_dosage()` | the **whole project's** accessions×markers dosage (subset at read). Looked up by pattern (ignoring `sz`); the VCF byte size invalidates a stale dosage if re-downloaded at a different size | ∞ |
| `calibration_lightweight.rds` | (diagnostics, when you save it) | the last calibration table | — |

**`state/`** — the single source of truth.

- `evals.sqlite` (`store.R`): one table `evals` with columns
  `id, config_hash, config_json, trial_id, scheme, score, n_test, status, reason,
  detail, seconds, ts`. Every evaluation — success **and** failure — is a row;
  the incumbent, surrogate training data, and "what to try next" are all recomputed
  from this table, so the process is killable/resumable. `status` ∈
  `ok | infeasible | suspect | error`; `reason` is the machine code; `detail` is the
  failure funnel.
- `report.md` (`report.R::write_report`): rewritten every `checkpoint_every`
  iterations — incumbent, subtask-method importance, the failure log
  (`failure_summary`), the "⚠ Suspected bugs" section, and the running-best curve.
- `STOP`: create it (`touch state/STOP`) to halt cleanly after the current eval.

**`logs/`** — `run_optimizer.R` writes a per-iteration heartbeat to stderr; the
launch command redirects it to `logs/run.out`. The startup canary result and any
`CANARY ALARM` / `step error` / `FATAL` messages appear here.

---

## 9. Validation & debugging tooling

### Failure classification (`R/conditions.R`)

Three deliberately distinct "didn't work" signals, so the loop knows whether to
continue or stop:

- **`infeasible(code, detail, funnel, suspect)`** — this *one* (trial, config,
  scheme) cannot be evaluated (too little genotyped overlap, no training trials,
  too few markers, …). Expected. `evaluate_config_on_trial` records it
  (`status = "infeasible"`, `reason = code`, `detail = funnel`) and the loop samples
  a fresh trial + config next iteration. If the funnel shows a **cliff** (lots of
  data in, ~zero out), it is tagged `suspect` instead — a likely bug, broken out in
  the report.
- **`fatal(message, code)`** — the whole run cannot proceed (bad settings, an
  unimplemented method). Re-raised past `evaluate.R` so `run_optimizer()` halts
  instead of spinning. (Implementation note: `evaluate_config_on_trial` uses a
  *single* `error` handler that branches on class — re-raising a fatal from a
  multi-handler `tryCatch` would be swallowed by the sibling `error` handler.)
- **`sample_failed(detail)`** — no usable trial could be sampled at all. Not a
  config's fault; the loop tolerates a few in a row then stops
  (`settings$max_sample_fail`).

`failure_summary()` (`report.R`) turns the stored failures into: counts by status,
dominant reasons, and the **failure rate of each subtask method** — i.e. which
configurations break most often.

### The canary bug-oracle (`R/diagnostics.R`)

The "skip and continue" design means a **bug that hides data looks identical to a
trial that is genuinely infeasible.** The oracle catches that, via a
**calibrate-then-freeze** bootstrap against ground truth that does **not** touch our
code — the five Predictathon teams' submission files.

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
independent re-derivation of the raw counts from T3** — including, per genotyping
project, whether the dosage matrix rownames actually intersect the trial's
accession names. A line like `overlap with accessions = 0` on a non-empty VCF is
the decisive synonym/name-mismatch signal. Reach for it whenever `subtaskC`
(§L8) or the `flow` funnel (§L11) looks wrong.

### Lessons baked in (don't re-learn the hard way)

- **Cache poisoning + big VCFs.** A partial/raced VCF download once cached a tiny
  dosage matrix *forever* (`get_project_dosage` used `max_age_days = Inf`). Now
  `.ensure_project_vcf()` validates completeness (`.vcf_complete`) and re-downloads
  truncated files. `get_project_dosage()` extracts and caches the **whole
  project's** dosage once (keyed by project+thin+VCF-size), reused across trials by
  subsetting at read, and **deletes the raw VCF once its full dosage is cached**. If
  genotype counts ever look low, clear `cache/dosage_*` (regenerable).
- **Interactive prompt.** `conn$vcf_archived()` can prompt to pick a file and hang a
  non-interactive run — always launch real-mode scripts with `</dev/null`.

---

## 10. Extending the search space

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
5. If you added a function worth stepping through, add its name to the appropriate
   group in `R/evaluation.R::EVAL_GROUPS` so `arm_evaluation()` reaches it.

Nothing else in the optimizer needs to change — the surrogate encoding, crossover,
mutation, store, and report all derive from `SUBTASKS`.
