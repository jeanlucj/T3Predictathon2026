# `dev/` — changing the system

For someone **modifying** the optimizer: choosing between algorithm designs, measuring where an
evaluation's time goes, or working out why something broke in a way `tools/` cannot answer.

Operating a run is `../tools/README.md`. Understanding the system is `../docs/`.

## What is in here

| file | when you reach for it |
|---|---|
| `surrogate_bakeoff.R` | You are choosing a surrogate design. Cross-validates the candidate arms (MERF, blocked, pooled) on the stored evaluations and draws a learning curve, so "which one predicts held-out configurations best, and does it improve as evaluations accumulate?" is answered with numbers rather than a preference. Single core, saturated, for a long time — ~780 fits at `ntree=200`, and `--no-curve` removes about three-quarters of that. |
| `profile_evaluation.R` | You are about to optimize something and want to know whether it matters. Runs a real evaluation with timers on each subtask. **~12 GB and hours, and it takes the dosage locks** — same hazard as `tools/diagnose_failures.R`, so not beside a live run. |
| `check_blas.R` | R is slower than it should be on a new machine. Says whether the linear algebra is actually threaded. Pins every core, briefly, so cap threads before running it beside anything. |
| `check_oom.sh` | Something in a SLURM allocation died and you suspect memory. Reads the job's cgroup counters and the `sacct` record. Must run **inside** the allocation being asked about: the cgroup high-water mark survives the process but disappears with the allocation. |
| `check_trial_renames.R` | The store and the catalogue disagree about a trial. Asks whether a trial's `study_name` changed between the two. This was the livelock of August 2026 (`docs/LESSONS.md`), and the reason the run universe is pinned to trial **ids** now — so this script should find nothing, and a hit means something upstream moved. |
| `EM_COMBINE_COMPARISON.md` | You are working on `em_combine`. What the sibling pipeline's combiner does differently, item by item, and which of those differences were deliberately declined. |

Everything here runs **from the optimizer root**, like `tools/`:

``` bash
cd <repo>/scripts/Analysis_Claude/optimizer
Rscript dev/surrogate_bakeoff.R --no-curve
```

`.Renviron` — and therefore the T3 credentials — is read from the working directory only.

## Where a new script goes

- **`tools/`** if an operator would run it during a normal run.
- **`dev/`** if it exists because something broke, or because the algorithm is being changed.
- **Nothing new at the root.** The root holds entry points (`run_optimizer.R`,
  `run_workers.sh`), the settings, and the two runbooks. That is all.

This rule is the point of the partition. Every file moved into `tools/` and `dev/` was
individually reasonable to create — each arrived in response to a bug, a server problem, or a
question about a run in flight — and 24 scripts in one namespace is what that adds up to
without a placement rule. The rule is what makes the next one land in the right place instead
of re-accreting the flat list.

A corollary for documentation: history goes in `docs/LESSONS.md`, with a one-line pointer from
the code. Not in the comments.

## Open threads

*Rewritten 2026-09-03 at build 0.8.7. Supersedes `NEXT_STEPS.md`, whose ranking was taken at
0.7.6 and predates the replicate livelock, the trial-id universe pinning, the MERF surrogate
switch and two OOMs.*

Two numbers still drive prioritisation, and they rank the work below very differently:

- **Evaluations are expensive.** Median ~29 min, mean ~122 min, max 33 hours.
- **Memory is the binding resource, not cores.** `peak_r_mb` median ~19 GB, max ~82 GB — and
  that is R's heap peak, an *under*-estimate of true RSS.

So each thread is scored on **throughput** (evaluations per unit time) and **quality** (how much
closer each evaluation gets to the right answer). Conflating the two is how this project has
mis-prioritised before: `docs/BACKGROUND.md` §4 records that the search is only nominally better
than random at realistic budgets, and `docs/LESSONS.md` #19 records the diagnosis — the
constraint was never measurement precision.

| # | thread | throughput | quality | effort |
|---|---|---|---|---|
| 1 | Honest cross-trial accuracy has never been measured | — | *it is the headline claim* | small |
| 2 | GRM / kernel caching | **high** | — | medium |
| 3 | Load only the panels actually used | **high** | — | medium |
| 4 | Estimate a method's memory before choosing it | medium | — | high |
| 5 | Extract `T3BrapiHelpers` | — | — | medium |
| 6 | Pedigree bridging for `em_combine` | — | medium | medium |
| 7 | `.acquire_lock` judges staleness by mtime, not liveness | — | correctness | small |
| 8 | Method modularization | — | indirect | medium |
| 9 | Multi-fidelity / early stopping; multi-objective CV0+CV00 | high | risk of bias | large |

### 1. Honest cross-trial accuracy has never been measured

`docs/BACKGROUND.md` §4: the canary scores *feasibility*, not tuned accuracy, and an incumbent's
reported score is an optimistic estimate — the maximum over many noisy means. The BLUP shrinkage
added in 0.7.6 reduces that winner's-curse bias but does not remove it. Since "we beat the
submissions" is the project's headline claim, a clean held-out comparison is small work with
outsized reporting value. **Do this before the manuscript.**

### 2. GRM / kernel caching

Caching stops at dosage matrices (`docs/DESIGN.md` §4). The kernel is rebuilt on **every**
evaluation, though it depends only on `(panel set, maf, max_missing, impute, method, ridge)` — a
serviceable cache key, and one many configurations share. Median `em_combine` 36.4 min against
`vanRaden_single` 22.4 puts kernel construction at roughly 14 minutes of a 36-minute evaluation,
but **that is a proxy, not a measurement**: the two kernels differ in more than construction
cost. Instrument `build_kernel()` and record its share in the store first, then decide.

Build it on `cached()` + `no_cache()`. A GRM key is a single deterministic key, so the memoizer
fits directly. Do not copy the dosage cache's filename-encoded scheme — that shape exists only
because marker density is negotiable and a GRM's is not.

### 3. Load only the panels actually used

Panel *selection* needs no dosage at all: which panel covers `need` best is answerable from
`.project_index` (accession → projects). `.group_by_panel` compares marker *sets*, so caching
per-project marker names — one small file per project, the same shape as the other maps — would
let grouping happen without loading either. Then load only the panels that survive. This is the
largest remaining reduction in the memory peak, which is the binding resource.

### 4. Estimate a method's memory before choosing it

If a random-forest estimate of an evaluation's peak memory exceeded a fraction of what is free,
the worker could take the next candidate down the list instead of firing up a method that will
OOM. Less efficient per evaluation, plausibly more efficient per compute-hour on a
memory-limited machine — and it would stop the pathology where every worker simultaneously tries
`em_combine`.

**Where the training data has to come from**, because this is the part that sinks the naive
version: not `evals`. `store_eval()` runs only on completion, so an OOM-killed evaluation —
precisely the one worth learning from — leaves no row at all, and the store is censored on
exactly the outcome being predicted. Three external sources, in descending usefulness:

- `tools/watch_memory.sh` already logs a per-pid `pids_rss` column each interval. The last
  sample before a pid disappears is the victim's RSS, and it is the only source that attributes
  a kill to a specific worker. `claims` now records `worker`, `host`, `pid` and `ts` per
  in-flight evaluation, so the join from a sample back to the `(config, trial)` it was working
  on is newly possible — and a killed worker's claim is the record of what it died doing.
  Nothing does that join yet. **This is the one to build on.**
- `dmesg` gives the kernel OOM killer's exact `anon-rss` per pid, but is usually unreadable as a
  non-root user on Ceres.
- `sacct -j <jobid> --format=JobID,MaxRSS,State` is real but reports the whole job step, which
  with N workers on one node is the aggregate high-water mark, not the culprit's.

### 5. Extract `T3BrapiHelpers`

`R/data_access.R` and `R/genotypes.R` hold BrAPI access, retry budgeting, the wizard indices and
the VCF parser — code that is not specific to this optimizer and that the wider T3 work would
use. Extracting it is still wanted.

A candidate list was drafted 2026-07-17 and has been dropped rather than carried: it was a
snapshot with unfilled `User notes:` lines, and a list of functions that were extraction
candidates six weeks and several refactors ago is not the document to work from. **Regenerate it
when the work is actually scheduled** — the current call sites are the authority, not a stale
inventory.

### 6. Pedigree bridging for `em_combine`

`EM_COMBINE_COMPARISON.md` item 3. The sibling pipeline can pass pedigree relationship matrices
as extra partials, stitching marker panels that share no accessions; here, disjoint panels
degrade to block-diagonal with no error (`docs/LESSONS.md` #13). A genuine capability gap, and it
bites hardest exactly where combining would be most valuable.

### 7. `.acquire_lock` judges staleness by mtime, not by whether the holder is alive

Every cache lock already carries `host=… pid=… worker=… time=…`, and `.acquire_lock` then breaks
it on `lock_stale_minutes` (default 90). So a download slower than that gets its lock stolen
while it is still running, and a lock whose holder actually died is held for the full 90 minutes
— both failure modes, from the same missing check. `R/store.R::.pid_alive()` is the check it
wants; the claims table already uses it. Small, and it is the mechanism behind
`docs/LESSONS.md` #24.

### 8. Method modularization

Adding a method today means editing `SUBTASKS` in `R/config_space.R` *and* adding a branch to
the dispatcher in `R/pipeline.R` — documented in the `config_space.R` header, and cleaner than
most research code, but it needs two files and the dispatcher idiom. What would make it genuinely
contributable: a `register_method()` call so a method arrives as one self-contained unit;
validation at registration rather than a failure deep in a run hours later; one worked example
plus the contract each subtask must honour (inputs, return shape, which `infeasible()` codes it
may raise); and a test template, so a contributed method arrives with an oracle. The value is
adoption, not accuracy — worth doing only if others are meant to use this.

### 9. Deferred by design

Multi-fidelity evaluation (`docs/BACKGROUND.md` §4) is the other large throughput idea — not
every configuration needs the full marker budget — but truncating evaluations biases scores in a
way that would need its own validation. Multi-objective CV0/CV00 changes the question being
asked. The remaining `EM_COMBINE` items (glmnet imputation, panel-size floor) are small.

### Closed, so it is not re-derived

- **Config and trial replication.** `docs/LESSONS.md` #19 said the trial effect was "not even
  estimable" because a fresh trial was sampled per evaluation. 0.7.5's `trial_replication` fixed
  that, 0.7.6 put the trial adjustment into the incumbent and elite pool via random-effects
  BLUPs, and 0.8.5's `config_replication` closed the other half. 0.8.6 ramped it.
- **Caching-layer consolidation and the VCF split** — 0.8.5. `cached()` gained `no_cache()`;
  `.neg_cache_hit()`/`.neg_cache_mark()` name the permanent negative cache; `R/genotypes.R` now
  holds the VCF parser, the download retry budget and the dosage cache. `get_project_dosage` was
  deliberately **not** folded into `cached()`: its key is `(project, thin, size)` encoded in a
  filename and resolved by globbing for the densest match, and generalising the memoizer to
  cover that would burden four simple callers with machinery one needs.
- **SciNet port** — 0.8.5. See `container/README.md`.
- **MERF as the surrogate** — promoted from a bake-off arm to the default config picker for
  SMAC. `surrogate_bakeoff.R` is how that decision is revisited on a larger store.
- **The replicate livelock** — fixed 2026-09-02. A mid-run trial rename moved rows out of the
  domain slice, so replication could never complete. The universe is pinned to trial **ids** at
  run start, there is a `runs` table recording it, and `max_consec_skip` bounds the retry.
- **EM_COMBINE item 2** (E-step regularization) — investigated 2026-07-31 and deliberately
  declined: pre-ridging the partials already bounds the E-step inversion, and the epsilon
  converts a convergent iteration into a drifting one.
- **Excessive comments and dead options** — 0.8.4. Rationale that exists to stop a bad idea
  being retried goes in `docs/LESSONS.md`, not inline; an option discovered to be a bad idea is
  removed rather than left available.
