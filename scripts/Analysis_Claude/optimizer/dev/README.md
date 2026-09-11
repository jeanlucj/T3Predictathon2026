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
| `analyze_failures.R` | The failure rate moved and you want to know whether it is the code or the domain. Answers from the store alone — failures by trial, `train_select` crossed with trial (the report's per-method table charges a failure to all six subtasks, so five of its six columns are the base rate), the `constant` bucket, an old/new domain diff keyed on `trial_id`, before/after rates restricted to the trials both stores share, and the true contender count the report caps at 8. Read-only, seconds, safe beside a live run — reach for `tools/diagnose_failures.R` only if one trial is still unexplained after this. |
| `holdout_test.R` | You want to know whether optimization actually helped. Runs the top-k optimized configurations **and** all five submissions on trials the optimizer never saw, then compares the two distributions — a mixed model with configuration nested in group as the primary test, plus Wilcoxon/t on the per-configuration means. Refuses any trial that is in the run's pinned universe, and writes to a store of its own so a held-out evaluation can never become training data. `--report` analyses what is already stored; `dev/run_holdout.sh` runs N copies. **`dev/HOLDOUT.md` is the walkthrough** — what calls what, worker sizing, running it beside the optimizer, and how to read the four report sections. This is open thread 1, and the manuscript's headline claim. |
| `check_domain_coverage.R` | A run aborted with `domain_not_covered`. Walks every name in `target_domain$trials` down the same ladder `.eligible_trials()` applies and names the step that dropped it: no focal-trait data, name not found in T3, or dropped by `programs`/`years`/`locations`. The distinction matters because the fixes differ — the first is a property of the data, the last is a self-inconsistent domain, since `trials` **intersects** with the other filters rather than overriding them. |
| `check_trial_renames.R` | The store and the catalogue disagree about a trial. Asks whether a trial's `study_name` changed between the two. This was the livelock of August 2026 (`docs/LESSONS.md`), and the reason the run universe is pinned to trial **ids** now — so this script should find nothing, and a hit means something upstream moved. |
| `EM_COMBINE_COMPARISON.md` | You are working on `em_combine`. What the sibling pipeline's combiner does differently, item by item, and which of those differences were deliberately declined. |

`run_holdout.sh` is the only launcher in here; it is `dev/` rather than the root because it
drives `holdout_test.R`, not the optimizer.

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

*Rewritten 2026-09-03 at build 0.8.7; robustness section added 2026-09-04 at 0.8.8.
Supersedes `NEXT_STEPS.md`, whose ranking was taken at 0.7.6 and predates the replicate
livelock, the trial-id universe pinning, the MERF surrogate switch and two OOMs. This section
is the project's live next-steps document — there is deliberately no separate file, because a
separate file is what let the last one go stale unnoticed.*

### Immediate — robustness, found launching the run of 2026-09-04

These are ahead of the ranked list rather than in it: they are not throughput or quality levers,
they are the difference between a failure you see and one you find hours later. They came out of
launches that died and looked, from outside, like nothing at all. **A is fixed (2026-09-08); B, C, D and E remain.**

**B. The message the surviving workers print points away from the cause.** Worker 1 fatals
immediately with the real reason; workers 2..N wait 600 s and then report `no run row after
10 min: the leader never recorded run <id>`. That reads like a store or leader-election problem
and says nothing about where to look. It should name `run_w1.out` explicitly, and better, read
the leader's failure back out and quote it. Nine tenths of the logs describe a symptom, and
they are the ones a person opens first.

**C. `target_domain$trials` matches by `study_name`.** `.apply_target_domain()` runs *before*
`pin_trial_universe()` resolves ids, so the domain spec is the last name-based coupling in the
run path — precisely the fragility that pinning the universe to trial ids was introduced to
remove, still present one step upstream. A rename on T3 can no longer move a run that is going;
it can still stop one from starting. On 2026-09-04, 15 of 32 named trials failed to resolve
because T3 had canonicalised them to `Big6_<year>_<LOC>` (`Big6_Fra_23` -> `Big6_2023_FRA`, and
so on). Accepting `study_db_id` in `target_domain$trials` alongside names would close it.

Two things that made this harder to see than it needed to be, both worth fixing with C:

- **A trial can be named in the domain and not measure the focal trait.** Trials selected on
  accession overlap with the focal trial are not thereby trials that record grain yield, and
  nothing checked. `trial_catalog()` searches on `observationVariableDbIds`, so such a trial is
  simply absent, and `pin_trial_universe()` reports it identically to a typo. Whatever picks
  candidate trials should filter on the trait at selection time.
- **The `domain_not_covered` message lists names without causes.** There are three distinct
  reasons with three different fixes — no focal-trait data, name not found, dropped by an
  earlier domain filter. `dev/check_domain_coverage.R` separates them after the fact; the
  fatal itself could do it at the point of failure.

**D. A stopping rule, and the resolution ceiling behind it.** The search has no termination
criterion but `max_iters`. The natural proposal — stop when the top-k contenders stop changing —
is the right family but too brittle as stated: membership moves as `se` shrinks even when nothing
is learned; 8 configurations at z = 1 against 280 at z = 2 means the boundary flips on noise; and
counting "new configurations tried" degenerates toward a wall clock, since with continuous
parameters nearly every mutation is a new configuration. Two conditions together are steadier:
nothing entering the top-k over the last N **scored** configurations, **and** the incumbent's mean
improving by less than epsilon.

**Do the cheap prerequisite first.** `choose_config()` already computes
`expected_improvement()` and returns `ei` with its pick, and it goes only to a log line
(`run_optimizer.R:97`) — not to a column on `evals`. A collapsing max-EI is the textbook stop for
this kind of search, and persisting that one number would let any plateau rule be tested against
history instead of argued about.

**The deeper point may be that no rule is needed.** With `sd_config` 0.015 and `sd_resid` 0.068
per evaluation, the fixed-effect SE of a configuration's mean is `0.068/sqrt(n)`: 0.0139 at the
24 evaluations a full domain sweep allows, and reaching 0.0075 would take 82. So configurations
closer than about one `sd_config` apart cannot be separated however long the search runs. The
ceiling is the domain size, not the budget — which argues for widening the domain over running
longer, and is `docs/LESSONS.md` #19 arriving from a third direction.

**But do not trust a plateau until E is done.** A stalled top-8 is not evidence of convergence
while the proposal operators structurally cannot refine a good configuration — see E. Try
refinement first; a stall that survives it means something.

**E. The proposal operators cannot refine a good configuration.** The search space is **648
categorical method combinations** (3 x 4 x 3 x 3 x 3 x 2 over the six subtasks) wrapped around
**14 continuous/integer parameters**, plus a few categorical sub-parameters. Of the three
operators in `R/config_space.R`, none can move within that interior at a fixed skeleton:

- `crossover(a, b)` copies a whole subtask block — method *and* its parameters — from one
  parent or the other. It moves between skeletons but **never creates a parameter value**: every
  number in the child already existed in a parent.
- `mutate_config(cfg, n_subtasks = 1:2)` deletes the chosen subtask's block entirely and calls
  `sample_block()` — a **fresh random draw of method and parameters together**. There is no
  parameter-only mutation.
- `sample_config()` draws everything at random.

So a new parameter value can only arrive by a random draw that simultaneously re-rolls the
method. Once the skeleton is right, the search has no way to walk downhill in 14 dimensions — it
can only re-roll and hope. The generation side compounds it: `n_elites = 8`, and
`propose_candidates()` draws all 60 crossovers and all 60 mutations from those 8, so **120 of 180
candidates come from 8 configurations** and a frozen leaderboard freezes two thirds of every
round.

Four pieces, cheapest first:

1. **A refinement operator.** Hold a contender's method skeleton fixed and perturb only its
   continuous/integer parameters — a local jitter scaled to each parameter's declared range,
   shrinking as the field matures. This is what converts "we found a good method combination"
   into "we found its best settings", and it is simply absent. `sample_block()` already reads
   each parameter's type and range from `SUBTASKS`, so a `jitter_params()` beside it reuses that
   metadata rather than restating it.

2. **The per-subtask argmax skeleton.** Take the best method in each subtask by the surrogate and
   combine them: one candidate, essentially free, and *no contender currently uses it*. Marginal
   bests need not combine well — interactions are exactly what a forest can represent and a
   marginal table cannot — which is the argument for evaluating it rather than assuming it.

3. **Coordinate ascent on the surrogate.** From each elite, sweep one subtask at a time to the
   value maximising predicted mean, iterating to a local optimum. A random forest has no
   gradient, but each coordinate is small enough to enumerate. `choose_config()` picks
   `which.max(ei)` — pure Expected Improvement, which rewards uncertainty — so nothing ever
   *constructs* the configuration the surrogate thinks is best; it is found only if it happens to
   be sampled, which across 648 skeletons is unlikely. Keep EI as the selector: these add
   candidates, they do not replace the acquisition, and picking on predicted mean alone would
   collapse the search.

4. **Make the 60/60/60 split adaptive.** Heavy random and crossover while the field is unformed,
   shifting toward refinement and exploitation once the leaders are stable. The stall is the
   natural trigger — the same signal D would use to stop, redirected to change the proposal mix
   first.

**How to tell whether any of it worked**: persist `ei`, which D also needs. With it on `evals`,
"did the new operators raise achievable EI, or confirm the plateau?" is answerable from the store
rather than by argument — which makes it the shared prerequisite for D and E both, and worth
doing before either.

------------------------------------------------------------------------

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

The same mtime-for-liveness defect was found and **fixed in `run_workers.sh` on 2026-09-08**: it
refused a relaunch because a cancelled job's `run_w5.out` was still warm, and `$LOG_DIR` is
durable, so the refusal followed the user to a different node. It now asks the OS instead. That
leaves `.acquire_lock` as the remaining instance, and it is the one where the consequence is
worse — a stolen lock corrupts a download rather than declining to start.

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
- **Panel selection** — 0.8.7 made `best_single_project` apply `.best_panel`'s floors (it had
  bypassed them, so a project covering no focal lines could win on training coverage alone), and
  0.8.8 replaced the max-union rule with the harmonic mean of focal and training coverage. Still
  open underneath it: the exactly-right criterion is mean GBLUP reliability over the test lines,
  which needs a GRM per candidate panel; and whether the harmonic mean or the geometric mean is
  right turns on whether accuracy has saturated in training size, which `n_test` in the store
  can settle empirically once the new run has rows.
- **A dead run reporting success to SLURM** — fixed 2026-09-08, after a second multi-day loss.
  `run_workers.sh` ended in a bare `wait`, which returns 0 whatever the children did, so every
  worker could die in the first minutes and `sacct` would still say `COMPLETED` —
  `--mail-type=FAIL` could not fire and the only signal was `t3opt` leaving `squeue`. It now
  waits on each pid and separates the two cases that matter: workers gone within the launch
  window with no STOP file is a failed *launch* and exits non-zero; a worker lost later is
  attrition, which the run survives by design, and exits 0 with a warning. Verified that the
  status propagates through `t3opt_ceres.sh`'s `set -e` and its `EXIT` trap.
- **The startup connect having no retry** — fixed 2026-09-08, the cause of that second loss.
  `t3_connect()` was the one BrAPI call outside `.brapi_try()`, so a 10 s connect timeout at
  relaunch killed all 22 workers, which the 20 s launch stagger put inside one outage window.
  It now retries on `brapi_connect_tries` (8, ~5 minutes) — deliberately longer than
  `brapi_tries`, since a worker that cannot connect has nothing else to do while a call inside
  the loop is holding a claim. Credential errors still fail fast: `.brapi_try()` documented
  that intent but only honoured it in its re-login branch, and now honours it when the thunk
  raises those classes too.
- **The replicate livelock** — fixed 2026-09-02. A mid-run trial rename moved rows out of the
  domain slice, so replication could never complete. The universe is pinned to trial **ids** at
  run start, there is a `runs` table recording it, and `max_consec_skip` bounds the retry.
- **EM_COMBINE item 2** (E-step regularization) — investigated 2026-07-31 and deliberately
  declined: pre-ridging the partials already bounds the E-step inversion, and the epsilon
  converts a convergent iteration into a drifting one.
- **Excessive comments and dead options** — 0.8.4. Rationale that exists to stop a bad idea
  being retried goes in `docs/LESSONS.md`, not inline; an option discovered to be a bad idea is
  removed rather than left available.
