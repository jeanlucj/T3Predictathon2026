# The held-out test

Do the optimized configurations actually beat the five submissions?

An incumbent's score in `state/report.md` is a maximum over noisy means **on the trials it was
selected from**, so it is optimistic by construction (`docs/BACKGROUND.md` §4). The only way to
know whether optimization helped is to run the optimized configurations *and* the submissions on
trials the optimizer never saw. That is what this does, and it is the measurement the project's
headline claim rests on.

## What calls what

The confusing part, so it is first. `holdout_test.R` is the **innermost** thing; the submit
script is the outermost:

```
container/submit_holdout.sh      you run this. Resolves ACCOUNT and OPTIMIZER_HOME, calls sbatch.
  └─ container/holdout_ceres.sh  the batch job. Loads apptainer, then inside the image:
     └─ dev/run_holdout.sh N T   launches N copies, waits, reports a non-zero exit if all died
        └─ dev/holdout_test.R    × N — one R process per copy
```

Each copy walks the same **(configuration × trial) grid** and takes cells via `claim_eval()`, the
same interlock the optimizer uses. There is no search and nothing to coordinate beyond that: the
grid is fixed before any evaluation starts.

## Three entry points

| you want to | run |
|---|---|
| the real thing, as a batch job | `cd container && ./submit_holdout.sh -- --trials=<ids>` |
| the same, inside an allocation you already hold | `./dev/run_holdout.sh 12 1 --trials=<ids>` |
| see the grid without evaluating anything | `Rscript dev/holdout_test.R --trials=<ids> --dry-run` |
| read the result | `Rscript dev/holdout_test.R --report` |
| any of the above, for one domain of several | add `--settings=settings.local.XYZ.R` (below) |

`--report` needs neither the optimizer's store nor T3 — the group of each row is recorded in the
holdout store — so the analysis runs anywhere, any time after the evaluations.

Everything runs **from the optimizer root**: `.Renviron`, and so the T3 credentials, is read from
the working directory only.

## Which domain's contenders — `--settings=`

**The contenders are not stored anywhere; every copy works them out again from the optimizer's
store when it starts.** To do that it must know which optimization run it is looking at, and the
store may hold several: `evals.sqlite` (and its backup) is an archive of every evaluation ever
made, across target domains.

The chain:

1. **Settings:** it builds the effective settings.
2. **Run id:** `run_id_for()` hashes the run-defining ones (`.RUN_SIGNATURE_KEYS` in
   `R/store.R`: `optimize_scheme`, `target_domain`, `focal_trait_db_id`, `min_trial_acc`, the
   replication settings, `replicate_every`, `contender_z`, `n_random_init`, `simulate`) plus the
   build.
3. **Lookup:** that id is looked up in the store's `runs` table, whose row holds the run's pinned
   trial universe.
4. **Filtering:** evaluations are filtered to that scheme, build and universe, and the pool is
   drawn from what is left.

**By default the settings come from `settings.R` + `settings.local.R`**, i.e. whatever domain
the checkout is currently set up to optimize. To test a different one without editing
`settings.R`, write a **domain file** and name it on the command line:

```r
# settings.local.big6.R  -- in the optimizer root, gitignored
settings_override <- list(
  target_domain = list(programs = NULL, years = NULL, locations = NULL,
                       trials = c("Big6_2023_FRA", "Big6_2023_URB", ...)),
  optimize_scheme = "CV00",
  build = "0.8.7"            # only if that run was made under an older build
)
```

```
./submit_holdout.sh -- --settings=settings.local.big6.R --trials=<ids>
Rscript dev/holdout_test.R --settings=settings.local.big6.R --trials=<ids> --dry-run
Rscript dev/holdout_test.R --settings=settings.local.big6.R --report
```

How it is applied:

- **Order:** tracked defaults, then `settings.local.R` (cluster paths, sizing — unchanged), then
  the domain file.
- **Each key in the domain file replaces the whole value.** Settings are not merged, so a
  `target_domain` must be written out in full, `NULL` fields included, exactly as it was when the
  run was made.
- **Only run-defining keys are allowed** (the signature keys above, plus `build`). A domain file
  that sets anything else — `db_path`, `cache_dir`, … — is rejected, so it can never move the
  store or cache that `settings.local.R` put on node-local disk.
- **The file is gitignored** (`settings.local.*.R`), so copy it to the cluster checkout yourself,
  as with `settings.local.R`. The path is relative to the optimizer root.

**Where the domain's trial list comes from.** Two tables are involved:
- **Evaluations (`evals`)** are the measurements: one row per configuration × trial × scheme.
- **Run records (`runs`)** hold one row per optimization run: its settings and its domain
  resolved into T3 trial ids.

The filter keeps any evaluation whose trial is in the list, whichever run produced it; the run
record only supplies the list. There are two ways to get it:

1. **From the run record (preferred).** Your settings are hashed into a run id and looked up. A
   single differing key or build gives a different id, so the file must reproduce the run's
   run-defining settings exactly. On success the script prints
   `domain: settings.local.big6.R -> run <id> (N trials)`.
2. **From today's catalogue (fallback).** If no run record matches, the script prints every run
   in the store and what differs, then resolves `target_domain` against the current T3
   catalogue, as the optimizer does when a run starts (`resolve_trial_universe()`):

   ```
   no run in the store matches these settings (run_id e6b5fc4412f0, build 0.8.8).
   Runs in the store:
     7d75f214584a  started 2026-09-14 ...  build 0.8.8  scheme CV00  9 trials
         differs in: target_domain
     Resolving target_domain against the current T3 catalogue instead.
   universe: 4 named, 4 samplable
   domain: settings.local.24Crk.R -> catalogue (no run record) (4 trials)
   ```

   Trials below `min_trial_acc` are dropped, exactly as for a run. If T3 has renamed a trial
   since the run, the name is not found and the script **stops** (`domain_not_covered`) rather
   than use a shorter list. Update the name in the domain file; `dev/check_domain_coverage.R`
   says why a name was dropped. The domain's *rule* is what counts here, so the other
   run-defining keys (`contender_z`, replication settings) need not match — but `optimize_scheme`
   and `build` still pick which evaluations are read.

The store restore copies `runs` along with `evals`, so a run record survives a job restarting
on a fresh node. Records lost before that fix (`docs/LESSONS.md` #32) are reached through the
fallback.

The same logic applies **without** `--settings`, using the domain in `settings.R`. It never
carries on without a trial list: that would pool contenders from every domain in the store and
leave the "not held out" guard with nothing to check against. A matching run with no trial list
stops the script, except in simulate mode; so does a missing run record in simulate mode.

## Choosing the worker count

**Memory, not cores.** Each copy runs a full `run_pipeline`, so the same figures apply as for the
optimizer: median ~19 GB, worst case ~82 GB. `--mem` must cover `N_WORKERS × 82 GB` if you want
the worst case to be safe, and `tools/report_memory.R` is how you check what it really costs.

The number appears in **three places and they must agree**:

| where | what |
|---|---|
| `container/submit_holdout.sh` | `--ntasks` and `--mem`, the sbatch flags |
| `container/holdout_ceres.sh` | `N_WORKERS` (default 12), passed to the launcher |
| `dev/run_holdout.sh` | its first argument |

Override together: `N_WORKERS=8 ./submit_holdout.sh --ntasks=8 --mem=700G -- --trials=<ids>`.
The job echoes all three and the memory per worker in its first lines; a mismatch is visible
there rather than as an OOM eight hours in.

Copies launch **20 s apart** (`OPTIMIZER_STAGGER` to change), matching `run_workers.sh`. That is
about 4 minutes to bring up 12, which is nothing against a day-long grid, and it avoids several
processes opening the store at once — `open_store()`'s `PRAGMA synchronous` contends under that,
and `.with_busy_retry` can exhaust its five attempts.

## Running it beside the optimizer

Yes, **on a different node**, which SLURM will do by default since `t3hold` is its own job name
and `--dependency=singleton` therefore neither blocks nor is blocked by `t3opt`.

- **Dosage locks do not collide.** `cache_dir` is `$TMPDIR/...`, node-local, so two jobs on
  different nodes share no locks. `docs/LESSONS.md` #24 is a same-node hazard.
- **Reading the optimizer's store is safe.** On another node `resolve_read_store()` finds no
  node-local `db_path` and falls back to `evals_backup.sqlite`, which is copied with its
  `-wal`/`-shm` sidecars before being read.
- **The durable cache is shared and that is fine** — `sync_cache_to_backup()` is an additive,
  self-throttled rsync.

**The caveat is not technical.** While the optimizer runs, the contender set is a moving target:
a holdout started now tests *the pool as of this moment*, which may not be the pool when the
optimizer stops. For the manuscript claim, run it after the optimizer finishes. If you run it
early anyway, record the timestamp — `--report` prints the configurations it used.

## Several domains at once

Yes: submit one job per domain, each with its own `--settings=`, and they run side by side.

```
./submit_holdout.sh -- --settings=settings.local.24Crk.R --trials=...
./submit_holdout.sh -- --settings=settings.local.big6.R  --trials=...
```

The domain tag (`XYZ` from `settings.local.XYZ.R`) is added to everything two jobs could
otherwise share:

| | without `--settings` | with `settings.local.XYZ.R` |
|---|---|---|
| job name | `t3hold` | `t3hold-XYZ` |
| holdout store and backup | `holdout.sqlite`, `holdout_backup.sqlite` | `holdout_XYZ.sqlite`, `holdout_XYZ_backup.sqlite` |
| per-copy logs | `holdout_w<N>.out` | `holdout_XYZ_w<N>.out` |
| cache-ready flag | `state/.holdout_cache_ready` | `state/.holdout_XYZ_cache_ready` |

- **The job name decides what queues.** `--dependency=singleton` holds a job back only while
  another job with the **same name** runs. So two jobs for one domain run one after the other,
  as they must since they share a store, and different domains don't wait for each other. A
  `--job-name=` before `--` still overrides the default.
- **Keep them on different nodes.** At 2000G each job normally gets a node to itself. Two jobs on
  one node would share its node-local cache directory, which is the same-node hazard in
  `docs/LESSONS.md` #24.
- **Shared, and safe to share:** the optimizer's store (each job reads its own copy) and the
  durable cache (additive rsync).

## Where things are written

| | |
|---|---|
| working store | `<dirname of the optimizer's db_path>/holdout.sqlite` — node-local on a cluster |
| durable backup | `<dirname of db_backup_path>/holdout_backup.sqlite`, written after every evaluation |
| with `--settings=settings.local.XYZ.R` | `holdout_XYZ.sqlite` and `holdout_XYZ_backup.sqlite` instead |
| logs | `$LOG_DIR/holdout_w<N>.out` (or `holdout_XYZ_w<N>.out`), one per copy |

**Never the optimizer's store.** Beyond keeping its slice counts clean, this makes it impossible
for a held-out evaluation to become training data if the target domain is ever widened. The
working store is node-local because SQLite's WAL cannot coordinate N writers over NFS —
`open_store()` warns if you get that wrong.

A new job restores the backup first, so resubmitting after a wall-clock kill resumes rather than
repeats.

**One holdout store per domain.** `--report` analyses everything in the store it opens, so results
for two domains in one file would be pooled into a single comparison. A domain file's name
therefore selects its own store. Pass the same `--settings=` to `--report` to read it, and keep
the file name fixed between the evaluation job and the report. `--out=<path>` still overrides
all of this.

## Choosing the trials

Any trial **not** in the pinned universe of the run chosen above; the script hard-stops on one that is, which
is the single error that would silently invalidate the result. Names are accepted and resolved to
ids against the catalogue, because T3 renames trials.

**How many?** The honest standard error of the group difference floors at
`sd_config × √(1/8 + 1/5) ≈ 0.0086` however many trials you run:

| held-out trials | 5 | 10 | 20 | 50 | ∞ |
|---|---|---|---|---|---|
| SE of the group difference | 0.0193 | 0.0149 | **0.0122** | 0.0102 | 0.0086 |

So the smallest detectable difference is **~0.017**, a little over one `sd_config`, and past ~20
trials more trials buy almost nothing — **the binding constraint is the number of configurations,
not trials.** That is also why the default is 8 optimized configurations: the five submissions are
fixed, so `√(1/n + 1/5)` cannot fall below `√(1/5)`, and doubling to 16 improves the SE by ~8%
for twice the compute.

Budget ~13 configurations × N trials at a ~2 h mean: 20 trials is ~260 evaluations, roughly a day
across 12 copies.

## Reading the output

`--report` prints four sections.

1. **Completion** — which cells scored, per configuration and per trial. Read it first: a
   configuration that failed half the grid has its mean computed on whichever half it survived.
2. **The two distributions** — each configuration's mean over the held-out trials, seeds and
   optimized side by side, with a one-sided Wilcoxon and Welch t on those 13 numbers.
3. **The mixed model — this is the result.** `z ~ group + (1|trial_id) + (1|config_hash)`.
   Configuration is *nested* within group, not confounded with it: a **fixed** `config_hash` term
   would alias `group`, a **random** one does not, and it supplies the correct error term.
   Without it, `group` is tested against evaluation-level noise — pseudoreplication that
   understates the SE by about 1.5×.
4. **How much to trust it** — the number of distinct **method skeletons** among the optimized
   configurations. SMAC builds them by crossover and mutation from shared elites, so they are not
   independent draws; if 8 configurations collapse to 3 skeletons, the effective sample size is
   nearer 3. Read sections 2 and 3 against the skeleton count.

The closing **"best optimized vs best seed"** line is descriptive and is *not* the claim: a
maximum over configurations is biased upward, which is the bias this whole exercise exists to
remove. Section 3 is the result.

## Which configurations get tested — a random sample, never the top k

**Ranking cannot break a tie, and on a small domain almost everything ties.** The BLUP shrinks
every configuration onto the grand mean, so hundreds can share one `mean_score`. `head(k)` after
`arrange()` is a *stable* sort and `agg` arrives ordered by `config_hash`, so "the top 8" is
literally the eight alphabetically-first hashes of the tied set — a sample of the hash function,
not of the optimized configurations. Testing those would invalidate the comparison.

So the script builds an **eligible pool** and samples it:

| | |
|---|---|
| `--select=mean` (default) | configurations tied with the best `mean_score`, within `--tie-tol` |
| `--select=ucb` | everything not ruled out at `contender_z` — the set the report counts |
| pool ≤ `--max-equal` (default 20) | **every one of them is tested** |
| pool > `--max-equal` | a **random** `--max-equal` of them, reproducible via `--seed` |

The pool size is printed either way, and it is diagnostic in its own right: above 50 the script
says so, because a pool that large means the optimizer has not separated configurations at all —
evidence about the domain, not a field of good candidates.

`mean` is the default rather than `ucb` because `.contenders()`'s UCB rule is an *elimination*
criterion — "which configurations can I not yet rule out?" — right for deciding where to spend
more evaluations, but it would admit a merely-uncertain configuration here. `--select=ucb`
reproduces the report's set.
