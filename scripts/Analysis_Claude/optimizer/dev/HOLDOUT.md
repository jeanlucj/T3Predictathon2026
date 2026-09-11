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

`--report` needs neither the optimizer's store nor T3 — the group of each row is recorded in the
holdout store — so the analysis runs anywhere, any time after the evaluations.

Everything runs **from the optimizer root**: `.Renviron`, and so the T3 credentials, is read from
the working directory only.

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
a holdout started now tests *"the top 8 as of this moment"*, which may not be the top 8 when the
optimizer stops. For the manuscript claim, run it after the optimizer finishes. If you run it
early anyway, record the timestamp — `--report` prints the configurations it used.

## Where things are written

| | |
|---|---|
| working store | `<dirname of the optimizer's db_path>/holdout.sqlite` — node-local on a cluster |
| durable backup | `<dirname of db_backup_path>/holdout_backup.sqlite`, written after every evaluation |
| logs | `$LOG_DIR/holdout_w<N>.out`, one per copy |

**Never the optimizer's store.** Beyond keeping its slice counts clean, this makes it impossible
for a held-out evaluation to become training data if the target domain is ever widened. The
working store is node-local because SQLite's WAL cannot coordinate N writers over NFS —
`open_store()` warns if you get that wrong.

A new job restores the backup first, so resubmitting after a wall-clock kill resumes rather than
repeats.

## Choosing the trials

Any trial **not** in the optimizer's pinned universe; the script hard-stops on one that is, which
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

Selection defaults to the top-k by **`mean_score`**, not by UCB. `.contenders()`'s UCB rule is an
*elimination* criterion — "which configurations can I not yet rule out?" — which is right for
deciding where to spend more evaluations but would admit a merely-uncertain configuration here.
`--select=ucb` reproduces the report's set.
