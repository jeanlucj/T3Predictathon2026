# Genomic-prediction pipeline optimizer — user guide

A background process that recombines the five Predictathon submissions -- and
literature methods they missed -- into a genomic-prediction pipeline that beats
any single submission at predicting grain yield in a **randomly chosen T3 trial**.

This file is for **running** the optimizer. If you want to understand, validate,
or modify the code, see **`EVALUATION.md`** (evaluation & validation runbook) and
**`DESIGN.md`** (architecture). `BACKGROUND.md` explains why it is built this way.

## Start here, then pick a runbook

Read this file for what the optimizer is, how it is configured, and the facts that hold
everywhere. Then follow **one** runbook end to end — each has a complete checklist, and you
should not need the other:

| your situation | runbook |
|---|---|
| You can `ssh` to a machine and start things on it yourself — a laptop, or a Cornell BioHPC server you have reserved | **`RUNBOOK_INTERACTIVE.md`** |
| A scheduler decides when and where your work runs — any SLURM cluster, including USDA-ARS SciNet | **`RUNBOOK_SLURM.md`** (SciNet/Ceres is §7) |

The difference that matters is not the size of the machine; it is **who decides when your job
runs**. If you type the command that starts the workers, you want the interactive runbook.

## What it does, in one paragraph

Every genomic-prediction algorithm performs six subtasks (select training trials,
preprocess phenotypes, select genotyping data, build a relationship matrix/kernel,
train a model, predict). The optimizer catalogues how each submission does each
subtask and treats the choice of method + parameters as a *configuration* to
optimize. A SMAC-style optimizer -- a random-forest surrogate that predicts a
configuration's accuracy, plus evolutionary crossover/mutation of the best
configurations -- searches that space. Each configuration is evaluated on a *fresh
random T3 trial* and scored by the correlation between predicted and observed
BLUEs, so the search is driven toward pipelines that generalize across the whole
database rather than toward the nine Predictathon trials specifically.

## Install once

```r
# Required (most already present in this repo's environment):
#   tidyverse, here, DBI, RSQLite, jsonlite, rpart, rrBLUP, MASS
#   BrAPI (github TriticeaeToolbox/BrAPI.R), T3BrapiHelpers
#   (VCF parsing is base-R streaming now -- vcfR is no longer needed.)
# Recommended for the G+E / RKHS model variants and shrinkage BLUEs:
install.packages(c("sommer", "lme4"))
```

No Bayesian-optimization or genetic-algorithm package is needed; the optimizer is
self-contained in `R/`.

## Configure (`settings.R`)

All knobs live in `optimizer_settings()` in `settings.R`. The ones you will set:

| Setting | What it controls |
|---|---|
| `simulate` | `TRUE` = fast offline synthetic world (for confirming the install); `FALSE` = the real T3 pipeline. |
| `target_domain` | Restrict the random focal trials to `programs` / `years` / `locations` (each `NULL` = no constraint). This tailors the pipeline to a subpopulation; leave all `NULL` for a global all-rounder. |
| `optimize_scheme` | The **one** CV scheme this run optimizes (`"CV0"` or `"CV00"`). CV0 and CV00 are distinct prediction tasks with potentially different optimal pipelines, so a run targets one; to optimize both, run twice (see below). Must be a member of `schemes`. |
| `schemes` | The CV schemes the **diagnostics** (`check_canaries`, `sweep_rich_trials`) sanity-check — kept at `c("CV0","CV00")` so bug checks cover both; this does **not** control what the optimizer targets (that's `optimize_scheme`). |
| `focal_trait`, `focal_trait_db_id` | The trait to optimize prediction of (defaults to grain yield, T3 variable id `84527`). |
| `max_iters`, `max_hours` | Stop after this many evaluations / this much wall-clock, whichever first. |
| `db_path`, `cache_dir`, `report_path`, `stop_file`, `log_dir` | Where state, cache, the report, the stop-file, and logs go. Splitting these across disks is environment-specific — see your runbook. |

**Machine-specific values never go in `settings.R`.** Put them in `settings.local.R`, which is
gitignored, so the tracked file stays identical everywhere and `git pull` never conflicts on
it. The same pattern covers credentials (`.Renviron`) and, on a cluster, the submission
settings (`container/submit.local.sh`). Each has a tracked `.example` to copy.

## Run

**1. Confirm the install (seconds, offline):**

```bash
cd scripts/Analysis_Claude/optimizer
Rscript tests/run_all.R        # expect "3/3 test files passed"
```

**2. Launch the real optimizer.** Set `simulate = FALSE` in `settings.R`, then:

```bash
nohup Rscript run_optimizer.R > logs/run.out 2>&1 &
```

In real mode the optimizer first runs a **startup self-check** (the "canary"
check) and prints the result near the top of `logs/run.out`. If you see
**`CANARY ALARM`** there, the build is failing to read data it should be able to
read -- **do not trust that run**; see `EVALUATION.md` ("Validation & debugging tooling")
before continuing. `canaries OK` means it is reading data correctly.

**3. Monitor:**

```bash
tail -f logs/run.out      # live iteration log
cat state/report.md       # current best pipeline + diagnostics (see below)
```

**4. Stop cleanly** (finishes the current evaluation, writes a final report):

```bash
touch state/STOP
```

It is fully resumable: kill it any time and re-launch -- it continues from the
SQLite store at `state/evals.sqlite`. Phenotypes and genotypes are cached under
`cache/`, so revisiting a trial is cheap.

**Optimizing both CV schemes.** A run optimizes one scheme (`optimize_scheme`), because
CV0 and CV00 are distinct tasks that can want different pipelines. To do both, run the
optimizer **twice** against the same store and cache -- e.g. set `optimize_scheme = "CV0"`
(or override it in `settings.local.R`) and run to completion, then set `"CV00"` and run
again. `state/evals.sqlite` is a shared archive: `choose_config()` reads only the rows for
the current `optimize_scheme` (and `target_domain`), so the two runs don't interfere, while
the second run reuses all genotypes the first cached under `cache/`. `state/report.md`
records which scheme it optimized. (`schemes` is separate -- it only lists which schemes the
diagnostics sanity-check.)

## Reading the output (`state/report.md`)

- **Incumbent** -- the current best pipeline (its six subtask choices) and its
  mean score, flagged when it beats the best submission.
- **Subtask-method importance** -- for each subtask, the mean score of every
  method, so you can see *which submitted ideas actually win* and which lose.
- **Failure log** -- how attempts broke down by status (`ok` / `infeasible` /
  `suspect` / `error`), the dominant reasons, and the failure rate of each method.
- **⚠ Suspected bugs** -- failures whose data funnel looks like a *bug* (real data
  not visible) rather than genuine infeasibility. A non-zero count here, or a
  `CANARY ALARM` in the log, means the build needs a developer's attention
  (`EVALUATION.md`).
- **Running best** -- the learning curve over evaluation order.

## Running several workers

One evaluation takes minutes to hours, so throughput comes from running several at
once.

**`run_workers.sh` does not replace `run_optimizer.R` — it launches N copies of it.**
That is the whole of the relationship. There is no separate parallel code path: each
worker runs the same sequential loop as always (choose a configuration, sample a
trial, evaluate, store), against one shared `evals.sqlite` and one shared `cache/`.
What the launcher adds is three environment settings per copy — `OPTIMIZER_WORKER=i`,
the BLAS thread count, and a separate log file — plus a stagger between starts. You
could do it by hand:

```bash
OPTIMIZER_WORKER=1 nohup Rscript run_optimizer.R </dev/null > logs/run_w1.out 2>&1 &
OPTIMIZER_WORKER=2 nohup Rscript run_optimizer.R </dev/null > logs/run_w2.out 2>&1 &
```

Because they share the store, they also share progress: each worker reads every
worker's results when it fits its surrogate, so N workers explore one search, not N
independent ones.

#### Setup, in order

The order matters — steps 1 and 2 must be done **before** any worker starts.

**1. Put the store on local disk** (`settings.local.R`, untracked):

```r
settings_override <- list(
  simulate       = FALSE,
  db_path        = "/workdir/<user>/optimizer/evals.sqlite",   # local: WAL works
  db_backup_path = "/home/<user>/t3_optimizer/state/evals_backup.sqlite"
)
```

This is the one step you cannot skip. `open_store` puts SQLite in WAL mode, which is
what lets workers write concurrently, and WAL *cannot* work over NFS — it coordinates
through an mmap'd `-shm` file. A home directory on a cluster is usually NFS. The
store is copied to `db_backup_path` every `db_backup_minutes` using `VACUUM INTO`
(safe on a live database) by whichever worker reaches the check first, so a wiped
`/workdir` costs at most one interval.

If you have an existing store on `/home`, copy it to the new path first — the workers
continue from whatever is there:

```bash
mkdir -p /workdir/<user>/optimizer
cp ~/t3_optimizer/state/evals.sqlite /workdir/<user>/optimizer/evals.sqlite
```

**2. Set the memory budget** for the number of workers you want — see "How many
workers?" below.

**3. Verify the configuration resolves as you expect**, before committing hours to it:

```r
Rscript -e 'source("settings.R"); s <- optimizer_settings(); str(s[c("db_path","simulate","dosage_budget_bytes","dosage_total_budget_bytes","stop_file")])'
```

**4. If the cache is COLD, run one worker on its own for a while first.**

"The cache" is `cache/` — every genotyping project's VCF, downloaded once and stored
as a parsed dosage matrix (`cache/dosage`), plus phenotypes, accessions and the trial
catalogue. Once a project is in there, no evaluation ever downloads it again. "Cold"
means those files do not exist yet; "warm" means they do.

Check before deciding — this step is usually unnecessary:

```bash
ls cache/dosage | wc -l      # projects already downloaded and parsed
du -sh cache/dosage
```

If that count covers most of the projects your `target_domain` touches, skip to step
5. If it is near zero (a fresh node, a `/workdir` that was purged, a new domain), then
run a single worker first. There is **no special command for this** — it is the
ordinary optimizer run:

```bash
Rscript run_optimizer.R </dev/null > logs/warm.out 2>&1 &
tail -f logs/warm.out
```

Let it work through a few evaluations — you are watching for
`Load project dosages` and download lines to stop dominating the log, and for
`cache/dosage` to stop growing quickly. Then stop it, clear the flag, and go to step 5:

```bash
source ./optimizer_paths.sh
touch "$STOP_FILE"      # it exits after the current evaluation
rm    "$STOP_FILE"      # or run_workers.sh will refuse to start
```

Nothing is lost by skipping this: the results go to the same store either way, and the
per-project lock means two workers never download the same VCF. The reason to do it is
that the first evaluations on a cold cache are the worst moment to have eight
processes running — eight first-time VCF downloads at once is both a lot of load on
the T3 server (see LESSONS.md #10) and the peak memory moment, since parsing a large
panel is itself memory-heavy. Paying that once, serially, is cheaper than tuning for
it.

**5. Launch:**

```bash
nohup ./run_workers.sh 8 > logs/workers.out 2>&1 &    # 8 workers, 2 BLAS threads each
nohup ./run_workers.sh 4 4 > logs/workers.out 2>&1 &  # 4 workers, 4 threads each
nohup ./monitor_memory.sh > /dev/null 2>&1 &          # watch what it actually costs
```

Launch both with `nohup ... &`. The workers are each `nohup`'d individually, so they
survive a hangup regardless — but `run_workers.sh` ends in `wait` (needed under SLURM, so
the batch job does not exit and tear down the allocation), which in an interactive shell
blocks the terminal for the whole run. `monitor_memory.sh` loops until stopped.

**6. Stop them** — all at once, with the stop file. Get the path from
`optimizer_paths.sh` rather than building it yourself: `OPTIMIZER_HOME` is set in
`.Renviron`, which **only R reads**, so in a shell `"$OPTIMIZER_HOME/state/STOP"` expands to
`"/state/STOP"` — `touch` fails and you believe you stopped a run that is still going.

```bash
source ./optimizer_paths.sh
touch "$STOP_FILE"
```

Each worker finishes its current evaluation and exits cleanly. Remove the file before
launching again — `run_workers.sh` refuses to start if it is still there.

#### What makes it safe

- **Worker 1 is the leader, and by now that means one thing only: it restores the cache at
  startup** and signals when it is ready, so the others do not read a half-populated tree.
  Everything else — `state/report.md`, the store backup, the cache backup — is done by *every*
  worker, deliberately. All three can only happen *between* evaluations, so assigning any of
  them to worker 1 alone made the real interval `max(setting, worker 1's current evaluation)`:
  observed at a day and a half against a 30-minute setting. See LESSONS #25.
- **The cache backup is still one-at-a-time**, just not one-*worker*. It is an rsync over
  thousands of files, so unlike the store's millisecond copy it would be wasteful in parallel;
  whichever worker finds the shared timestamp stale claims it first, then transfers.
  It also restores the cache from backup at startup and signals when it is done, which
  the other workers wait for. So worker 1 should be running — the others tolerate its
  absence but nothing will write a report.
- **The shared cache is locked per project.** Downloading a project's VCF clears any
  existing copy first, so two workers on one project would delete each other's
  in-flight multi-GB file. A worker that loses the lock waits for the winner's result
  rather than repeating the download.
- **Cache writes are atomic**, so a reader never sees a half-written file.

A worker that dies takes only its current evaluation with it: everything else is
already in the store. Restart it with the same `OPTIMIZER_WORKER` number.

#### How many workers?

Memory, not cores, is the binding constraint: R's heap is per-process, so N workers
cost N times the memory of one. Get the number from measurement, not arithmetic:

```bash
nohup ./monitor_memory.sh > /dev/null 2>&1 &   # attaches to a RUNNING job; no restart needed
Rscript report_memory.R                        # per-evaluation peaks + "how many workers fit"
```

The setting that bounds a worker's peak is **`dosage_total_budget_bytes`**, not
`dosage_budget_bytes` — the latter caps one project at parse time, while the pipeline
holds *every* project covering a trial at once.

**Measured anchor (2026-07-29):** one worker at `dosage_budget_bytes = 16e9`, with **no total
cap set**, held **63 GB RSS** on a 128 GB node. That is ~2.3 GB of process memory per GB of
resident dosage, which the table below is derived from. *Uncapped*, 16e9 therefore supports
exactly one worker — but with `dosage_total_budget_bytes` set it is fine at eight, because
the cap is what bounds a worker, not the per-project budget.

| RAM | workers | `dosage_budget_bytes` | `dosage_total_budget_bytes` |
|---|---|---|---|
| 256 GB | 4 | 8e9 | 19e9 |
| 256 GB | 8 | 4e9 | 9e9 |
| 512 GB | 4 | 16e9 | 35e9 |
| 512 GB | 8 | 8e9 | 18e9 |

`report_memory.R` sizes from **`peak_rss_mb`** — the evaluation's true peak RSS. Do not size
from `peak_r_mb`: that is R's own heap peak, and `gc()` cannot see BLAS, `sommer` or `dist`
allocations, so it understated the real requirement by **2.7x** on the run measured above (33.6
GB of heap against 89.4 GB of RSS). `peak_rss_mb` needs Linux (it resets the kernel's `VmHWM`
via `/proc/self/clear_refs`); elsewhere it is `NA` and the report says loudly that its figures
understate.

A large `peak_rss_mb / peak_r_mb` ratio means the cost is in compiled code — the GRM copies in
`em_combine`, the RKHS distance matrices, a `sommer` fit — and **`dosage_total_budget_bytes`
does not bound any of that**; it caps dosage bytes only.

**Start at half the worker count you want**, confirm with `report_memory.R`, then scale up:
the 2.3x multiplier comes from a single node-hour and is more likely too low than too high.
Keep `workers x threads` at or under the core count; `run_workers.sh` sets the thread
count for every BLAS R might be linked against.

The `dosage_budget_bytes` column applies when you are building a **fresh** cache. If you
already have one, see the next section — you probably should not change it.

## Reusing a cache built at a different budget

`dosage_budget_bytes` is a **parse-time** setting: it fixes a project's marker density when
that project is first downloaded. Lowering it later does **not** re-thin an existing cache —
`get_project_dosage` only ever re-parses to make a cache *denser*, never coarser.

So if your cache was built at 16e9, **leave `dosage_budget_bytes = 16e9` and control memory
with `dosage_total_budget_bytes` alone.** The cap thins at serve time by column-subsetting
the cached matrix: no re-download, and you keep the denser markers. Lowering the budget would
take effect only after deleting `cache/dosage` and re-downloading everything — many GB of T3
traffic, for *worse* density.

What the higher budget does still cost is a **transient floor**. `serve_cached()` calls
`readRDS()` on the whole cached file and only then subsets it, so the largest single cached
project is briefly resident in full. Across the T3 wheat archive, under an 18e9 cap:

| cache built at | retained after the cap | + largest file (transient) | dosage peak per worker |
|---|---|---|---|
| 16e9 | 14.0 GB | 10.2 GB | **24.2 GB** |
| 8e9 | 9.8 GB | 6.8 GB | **16.6 GB** |

**Worked example — 512 GB, 8 workers, an existing 16e9 cache:** keep
`dosage_budget_bytes = 16e9`, set `dosage_total_budget_bytes = 18e9`. That is ~24 GB of
dosage per worker (~193 GB across all eight, 38% of the machine), or ~330 GB of RSS at the
measured 2.3x — within the 70% guideline either way. Rebuilding the cache at 8e9 would save
about 3.4 GB per worker, which is not worth the re-download.

Scaling up needs no restart — add workers *above* the running ids:

```bash
./run_workers.sh 4                            # workers 1-4
OPTIMIZER_FIRST_WORKER=5 ./run_workers.sh 4   # later: adds workers 5-8
```

Running `./run_workers.sh 4` twice instead would start a second worker 1 — two leaders, and
a truncated `logs/run_w1.out` under the first batch. The script refuses and names the id to
use.

Note that changing `dosage_budget_bytes` changes marker density and therefore scores,
and density is not a configuration parameter — so rows made under different budgets
are not comparable. Every row records the budget it ran under (`dosage_budget`), and
`report_memory.R` warns when a store mixes them.

## Diagnostics — what each script answers

Seven standalone scripts. What each one tells you is the same everywhere; **how you invoke it
differs**, so the commands live in your runbook rather than here.

| script | question it answers | invocation |
|---|---|---|
| `profile_evaluation.R` | where does an evaluation's time actually go? | interactive §7 · SLURM §6 |
| `surrogate_bakeoff.R` | which surrogate design predicts held-out configurations best, and does it improve as evaluations accumulate? | interactive §7 · SLURM §6 |
| `peek_failures.R` | why did the non-`ok` evaluations fail? | interactive §7 · SLURM §6 |
| `peek_config.R` | what configuration did a given eval run? | interactive §7 · SLURM §6 |
| `prewarm_indices.R` | fill the local wizards so the pipeline needs no BrAPI call | interactive §7 · SLURM §6 |
| `blas_check.R` | is R's linear algebra actually threaded? | interactive §1 |
| `report_memory.R` | how much memory did evaluations really use, and how many workers fit? | interactive §5 · SLURM §2 |

**On a cluster these run inside the container**, which is where the packages are:

``` bash
./container/run_in_container.sh exec peek_failures.R     # any script, plus its arguments
./container/run_in_container.sh shell                    # interactive R inside the image
```

A bare `Rscript` uses the cluster module's R, which has **no packages installed** — the symptom
is `there is no package called 'here'`. Either go through the container as above, or populate a
personal library once with `Rscript setup_fallback_libs.R` (see
`container/FALLBACK_modules.md`); the container is the better answer unless you need RStudio.

Two properties hold in both environments:

- **All are read-only against the store except `prewarm_indices.R`**, which writes cache files.
  The read-only ones copy the database *and* its `-wal`/`-shm` sidecars before reading, so they
  are safe to run against a live store while workers are going.
- **`peek_failures.R` and `surrogate_bakeoff.R` need a store that already has rows.** On a
  fresh install they stop with `no store at <db_path>`, which is correct behaviour rather than
  a fault.

## Where to go next

- **`EVALUATION.md`** -- the evaluation runbook: step through each module one at a
  time (`arm_evaluation`), run/extend the tests, trace one evaluation, and use the
  validation tooling (the canary oracle, `diagnose_trial`).
- **`DESIGN.md`** -- the architecture: every function, and where outputs go.
- **`BACKGROUND.md`** -- the statistical and data-management challenges and how
  each is addressed.
- **`RUNBOOK_INTERACTIVE.md`** / **`RUNBOOK_SLURM.md`** -- the do-lists for actually launching
  a run, one per environment.
- **`EVALUATION_CHECKLIST.md`** -- whether the pipeline is *correct*. A slower question than
  whether it launches; do this before trusting a run's results.
