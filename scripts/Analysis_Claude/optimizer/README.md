# Genomic-prediction pipeline optimizer — user guide

A background process that recombines the five Predictathon submissions -- and
literature methods they missed -- into a genomic-prediction pipeline that beats
any single submission at predicting grain yield in a **randomly chosen T3 trial**.

This file is for **running** the optimizer. If you want to understand, validate,
or modify the code, see **`EVALUATION.md`** (evaluation & validation runbook) and
**`DESIGN.md`** (architecture). `BACKGROUND.md` explains why it is built this way.

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
| `focal_trait`, `focal_trait_db_id` | The trait to optimize prediction of (defaults to grain yield, T3 variable id `84527`). |
| `max_iters`, `max_hours` | Stop after this many evaluations / this much wall-clock, whichever first. |
| `db_path`, `cache_dir`, `report_path`, `stop_file`, `log_dir` | Where state, cache, the report, the stop-file, and logs go (see the BioHPC section for splitting these across disks). |

## Run

**1. Confirm the install (seconds, offline):**

```bash
cd scripts/Analysis_Claude/optimizer
Rscript tests/run_all.R        # expect "2/2 test files passed"
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

## Running on a remote server (e.g. Cornell BioHPC)

The optimizer was built to survive being stopped and resumed, which is exactly
what you need on a rented server. The design feature that makes remote operation
easy: **all durable state is two directories** -- `state/` (the SQLite store and
report) and `cache/` (downloaded phenotypes/genotypes). Nothing else needs to be
preserved.

### What is different about the setup

1. **Install the packages once, to a user library.** A server R has no packages
   in your account yet. Load/choose an R (on BioHPC, e.g. `module avail R` then
   `module load R/4.x`, or use an R under `/programs`), then:
   ```r
   # one-time, in an R session on the server
   install.packages(c("tidyverse","DBI","RSQLite","jsonlite","rpart",
                       "rrBLUP","MASS","sommer","lme4","remotes"))
   remotes::install_github("TriticeaeToolbox/BrAPI.R")   # BrAPI
   # T3BrapiHelpers (snake_case fns) from github jeanlucj/T3_brapi_helpers
   ```
   Set a persistent library so reinstalls survive reservations:
   `export R_LIBS_USER=$HOME/Rlibs` (and `mkdir -p $HOME/Rlibs`) in your
   `~/.bashrc`.

2. **Outbound internet.** Real mode downloads from T3/Wheat (BrAPI + archived
   VCFs). Interactive reserved machines normally have internet; some clusters
   block it on batch/compute nodes. If your compute node has no internet,
   **pre-populate `cache/` on an internet-enabled (login/reservation) node first**
   -- run a handful of iterations there so the trials you'll use are cached -- then
   run the heavy loop offline against that cache.

3. **Put the working directories on the right disk -- via the environment, not by
   editing `settings.R`.** On BioHPC you typically reserve a machine and `ssh` in;
   `/workdir/<user>` is fast *local* scratch but is **not guaranteed to persist after
   your reservation**, while `/home` is persistent network storage (smaller, slower).
   The right split is the large, regenerable **cache on local scratch** and the small,
   irreplaceable **state on persistent `/home`**.

   You get that split by setting **one environment variable** -- no code edit:
   ```bash
   # in the optimizer folder's .Renviron (alongside T3_USERNAME / T3_PASSWORD):
   OPTIMIZER_PATH=/home/<user>/t3_optimizer
   ```
   `settings.R` reads this and switches to **remote mode automatically**: `state/`,
   `logs/`, `report.md`, and the `STOP` flag go under `$OPTIMIZER_PATH` (persistent),
   while `cache/` stays on the work disk and is backed up to `$OPTIMIZER_PATH/cache`
   (see "Automatic cache backup" below). A laptop with no `OPTIMIZER_PATH` stays in
   local mode. **`remote_server` is derived from the environment** (`.detect_remote_server()`
   in `settings.R`):
   - it is **TRUE when `OPTIMIZER_PATH` is set**, else FALSE (auto-detect), and
   - `OPTIMIZER_REMOTE=true|false` (also in `.Renviron` or the shell) forces it either
     way if you ever need to.

   **Why env-driven matters:** the tracked `settings.R` is *identical* on your laptop
   and the server -- each machine picks its mode from its own environment -- so
   `git pull` on the server **never conflicts on `settings.R`**. (Before this, you had
   to edit `remote_server <- TRUE` on the server, and every upstream change to
   `settings.R` collided on pull.) `.Renviron` is gitignored and is read only at R
   **startup**, so after creating/editing it, **restart R**; verify with
   `Sys.getenv("OPTIMIZER_PATH")`.

4. **Change any OTHER setting per machine via an untracked `settings.local.R`.** The
   environment handles only the directory split; for anything else you want different
   on the server -- `simulate = FALSE`, a bigger `dosage_budget_bytes`, a `max_hours`
   under your reservation limit, `run_startup_canary = FALSE`, custom paths -- **do not
   edit the tracked `settings.R`** (that reintroduces the `git pull` conflict). Instead
   copy `settings.local.R.example` to `settings.local.R` (gitignored) and set just your
   deltas:
   ```r
   # settings.local.R  (untracked)
   settings_override <- list(
     simulate            = FALSE,
     dosage_budget_bytes = 16e9,
     max_hours           = 23.5,
     run_startup_canary  = FALSE
   )
   ```
   `optimizer_settings()` layers these on top of its defaults, so upstream owns the
   tracked defaults and your server owns only its overrides. An unknown key warns
   (typo guard). Overrides are top-level settings; to change a list-valued one (e.g.
   `target_domain`) give the whole list. Unlike `.Renviron`, `settings.local.R` is read
   every time `optimizer_settings()` is called, so no restart is needed after editing it.

### What to save when the reservation ends, and how to resume

- **Must save:** `state/evals.sqlite` (or wherever you pointed `db_path`). This
  single file *is* the optimizer's state -- the incumbent, the surrogate's
  training data, and what to try next are all recomputed from it on startup. Put
  it on durable storage by setting `OPTIMIZER_PATH=/home/<user>/t3_optimizer` in
  `.Renviron` (see item 3 above) -- remote mode turns on automatically and sends
  `state/` and `logs/` to `$OPTIMIZER_PATH` while the cache stays on the work disk.
- **Worth saving:** `cache/` -- large but regenerable; keeping it lets a resumed
  run skip re-downloading the trials it already touched. `state/report.md`/`logs/`
  are convenience only.

**Automatic cache backup.** In remote mode (i.e. `OPTIMIZER_PATH` set), `run_optimizer()`
backs the cache up to `$OPTIMIZER_PATH/cache` on its own: it restores from there at startup if
the work cache is empty (fresh node), rsyncs additively every `cache_sync_minutes`
(default 30) during the run, and flushes once more on a clean stop or an R-level error.
So a graceful stop, an error, or a hit budget loses nothing.

**Surviving an abrupt kill.** An in-process backup cannot run if R itself is
`SIGKILL`ed (OOM-killer, node reboot, `scancel`). For that, run an external rsync
loop in a separate `tmux`/`nohup` shell so at most one interval is ever at risk:
```bash
while true; do
  rsync -a --exclude 'raw_project/' /workdir/<user>/.../cache/ "$OPTIMIZER_PATH/cache/"
  sleep 1800
done
```
Note this only protects the *cache*; `evals.sqlite` is already on durable `$OPTIMIZER_PATH`
storage, so optimizer **progress** is never lost to a kill -- only some re-downloadable cache.

Manual copy off the machine before your time is up (belt-and-suspenders):
```bash
rsync -av /workdir/<user>/optimizer_cache/  $HOME/optimizer_cache/      # optional, large
cp /workdir/<user>/.../state/evals.sqlite    $HOME/optimizer_state/      # essential
# or pull to your laptop:  rsync -av biohpc:/home/<user>/optimizer_state/ ./state/
```

To **resume** -- on the same server, a new reservation, or your laptop -- put
`evals.sqlite` (and `cache/` if you kept it) back at the paths `settings.R` names,
then launch `run_optimizer.R` again. It reads the store and continues; already-run
seeds and configurations are skipped by hash, so no work is repeated.

### Using SLURM

Yes. Because the loop is checkpointed and resumable, it fits SLURM cleanly. Two
patterns:

- **One long job.** Request a node for N hours and run the loop, but set
  `max_hours` in `settings.R` to a bit **less** than the SLURM `--time` so the
  loop exits and writes a final report *before* SLURM kills the job:
  ```bash
  #!/bin/bash
  #SBATCH --job-name=gp_optimizer
  #SBATCH --time=24:00:00
  #SBATCH --cpus-per-task=8          # threads for BLAS / sommer
  #SBATCH --mem=32G
  module load R/4.x
  export R_LIBS_USER=$HOME/Rlibs
  cd /path/to/scripts/Analysis_Claude/optimizer
  Rscript run_optimizer.R </dev/null   # set max_hours = 23.5 in settings.R
  ```
  (Redirect stdin from `/dev/null` so an archived-VCF download can never block on
  an interactive prompt.)
- **Chained / requeued jobs.** Submit short jobs that each run for a while and
  exit; resubmit (`sbatch` again, a job-dependency chain, or `--requeue` on
  preemption) to keep extending the *same* store. This rides under tight
  wall-time limits and survives preemption, since each job just resumes from
  `evals.sqlite`.

SLURM caveats specific to this optimizer:
- **One writer per store.** SQLite allows a single writer; do **not** point a
  SLURM *array* of parallel jobs at the same `evals.sqlite` (lock contention /
  corruption). The loop is currently sequential, so use one job per store. To
  use many cores, give one job several `--cpus-per-task` for threaded linear
  algebra rather than launching concurrent writers. (Parallel evaluation with
  separate per-worker stores merged afterward is a noted future enhancement.)
- **Stop cleanly with `max_hours`, not just the STOP file.** The STOP file is for
  interactive use; for batch, `max_hours < --time` is what guarantees a graceful
  final checkpoint.
- **Confirm internet on the partition you use**, or pre-cache (point 2 above).

The exact module names, queue/partition names, and `/workdir` vs `/home` policy
vary -- check the current Cornell BioHPC documentation for specifics.

## Where to go next

- **`EVALUATION.md`** -- the evaluation runbook: step through each module one at a
  time (`arm_evaluation`), run/extend the tests, trace one evaluation, and use the
  validation tooling (the canary oracle, `diagnose_trial`).
- **`DESIGN.md`** -- the architecture: every function, and where outputs go.
- **`BACKGROUND.md`** -- the statistical and data-management challenges and how
  each is addressed.
