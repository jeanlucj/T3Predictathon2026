---
editor_options: 
  markdown: 
    wrap: 72
---

# Runbook: SLURM clusters (and USDA-ARS SciNet)

**Use this when a SLURM decides where and when your work runs.** If you
`ssh` to a machine and use the command line there — use
**`RUNBOOK_INTERACTIVE.md`** instead.

The *why* behind these commands is in **`README.md`**. Sections 1–6
apply to any SLURM cluster; **§7 applies to SciNet (Ceres)**, the
cluster this project targets.

------------------------------------------------------------------------

## A checklist

From no account to a running search. Concrete commands are for Ceres;
§§1–6 give the generic form for any cluster.

**Access**

- [ ] `step ca bootstrap` / `step ssh config` done on the **laptop**
  (§7.6)
- [ ] `ssh <first.last>@ceres.scinet.usda.gov` opens a browser and
  succeeds
- [ ] `my_quotas` — project storage exists and has room
- [ ] `sacctmgr -Pns show user format=account,defaultaccount` — note the
  account

**Code and image**

- [ ] `git clone` the repo into `/project/<account>` (brings
  `container/` too)
- [ ] `salloc -N1 -n8 --mem=32G -t 2:00:00 -A <account>`
- [ ] `module load apptainer` — **not** on `PATH` by default
- [ ] `cd .../optimizer/container && ./build.sh` (20–45 min)
- [ ] `apptainer test optimizer.sif` passes
- [ ] `./run_in_container.sh test`

**The three `.local` files** — these are untracked (§7.5) NOTE: make
your life easier editting these in RStudio with Ceres ondemand

- [ ] `.Renviron` from `.Renviron.example`: `T3_USERNAME`,
  `T3_PASSWORD`, `OPTIMIZER_HOME` (project storage, not node drive)
- [ ] `settings.local.R` optimizer settings tailored to the cluster
  node. From `container/settings.local.R.scinet`
- [ ] `container/submit.local.sh` what SLURM runs. From
  `submit.local.sh.example`

**Prove the environment**

- [ ] `module load r/4.5.3` — check current version of R; only for Rscript checks not in container
- [ ] Settings resolve, inside an allocation:
  `Rscript -e 'source("settings.R"); s <- optimizer_settings(); cat(s$db_path, "\n")'`
  db_path must be in `$TMPDIR` (node-local), not `$OPTIMIZER_HOME` (networked, permanent)  
  check out `s$db_backup_path` and `s$cache_dir` which should be under `$OPTIMIZER_HOME`

- [ ] `container/run_in_container.sh exec prewarm_indices.R --only=projects --limit=5`
  — credentials and outbound HTTPS together

- [ ] `container/run_in_container.sh exec peek_backup.R` — the pre-flight
  report. Read three things off it before submitting:
  - **backup**: rows, and when it was written. `MISSING` means this run
    starts from scratch.
  - **what this run will see**: the row count *after* the domain, scheme
    and build filters. This is the resume-or-restart number, and it can
    be far below the total.
  - **iteration 1**: which phase the first evaluation lands in —
    `seed:PredictionN`, `replicate`, `random_init` or `acquisition`.

  You do **not** restore the store by hand. The leader merges the backup
  into node-local scratch at startup, keyed on config × trial × scheme,
  so nothing already there is lost. `peek_backup.R` reports exactly what
  that merge will insert.

- [ ] The same report's **cache** section lists files per category.
  `proj_acc` / `acc` / `attempted` are the pre-warmed index maps;
  `dosage` is the bulk.

**Submit**

- [ ] `./submit.local.sh --test-only` — shows the allocation without
  queueing
- [ ] Shakeout first:
  `./submit.local.sh --qos=debug --time=00:30:00 --ntasks=4 --mem=200G`
- [ ] `squeue -u $USER` shows it running; `$OPTIMIZER_HOME/logs/slurm-<jid>.out` has no
  WAL warning
- [ ] Only then the real submission: `./submit.local.sh`
- [ ] Queue more: just submit again. `--dependency=singleton` makes them
  run in sequence (`./submit.local.sh; ./submit.local.sh` → a second
  three-week block)

**While it runs** — §6

- [ ] `squeue -u $USER` says `RUNNING`, and
  `find "$OPTIMIZER_HOME/logs"/run_w*.out -mmin -30` lists every worker.
  A green `squeue` alone is not enough: one OOM-killed worker leaves the
  job running with fewer
- [ ] `peek_workers.R`, from a shell on the job's node
  (`srun --overlap --jobid=<jid> --pty bash`) — who is alive, on what,
  and its WORKERS LOST section
- [ ] `report.md`'s `- store backup:` line is fresh and not `** STALE **`

**After the first job finishes**

- [ ] `sacct -j <jid> --format=JobID,State,Elapsed,MaxRSS,ReqMem` — was
  the memory request anywhere near right?
- [ ] `report.md` shows a rising best score and a fresh
  `- store backup:` line
- [ ] `./run_in_container.sh exec report_memory.R` before raising the
  worker count (a bare `Rscript` has no packages on Ceres)

------------------------------------------------------------------------

## 1. SLURM in five commands

Three kinds of machine:

- **Login node** — where you land on `ssh`. Shared by everyone. Edit
  files, submit jobs, run 30-second commands. **Do not compute here**; a
  long R job on a login node is the classic way to annoy an
  administrator.
- **Compute nodes** — where work actually runs. You never `ssh` to one;
  you ask for one.
- **The scheduler** — SLURM. You describe the resources you want and it
  queues you until they are free.

``` bash
sbatch  myjob.sh          # submit a batch job; prints a job id
squeue  -u $USER          # what is queued/running for me, and why it is waiting
scancel <jobid>           # kill it
sacct   -j <jobid> --format=JobID,State,Elapsed,MaxRSS,ReqMem
                          # what a FINISHED job actually used -- the one to learn from
sinfo   -s                # partitions and how busy they are
```

`sacct`'s `MaxRSS` is the number to check after your first run: it tells
you whether the memory you requested was anywhere near what you needed.

## 2. Mapping this project onto SLURM

**One job = one node running N workers. Not one job per worker.**

The reason is the store. SQLite in WAL mode coordinates through an
mmap'd `-shm` file, which network filesystems do not provide —
`open_store()` warns about exactly this, and the consequence of ignoring
it is a corrupted store. So every worker sharing an `evals.sqlite` must
be on the same physical node, with `db_path` on **node-local** disk.

That gives the layout:

| setting | where it must point |
|------------------------------------|------------------------------------|
| `db_path` | node-local scratch (often `$TMPDIR`) |
| `db_backup_path` | durable project storage — `backup_store()` `VACUUM INTO`s here every 30 min |
| `cache_dir` | node-local scratch if it fits, else project storage |
| `cache_backup_dir` | durable project storage |
| `OPTIMIZER_HOME` | durable project storage (drives the three above) |

`run_workers.sh` ends in `wait`, which is exactly what a batch script
needs — without it the script would exit immediately and SLURM would
tear down the allocation with the workers still running.

**Sizing: this job is memory-bound, not core-bound.** Measured over 121
real evaluations, peak R heap was a median of **19 GB** and a maximum of
**82 GB** per worker — and that is R's heap peak, an under-estimate of
true RSS. So:

```         
workers = usable node RAM / 80 GB      (conservative, start here)
workers = usable node RAM / 25 GB      (after report_memory.R confirms headroom)
```

Not `cores / threads`. Start conservative, run
`./run_in_container.sh exec report_memory.R` after a few hours, and
raise it. Because SLURM kills rather than swaps,
err low.

### Wall-clock limits and resumability

Clusters cap job length (commonly 1–3 weeks); the optimizer is designed
to run indefinitely. This is not a problem, because the store is
resumable. Queue several jobs instead, each picking up where the last
left off. Configurations already run are skipped by hash.

**The job is *meant* to be killed at its wall clock.** There is no clean
shutdown and no wall-clock setting, deliberately. What a kill costs:

- **No stored results.** The store is backed up after every evaluation,
  by whichever worker just finished — not by a designated leader that
  might be mid-evaluation. (The cache still syncs on a 30-minute
  interval; an rsync of that tree is not free the way a VACUUM is.)
- **The evaluations in flight**, about 29 minutes each at the median.
  No shutdown scheme saves those: a wall-clock stop is only ever tested
  *between* evaluations, so whatever is running is lost either way.
- **Nothing else.** `backup_store()` writes to a temp file and renames,
  so a killed backup leaves the previous good one; and SQLite in WAL
  mode is crash-safe by design.

A self-imposed stop below `--time` would buy back only that ≤30-minute
window, while making *every* worker idle for the margin — around 11
worker-hours per job at 22 workers. It also cannot tell whether the next
evaluation needs 29 minutes or 33 hours, so it cannot avoid starting one
that will not finish.

### Swapping a running job onto a new build

`scancel` alone forfeits up to 30 minutes of results, because `$TMPDIR`
is wiped at job exit and anything since the last `backup_store()` goes
with it. Stopping through the stop-file first makes every worker take a
final backup, report and cache flush on the way out:

``` bash
source ./optimizer_paths.sh
touch "$STOP_FILE"                                   # stops between evaluations
tail -f "$OPTIMIZER_HOME/logs/run_w1.out"            # wait for "store backed up to"
scancel <jobid>                                      # only for workers still mid-evaluation
git pull
container/run_in_container.sh exec peek_backup.R     # what the next job will start from
./submit.local.sh
```

You do **not** need `rm "$STOP_FILE"` — `run_workers.sh` clears a stale
one before it spawns anything, since the file is on durable storage and
would otherwise halt every subsequent job before its first iteration. A
STOP touched while a run is going still works; the loop tests it every
iteration.

You do **not** need to rebuild the container. `container/*.sif` is
gitignored and the optimizer's R code is bind-mounted, which is the
whole point of that split — only a change to `optimizer.def` itself
calls for `./build.sh`.

The store survives a build bump unless `BUILD_CHANGES` (`R/optimizer.R`)
names it, and it only names changes that alter what a configuration
*computes*. A bump that changes only *selection* leaves every stored row
usable.

> **Expect a burst of revisits after a bump that changes replication.**
> `.trial_target()` reads `n_scored` from the store you are resuming
> from, so a schedule arrives already ramped rather than ramping: around
> 400 scored configurations it asks for 6 configs per trial, around 900
> for 8. Every existing trial below that re-enters the trial backlog, so
> the first hours revisit old trials instead of drawing new ones. That is
> the schedule working, not a stall.

The simplest way is **`--dependency=singleton`** with a fixed
`--job-name`, which is what `submit.local.sh` sets. SLURM then refuses
to start a job while an earlier one of the same name and user is still
going, so submitting three times gives three consecutive blocks:

``` bash
./submit.local.sh; ./submit.local.sh; ./submit.local.sh
```

No job ids to capture, and it doubles as a guarantee that **only one
optimization runs at a time** — worth having, since several concurrent
runs would each back up to the same `db_backup_path` and overwrite one
another.

<details>

<summary>The other dependency types, and why <code>afterany</code> is
the one you would reach for</summary>

| type | releases the successor when the predecessor… |
|------------------------------------|------------------------------------|
| `after` | merely **starts** |
| `afterok` | finishes **successfully** (exit 0) |
| `afternotok` | **fails** |
| `afterany` | ends **for any reason** — success, failure, timeout, cancellation |
| `singleton` | …and every earlier job of the same name and user has ended |

If you ever chain explicitly —
`sbatch --dependency=afterany:<jobid> myjob.sh` — it must be `afterany`,
not `afterok`. **A job that runs to its wall clock ends in SLURM state
`TIMEOUT`, which counts as a failure**, and that is the *normal* ending
here. Under `afterok` the successor would never be released and the
chain would stop silently after one link.

A held job shows in `squeue` as `PD` with reason `Dependency`.

</details>

**The store is put back automatically.** Node-local scratch is wiped
between jobs, so the leader merges `db_backup_path` into the fresh store
at startup — `restore_store_from_backup()` (`R/store.R`), the store's
counterpart to `restore_cache_from_backup()`. Both are **additive**:
whatever is already on the work disk survives. A cell is redundant when
its (config, trial, scheme) is already present, since the pipeline is
deterministic in those three, so re-running the merge inserts nothing.
Its `+N row(s)` line is in `logs/run_w1.out`.

Two things that still work from the **login node** while a job runs:

``` bash
tail -f "$OPTIMIZER_HOME/logs/run_w1.out"           # watch a worker
source ./optimizer_paths.sh && touch "$STOP_FILE"   # clean stop; workers finish and exit
```

## 3. A minimal job script

For a cluster other than Ceres. (On Ceres, use
`container/t3opt_ceres.sh` with `submit.local.sh` — §7.5.)

``` bash
#!/bin/bash
#SBATCH --job-name=t3opt        # shows in squeue
#SBATCH --partition=<partition> # which queue (see `sinfo -s`)
#SBATCH --nodes=1               # ONE node -- see §2
#SBATCH --ntasks=16             # logical cores on that node
#SBATCH --mem=600G              # TOTAL memory for the job (or --mem-per-cpu)
#SBATCH --time=3-00:00:00       # D-HH:MM:SS. Defaults are usually SHORT -- always set it
#SBATCH --account=<account>     # billing/allocation; required on many clusters
#SBATCH --output=logs/slurm-%j.out   # %j = job id
#SBATCH --error=logs/slurm-%j.err

cd /path/to/optimizer           # .Renviron is read from the WORKING DIRECTORY
module load r/<version>         # or use a container -- see §7.4

./run_workers.sh 8 2            # 8 workers, 2 BLAS threads each
```

Two directives people get wrong:

- **`--time`** — the default is often an hour or two, and the job is
  killed the moment it expires. Always set it explicitly.
- **`--mem` vs `--mem-per-cpu`** — these are alternatives, not
  additions. Exceeding either is an instant kill, not a slowdown.

## 4. Interactive work on a cluster

For the short read-only scripts, do not write a batch script. Grab a
node:

``` bash
salloc --nodes=1 --ntasks=4 --mem=32G --time=1:00:00 --account=<account>
# on some clusters salloc drops you onto the node; on others you then need:
#   srun --pty --preserve-env bash
cd /path/to/optimizer && Rscript peek_failures.R
exit                            # releases the allocation -- do not forget
```

On Ceres `salloc` drops you straight onto the node, and everything runs
through the container — see §6.

## 5. Troubleshooting

| symptom | cause |
|------------------------------------|------------------------------------|
| Job dies immediately, `sacct` shows `OUT_OF_MEMORY` | `--mem` too low. Check `MaxRSS`, raise it or cut worker count. |
| Job killed at a round time (1:00:00, 2:00:00) | `--time` default. Set it explicitly. |
| `squeue` shows `PD` forever, reason `PartitionConfig` | Asked for more cores/memory/time than the partition allows. |
| `T3 login was REJECTED` / `could not build descriptor` | `.Renviron` not read — the job did not `cd` to the optimizer directory. |
| `could not put the store in WAL mode` warning | `db_path` is on a network filesystem. Move it to node-local disk. |
| `apptainer: command not found` | no `module load apptainer`. The `container/` scripts do this themselves; an interactive shell does not. |
| Job starts from zero rows | no backup at `db_backup_path`, or the leader's merge failed. `logs/run_w1.out` has the `store restore:` line; `peek_backup.R` says what the backup holds. |
| Workers all evaluate the same configuration | Pre-0.7.4 code; `git pull`. |

## 6. Monitoring and diagnostics

**Monitoring** is "is it healthy right now, and has anything broken since
I submitted" — the first half of this section. **Diagnostics** is "why
did that happen" — the scripts in the second half. The one distinction
that governs everything here:

> **A dead worker is not a dead job.** SLURM kills the whole job only
> when the cgroup memory limit is exceeded. One worker OOM-killed inside
> a job that is otherwise within its limit just leaves you running with
> fewer workers, and `squeue` still says `RUNNING`. So a green `squeue`
> is not evidence that all N workers are alive.

### From the login node, while the job runs

These are the 30-second commands a login node is for. Logs live on
project storage, so they are readable from anywhere:

``` bash
squeue -u $USER                                  # is it running, or still queued and why
squeue -j <jid> -o '%.18i %.9P %.8T %.10M %.6D %R'
sacct  -j <jid> --format=JobID,State,Elapsed,MaxRSS,ReqMem,ExitCode
sstat  -j <jid>.batch --format=JobID,MaxRSS,AveRSS   # a RUNNING job (sacct is for finished ones)
seff   <jid>                                     # efficiency summary; most useful after it ends
scontrol show job <jid>                          # everything, including the node it landed on

tail -f "$OPTIMIZER_HOME/logs/run_w1.out"        # one worker
tail -n 40 "$OPTIMIZER_HOME/logs/slurm-<jid>.out"  # the job banner + worker launch lines
cat "$OPTIMIZER_HOME/state/report.md"            # best pipeline, backlog, backup freshness
```

Sweep every worker at once rather than tailing one:

``` bash
cd "$OPTIMIZER_HOME/logs"
grep -l -E 'FATAL|step error:|NOT being backed up|WAL mode' run_w*.out   # who hit trouble
find run_w*.out -mmin -30                        # whose log is still moving
grep -c 'gave up waiting' run_w*.out             # cache-lock contention
```

A log **missing** from that `find` has been silent for half an hour —
normal for a worker inside a long evaluation (33 h is the observed
maximum, LESSONS #28), suspicious for one that is not. The next section
tells them apart.

### On the compute node — the only place in-flight work is visible

`peek_workers.R` reads the `claims` table, which lives in the **live**
store on node-local `$TMPDIR`. A separate `salloc` gets a different node
and cannot see it, so attach to the running job instead:

``` bash
squeue -u $USER                                  # get <jid>
srun --overlap --jobid=<jid> --pty bash          # a shell ON that node
cd /project/<account>/T3Predictathon2026/scripts/Analysis_Claude/optimizer
./run_in_container.sh exec peek_workers.R
exit                                             # leaves the job running
```

It prints, in one pass: every worker with the trial it is on and how long
it has held it; how many `run_optimizer.R` processes are alive; how long
since each worker's log moved; and **WORKERS LOST MID-EVALUATION**, the
record of workers that died. Run it with no live store and it says so,
with this recipe.

### Was a worker OOM-killed?

Job-level, from anywhere:

``` bash
sacct -j <jid> --format=JobID,JobName,State,ExitCode,MaxRSS,ReqMem
```

`State=OUT_OF_MEMORY`, or `ExitCode 0:125`, means SLURM killed the job.
`MaxRSS` is the whole step — with N workers on one node it is the
aggregate, not the culprit.

Worker-level, which `sacct` cannot see:

- `peek_workers.R`'s **WORKERS LOST** rows — a claim whose owner is gone,
  kept before the claim is reclaimed. This is the only durable trace: a
  killed evaluation writes no store row, so `evals` is censored on
  exactly the outcome you are looking for.
- Fewer processes than workers you asked for.
- A `run_w<i>.out` that stops mid-evaluation with no `optimizer stop` line.
- What it was using when it went: find the pid in the `pids_rss` column
  of `$OPTIMIZER_HOME/logs/memory_<node>.tsv`, sampled every 60 s by
  `monitor_memory.sh`, which the job script starts on the node.

``` bash
grep -o "$PID:[0-9]*" "$OPTIMIZER_HOME/logs/memory_"*.tsv | tail -5
```

`dmesg` carries the kernel's exact `anon-rss` for the victim but is
normally unreadable as a non-root user on Ceres — do not plan on it.

### Is it still making progress?

`report.md`'s best-score line cannot answer this: a run that has stopped
*improving* looks identical to one that has stopped *running*. Count
evaluations per hour instead:

``` bash
./run_in_container.sh exec peek_backup.R         # rows in the durable backup, and its age
sqlite3 "$OPTIMIZER_HOME/state/evals_backup.sqlite" \
  "SELECT substr(ts,1,13) AS hour, COUNT(*) FROM evals GROUP BY hour ORDER BY hour DESC LIMIT 24;"
```

The backup is written after **every** loop iteration, including ones that
store nothing, so its mtime says the workers are alive and nothing about
whether they are producing. The query above is the honest signal, and
`report.md`'s `- newest row:` line is the same measure in words (`- store
backup:` reports only whether backups themselves are keeping up).

### Are the container and the code the same build?

The job banner echoes the build `sed`-scraped from `settings.R` **on the
host**; each worker's first log line prints the build the **container's
R** actually loaded. If they disagree, the bind mount is not delivering
the code you edited:

``` bash
grep '^build     :' "$OPTIMIZER_HOME/logs/slurm-<jid>.out"   # host's settings.R
head -5 "$OPTIMIZER_HOME/logs/run_w1.out" | grep 'optimizer build'  # container's R
```

### The diagnostic scripts

What each script answers is in `README.md`. On a cluster they run inside
the container, from an interactive allocation.

> **A bare `Rscript` will not work.** The module's R has no packages, so
> you get `there is no package called 'here'`. Everything below goes
> through `run_in_container.sh`. The one exception is RStudio via Ceres
> OnDemand, which hands you the module's R and cannot easily reach the
> container — for that, populate a personal library once with
> `Rscript setup_fallback_libs.R` — see `container/FALLBACK_modules.md`.

``` bash
salloc -N1 -n4 --mem=32G -t 1:00:00 -A <account>
module load apptainer
cd /project/<account>/T3Predictathon2026/scripts/Analysis_Claude/optimizer

./run_in_container.sh exec peek_backup.R                      # what a run would start from
./run_in_container.sh exec peek_workers.R                     # who is alive (ON the job's node)
./run_in_container.sh exec peek_failures.R                    # why evals failed
./run_in_container.sh exec peek_config.R --ids=4,23           # what those evals ran
./run_in_container.sh exec report_memory.R                    # how many workers fit
./run_in_container.sh exec report_timing.R                    # where the wall time goes
./run_in_container.sh exec surrogate_bakeoff.R                # surrogate comparison + curve
./run_in_container.sh exec profile_evaluation.R --trial=<id>  # where the time goes
./run_in_container.sh exec prewarm_indices.R --only=projects --limit=5
./run_in_container.sh test                                    # the whole test suite
exit
```

**All are read-only against the store except `prewarm_indices.R`**,
which writes cache files. The read-only ones copy the database *and* its
`-wal`/`-shm` sidecars first, so none of them can corrupt a live store.

**On a cluster the store is the easy half.** Which of these you can run
beside 22 working workers is decided by memory, and by the dosage cache
locks — see the table in `README.md`, *Which of these can run while
workers are evaluating*. Three points specific to here:

- **An `srun --overlap` step joins the job's cgroup.** Its memory counts
  against the same `--mem` the workers are using, and that is sized
  tight: `--mem=1800G` against 22 workers × 82 GB measured worst case =
  1804 GB. SLURM **kills rather than throttles**. Typical headroom is
  wide (22 × the *median* 19 GB is ~418 GB), so this is a judgement about
  free memory now, not a blanket no. Check first:

    ``` bash
    sstat -j <jid>.batch --format=JobID,MaxRSS,AveRSS
    tail -3 "$OPTIMIZER_HOME/logs/memory_"*.tsv   # rss_total_mb, mem_avail_mb
    ```

- **Cap threads yourself.** Only `run_workers.sh` exports the
  `*_NUM_THREADS` variables, and `run_in_container.sh` deliberately does
  not use `--cleanenv`, so a script started from an `srun` shell inherits
  no cap while the workers sit at one thread each:

    ``` bash
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
      ./run_in_container.sh exec surrogate_bakeoff.R --no-curve
    ```

- **`profile_evaluation.R` and `diagnose_failures.R` do not belong on a
  busy node.** Both run a real `run_pipeline` — tens of GB and hours —
  and both take the per-project dosage locks, which `.acquire_lock`
  breaks after `lock_stale_minutes` (90) on the assumption the holder
  died. A worker deep in a large VCF download can lose its lock that way,
  and the raw file is deleted before re-download, so both sides then
  fight over it (LESSONS #24). `diagnose_failures.R` also has **replay on
  by default**; `peek_failures.R` answers the same question in seconds
  with no network.

**Here they read the BACKUP, and they say so.** `db_path` is
`$TMPDIR/t3opt_<user>/evals.sqlite` — node-local scratch belonging to the
job that wrote it — so a diagnostic run from a *separate* allocation
cannot see it. Each script falls back to
`$OPTIMIZER_HOME/state/evals_backup.sqlite`, printing which file it
opened, how old it is, and that anything evaluated since that backup is
not in it. That is expected, not a fault. To read some other file:

``` bash
./run_in_container.sh exec report_memory.R --store=/path/to/evals.sqlite
```

On a **login node** all of them refuse outright with the `salloc` recipe
— `run_in_container.sh` checks before it even looks for the image, so
you get that message rather than "no image, run ./build.sh first".

`peek_failures.R` needs a store to exist. On a fresh cluster install
there is neither a live one nor a backup until a job has run — until
then it stops naming both paths, which is correct behaviour, not a
fault. `peek_backup.R` is the one to use beforehand: comparing the two
*is* its job, so it reports cleanly when there is neither.

`surrogate_bakeoff.R` prints, in its footer, **the smallest difference
the current store can resolve**. A gap smaller than that is not evidence
of no effect; it is too little data.

------------------------------------------------------------------------

# 7. SciNet (Ceres)

*Facts carry their source. Anything not confirmed is in §7.8 rather than
asserted.*

The port is **cheaper than it looks**, for one reason: the optimizer's
state lives entirely in a resumable SQLite store, so a cluster's
wall-clock limit is an inconvenience rather than an obstacle. The work
is configuration and packaging, not redesign.

**The question that could have sunk it is answered: compute nodes DO
have outbound internet.** Confirmed 2026-08-06 by building
`optimizer.def` on a Ceres compute node — that build pulls from Docker
Hub, Posit Package Manager and GitHub, so a successful build is direct
evidence of outbound HTTPS. The pipeline's own calls to
`wheat.triticeaetoolbox.org` should be equally fine; confirm that
specific host with `curl -sI` in the first interactive allocation, since
a proxy could in principle allow package repositories and not arbitrary
hosts.

## 7.1 The machine

**Ceres** — Ames, IA; SLURM; over 13,000 physical cores (26,000
logical), 240 TB total RAM, 5.4 PB VAST storage
([overview](https://scinet.usda.gov/about/compute)).

> SciNet also runs **Atlas** (Starkville, MS). **We are on Ceres**, so
> Atlas is not covered here. Data is not synced between the two, so this
> is a decision rather than a preference.

Node types ([technical
overview](https://scinet.usda.gov/guides/resources/ceres)):

| node             | count | logical cores | RAM      |
|------------------|-------|---------------|----------|
| Intel Xeon 6240R | 72    | 96            | 768 GB   |
| Intel Xeon 6248R | 15    | 96            | 1,536 GB |
| AMD Epyc 9754    | 20    | 256           | 2,304 GB |
| AMD Epyc 9745    | 50    | 256           | 2,304 GB |

Partitions and limits ([partitions &
queues](https://scinet.usda.gov/guides/use/partitions-queues)):

- **`ceres`** — the single community partition since mid-2025, replacing
  short/medium/long/mem. **Default wall time 2 hours; maximum 3 weeks.**
  Per-user ceiling 2,000 cores and 13 TB across all running jobs.
- **QOS `debug`** — 30 minutes, higher priority. Useful for shaking out
  a job script.
- **QOS `long`** — up to 60 days, but capped at 144 cores across your
  running jobs.
- **`scavenger`** — 51 nodes × 96 cores, up to 21 days, lower priority,
  3,000 MB/core default.

Resources ([resource
allocation](https://scinet.usda.gov/guides/use/resource-allocation)):
hyper-threading is on, the minimum allocation is 2 logical cores (odd
requests round up), memory is requested with `--mem-per-cpu`, and **a
job exceeding its memory allocation is killed**, not throttled.

Working habits on Ceres:

- **Your Slurm account** — needed for `-A` on every `sbatch` and
  `salloc`. Look it up yourself, no support ticket required; the second
  column marks your default:

  ``` bash
  sacctmgr -Pns show user format=account,defaultaccount
  ```

  Ceres does provide a default account, so `-A` can often be omitted —
  but it goes in `container/submit.local.sh` explicitly rather than
  inherited, because a default account is the kind of thing that changes
  without notice.

- **`salloc` drops you straight onto the node** — no follow-up
  `srun --pty` needed.

- **`module load apptainer` is required** — it is *not* on `PATH` by
  default. You need this yourself whenever you work interactively; the
  `container/` scripts do it for you via `lib_apptainer.sh`.

Storage: **home is 30 GB per user**; project space is requested
separately ([compute page](https://scinet.usda.gov/about/compute)).
Check with `my_quotas` on Ceres.

## 7.2 Storage, and the one hard constraint

### The store cannot live on network storage

SQLite's WAL mode coordinates through an mmap'd `-shm` file, which
network filesystems do not provide. `open_store()` (`R/store.R`) already
warns when the pragma fails to take, and the consequence of proceeding
is a corrupted store shared between workers.

**Ceres gives us a good resolution.** Every compute node has **1.5 TB of
local NVMe SSD**, reached through `$TMPDIR`. That is genuine node-local
disk, so WAL works, and it is faster than `/project` besides.

> **Where it actually is.** The [scratch-space
> docs](https://scinet.usda.gov/guides/use/scratch) describe a per-job
> `/local/bgfs/<jobid>`. **That path does not exist on the nodes
> checked** (2026-08-10): `$TMPDIR` is `/tmp`, `mount` shows
> `/dev/nvme0n1p1 on /tmp type xfs`, and `df -h /tmp` reports the full
> 1.5 TB. So `/tmp` *is* the local SSD there. Trust `$TMPDIR` rather
> than any documented literal — same lesson as the R version, where the
> published page was stale.
>
> A consequence: `/tmp` is a plain partition, not visibly a per-job
> namespace, so it may be shared with anything else running on the node.
> Hence the directory layout below.

`container/settings.local.R.scinet` is the ready-made config:

``` r
.base <- file.path(Sys.getenv("TMPDIR"), paste0("t3opt_", Sys.info()[["user"]]))

db_path   = file.path(.base, "evals.sqlite")
cache_dir = file.path(.base, "cache")
db_backup_path = file.path(Sys.getenv("OPTIMIZER_HOME"), "state", "evals_backup.sqlite")
```

**It is a template, and must be copied into place to take effect:**

``` bash
cd <repo>/scripts/Analysis_Claude/optimizer
cp container/settings.local.R.scinet settings.local.R
```

`settings.R:38` looks in exactly one place —

``` r
.local_overrides <- function(path = here::here("settings.local.R")) {
  if (!file.exists(path)) return(list())
```

— so a file left in `container/`, or left under its `.scinet` name, is
**never read**, and `.local_overrides()` quietly returns an empty list.
There is no error and no warning: the run simply falls back to the
default `db_path` under `OPTIMIZER_HOME`, which is network storage,
which is the one place the store must not be. Confirm it took, inside an
allocation:

``` bash
Rscript -e 'source("settings.R"); s <- optimizer_settings(); cat(s$db_path, "\n")'
```

That must print a path under `$TMPDIR`. A path under `/project` means
the copy did not land where `here::here()` looks.

**Nothing is scoped by job id, and that is deliberate.** Only one
optimization runs at a time (`--dependency=singleton`, §2), so the paths
only have to be *stable*, not unique. The store is restored from the
durable backup at the start of every job regardless.

To be clear about what is and is not a hazard here, because it is easy
to over-worry:

- **Several workers sharing one store is fully supported.** That is what
  WAL on node-local disk is for, and it is the ordinary N-worker
  configuration. Nor is the *archive* single-purpose:
  `filter_evals_to_universe` / `filter_evals_to_scheme` /
  `filter_evals_to_build` (`R/optimizer.R:32-45`) exist precisely so
  runs with different domains and schemes coexist in one file.
- **What genuinely does not work is two concurrent runs**, because both
  would `VACUUM INTO` the same `db_backup_path` and the later one would
  win, dropping the other's evaluations from the only durable copy.
  `singleton` is what prevents that.

The **cache** is shared across jobs on a node: if `/tmp` survives, the
next job reuses multi-GB downloads rather than re-fetching them; if it
does not, `restore_cache_from_backup()` repopulates it from
`cache_backup_dir` at startup.

Paths are computed at runtime rather than hardcoded, and the file
refuses to load when `SLURM_JOB_ID` is unset — outside a job `$TMPDIR`
is shared storage, where WAL cannot work.

`backup_store()` uses `VACUUM INTO`, which emits a self-contained file
with no WAL sidecar, so the durable copy is directly usable. It runs
after every evaluation, from whichever worker just finished — it costs
about 0.01 s, so there is no interval worth configuring.

**`$TMPDIR` IS ERASED WHEN THE JOB EXITS** — which is why the leader
merges the backup back in at startup (§4).

> Related trap, learned the hard way: never copy `evals.sqlite` without
> its `-wal` and `-shm` sidecars. The main file can legitimately be 0
> bytes with all the data in the WAL. This is also why the restore is a
> SQL merge rather than a file copy, and why the read-only tools go
> through `.copy_store_with_sidecars()`.

### The cache

**Largely dissolved by the 1.5 TB of node-local SSD**: the working cache
lives on `$TMPDIR` alongside the store, where size is a non-issue and
reads are fast. What still has to go somewhere durable is the cache
*backup* (`cache_backup_dir`), since `$TMPDIR` is erased at job end and
`restore_cache_from_backup()` repopulates a fresh node from it.

The cache is dominated by dosage matrices — 601 MB on a laptop with a
small target domain, and the server runs `dosage_budget_bytes = 16e9`.
Point `OPTIMIZER_HOME` at project storage, never home. The same applies
to Apptainer's image cache (§7.7), which is the more common way people
fill their home directory without noticing.

## 7.3 Sizing

Measured over 121 real evaluations: `peak_r_mb` median **19 GB**, max
**82 GB** per worker (R heap peak — an *under*-estimate of RSS).

| node | RAM | workers \@ 82 GB (start here) | workers \@ 25 GB (after measuring) |
|------------------|------------------|------------------|------------------|
| Xeon 6240R | 768 GB | \~9 | \~30 |
| Epyc 9745 | 2,304 GB | \~28 | \~90 |

Start at the conservative column, run `./run_in_container.sh exec
report_memory.R` after a few hours, then raise the count. A single Epyc node at \~28 workers is
already a large multiple of the BioHPC server's throughput. Request
memory to match: `--mem-per-cpu` × cores ≥ workers × 82 GB.

## 7.4 The container

### The one idea to hold onto

**The container holds the environment, not the analysis.**

Inside the image: R 4.5.3, the sixteen CRAN packages, `BrAPI`,
`T3BrapiHelpers`, `rsync`, and the C libraries those need. That is the
whole contents.

Outside the image, on ordinary disk: the optimizer's own R code,
`.Renviron`, `settings.local.R`, the store, the cache, the logs.

```         
    HOST DISK (real directories)              INSIDE THE CONTAINER
    ────────────────────────────              ────────────────────
    /project/.../optimizer/                   R 4.5.3
        R/, tests/, settings.R      ──┐       tidyverse, lme4, sommer,
        run_workers.sh                │       rrBLUP, DBI, RSQLite, ...
        .Renviron                     │       BrAPI, T3BrapiHelpers
        settings.local.R              ├──►    rsync
    /project/.../t3_optimizer/        │
        state/  cache/  logs/       ──┘       (and nothing else)
```

At run time Apptainer makes the host directories on the left *visible*
to the R on the right. Nothing is copied; the container reads and writes
the real files.

Three consequences, and they are the answer to most "do I need to
rebuild?" questions:

- **`git pull` updates the optimizer with no rebuild.** The code was
  never in the image.
- **`.Renviron` and `settings.local.R` can be written or edited at any
  time**, before or after the image exists. They are read from disk at
  each run.
- **You rebuild only to change R or a package version.** Nothing else
  lives in there.

### Glossary

| term | meaning |
|------------------------------------|------------------------------------|
| **image** / **`.sif`** | one big read-only file holding an entire Linux software environment. `optimizer.sif` is \~1–2 GB. Copy it and you have copied R and every package, exactly. |
| **definition file** (`.def`) | the recipe the image is built from. `optimizer.def` is the source; the `.sif` is the compiled output. |
| **`%post`** | the section of the recipe that runs *while building* — apt installs, `install.packages()`. |
| **`%test`** | commands run by `apptainer test optimizer.sif` to prove the image is sound. |
| **`%labels`** | free-text metadata stored inside the image, readable later with `apptainer inspect`. |
| **bind mount** (`--bind`) | making a host directory visible inside the container. Without one, the container cannot see your files at all. |
| **pin** | fixing a version so a later rebuild produces the same thing — an R version, a dated CRAN snapshot, a git commit SHA. |
| **fakeroot** | the mechanism letting an ordinary user build an image without being root. Automatic; nothing to type. |
| **`$PWD`** / **working directory** | the directory a command is run *from*. R reads `.Renviron` from here — see §7.7. |

### Why bother

The dependency list is the awkward part of the port: `BrAPI` and
`T3BrapiHelpers` are GitHub-only, alongside `lme4`, `sommer`, `rrBLUP`,
`RSQLite`, `rpart` and `tidyverse`. Installing that against whatever R
the cluster modules offer, without root, is a multi-day yak-shave. A
container turns it into one build step, done once.

SciNet runs **Apptainer** (formerly Singularity), not Docker — a
root-owned Docker daemon is not acceptable on a shared cluster.
Apptainer imports Docker images directly, and users can build their own
from a definition file on compute nodes without admin privileges
([container
guide](https://scinet.usda.gov/guides/software/singularity)).

### What "pinning" buys, concretely

A score is only reproducible if you can say what produced it. Two
halves:

- `OPTIMIZER_BUILD` (`settings.R`) pins the **code**, and is already
  stamped into every stored evaluation.
- The image pins **everything the code runs on** — R 4.5.3, a dated CRAN
  snapshot, and two exact git commits.

`optimizer.def`'s `%labels` block writes the second half *into* the
image, so months later you can ask a `.sif` what it is without having
the recipe to hand:

``` bash
$ apptainer inspect optimizer.sif
    R.Version             4.5.3
    CRAN.Snapshot         2026-08-01
    T3BrapiHelpers.Sha    6c756462b5a315a992bdd7a26585d912a5452013
    BrAPI.Sha             51d8d450d8ec4f9fc13248165b6382f4a24030b0
```

The build number is **not** a label. Labels have to do with the environment. 
Build numbers have to do with the code that runs in the environment. 
`t3opt_ceres.sh` reads the build out of `settings.R` 
and prints it as the `build :` line of `slurm-<jid>.out`.

Without the snapshot date, "rebuild the image" quietly means "install
whatever CRAN holds today", and the image stops being a pin at all.

**Why commits and not branches.**
`remotes::install_github("owner/repo")` installs whatever is on the
default branch *today*, so two builds a month apart can differ.
Appending `@<sha>` fixes an exact revision. That matters more than usual
here: `R/data_access.R:372` reaches into a private function of
`T3BrapiHelpers` —

``` r
make_row <- getFromNamespace("make_row_from_trial_result", "T3BrapiHelpers")
```

— and a private function carries no promise of staying put. Renamed
upstream, the failure would appear hours into a run, inside a catalogue
fetch. `optimizer.def`'s `%test` asserts that function exists, which
turns a silent time bomb into a build-time error.

**Pin the container's R to match the module.** SciNet has `r/4.5.3`, so
`optimizer.def` uses `rocker/r-ver:4.5.3`. The container route and the
fallback then run the same R, and switching between them cannot by
itself move a result.

### The files

`container/` holds:

| file | what it is |
|------------------------------------|------------------------------------|
| `optimizer.def` | the recipe — R 4.5.3 plus the pinned packages |
| `build.sh` | turns the recipe into `optimizer.sif`; sets the scratch dirs, refuses to run on a login node, prints the size and sha256 checksum |
| `run_in_container.sh` | runs the optimizer, the tests, or any script inside the image, with the binds set correctly |
| `t3opt_ceres.sh` | the SLURM batch script — a template; **do not edit it** |
| `submit.local.sh.example` | copy to `submit.local.sh` (gitignored) and put YOUR account, paths and sizing there |
| `lib_apptainer.sh` | `module load apptainer`, which Ceres needs and other hosts do not |
| `settings.local.R.scinet` | copy to `../settings.local.R` — the `$TMPDIR` store/cache paths |
| `FALLBACK_modules.md` | the no-container route, if this stalls |

Normally `container/` sits **inside** the optimizer directory —
`run_in_container.sh` finds the code by looking at its own parent. If
you keep the image somewhere else, set `OPTIMIZER_REPO` to the optimizer
directory and it will say so clearly rather than failing further
downstream.

**`container/` is tracked**, so `git clone` or `git pull` delivers the
whole recipe to any build host — nothing has to be copied by hand. The
built image is a different matter: `build.sh` writes `optimizer.sif`
into `container/`, and `.gitignore` excludes `container/*.sif` so a 1–2
GB binary is never swept into a commit.

### Building it

**Build on a Ceres compute node.** Nothing gets transferred:
`optimizer.def` bootstraps from `rocker/r-ver:4.5.3` on Docker Hub, and
the recipe itself arrives with `git clone`. Only text crosses the
network; the 1–2 GB image is created where it will be used. `--fakeroot`
is implied automatically for an unprivileged `apptainer build`, so no
admin involvement is needed.

Not on a login node — SciNet documents builds on compute nodes, and
`build.sh` refuses anyway.

``` bash
# 1. Code and recipe together -- container/ is tracked, so one clone brings both. The clone is
#    needed regardless because step 3 runs the test suite, which lives in the repo and is
#    deliberately NOT in the image.
ssh <first.last>@ceres.scinet.usda.gov
cd /project/<account>
git clone https://github.com/jeanlucj/T3Predictathon2026.git

# 2. Build, inside an allocation
salloc -N1 -n8 --mem=32G -t 2:00:00 -A <account>
module load apptainer                      # NOT loaded by default. build.sh does this for
                                           # itself, but in a SUBSHELL -- so your interactive
                                           # shell still needs it for the `apptainer test`
                                           # below and for anything else you run by hand.
cd /project/<account>/T3Predictathon2026/scripts/Analysis_Claude/optimizer/container
./build.sh                                 # 20-45 min
apptainer test optimizer.sif               # asserts R 4.5.3, every package, rsync, and the
                                           # private T3BrapiHelpers function

# 3. Prove it. No credentials or network needed -- the test suite stubs both.
./run_in_container.sh test                 # tests/run_all.R INSIDE the image: 32/8007/317

# 4. Needs a real .Renviron: proves BrAPI authenticates and fetches THROUGH the container,
#    which the test suite cannot show because it stubs the network.
./run_in_container.sh exec prewarm_indices.R --only=projects --limit=5
```

Note the order: step 4 is the first thing that needs credentials, and
the first thing that needs outbound HTTPS. Everything before it works on
an empty account.

Nothing about `.Renviron` or `settings.local.R` has to be settled before
building — they are never in the image, and on SciNet they need
different values anyway (`OPTIMIZER_HOME` on project storage, `db_path`
on `$TMPDIR`). Write them once you get there.

## 7.5 Never edit a tracked file on the server

Anything specific to your account or machine — the Slurm account,
`/project` paths, worker count, `--mem`, `--time` — goes in a
**`.local`** file, which git ignores. Editing the tracked file instead
means a merge conflict on every `git pull`, forever.

| you edit (gitignored) | copied from (tracked) | carries |
|------------------------|------------------------|------------------------|
| `container/submit.local.sh` | `submit.local.sh.example` | Slurm account, node sizing |
| `settings.local.R` | `container/settings.local.R.scinet` | `dosage_budget_bytes` |
| `.Renviron` | `.Renviron.example` | T3 credentials, **`OPTIMIZER_HOME`** |

This is the arrangement `settings.R:31-32` already describes for
`settings.local.R`: the tracked file stays identical on every machine,
so `git pull` can never conflict on it.

### `.local` files hold values. Tracked files hold logic.

The rule has a sharp edge worth stating, because it has already cost a
debugging session. **A `.local` file is gitignored, so `git pull` never
updates it.** Any *mechanism* copied into one is frozen at the moment you
copied it, and later fixes silently do not reach you — with no conflict,
no warning, and a run that works but behaves like an older version.

That is exactly what happened on the first shakeout: `run_workers.sh` is
tracked, so its new log location arrived, while the `submit.local.sh` in
use predated the `--output`/`--chdir` flags and quietly wrote
`slurm-<jid>.out` to the repo instead.

So both templates now carry values only, and source their mechanism from
tracked files:

| your file sets | mechanism lives in |
|---|---|
| `submit.local.sh` — account, sizes | `container/lib_submit.sh` → `submit_optimizer()` |
| `settings.local.R` — `dosage_budget_bytes` | `settings.R` → `cluster_scratch_paths()` |

Fix a flag or a path rule in the tracked file, `git pull`, and your
unchanged local file picks it up. If you find yourself pasting logic into
a `.local` file, that is the signal it belongs in a tracked one.

As a backstop, `t3opt_ceres.sh` warns when a job was submitted by a
`submit.local.sh` predating `lib_submit.sh` — the check that would have
caught the log problem at submission rather than afterwards.

**`.Renviron` is the single source for `OPTIMIZER_HOME`** —
`submit.local.sh` reads it from there rather than restating it, and
derives `REPO` from its own location. So the only things you edit in
`submit.local.sh` are the account and the sizing.

That direction is forced, not a preference: **R overrides the inherited
environment from `.Renviron` at startup** (see the note in
`optimizer_paths.sh`). A second value set in the shell would lose to it
inside the container, so the shell would bind and restore one directory
while R wrote to another — and the symptom looks like a permissions
fault, not a config mismatch. `t3opt_ceres.sh` therefore *requires*
`OPTIMIZER_HOME` and `REPO` in its environment rather than defaulting
them; a default there could contradict `.Renviron` just as easily.

One constraint follows: the value in `.Renviron` must be a **literal
absolute path**. R would expand `${HOME}` or `~` and the shell read does
not, so `submit.local.sh` refuses those rather than letting the two
disagree.

It works for the batch script because **`sbatch` command-line flags
override `#SBATCH` directives**. The directives inside `t3opt_ceres.sh`
are documented defaults; `submit.local.sh` passes the real values on the
command line and the tracked file is never touched:

``` bash
cd container
cp submit.local.sh.example submit.local.sh && chmod +x submit.local.sh
# edit submit.local.sh, then
./submit.local.sh --test-only     # shows the allocation it would get, queues nothing
./submit.local.sh
```

Extra arguments pass straight through to `sbatch`, so queueing more work
and one-off overrides need no new files either:

``` bash
./submit.local.sh; ./submit.local.sh          # two blocks, run in sequence via singleton
./submit.local.sh --job-name=t3opt-debug \
  --qos=debug --time=00:30:00 --ntasks=4 --mem=200G   # 30-min shakeout, alongside
```

Submitting `t3opt_ceres.sh` directly, with its placeholders intact,
stops immediately and points you here rather than failing hours later on
a path that does not exist.

## 7.6 Logging in — SmallStep, and why `Permission denied (publickey)` is normal

SciNet does **not** use ordinary SSH keys. It uses **SmallStepCLI**,
which mints a **short-lived SSH certificate (\~16 hours)** after you
authenticate through USDA eAuthentication **in a web browser** ([SSH
access](https://scinet.usda.gov/guides/access/ssh-login)). Nothing else
is accepted, so `Permission denied (publickey)` is the standard symptom
of *no valid certificate* — never set up, or simply expired since
yesterday.

Two things that produce that identical message and send people hunting
in the wrong place:

- **Usernames are `first.last`**, not a BioHPC-style short name.
- **The browser step cannot happen on a headless machine.** So you
  cannot `ssh` or `scp` into SciNet *from* another cluster — it has to
  originate somewhere with a browser, i.e. the laptop. This is the main
  reason the container is built on Ceres rather than moved there.

Two hostnames, and they are different:

| purpose                        | host                        |
|--------------------------------|-----------------------------|
| log in (mints the certificate) | `ceres.scinet.usda.gov`     |
| bulk data transfer             | `ceres-dtn.scinet.usda.gov` |

Log in to the first once; the certificate then covers the second for
\~16 hours. For anything large,
[Globus](https://scinet.usda.gov/guides/data/transfer/globus) is
web-authenticated and avoids SSH entirely — worth setting up if you ever
need to move the cache or a batch of results, though the current design
keeps everything either generated in place or small.

If plain `ssh` from the laptop still fails after `step ca bootstrap` and
`step ssh config`, that is account provisioning rather than local
configuration: mail
[scinet_vrsc\@usda.gov](mailto:scinet_vrsc@usda.gov){.email}.

## 7.7 Container traps for this project

`container/run_in_container.sh` and `container/t3opt_ceres.sh` already
handle all of these. They are written out because when one goes wrong
the symptom is usually silence, not an error.

**1. The container cannot see a directory unless you bind it.**

A container starts with almost nothing visible. `--bind /a/b:/a/b` makes
the host's `/a/b` appear at `/a/b` inside. Your home directory is
mounted automatically; **`/project` is not**, and neither is `/workdir`.
Forget the bind and R reports a missing file for a file that plainly
exists — which is confusing until you know the rule.

**2. Bind each path *onto itself*.**

`--bind` can rename: `--bind /project/me/opt:/opt/optimizer` puts the
host's `/project/me/opt` at `/opt/optimizer` inside. Do not do that
here.

`settings.R` reads `OPTIMIZER_HOME` and *computes* `db_path`,
`cache_dir`, `log_dir` and the rest from it. Those computed paths are
host paths. If `OPTIMIZER_HOME` says `/project/me/t3_optimizer` but you
mounted that at `/opt/whatever`, the computed paths point at nothing
inside the container — and `settings.R` does not stop. It falls back to
local mode and writes a fresh, empty store in the working directory. The
run looks fine and produces nothing you want. Same path in and out, and
the problem cannot arise.

**3. R reads `.Renviron` from the working directory — and Apptainer no
longer sets one.**

R looks for `.Renviron` in the directory the command was launched
*from*, and if there is none, in your home directory. It does not search
parent directories, and it does not merge the two. The T3 credentials
live in `optimizer/.Renviron`, so **the optimizer must be launched from
the optimizer directory** or `T3_USERNAME` and `T3_PASSWORD` are simply
absent.

Older Singularity started you in whatever directory you were in on the
host. **Apptainer 1.1 and later no longer do that** — you land somewhere
else entirely, so the working directory must be set explicitly:

``` bash
# REPO is just a shorthand for the optimizer directory -- the folder holding settings.R,
# run_workers.sh and .Renviron. Set it to wherever you cloned the repo:
REPO=/project/<account>/T3Predictathon2026/scripts/Analysis_Claude/optimizer

apptainer exec --bind "$REPO:$REPO" optimizer.sif \
  bash -c "cd '$REPO' && ./run_workers.sh 8 2"
#          ^^^^^^^^^^ this is the part that matters
```

You do not normally type this. `container/run_in_container.sh` works
`REPO` out for itself (`REPO="${OPTIMIZER_REPO:-$(dirname "$HERE")}"` —
the parent of `container/`, since that folder normally lives inside the
optimizer directory), and `container/t3opt_ceres.sh` reads it from
`submit.local.sh`. The form is spelled out here because it is what those
scripts are doing on your behalf.

`run_workers.sh` opens with `cd "$(dirname "$0")"` — shell for "change
into the directory this script lives in". So launching *through that
script* puts you in the optimizer directory no matter where it was
called from, and `.Renviron` is found. That is why the invocation goes
through `run_workers.sh` rather than calling `Rscript` directly.

Getting this wrong produces **no error**. The BrAPI connection falls
back to anonymous, and every catalogue fetch returns 401 while the run
continues. `diagnose_failures.R` refuses to start outside the optimizer
directory for exactly this reason.

**4. Never put `.Renviron` inside the image.**

It would work, and it would mean your T3 password is inside a 1–2 GB
file you copy between machines and might share. Credentials stay on disk
and get bound in.

**5. `rsync` has to be inside the image.**

`sync_cache_to_backup()` shells out to the `rsync` command, and
`rocker/r-ver` does not ship one. Missing, it prints
`rsync not on PATH -- skipping` and continues with **no cache backup at
all** — days of downloads left on a wipeable disk. `optimizer.def`
installs it and `%test` asserts it is there.

**6. Point Apptainer's own scratch at `$TMPDIR` before building.**

Building unpacks gigabytes of intermediate layers, into `~/.apptainer`
by default — which fills a 30 GB home quickly. `APPTAINER_CACHEDIR` and
`APPTAINER_TMPDIR` redirect that. Ceres sets them for SLURM jobs; BioHPC
does not. `build.sh` sets them itself, so this only bites if you run
`apptainer build` by hand.

**What the container does not solve.** It fixes *what software is
available*, and nothing else. The store still cannot live on network
storage, and it still has to be merged back into `$TMPDIR` at the start
of every job. Neither has anything to do with packaging, and neither
gets easier by containerising.

## 7.8 Open questions and the support email

*Numbered as originally asked. Items 1, 2, 5 and 7 were answered without
support, so the email below asks only the remaining three and its
numbering does not match this list.*

- [x] **1. Outbound HTTPS from compute nodes** — **YES.**
  `optimizer.def` built on a Ceres compute node on 2026-08-06, reaching
  Docker Hub, PPM and GitHub. The cache-staging redesign is not needed.
  Still confirm `wheat.triticeaetoolbox.org` specifically.
- [x] **2. `$TMPDIR` size / lifetime** — **1.5 TB of local NVMe**, and
  on the nodes checked (2026-08-10) `$TMPDIR` is `/tmp`, *not* the
  documented `/local/bgfs/<jobid>`. Cache *and* store both go there,
  under a per-user directory. Whether it survives between jobs is still
  unconfirmed, so the store is job-scoped and the leader merges it back
  from backup at job start either way.
- [ ] **3. `/project` filesystem type** — confirms the WAL constraint
  and the backup pattern. Lower stakes now that both hot paths sit on
  `$TMPDIR`; it only governs the backups.
- [ ] **4. Project storage granted** — quota and path.
- [x] **5. Account string** — **self-serve, no support needed:**
  `sacctmgr -Pns show user format=account,defaultaccount` lists your
  accounts and marks the default. Put it in `container/submit.local.sh`,
  not in the tracked batch script.
- [ ] **6. Long-job policy** — chained jobs vs QOS `long`.
- [x] **7. R module** — **answered: SciNet has `r/4.5.3`.** (The public
  preinstalled-software page lists only up to 4.4.1 and is stale;
  `module spider r` is the authority.) The container is pinned to
  `rocker/r-ver:4.5.3` to match, so the two routes run the same R.
  Fallback procedure in `container/FALLBACK_modules.md`.

**1, 2, 5 and 7 are closed.** Only 3, 4 and 6 remain — filesystem type,
storage quota and long-job policy — and none of them blocks starting.
The technical risk in this port is retired.

### Email to SciNet support — ready to send

Send to [**scinet_vrsc\@usda.gov**](mailto:scinet_vrsc@usda.gov){.email}
— the VRSC (Virtual Research Support Core), which handles batch scripts,
storage, software and performance questions ([contact
page](https://scinet.usda.gov/support/contact)), including login issues.
Questions of general interest can also go to the [SCINet
Forum](https://forum.scinet.usda.gov/); a formal software-installation
request has its own [request
process](https://scinet.usda.gov/support/request).

------------------------------------------------------------------------

**Subject:** New workload — long-running memory-heavy R job: storage and
long-job policy

Hello,

I am preparing to move a long-running analysis onto Ceres and would like
to check a few things before I request an allocation. The workload
profile:

- **R**, single node, 8–30 concurrent worker processes.
- **Memory-heavy**: measured 19 GB median and 80 GB peak *per worker*,
  so I expect to reserve most of one large-memory node.
- **Long-running**: weeks to months, but fully checkpointed — it can be
  stopped and resumed from a SQLite file with no loss.
- It queries an **external web API** (a public plant-breeding database)
  throughout, and caches the downloads to disk.

I have already built my Apptainer container on a compute node and
confirmed it works, so questions about software and outbound network
access are settled. What remains:

1.  **Is `/project` storage a network filesystem?** I ask because SQLite
    in WAL mode coordinates through an mmap'd shared-memory file and is
    documented as unsafe on NFS-style filesystems. My plan is to keep
    the live database on `$TMPDIR` (node-local SSD) and back it up to
    `/project` every 30 minutes, restoring from that backup at the start
    of each job. If you have a recommended pattern for SQLite or other
    file-locking workloads on SCINet storage, I would rather follow it
    than invent one.

2.  **Project storage** — what is the process to request it, and what
    quota is reasonable to ask for? I expect 100–500 GB of cached
    genomic data. (I understand home is 30 GB, which is not enough.)

3.  **Long jobs.** The `ceres` partition allows 3 weeks and QOS `long`
    allows 60 days but caps at 144 cores. For a months-long job my plan
    is to submit successive 3-week jobs sharing a job name with
    `--dependency=singleton`, so they run one after another; the work is
    checkpointed, so a job ending at its wall clock loses almost
    nothing. Is that the preferred approach, or would QOS `long` be more
    appropriate? Is there any guidance on repeatedly occupying most of a
    large-memory node for an extended period?

Thank you, Jean-Luc Jannink

------------------------------------------------------------------------
