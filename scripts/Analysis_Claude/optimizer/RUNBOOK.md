# Runbook: launching the optimizer

The do-list, in the order the work actually happens. The *why* behind any of it is in
**`README.md`**; this file is what you tick off. `EVALUATION_CHECKLIST.md` is the separate,
slower question of whether the pipeline is *correct* — do that before trusting a run's
results, not before launching one.

The default path here runs **several workers in parallel** against one store. Running a
single worker is the same list minus the parallel-only items, marked **⑂** below.

---

## Quick reference

Once the drill is familiar, the whole thing is:

```bash
export OPTIMIZER_WORK=/workdir/<user>/T3Predictathon2026/scripts/Analysis_Claude/optimizer
cd $OPTIMIZER_WORK
git pull && loadR
source ./optimizer_paths.sh                  # $STOP_FILE, $DB_PATH, $REPORT_PATH, ...
Rscript tests/run_all.R                      # expect 3/3
find $OPTIMIZER_WORK/cache -type f | wc -l   # how many files in the workdir cache
find $OPTIMIZER_HOME/cache -type f | wc -l   # how many files in the home cache
# If the home cache has more then:
rsync -a $OPTIMIZER_HOME/cache/ cache/       # trailing slashes!
# If the workdir cache has more then:
rsync -a cache/ $OPTIMIZER_HOME/cache/       # trailing slashes!
rm -f "$STOP_FILE"
nohup ./run_workers.sh 8 > logs/workers.out 2>&1 &   # 8 workers, 2 BLAS threads each
nohup ./monitor_memory.sh > /dev/null 2>&1 &         # AFTER the workers
tail -f logs/run_w1.out
```

Stopping:

```bash
touch "$STOP_FILE"                           # halts all workers, and the monitor
```

Two things that matter here:

**`source ./optimizer_paths.sh` first.** `OPTIMIZER_HOME` is set in `.Renviron`, which **only
R reads** — in a shell it is empty, so `"$OPTIMIZER_HOME/state/STOP"` expands to
`"/state/STOP"`, `touch` fails with permission denied, and you conclude you stopped a run
that is still going. `optimizer_paths.sh` asks R for the real paths (and so also honours any
`settings.local.R` override of `db_path` / `stop_file`). The scripts do this internally
already; sourcing it is for *your* shell.

**Both launches need `nohup ... &`.** `run_workers.sh` ends in `wait` (so it blocks until the
workers finish, which is what a SLURM batch job needs) and `monitor_memory.sh` loops forever.
Without it the first command occupies the terminal and the rest never run.

---

## 1. Code and environment

- [ ] **Latest optimizer.** `git pull`
- [ ] **Right R.** `which R` — if not the version you want, `loadR` (alias in `~/.bashrc`)
      or `module load R/4.6.1`
- [ ] **`.Renviron` exists and is filled in.** It is gitignored, so it never arrives from
      GitHub — you create it on each machine. R reads it **only at startup**, so restart R
      after editing.
      ```
      T3_USERNAME=<your-t3-username>
      T3_PASSWORD=<your-t3-password>
      OPTIMIZER_HOME=/home/<user>/T3optimizer
      ```
      `OPTIMIZER_HOME` is what turns on remote mode — there is no `remote_server` flag to
      edit any more, it is derived from this variable. Confirm with
      `Rscript -e 'Sys.getenv("OPTIMIZER_HOME")'`.

      > **Never run R here with `--vanilla`.** It implies `--no-environ`, which skips
      > `.Renviron` — you lose `OPTIMIZER_HOME` *and* the T3 credentials, and the run
      > dies at login looking like a missing `.Renviron`.
- [ ] **Packages installed** — this takes a while on a fresh machine. Includes the two
      non-CRAN ones: `BrAPI.R` and `T3BrapiHelpers`. Confirm with
      `Rscript -e 'ip <- installed.packages(); which(ip[,1] == "BrAPI")'`.

### Are R's linear-algebra threads actually being used? (BioHPC)

- [ ] **Measure it, do not infer it:**
      ```bash
      Rscript blas_check.R
      ```
      It prints the BLAS R is linked against, the thread environment, and then times a
      compute-bound matrix multiply, reporting **"cores busy"** = user CPU / elapsed. ~1.0 is
      single-threaded; ~N means N threads.

If it reports single-threaded, there are four causes. **Measured on the BioHPC server
2026-08-01, the answer was #1** — and the diagnostic output looked healthy right up to the
timing line.

1. **R is linked against a SERIAL OpenBLAS build.** The BLAS path says `openblas`, so it
   looks right, but the library has no threading code and **no environment variable can
   change that** — `OPENBLAS_NUM_THREADS=8` and `OMP_NUM_THREADS=8` gave 10.30 s and 10.33 s
   against 10.25 s unset, all at 1.0 cores on a 64-core machine. RHEL ships three builds side
   by side; `ls -la /usr/lib64/libopenblas*` shows them:

   | library | threading |
   |---|---|
   | `libopenblas-r*.so` | **serial** — R links here by default |
   | `libopenblasp-r*.so` | pthread → honours `OPENBLAS_NUM_THREADS` |
   | `libopenblaso-r*.so` | OpenMP → honours `OMP_NUM_THREADS` |

   Switch without root:
   ```bash
   LD_PRELOAD=/usr/lib64/libopenblasp.so.0 OPENBLAS_NUM_THREADS=8 Rscript blas_check.R
   ```

2. **Reference BLAS** — a path containing `libRblas`. The *other* single-threaded case, and it
   looks different in the diagnostic. On BioHPC, R builds from 4.4.3 onward are compiled with
   OpenBLAS; 4.0.5–4.4.2 are not, and the unmodified default is 4.2.3. `module load R/4.6.1`
   puts you past this. Note that **switching R version means reinstalling your packages** —
   they live in `$HOME/R` under a per-version path.

3. **`run_workers.sh` capping it — on purpose.** Its second argument is BLAS threads per
   worker and **defaults to 2**. `./run_workers.sh 8` therefore uses 8 × 2 = 16 cores and
   leaves the rest idle by design. Since this workload is memory-limited, spend the surplus on
   threads: `./run_workers.sh 6 8`. Keep `workers × threads` at or under the core count.
   **Caveat**: the launcher exports `OPENBLAS_NUM_THREADS`, which a *serial* build (#1) ignores
   entirely — it is not a safety net.

4. **There is very little linear algebra to thread.** See below.

> **Temper expectations — this pipeline is not BLAS-bound.** Benchmarked at realistic sizes:
> the VanRaden GEMM (`tcrossprod`, 500 accessions x 130k markers) **0.13 s**;
> `rrBLUP::mixed.solve` on 600 accessions **0.11 s**; whole-population marker QC **~0.5 s**;
> a 130k-marker VCF parse **~9 s**; and the largest single item, `eigen()` on 3000x3000,
> **25 s**. That is well under a minute against a **median evaluation of 29 minutes**, so
> upwards of 95% of the wall time is elsewhere — almost certainly VCF downloads and BrAPI
> calls, which no BLAS setting touches. The pipeline's products are also tall-and-thin and
> memory-bandwidth-bound rather than FLOP-bound, so they thread poorly even when threads are
> available.
>
> Fixing the thread cap is worth doing for correctness, but **it will not convert spare cores
> into throughput**. A profiled real evaluation on the server (`profile_evaluation.R`, trial
> 10938, 803 s) put **`build_kernel` at 6.4 s — 0.8% of the run**, against `.find_related` at
> 54% and `.group_by_panel` at 13%, with CPU/wall = 0.98: serial R computation, not waiting
> and not linear algebra. One busy core per worker is what that looks like in `top`, and it is
> correct behaviour, not a broken BLAS.

- [ ] **Find where the time actually goes** before optimising anything:
      ```bash
      Rscript profile_evaluation.R --trial=<id>
      ```
      It runs one real evaluation with per-subtask timing and reports the split plus CPU/wall.
      Attack the stage it names — on the run above that was `choose_geno_sources`, not
      anything BLAS touches.

## 2. Settings

Put these in **`settings.local.R`** (untracked, gitignored), never in the tracked
`settings.R` — that is what stops `git pull` conflicting on every machine. Copy
`settings.local.R.example` to start.

- [ ] `simulate = FALSE` — otherwise you are optimizing the synthetic world
- [ ] `optimize_scheme` — **one** scheme per run (`"CV0"` or `"CV00"`). They are different
      tasks with different best pipelines; to do both, run twice against the same store.
      (`schemes` is a different setting, used only by the diagnostics.)
- [ ] `max_hours` — a little **under** your reservation, so the loop exits and writes a final
      report before the scheduler kills it
- [ ] `target_domain` — the trial list appropriate for your target trial
- [ ] **Memory budget** — `dosage_budget_bytes` and `dosage_total_budget_bytes`, from the
      table just below
- [ ] **⑂ `db_path` on `/workdir`** — **required** for parallel workers. SQLite's WAL mode is
      what lets them write concurrently, and it cannot work over NFS, which is what `/home`
      usually is.
- [ ] **⑂ `db_backup_path` on `/home`** — the leader copies the store there every
      `db_backup_minutes` (default 30) with `VACUUM INTO`, so a purged `/workdir` costs at
      most one interval.
- [ ] **⑂ Move an existing store to the new `db_path`** — the workers continue from whatever
      is at that path, so it must be the store you care about:
      ```bash
      mkdir -p /workdir/<user>/T3optimizer
      cp $HOME/T3optimizer/state/evals.sqlite /workdir/<user>/T3optimizer/evals.sqlite
      ```

```r
# settings.local.R
settings_override <- list(
  simulate                  = FALSE,
  optimize_scheme           = "CV00",
  max_hours                 = 23.5,
  dosage_budget_bytes       = 16e9,     # keep at whatever the CACHE was built with
  dosage_total_budget_bytes = 18e9,     # 512 GB, 8 workers -- this is what bounds memory
  db_path        = "/workdir/<user>/T3optimizer/evals.sqlite",
  db_backup_path = "/home/<user>/T3optimizer/state/evals_backup.sqlite"
)
```

### The memory budget

`dosage_budget_bytes` caps **one project** at parse time and sets marker density.
`dosage_total_budget_bytes` caps the **sum over every project covering a trial**, which is
what actually bounds a worker's peak.

**Measured anchor (2026-07-29):** one worker at `dosage_budget_bytes = 16e9`, with **no total
cap set**, held **63 GB RSS** on a 128 GB node — half the machine for a single process. That
works out to about **2.3 GB of process memory per GB of resident dosage**, which is what the
table below is derived from. *Uncapped*, 16e9 supports exactly one worker; **with**
`dosage_total_budget_bytes` set it runs eight fine, because the cap is what bounds a worker.

| RAM | workers | `dosage_budget_bytes` | `dosage_total_budget_bytes` |
|---|---|---|---|
| 256 GB | 4 | 8e9 | 19e9 |
| 256 GB | 8 | 4e9 | 9e9 |
| 512 GB | 4 | 16e9 | 35e9 |
| 512 GB | 8 | 8e9 | 18e9 |

The `dosage_budget_bytes` column is for building a **fresh** cache. **If you already have a
cache, do not lower it** — lowering it re-thins nothing (`get_project_dosage` only ever
re-parses to make a cache *denser*), so you would pay for a full re-download to get *worse*
density. Leave it at whatever the cache was built with and let
`dosage_total_budget_bytes` bound the memory; the cap thins at serve time by
column-subsetting the cached matrix, which is free.

> **Existing 16e9 cache, 512 GB, 8 workers:** keep `dosage_budget_bytes = 16e9` and set
> `dosage_total_budget_bytes = 18e9`. ~24 GB of dosage per worker, ~193 GB across all
> eight (38% of the machine). Fine. Rebuilding at 8e9 would save ~3.4 GB per worker —
> not worth the re-download.

The one residual cost of the denser cache is a **transient floor**: `readRDS` loads a cached
file whole before the cap subsets it, so the largest single project (10.2 GB at 16e9, 6.8 GB
at 8e9) is briefly resident on top of everything retained. README has the full table.

`report_memory.R` sizes from **`peak_rss_mb`** (true peak RSS), not `peak_r_mb` (R's heap
peak, which cannot see BLAS/`sommer`/`dist` and understated the real figure by 2.7x). The RSS
peak needs Linux; elsewhere it is `NA` and the report says its numbers understate.

**Start at half the worker count you want**, confirm with `report_memory.R` and
`monitor_memory.sh`, then scale up. The 2.3× multiplier is derived from a single node-hour
and is more likely too low than too high — it assumed every cached project covered that
trial, and if fewer did, the real multiplier is larger.

**Scaling up does not need a restart.** Add workers above the running ids:

```bash
./run_workers.sh 4                              # workers 1-4 (1 is the leader)
# ... a day later, once report_memory.R says there is room ...
OPTIMIZER_FIRST_WORKER=5 ./run_workers.sh 4     # adds workers 5-8
```

Do **not** just run `./run_workers.sh 4` a second time: it always starts at id 1, so you
would get two leaders both rsyncing the cache and backing up the store, ambiguous `worker`
values in the store, and — because the log redirect truncates — the second batch wiping the
first batch's logs. The script now refuses to do this and tells you the right id to use.

Changing `dosage_budget_bytes` changes marker density and therefore scores. Density is not a
config parameter, so old and new rows are not strictly comparable; each row records its
`dosage_budget` so at least the mixing stays visible.

### Checks

These come **after** the settings, not before. `optimizer_settings()` *throws* on an
`optimize_scheme` that is not in `schemes`, and on remote mode with no `OPTIMIZER_HOME`; a
typo'd key warns rather than throws, but surfaces in the same run. `test_subtasks.R` and
`test_sim_loop.R` both call it, so the suite is a settings check as much as a code check —
running it before you have set the settings only tells you the code was fine yesterday.

- [ ] **Verify it all resolved** before committing hours to it:
      ```bash
      Rscript -e 'source("settings.R"); str(optimizer_settings()[c("db_path","simulate",
        "optimize_scheme","max_hours","dosage_budget_bytes","dosage_total_budget_bytes","stop_file")])'
      ```
- [ ] **Tests pass.**
      ```bash
      Rscript tests/run_all.R          # expect "3/3 test files passed"
      Rscript tests/run_all.R --all    # adds the slow sim loop: "4/4"
      ```
      A missing or wrong `.Renviron`, or a typo in `settings.local.R`, usually surfaces here.

## 3. Cache

The cache is every genotyping project's VCF, downloaded once and stored parsed
(`cache/dosage`), plus phenotypes, accessions and the trial catalogue. It lives on the work
disk and is backed up to `$HOME/T3optimizer/cache`. It is large but regenerable.

- [ ] **How full is it?**
      ```bash
      find cache -type f | wc -l
      ls cache/dosage | wc -l          # projects already downloaded + parsed
      du -sh cache/dosage
      ```
- [ ] **Restore it** if you have just copied the repo to a fresh node — the checkout brings
      no cache, because it lives under `$HOME/T3optimizer/cache/`.
      ```bash
      rsync -a $HOME/T3optimizer/cache/ /workdir/<user>/T3Predictathon2026/scripts/Analysis_Claude/optimizer/cache/
      ```
      **The trailing slashes matter.** Without the one on the source, rsync
      *creates a nested directory* instead of merging. This takes a while; drop `-v` unless
      you want the file list. (The leader worker also does this automatically at startup, and
      the other workers wait for it to finish.)
- [ ] **Check `cache/unparseable/`.** A project in there is negative-cached and will never be
      retried. T3 may have fixed the VCF since — delete the entry to force a retry.
- [ ] **⑂ If the cache is cold, run ONE worker first.** Not required, but the first
      evaluations on a cold cache are nearly all download, and eight simultaneous first-time
      VCF downloads is both a lot of load on T3 and the peak memory moment. There is no
      special command — it is the ordinary run:
      ```bash
      Rscript run_optimizer.R </dev/null > logs/warm.out 2>&1 &
      tail -f logs/warm.out                  # watch for downloads to stop dominating
      touch "$STOP_FILE"                     # then stop it
      rm    "$STOP_FILE"                     # and clear the flag
      ```

## 4. Launch

- [ ] **Clear a leftover stop file** — `run_workers.sh` refuses to start while it exists:
      ```bash
      rm -f "$STOP_FILE"
      ```
- [ ] **Start the workers.**
      ```bash
      nohup ./run_workers.sh 8 > logs/workers.out 2>&1 &     # 8 workers, 2 BLAS threads each
      nohup ./run_workers.sh 4 4 > logs/workers.out 2>&1 &   # 4 workers, 4 threads each
      ```
      `run_workers.sh` does not replace `run_optimizer.R` — it launches N copies of it, each
      with `OPTIMIZER_WORKER=i`, a thread limit, and its own `logs/run_w<i>.out`. Keep
      `workers × threads` at or under the core count.

      **Launch it with `nohup ... &`.** The workers themselves are each already `nohup`'d, so
      they survive a hangup either way — but the launcher ends in `wait` (which a SLURM batch
      script needs, so the allocation is not torn down early), and in an interactive shell
      that just blocks your terminal for the whole run.

      *Single worker instead:*
      ```bash
      nohup Rscript run_optimizer.R </dev/null > logs/run.out 2>&1 &
      ```
- [ ] **Start the memory monitor — after the workers**, since it exits when the STOP file
      appears:
      ```bash
      nohup ./monitor_memory.sh > /dev/null 2>&1 &
      ```
      It writes to `logs/memory_<host>.tsv`, so its stdout is only the one startup line.
      Inside tmux/screen a plain `./monitor_memory.sh &` is fine.
- [ ] **Confirm it came up.**
      - No WAL warning in `logs/run_w1.out` — a warning means the store is still on NFS, so
        stop, move `db_path`, and restart
      - `logs/run_w1.out` says `worker=1, leader`
      - `ps -ef | grep run_optimizer` shows N workers
      - rows accruing:
        ```bash
        Rscript -e 'source("settings.R"); source("R/store.R")
                    con <- open_store(optimizer_settings()$db_path); print(n_evals(con))'
        ```

## 5. Monitor

```bash
cat  "$REPORT_PATH"                       # best pipeline so far + learning curve
tail -f logs/run_w1.out                   # live log (one per worker)
Rscript report_memory.R                   # peak per evaluation, workers-that-fit table
Rscript report_timing.R                   # where the wall time goes
sort -k3 -n logs/memory_*.tsv | tail -5   # most memory-intensive moments
find cache -type f | wc -l                # number of files in cache
```

Red flags:

| Symptom | Means | Do |
|---|---|---|
| WAL warning at startup | store is on NFS; concurrent writes unsafe | stop, move `db_path` to `/workdir`, restart |
| `rss_total_mb` approaching RAM | too many workers for the budget | fewer workers, or lower `dosage_total_budget_bytes` |
| only worker 1's rows in the store | the others died | check `logs/run_w*.out` |
| `report.md` not updating | no worker has finished an evaluation since the last write | normal if evaluations are long; check `logs/run_w*.out` for progress |
| every evaluation `infeasible` | usually a data/name bug, not the configs | run the `EVALUATION_CHECKLIST.md` diagnostics |
| `serving 1 marker in N` messages | the aggregate cap is binding | fine if deliberate; else raise `dosage_total_budget_bytes` |
| `peak_rss_mb` ≫ `peak_r_mb` in `report_memory.R` | the cost is in compiled code (GRM/BLAS/sommer), which `dosage_total_budget_bytes` does **not** bound | size from `peak_rss_mb`; to cap it you must bound the kernel stage, not the dosage budget |

## 6. Stop and save

- [ ] **Stop.** One file halts every worker after its current evaluation, and the monitor
      with them:
      ```bash
      touch "$STOP_FILE"
      ```
      Remove it before relaunching.
- [ ] **Save the store — this is the essential one.** The run *is* that file: the incumbent,
      the surrogate's training data, and what to try next are all recomputed from it.
      ```bash
      cp /workdir/<user>/T3optimizer/evals.sqlite $HOME/T3optimizer/state/
      ```
      (With `db_backup_path` set, the leader has been doing this every 30 min anyway.)
- [ ] **Save the cache — optional, large, regenerable.** Keeping it lets the next run skip
      re-downloading. Note this is the **opposite direction** from the restore in step 3:
      ```bash
      rsync -a /workdir/<user>/T3Predictathon2026/scripts/Analysis_Claude/optimizer/cache/ $HOME/T3optimizer/cache/
      ```
- [ ] **Pull to the laptop** if the reservation is ending — FileZilla, or
      `rsync -av biohpc:/home/<user>/T3optimizer/state/ ./state/`.

To **resume** later — same server, new reservation, or laptop — put `evals.sqlite` (and
`cache/` if you kept it) back at the paths the settings name, and launch again.
Configurations already run are skipped by hash, so nothing is repeated.

---

## 7. Running under SLURM (clusters)

Steps 1–6 assume a machine you can `ssh` into and run things on directly. On a cluster you
instead *ask a scheduler* for a machine, and it runs your script when one is free. This
section starts from the beginning; if you already use SLURM daily, skip to "Mapping this
project onto it".

### The mental model

Three kinds of machine, and the distinction matters:

- **Login node** — where you land on `ssh`. Shared by everyone. Edit files, submit jobs, run
  30-second commands. **Do not compute here**; a long R job on a login node is the classic way
  to annoy an administrator.
- **Compute nodes** — where work actually runs. You never `ssh` to one; you ask for one.
- **The scheduler** — SLURM. You describe the resources you want and it queues you until they
  are free.

Five commands cover nearly everything:

```bash
sbatch  myjob.sh          # submit a batch job; prints a job id
squeue  -u $USER          # what is queued/running for me, and why it is waiting
scancel <jobid>           # kill it
sacct   -j <jobid> --format=JobID,State,Elapsed,MaxRSS,ReqMem
                          # what a FINISHED job actually used -- the one to learn from
sinfo   -s                # partitions and how busy they are
```

`sacct`'s `MaxRSS` is the number to check after your first run: it tells you whether the
memory you requested was anywhere near what you needed.

### A minimal job script

```bash
#!/bin/bash
#SBATCH --job-name=t3opt        # shows in squeue
#SBATCH --partition=<partition> # which queue (see `sinfo -s`)
#SBATCH --nodes=1               # ONE node -- see the store constraint below
#SBATCH --ntasks=16             # logical cores on that node
#SBATCH --mem=600G              # TOTAL memory for the job (or --mem-per-cpu)
#SBATCH --time=3-00:00:00       # D-HH:MM:SS. Defaults are usually SHORT -- always set it
#SBATCH --account=<account>     # billing/allocation; required on many clusters
#SBATCH --output=logs/slurm-%j.out   # %j = job id
#SBATCH --error=logs/slurm-%j.err

cd /path/to/optimizer           # .Renviron is read from the WORKING DIRECTORY
module load r/<version>         # or use a container -- see USDA_ARS_SCINET.md

./run_workers.sh 8 2            # 8 workers, 2 BLAS threads each
```

Submit with `sbatch myjob.sh`. Two directives people get wrong:

- **`--time`** — the default is often an hour or two, and the job is killed the moment it
  expires. Always set it explicitly.
- **`--mem` vs `--mem-per-cpu`** — these are alternatives, not additions. Exceeding either is
  an instant kill, not a slowdown.

### Interactive work

For the short read-only scripts here — `peek_failures.R`, `peek_config.R`,
`surrogate_bakeoff.R` — do not write a batch script. Grab a node interactively:

```bash
salloc --nodes=1 --ntasks=4 --mem=32G --time=1:00:00 --account=<account>
# on some clusters salloc drops you onto the node; on others you then need:
#   srun --pty --preserve-env bash
cd /path/to/optimizer && Rscript peek_failures.R
exit                            # releases the allocation -- do not forget
```

### Mapping this project onto it

**One job = one node running N workers. Not one job per worker.**

The reason is the store. SQLite in WAL mode coordinates through an mmap'd `-shm` file, which
network filesystems do not provide — `open_store()` warns about exactly this, and the
consequence of ignoring it is a corrupted store. So every worker sharing an `evals.sqlite`
must be on the same physical node, with `db_path` on **node-local** disk.

That gives the layout:

| setting | where it must point |
|---|---|
| `db_path` | node-local scratch (often `$TMPDIR`) |
| `db_backup_path` | durable project storage — `backup_store()` `VACUUM INTO`s here every 30 min |
| `cache_dir` | node-local scratch if it fits, else project storage |
| `cache_backup_dir` | durable project storage |
| `OPTIMIZER_HOME` | durable project storage (drives the three above) |

`run_workers.sh` ends in `wait`, which is exactly what a batch script needs — without it the
script would exit immediately and SLURM would tear down the allocation with the workers still
running.

**Sizing: this job is memory-bound, not core-bound.** Measured over 121 real evaluations,
peak R heap was a median of **19 GB** and a maximum of **82 GB** per worker — and that is R's
heap peak, an under-estimate of true RSS. So:

```
workers = usable node RAM / 80 GB      (conservative, start here)
workers = usable node RAM / 25 GB      (after report_memory.R confirms headroom)
```

Not `cores / threads`. Start conservative, run `Rscript report_memory.R` after a few hours,
and raise it — the same procedure as §2's memory budget.

### Wall-clock limits and resumability

Clusters cap job length (commonly 1–3 weeks); the optimizer is designed to run indefinitely.
This is not a problem, because the store is resumable: a job that hits its limit loses only
the evaluations in flight. Chain jobs instead —

```bash
sbatch --dependency=afterany:<previous_jobid> myjob.sh
```

— each one picking up from the store. Configurations already run are skipped by hash.

Two things that still work from the **login node** while a job runs:

```bash
tail -f logs/run_w1.out                       # watch a worker
source ./optimizer_paths.sh && touch "$STOP_FILE"   # clean stop; workers finish and exit
```

### Troubleshooting

| symptom | cause |
|---|---|
| Job dies immediately, `sacct` shows `OUT_OF_MEMORY` | `--mem` too low. Check `MaxRSS`, raise it or cut worker count. |
| Job killed at a round time (1:00:00, 2:00:00) | `--time` default. Set it explicitly. |
| `squeue` shows `PD` forever, reason `PartitionConfig` | Asked for more cores/memory/time than the partition allows. |
| `T3 login was REJECTED` / `could not build descriptor` | `.Renviron` not read — the job did not `cd` to the optimizer directory. |
| `could not put the store in WAL mode` warning | `db_path` is on a network filesystem. Move it to node-local disk. |
| Workers all evaluate the same configuration | Pre-0.7.4 code; `git pull`. |

For SciNet specifically — partitions, accounts, containers, storage — see
`USDA_ARS_SCINET.md`.
