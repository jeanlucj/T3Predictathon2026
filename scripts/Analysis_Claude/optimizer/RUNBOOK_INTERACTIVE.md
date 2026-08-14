---
editor_options: 
  markdown: 
    wrap: 72
---

# Runbook: interactive (laptop, or a machine you rent)

**Use this when you can `ssh` to a machine and run things on it directly** — your laptop, or a
reserved Cornell BioHPC server. If a scheduler decides when your work runs, use
**`RUNBOOK_SLURM.md`** instead.

The *why* behind any of it is in **`README.md`**; this file is what you tick off.
`EVALUATION_CHECKLIST.md` is the separate, slower question of whether the pipeline is
*correct* — do that before trusting a run's results, not before launching one.

The default path here runs **several workers in parallel** against one store. A single worker
is the same list minus the parallel-only items, marked **⑂**.

------------------------------------------------------------------------

## The complete checklist

Everything, in order, from a fresh machine to a running search. Details for each are in the
numbered section that follows.

**Environment**

- [ ] `git pull`
- [ ] Right R on `PATH` (`which R`; `loadR` or `module load R/4.6.1`)
- [ ] `.Renviron` exists, with `T3_USERNAME`, `T3_PASSWORD`, `OPTIMIZER_HOME`
- [ ] Packages installed, including `BrAPI` and `T3BrapiHelpers`
- [ ] `source ./optimizer_paths.sh` so your shell knows `$STOP_FILE`, `$DB_PATH`, …

**Settings**

- [ ] `settings.local.R` created from `settings.local.R.example` (never edit `settings.R`)
- [ ] `simulate = FALSE`
- [ ] `optimize_scheme` set to one scheme
- [ ] `target_domain` matches your target trial
- [ ] Memory budget: `dosage_budget_bytes`, `dosage_total_budget_bytes`
- [ ] **⑂** `db_path` on `/workdir` (local disk — WAL cannot work over NFS)
- [ ] **⑂** `db_backup_path` on `/home`
- [ ] Settings resolve as expected (`str(optimizer_settings()[...])`)
- [ ] `Rscript tests/run_all.R` → 3/3

**Cache and store**

- [ ] Existing store copied to the new `db_path`, if you are resuming
- [ ] Cache restored from `$HOME/T3optimizer/cache/` (mind the trailing slashes)
- [ ] `cache/unparseable/` reviewed
- [ ] **⑂** One worker run first if the cache is cold

**Launch**

- [ ] `nohup ./run_workers.sh 8 > logs/workers.out 2>&1 &`
- [ ] `nohup ./monitor_memory.sh > /dev/null 2>&1 &` — *after* the workers
- [ ] No WAL warning in `logs/run_w1.out`; it says `worker=1, leader`
- [ ] Rows accruing in the store

**While it runs**

- [ ] `cat "$REPORT_PATH"` shows a rising best score
- [ ] `Rscript report_memory.R` says there is headroom before you add workers

**Stopping**

- [ ] `touch "$STOP_FILE"`
- [ ] Store copied to `$HOME/T3optimizer/state/`
- [ ] Cache rsynced back to `$HOME/T3optimizer/cache/` (optional)
- [ ] Pulled to the laptop if the reservation is ending

------------------------------------------------------------------------

## Quick reference

Once the drill is familiar, the whole thing is:

``` bash
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
nohup ./run_workers.sh 8 > logs/workers.out 2>&1 &   # 8 workers, 2 BLAS threads each
nohup ./monitor_memory.sh > /dev/null 2>&1 &         # AFTER the workers
tail -f logs/run_w1.out
```

Stopping:

``` bash
touch "$STOP_FILE"                           # halts all workers, and the monitor
```

Which pid is which worker, what each is working on, and what has died:

``` bash
Rscript peek_workers.R
```

(The old `/proc/$p/environ` loop only worked on Linux and only told you the mapping.
`peek_workers.R` also names the trial each worker is on and keeps a record of workers lost
mid-evaluation. Worker 1 is no longer a special case — see §5.)

Two things that matter here:

**`source ./optimizer_paths.sh` first.** `OPTIMIZER_HOME` is set in `.Renviron`, which **only R
reads** — in a shell it is empty, so `"$OPTIMIZER_HOME/state/STOP"` expands to `"/state/STOP"`,
`touch` fails with permission denied, and you conclude you stopped a run that is still going.
`optimizer_paths.sh` asks R for the real paths (and so also honours any `settings.local.R`
override of `db_path` / `stop_file`). The scripts do this internally already; sourcing it is
for *your* shell.

**Both launches need `nohup ... &`.** `run_workers.sh` ends in `wait` (so it blocks until the
workers finish, which is what a SLURM batch job needs) and `monitor_memory.sh` loops forever.
Without it the first command occupies the terminal and the rest never run.

------------------------------------------------------------------------

## 1. Code and environment

- [ ] **Latest optimizer.** `git pull`

- [ ] **Right R.** `which R` — if not the version you want, `loadR` (alias in `~/.bashrc`) or
  `module load R/4.6.1`

- [ ] **`.Renviron` exists and is filled in.** It is gitignored, so it never arrives from
  GitHub — you create it on each machine. R reads it **only at startup**, so restart R after
  editing.

    ``` bash
    T3_USERNAME=<your-t3-username>
    T3_PASSWORD=<your-t3-password>
    OPTIMIZER_HOME=/home/<user>/T3optimizer
    ```

    `OPTIMIZER_HOME` is what turns on remote mode — there is no `remote_server` flag to edit
    any more, it is derived from this variable. Confirm with
    `Rscript -e 'Sys.getenv("OPTIMIZER_HOME")'`.

    > **Never run R here with `--vanilla`.** It implies `--no-environ`, which skips
    > `.Renviron` — you lose `OPTIMIZER_HOME` *and* the T3 credentials, and the run dies at
    > login looking like a missing `.Renviron`.

- [ ] **Packages installed** — this takes a while on a fresh machine. Includes the two
  non-CRAN ones: `BrAPI.R` and `T3BrapiHelpers`. Confirm with:

    ``` bash
    Rscript -e 'ip <- installed.packages(); which(ip[,1] == "BrAPI")'
    ```

### Are R's linear-algebra threads actually being used? (BioHPC)

- [ ] **Measure it, do not infer it:**

    ``` bash
    Rscript blas_check.R
    ```

    It prints the BLAS R is linked against, the thread environment, and then times a
    compute-bound matrix multiply, reporting **"cores busy"** = user CPU / elapsed. ~1.0 is
    single-threaded; ~N means N threads.

If it reports single-threaded, there are four causes. **Measured on the BioHPC server
2026-08-01, the answer was #1** — and the diagnostic output looked healthy right up to the
timing line.

1.  **R is linked against a SERIAL OpenBLAS build.** The BLAS path says `openblas`, so it looks
    right, but the library has no threading code and **no environment variable can change
    that** — `OPENBLAS_NUM_THREADS=8` and `OMP_NUM_THREADS=8` gave 10.30 s and 10.33 s against
    10.25 s unset, all at 1.0 cores on a 64-core machine. RHEL ships three builds side by side;
    `ls -la /usr/lib64/libopenblas*` shows them:

    | library              | threading                                |
    |----------------------|------------------------------------------|
    | `libopenblas-r*.so`  | **serial** — R links here by default     |
    | `libopenblasp-r*.so` | pthread → honours `OPENBLAS_NUM_THREADS` |
    | `libopenblaso-r*.so` | OpenMP → honours `OMP_NUM_THREADS`       |

    Switch without root:

    ``` bash
    LD_PRELOAD=/usr/lib64/libopenblasp.so.0 OPENBLAS_NUM_THREADS=8 Rscript blas_check.R
    ```

2.  **Reference BLAS** — a path containing `libRblas`. The *other* single-threaded case, and it
    looks different in the diagnostic. On BioHPC, R builds from 4.4.3 onward are compiled with
    OpenBLAS; 4.0.5–4.4.2 are not, and the unmodified default is 4.2.3. `module load R/4.6.1`
    puts you past this. Note that **switching R version means reinstalling your packages** —
    they live in `$HOME/R` under a per-version path.

3.  **`run_workers.sh` capping it — on purpose.** Its second argument is BLAS threads per
    worker and **defaults to 2**. `./run_workers.sh 8` therefore uses 8 × 2 = 16 cores and
    leaves the rest idle by design. Since this workload is memory-limited, spend the surplus on
    threads: `./run_workers.sh 6 8`. Keep `workers × threads` at or under the core count.
    **Caveat**: the launcher exports `OPENBLAS_NUM_THREADS`, which a *serial* build (#1) ignores
    entirely — it is not a safety net.

4.  **There is very little linear algebra to thread.** See below.

> **Temper expectations — this pipeline is not BLAS-bound.** Benchmarked at realistic sizes:
> the VanRaden GEMM (`tcrossprod`, 500 accessions × 130k markers) **0.13 s**;
> `rrBLUP::mixed.solve` on 600 accessions **0.11 s**; whole-population marker QC **~0.5 s**; a
> 130k-marker VCF parse **~9 s**; and the largest single item, `eigen()` on 3000×3000, **25 s**.
> That is well under a minute against a **median evaluation of 29 minutes**, so upwards of 95%
> of the wall time is elsewhere — almost certainly VCF downloads and BrAPI calls, which no BLAS
> setting touches. The pipeline's products are also tall-and-thin and memory-bandwidth-bound
> rather than FLOP-bound, so they thread poorly even when threads are available.
>
> Fixing the thread cap is worth doing for correctness, but **it will not convert spare cores
> into throughput**. A profiled real evaluation on the server (`profile_evaluation.R`, trial
> 10938, 803 s) put **`build_kernel` at 6.4 s — 0.8% of the run**, against `.find_related` at
> 54% and `.group_by_panel` at 13%, with CPU/wall = 0.98: serial R computation, not waiting and
> not linear algebra. One busy core per worker is what that looks like in `top`, and it is
> correct behaviour, not a broken BLAS.

- [ ] **Find where the time actually goes** before optimising anything:

    ``` bash
    Rscript profile_evaluation.R --trial=<id>
    ```

    It runs one real evaluation with per-subtask timing and reports the split plus CPU/wall.
    Attack the stage it names — on the run above that was `choose_geno_sources`, not anything
    BLAS touches.

## 2. Settings

Put these in **`settings.local.R`** (untracked, gitignored), never in the tracked `settings.R`
— that is what stops `git pull` conflicting on every machine. Copy `settings.local.R.example`
to start.

- [ ] `simulate = FALSE` — otherwise you are optimizing the synthetic world
- [ ] `optimize_scheme` — **one** scheme per run (`"CV0"` or `"CV00"`). They are different
  tasks with different best pipelines; to do both, run twice against the same store.
  (`schemes` is a different setting, used only by the diagnostics.)
- [ ] `max_iters` — leave unset (unbounded) unless you want a short bounded run; end a
  long one with `touch state/STOP`, which stops it cleanly between evaluations
- [ ] `target_domain` — the trial list appropriate for your target trial
- [ ] **Memory budget** — `dosage_budget_bytes` and `dosage_total_budget_bytes`, from the table
  just below
- [ ] **⑂ `db_path` on `/workdir`** — **required** for parallel workers. SQLite's WAL mode is
  what lets them write concurrently, and it cannot work over NFS, which is what `/home` usually
  is.
- [ ] **⑂ `db_backup_path` on `/home`** — the store is copied there after every evaluation
  with `VACUUM INTO` (~0.01 s), so a purged `/workdir` costs only the evaluations in flight.
- [ ] **⑂ Move an existing store to the new `db_path`** — the workers continue from whatever is
  at that path, so it must be the store you care about:

    ``` bash
    mkdir -p /workdir/<user>/T3optimizer
    cp $HOME/T3optimizer/state/evals.sqlite /workdir/<user>/T3optimizer/evals.sqlite
    ```

``` r
# settings.local.R
settings_override <- list(
  simulate                  = FALSE,
  optimize_scheme           = "CV00",
  dosage_budget_bytes       = 16e9,     # keep at whatever the CACHE was built with
  dosage_total_budget_bytes = 18e9,     # 512 GB, 8 workers -- this is what bounds memory
  db_path        = "/workdir/<user>/T3optimizer/evals.sqlite",
  db_backup_path = "/home/<user>/T3optimizer/state/evals_backup.sqlite"
)
```

### The memory budget

`dosage_budget_bytes` caps **one project** at parse time and sets marker density.
`dosage_total_budget_bytes` caps the **sum over every project covering a trial**, which is what
actually bounds a worker's peak.

**Measured anchor (2026-07-29):** one worker at `dosage_budget_bytes = 16e9`, with **no total
cap set**, held **63 GB RSS** on a 128 GB node — half the machine for a single process. That
works out to about **2.3 GB of process memory per GB of resident dosage**, which is what the
table below is derived from. *Uncapped*, 16e9 supports exactly one worker; **with**
`dosage_total_budget_bytes` set it runs eight fine, because the cap is what bounds a worker.

| RAM    | workers | `dosage_budget_bytes` | `dosage_total_budget_bytes` |
|--------|---------|-----------------------|-----------------------------|
| 256 GB | 4       | 8e9                   | 19e9                        |
| 256 GB | 8       | 4e9                   | 9e9                         |
| 512 GB | 4       | 16e9                  | 35e9                        |
| 512 GB | 8       | 8e9                   | 18e9                        |

The `dosage_budget_bytes` column is for building a **fresh** cache. **If you already have a
cache, do not lower it** — lowering it re-thins nothing (`get_project_dosage` only ever
re-parses to make a cache *denser*), so you would pay for a full re-download to get *worse*
density. Leave it at whatever the cache was built with and let `dosage_total_budget_bytes`
bound the memory; the cap thins at serve time by column-subsetting the cached matrix, which is
free.

> **Existing 16e9 cache, 512 GB, 8 workers:** keep `dosage_budget_bytes = 16e9` and set
> `dosage_total_budget_bytes = 18e9`. ~24 GB of dosage per worker, ~193 GB across all eight
> (38% of the machine). Fine. Rebuilding at 8e9 would save ~3.4 GB per worker — not worth the
> re-download.

The one residual cost of the denser cache is a **transient floor**: `readRDS` loads a cached
file whole before the cap subsets it, so the largest single project (10.2 GB at 16e9, 6.8 GB at
8e9) is briefly resident on top of everything retained. `README.md` has the full table.

`report_memory.R` sizes from **`peak_rss_mb`** (true peak RSS), not `peak_r_mb` (R's heap peak,
which cannot see BLAS/`sommer`/`dist` and understated the real figure by 2.7×). The RSS peak
needs Linux; elsewhere it is `NA` and the report says its numbers understate.

**Start at half the worker count you want**, confirm with `report_memory.R` and
`monitor_memory.sh`, then scale up. The 2.3× multiplier is derived from a single node-hour and
is more likely too low than too high — it assumed every cached project covered that trial, and
if fewer did, the real multiplier is larger.

**Scaling up does not need a restart.** Add workers above the running ids:

``` bash
./run_workers.sh 4                              # workers 1-4 (1 is the leader)
# ... a day later, once report_memory.R says there is room ...
OPTIMIZER_FIRST_WORKER=5 ./run_workers.sh 4     # adds workers 5-8
```

Do **not** just run `./run_workers.sh 4` a second time: it always starts at id 1, so you would
get two leaders both rsyncing the cache, ambiguous `worker` values in the store, and — because
the log redirect truncates — the second batch wiping the first batch's logs. The script now
refuses to do this and tells you the right id to use.

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

    ``` bash
    Rscript -e 'source("settings.R"); str(optimizer_settings()[c("db_path","simulate",
      "optimize_scheme","dosage_budget_bytes","dosage_total_budget_bytes","stop_file")])'
    ```

- [ ] **Tests pass.**

    ``` bash
    Rscript tests/run_all.R          # expect "3/3 test files passed"
    Rscript tests/run_all.R --all    # adds the slow sim loop: "4/4"
    ```

    A missing or wrong `.Renviron`, or a typo in `settings.local.R`, usually surfaces here.

## 3. Cache

The cache is every genotyping project's VCF, downloaded once and stored parsed
(`cache/dosage`), plus phenotypes, accessions and the trial catalogue. It lives on the work
disk and is backed up to `$HOME/T3optimizer/cache`. It is large but regenerable.

- [ ] **How full is it?**

    ``` bash
    find cache -type f | wc -l
    ls cache/dosage | wc -l          # projects already downloaded + parsed
    du -sh cache/dosage
    ```

- [ ] **Restore it** if you have just copied the repo to a fresh node — the checkout brings no
  cache, because it lives under `$HOME/T3optimizer/cache/`.

    ``` bash
    rsync -a $HOME/T3optimizer/cache/ \
      /workdir/<user>/T3Predictathon2026/scripts/Analysis_Claude/optimizer/cache/
    ```

    **The trailing slashes matter.** Without the one on the source, rsync *creates a nested
    directory* instead of merging. This takes a while; drop `-v` unless you want the file list.
    (The leader worker also does this automatically at startup, and the other workers wait for
    it to finish.)

- [ ] **Check `cache/unparseable/`.** A project in there is negative-cached and will never be
  retried. T3 may have fixed the VCF since — delete the entry to force a retry.

- [ ] **⑂ If the cache is cold, run ONE worker first.** Not required, but the first evaluations
  on a cold cache are nearly all download, and eight simultaneous first-time VCF downloads is
  both a lot of load on T3 and the peak memory moment. There is no special command — it is the
  ordinary run:

    ``` bash
    Rscript run_optimizer.R </dev/null > logs/warm.out 2>&1 &
    tail -f logs/warm.out                  # watch for downloads to stop dominating
    touch "$STOP_FILE"                     # then stop it
    ```

  Leave the stop file in place — `run_workers.sh` clears a leftover one itself.

## 4. Launch

- [ ] **Start the workers.**

    ``` bash
    nohup ./run_workers.sh 8 > logs/workers.out 2>&1 &     # 8 workers, 2 BLAS threads each
    nohup ./run_workers.sh 4 4 > logs/workers.out 2>&1 &   # 4 workers, 4 threads each
    ```

    `run_workers.sh` does not replace `run_optimizer.R` — it launches N copies of it, each with
    `OPTIMIZER_WORKER=i`, a thread limit, and its own `logs/run_w<i>.out`. Keep
    `workers × threads` at or under the core count.

    **Launch it with `nohup ... &`.** The workers themselves are each already `nohup`'d, so they
    survive a hangup either way — but the launcher ends in `wait` (which a SLURM batch script
    needs, so the allocation is not torn down early), and in an interactive shell that just
    blocks your terminal for the whole run.

    *Single worker instead:*

    ``` bash
    nohup Rscript run_optimizer.R </dev/null > logs/run.out 2>&1 &
    ```

- [ ] **Start the memory monitor — after the workers**, since it exits when the STOP file
  appears:

    ``` bash
    nohup ./monitor_memory.sh > /dev/null 2>&1 &
    ```

    It writes to `logs/memory_<host>.tsv`, so its stdout is only the one startup line. Inside
    tmux/screen a plain `./monitor_memory.sh &` is fine.

- [ ] **Confirm it came up.**
    - No WAL warning in `logs/run_w1.out` — a warning means the store is still on NFS, so stop,
      move `db_path`, and restart
    - `logs/run_w1.out` says `worker=1, leader` — it did the cache restore. Beyond that it is
      an ordinary worker, so this is a startup check, not something to keep watching
    - all N workers present, and rows accruing:

    ``` bash
    Rscript peek_workers.R      # process count, who holds what, anything lost
    ```

## 5. Monitoring

**Where the logs are depends on the mode.** `run_workers.sh` writes to `settings$log_dir`,
which `settings.R` derives from `perm_dir`: `./logs` on a laptop, but **`$OPTIMIZER_HOME/logs`
whenever `OPTIMIZER_HOME` is set** — which includes a rented BioHPC server. `optimizer_paths.sh`
exports it as `$LOG_DIR`, so `source ./optimizer_paths.sh` first and the commands below work in
either case. (Logs living on durable storage is the point: they outlive a purged work disk.)

``` bash
source ./optimizer_paths.sh               # $LOG_DIR, $REPORT_PATH, $STOP_FILE, ...
Rscript peek_workers.R                    # START HERE: who is alive, on what, what died
cat  "$REPORT_PATH"                       # best pipeline so far + learning curve
tail -f "$LOG_DIR/run_w1.out"             # live log (one per worker)
Rscript report_memory.R                   # peak per evaluation, workers-that-fit table
Rscript report_timing.R                   # where the wall time goes
sort -k3 -n "$LOG_DIR"/memory_*.tsv | tail -5   # most memory-intensive moments
find cache -type f | wc -l                # number of files in cache
# Out of memory causing processes to be killed
dmesg -T 2>/dev/null | grep -iE 'killed process|out of memory' | tail -20
journalctl -k --since "3 days ago" 2>/dev/null | grep -i -B2 -A6 'out of memory' | tail -40
```

**`peek_workers.R` is the one to run first.** It answers "what is each worker doing, and has
anything died since I launched" from three sources at once, because no single one is enough:
the claims table (a worker killed *mid*-evaluation), the process count, and how long since each
worker's log was touched. Its **WORKERS LOST MID-EVALUATION** section is the history — a claim
whose owner is gone is recorded before it is reclaimed, which is the only durable trace of an
OOM kill, since a killed evaluation writes no store row.

**Worker 1 is no longer special.** Its only exclusive jobs are restoring the cache and clearing
stale claims at startup; backups are every worker's job (LESSONS #25). Once the run is going,
`worker=1, leader` in the log means it did the restore, not that it must stay alive — what
matters is how many workers are alive, which is the first line of `peek_workers.R`'s output.

**Sweep every worker's log at once**, rather than tailing one:

``` bash
# always a problem
grep -l -E 'FATAL|step error:|could not put the store in WAL mode' "$LOG_DIR"/run_w*.out
# what kind, and how much of it
grep -h -o -E 'FATAL|step error:|could not put the store in WAL mode' "$LOG_DIR"/run_w*.out \
  | sort | uniq -c | sort -rn
grep -c 'gave up waiting' "$LOG_DIR"/run_w*.out # cache-lock contention between workers
find "$LOG_DIR"/run_w*.out -mmin -30            # which logs are still moving
```

`the store is NOT being backed up` is deliberately **not** in that list: with
`db_backup_path = NULL` — the default on a laptop — every worker logs it on every iteration, so
it would match always and mean nothing. On a machine where you *did* configure a backup, add it
to the pattern; there it is a real alarm.

The last one inverts usefully: a file **missing** from that list has logged nothing for half an
hour. That is normal for a worker inside a long evaluation — 33 h is the observed maximum
(LESSONS #28) — and suspicious for one holding no claim. `peek_workers.R` makes exactly that
distinction for you.

**Throughput — is it still making progress, and how fast?** `report.md`'s best-score line
cannot answer this, because a run that has stopped improving looks identical to one that has
stopped running:

``` bash
Rscript -e 'source("settings.R"); source("R/store.R")
  con <- open_store(optimizer_settings()$db_path)
  e <- read_evals(con); close_store(con)
  e$hour <- substr(e$ts, 1, 13)
  print(tail(as.data.frame(table(hour = e$hour)), 24))'
```

Red flags:

| Symptom | Means | Do |
|------------------------|------------------------|------------------------|
| WAL warning at startup | store is on NFS; concurrent writes unsafe | stop, move `db_path` to `/workdir`, restart |
| `rss_total_mb` approaching RAM | too many workers for the budget | fewer workers, or lower `dosage_total_budget_bytes` |
| fewer workers than you launched | some died — OOM is the usual cause | `Rscript peek_workers.R`: its WORKERS LOST section names them and the trial each was on |
| `report.md` not updating | no worker has finished an evaluation since the last write | normal if evaluations are long; check `logs/run_w*.out` for progress |
| every evaluation `infeasible` | usually a data/name bug, not the configs | run the `EVALUATION_CHECKLIST.md` diagnostics |
| `serving 1 marker in N` messages | the aggregate cap is binding | fine if deliberate; else raise `dosage_total_budget_bytes` |
| `peak_rss_mb` ≫ `peak_r_mb` in `report_memory.R` | the cost is in compiled code (GRM/BLAS/sommer), which `dosage_total_budget_bytes` does **not** bound | size from `peak_rss_mb`; to cap it you must bound the kernel stage, not the dosage budget |

The store backup appears in `report.md` as `- store backup: N min ago, up to date` or
`… , 37 row(s) behind`, flagged `** STALE **` when it is behind and more than 5 minutes old.
The backup fires after every evaluation, so a persistent lag means backups are failing — read
the `store backup -> … FAILED` line for the reason. An *age* of hours with `up to date` is
normal: it just means no worker has finished an evaluation recently.

## 6. Stop and save

- [ ] **Stop.** One file halts every worker after its current evaluation, and the monitor with
  them:

    ``` bash
    touch "$STOP_FILE"
    ```

    Leave it there — the next launch clears it.

- [ ] **Save the store — this is the essential one.** The run *is* that file: the incumbent,
  the surrogate's training data, and what to try next are all recomputed from it.

    ``` bash
    cp /workdir/<user>/T3optimizer/evals.sqlite $HOME/T3optimizer/state/
    ```

    With `db_backup_path` set this has been happening every 30 minutes anyway. If you want to
    be certain the copy is complete, take it with `VACUUM INTO` instead — that reads the WAL
    too and emits a single self-contained file:

    ``` bash
    sqlite3 /workdir/<user>/T3optimizer/evals.sqlite \
      "VACUUM INTO '$HOME/T3optimizer/state/evals_manual_$(date +%Y%m%d_%H%M).sqlite'"
    ```

    > **Never copy `evals.sqlite` without its `-wal` and `-shm` sidecars.** In WAL mode the main
    > file only changes at a checkpoint, so it can be days old while all the recent rows sit in
    > the `-wal`. A plain `cp` of just the one file silently loses everything since the last
    > checkpoint.

- [ ] **Save the cache — optional, large, regenerable.** Keeping it lets the next run skip
  re-downloading. Note this is the **opposite direction** from the restore in step 3:

    ``` bash
    rsync -a /workdir/<user>/T3Predictathon2026/scripts/Analysis_Claude/optimizer/cache/ \
      $HOME/T3optimizer/cache/
    ```

- [ ] **Pull to the laptop** if the reservation is ending — FileZilla, or
  `rsync -av biohpc:/home/<user>/T3optimizer/state/ ./state/`.

- [ ] **Rescue anything gitignored** that lives only on that machine: `settings.local.R`, and
  any local notes. `git pull` will not bring them back on the next machine.

To **resume** later — same server, new reservation, or laptop — put `evals.sqlite` (and
`cache/` if you kept it) back at the paths the settings name, and launch again. Configurations
already run are skipped by hash, so nothing is repeated.

------------------------------------------------------------------------

## 7. Diagnostics

What each script answers is in `README.md`. Here is how to run them on this kind of machine —
plainly, from this directory:

``` bash
Rscript profile_evaluation.R --trial=<id>          # one real evaluation, per-subtask timing
Rscript surrogate_bakeoff.R                        # surrogate comparison + learning curve
Rscript peek_failures.R                            # every non-ok eval, funnel unpacked
Rscript peek_config.R --ids=4,23                   # the configurations those evals ran
Rscript blas_check.R                               # is linear algebra threaded (see §1)
nohup Rscript prewarm_indices.R > logs/prewarm.out 2>&1 &   # fill the local wizards
```

**All are read-only against the store except `prewarm_indices.R`**, which writes cache files.
The read-only ones copy the database *and* its `-wal`/`-shm` sidecars before reading, so they
are safe to run against a live store while workers are going.

### Where the time goes

Pick a trial with several stored evaluations so the cache is warm; otherwise you are timing
downloads. `--kernel=em_combine` profiles the heavy kernel instead of the default. The report
gives the six subtasks, the expensive internals, and **CPU/wall** — which separates "waiting on
the network" from "compute-bound on one core".

### Which surrogate is best

``` bash
Rscript surrogate_bakeoff.R --no-curve --reps=12   # just the full-store table, tighter bars
Rscript surrogate_bakeoff.R --curve-reps=4 --n-grid=4    # faster curve
```

Three arms — `pooled` (config means, no trial), `blocked` (`trial_id` in the forest,
marginalised), `merf` (trial as a random effect outside the tree). Repeated cross-validation
over configurations, so every number carries an error bar, and the footer states **the smallest
difference the current store can resolve**. A gap smaller than that is not evidence of no
effect; it is too little data. That line exists because a single 5-fold split once suggested an
improvement that twelve splits showed was zero.

The learning curve re-runs the comparison on the first *n* evaluations by timestamp and writes
`logs/surrogate_bakeoff.tsv` (tidy, so the figure can be redrawn) and
`logs/surrogate_bakeoff.png`.

> **Read the curve carefully.** Rows are in time order, so it mixes *more data* with
> *differently distributed data*: early evaluations are seeds and random draws, later ones are
> acquisition picks concentrated in one region. A falling segment can mean the search narrowed,
> not that the surrogate got worse. The script prints this above the table.

### Why evaluations failed

`peek_failures.R` reads the funnel stored in `detail` and prints a verdict per row — whether
the kernel covered training lines but no focal lines, covered neither, or was merely thin. Rows
whose status comes from scoring carry a reason string instead of a funnel; those are printed
separately, because the reason *is* the explanation.

### Fill the local wizards

`prewarm_indices.R` fetches `project -> accessions` (110 projects, ~1 min) and
`trial -> accessions` (the catalogue, ~25 min) so candidate discovery needs no BrAPI call.
Resumable, rate-limited, and safe to run beside working workers — every fetch is cached
independently and workers contribute to the same files. Finish when both lines read
`covers universe: TRUE`.
