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
cd /workdir/jj332/T3Predictathon2026/scripts/Analysis_Claude/optimizer
git pull && loadR
source ./optimizer_paths.sh                  # $STOP_FILE, $DB_PATH, $REPORT_PATH, ...
Rscript tests/run_all.R                      # expect 3/3
rsync -a $HOME/T3optimizer/cache/ cache/     # trailing slashes!
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

**`source ./optimizer_paths.sh` first.** `OPTIMIZER_PATH` is set in `.Renviron`, which **only
R reads** — in a shell it is empty, so `"$OPTIMIZER_PATH/state/STOP"` expands to
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
      OPTIMIZER_PATH=/home/jj332/T3optimizer
      ```
      `OPTIMIZER_PATH` is what turns on remote mode — there is no `remote_server` flag to
      edit any more, it is derived from this variable. Confirm with
      `Rscript -e 'Sys.getenv("OPTIMIZER_PATH")'`.

      > **Never run R here with `--vanilla`.** It implies `--no-environ`, which skips
      > `.Renviron` — you lose `OPTIMIZER_PATH` *and* the T3 credentials, and the run
      > dies at login looking like a missing `.Renviron`.
- [ ] **Packages installed** — this takes a while on a fresh machine. Includes the two
      non-CRAN ones: `BrAPI.R` and `T3BrapiHelpers`. Confirm with
      `Rscript -e 'ip <- installed.packages(); which(ip[,1] == "BrAPI")'`.
      
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

```r
# settings.local.R
settings_override <- list(
  simulate                  = FALSE,
  optimize_scheme           = "CV00",
  max_hours                 = 23.5,
  dosage_budget_bytes       = 16e9,     # keep at whatever the CACHE was built with
  dosage_total_budget_bytes = 18e9,     # 512 GB, 8 workers -- this is what bounds memory
  db_path        = "/workdir/jj332/T3optimizer/evals.sqlite",
  db_backup_path = "/home/jj332/T3optimizer/state/evals_backup.sqlite"
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
would get two leaders both writing the report and rsyncing the cache, ambiguous `worker`
values in the store, and — because the log redirect truncates — the second batch wiping the
first batch's logs. The script now refuses to do this and tells you the right id to use.

Changing `dosage_budget_bytes` changes marker density and therefore scores. Density is not a
config parameter, so old and new rows are not strictly comparable; each row records its
`dosage_budget` so at least the mixing stays visible.

### Checks

These come **after** the settings, not before. `optimizer_settings()` *throws* on an
`optimize_scheme` that is not in `schemes`, and on remote mode with no `OPTIMIZER_PATH`; a
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
- [ ] **⑂ Move an existing store to the new `db_path`** — the workers continue from whatever
      is at that path, so it must be the store you care about:
      ```bash
      mkdir -p /workdir/jj332/T3optimizer
      cp $HOME/T3optimizer/state/evals.sqlite /workdir/jj332/T3optimizer/evals.sqlite
      ```

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
      rsync -a $HOME/T3optimizer/cache/ /workdir/jj332/T3Predictathon2026/scripts/Analysis_Claude/optimizer/cache/
      ```
      **The trailing slashes are load-bearing.** Without the one on the source, rsync
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
tail logs/memory_*.tsv                    # whole-machine RSS over time
```

Red flags:

| Symptom | Means | Do |
|---|---|---|
| WAL warning at startup | store is on NFS; concurrent writes unsafe | stop, move `db_path` to `/workdir`, restart |
| `rss_total_mb` approaching RAM | too many workers for the budget | fewer workers, or lower `dosage_total_budget_bytes` |
| only worker 1's rows in the store | the others died | check `logs/run_w*.out` |
| `report.md` not updating | worker 1 is not running | restart it — it is the leader |
| every evaluation `infeasible` | usually a data/name bug, not the configs | run the `EVALUATION_CHECKLIST.md` diagnostics |
| `serving 1 marker in N` messages | the aggregate cap is binding | fine if deliberate; else raise `dosage_total_budget_bytes` |

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
      cp /workdir/jj332/T3optimizer/evals.sqlite $HOME/T3optimizer/state/
      ```
      (With `db_backup_path` set, the leader has been doing this every 30 min anyway.)
- [ ] **Save the cache — optional, large, regenerable.** Keeping it lets the next run skip
      re-downloading. Note this is the **opposite direction** from the restore in step 3:
      ```bash
      rsync -a /workdir/jj332/T3Predictathon2026/scripts/Analysis_Claude/optimizer/cache/ $HOME/T3optimizer/cache/
      ```
- [ ] **Pull to the laptop** if the reservation is ending — FileZilla, or
      `rsync -av biohpc:/home/jj332/T3optimizer/state/ ./state/`.

To **resume** later — same server, new reservation, or laptop — put `evals.sqlite` (and
`cache/` if you kept it) back at the paths the settings name, and launch again.
Configurations already run are skipped by hash, so nothing is repeated.
