# `tools/` — operating a run

Everything here is for someone **running** the optimizer: getting it ready, watching it, sizing
it, and finding out why an evaluation failed. Nothing here changes the algorithm — that is
`dev/`.

The verb in each name says **when** you reach for it:

- **`prepare_`** / **`check_`** — before a run: fill a cache, confirm what the run will start from.
- **`watch_`** — during a run, on the machine the workers are on.
- **`report_`** — after evaluations exist: what did they cost?
- **`inspect_`** — reads what the store already holds. Seconds, no network.
- **`diagnose_`** — re-runs the pipeline to find out something the store cannot tell you. Hours.

`inspect_` versus `diagnose_` is the distinction worth internalising: `inspect_failures.R`
answers "what did the funnel record?" from stored rows, and `diagnose_failures.R` answers "why,
and would any configuration have worked?" by evaluating the trial again.

## The inventory

Every row is a file in this directory, and every file in this directory has a row.

| when | script | question it answers | what it costs | safe while workers evaluate |
|---|---|---|---|---|
| before | `prepare_indices.R` | fill the local wizards so the pipeline needs no BrAPI call | ~7,700 BrAPI calls; writes the cache | `--only=projects` only |
| before | `check_backup.R` | what would a run start from — what is in the durable backup, and how much of it is in the target domain? | temp copy of the store, one `lmer` + one surrogate fit | yes, not free |
| during | `watch_workers.R` | which workers are alive, what is each one evaluating, and has anything died? | trivial; read-only handle on the **live** store | **yes** |
| during | `watch_memory.sh` | what is the machine as a whole doing right now, per worker pid? | trivial; `ps` plus the OS counters, appends a TSV | **yes** |
| after | `report_memory.R` | how much memory did evaluations really use, and how many workers fit? | trivial; live store, read | **yes** |
| after | `report_timing.R` | where does the wall time go, by kernel and `geno_select`? | trivial; live store, read | **yes** |
| wrong | `inspect_failures.R` | why did the non-`ok` evaluations fail, according to the stored funnel? | trivial; sidecar copy of the store | **yes** |
| wrong | `inspect_config.R` | what configuration did a given eval run? | trivial; sidecar copy of the store | **yes** |
| wrong | `diagnose_failures.R` | re-run a failing trial: is the data missing, or is the configuration wrong? | **~10 GB per project, hours**; full pipeline, **takes dosage locks** | **no** |
| any | `optimizer_paths.sh` | `source` it to get `$STOP_FILE`, `$DB_PATH`, `$LOG_DIR`, `$CACHE_DIR` into your shell | one `Rscript -e` | **yes** |

`watch_memory.sh` and `optimizer_paths.sh` are shell, the rest is `Rscript`. All of them run
**from the optimizer root**, not from this directory:

``` bash
cd <repo>/scripts/Analysis_Claude/optimizer
Rscript tools/inspect_failures.R
source ./tools/optimizer_paths.sh
```

That is not a convention — `.Renviron` carries `T3_USERNAME`, `T3_PASSWORD` and
`OPTIMIZER_HOME`, and **R reads `.Renviron` from the working directory only** (no parent walk).
Run from anywhere else and every catalogue fetch 401s while the cached counts still look
healthy. Each script tests for `settings.R` in the cwd and stops with the right invocation
rather than letting that happen.

## Monitoring is not diagnostics

- **Monitoring** — *is it healthy right now, and has anything broken since I launched it?*
  Worker liveness, throughput, OOM kills, backup freshness. Start with `watch_workers.R`.
- **Diagnostics** — *why did that happen?* Start with `inspect_failures.R`, and reach for
  `diagnose_failures.R` only when the stored funnel is not enough.

The distinction that trips people up on a cluster: **a dead worker is not a dead job.** SLURM
kills the job only on the cgroup limit, so one OOM-killed worker leaves `squeue` saying
`RUNNING` while you quietly run with fewer. `watch_workers.R` is what sees that.

The invocations themselves are environment-specific and live in the runbooks —
`RUNBOOK_INTERACTIVE.md` §5 and §7, `RUNBOOK_SLURM.md` §6 — which cite the rows above rather
than restating them.

## Which store they read

`resolve_read_store()` picks it: the live `db_path` when it holds rows, otherwise
`db_backup_path`. On a cluster the second is the usual case — `db_path` is node-local scratch
belonging to the job that wrote it, so a diagnostic run from any other allocation cannot see
it. Each script prints which file it opened and how stale it is; anything evaluated since that
backup is not in it. Override with `--store=<path>`. On a login node they all refuse, with the
`salloc` recipe.

`watch_workers.R` is the exception, and has to be: in-flight work is recorded only in the
**live** store, so it must run on the node the workers are on — a backup's copy of that table
is whatever happened to be running when the last `VACUUM INTO` fired. Under SLURM that means
attaching to the running job (`srun --overlap --jobid=<jid> --pty bash`), not a fresh `salloc`.
It says so when the live store is absent.

`inspect_failures.R` and `check_backup.R` need a store that already has rows. On a fresh
install they stop with `no store at <db_path>`, which is correct behaviour rather than a fault.

## On a cluster, run them inside the container

``` bash
./container/run_in_container.sh exec tools/inspect_failures.R    # any script, plus its arguments
./container/run_in_container.sh shell                            # interactive R inside the image
```

A bare `Rscript` uses the cluster module's R, which has **no packages installed** — the symptom
is `there is no package called 'here'`. Either go through the container as above, or populate a
personal library once with `Rscript container/setup_fallback_libs.R` (see
`container/FALLBACK_modules.md`); the container is the better answer unless you need RStudio.

**An extra process on a cluster is not free.** A step started with `srun --overlap --jobid=<jid>`
**joins the job's cgroup**, so its memory counts against the same `--mem` the workers are using —
and that is sized tight (`container/t3opt_ceres.sh` asks `--mem=1800G` for 22 workers whose
measured worst case is 82 GB each, i.e. 1804 GB). SLURM **kills rather than throttles**. Typical
headroom is large, because 22 × the *median* 19 GB is ~418 GB — so the "yes" column above is a
judgement about free memory right now, not a blanket promise. Check first:

``` bash
sstat -j <jid>.batch --format=JobID,MaxRSS,AveRSS
tail -3 "$OPTIMIZER_HOME/logs/memory_"*.tsv     # rss_total_mb and mem_avail_mb columns
```

**Cap threads on anything you launch by hand.** Only `run_workers.sh` exports
`OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS`/…, and `container/run_in_container.sh` deliberately avoids
`--cleanenv` so your shell's environment passes straight through. A script started from an
`srun` shell therefore inherits *no* cap while the workers sit at one thread each.

## Why `diagnose_failures.R` is the one "no"

It runs a real `run_pipeline`, so beyond the hours and the tens of GB it takes the per-project
dosage locks — and `.acquire_lock` breaks any lock older than `lock_stale_minutes` (90) on the
assumption its holder died. A worker legitimately 90+ minutes into a large VCF download can have
its lock stolen, and `.fetch_and_cache_dosage` deletes the raw file before re-downloading, so
both sides then fight over the same file. That is `docs/LESSONS.md` #24 — the exact hazard the
lock exists to prevent. `--cached-only` avoids the downloads (and so the locks) at the cost of
reporting coverage as a lower bound; the closing `run_pipeline` replay is off unless you pass
`--replay`. `inspect_failures.R` answers the same first question in seconds from the stored
funnel, with no network and no re-running. `dev/profile_evaluation.R` is red for the same
reasons.

## Where a new script goes

`tools/` if an operator would run it during a normal run; `dev/` if it exists because something
broke or because the algorithm is being changed; nothing new at the root. See `dev/README.md`.
