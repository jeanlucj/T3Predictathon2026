# `container/` — the Apptainer image, and submitting jobs that use it

The files here group into six jobs, below. Nothing is optional scaffolding: every one of them
is referenced from somewhere else in the tree.

**Build the image.**

| file | what it is |
|---|---|
| `optimizer.def` | the recipe. R pinned to 4.5.3 to match SciNet's `r/4.5.3` module, so the container route and the module fallback run the same R and switching between them cannot itself move a result. |
| `build.sh` | builds `optimizer.sif` from it. Needs a Linux host with Apptainer — a BioHPC Rocky 9 machine or a SciNet **compute** node. There is no native macOS build. |

The `.sif` is deliberately untracked (`.gitignore`): the recipe travels by `git pull`, the 1–2 GB
image does not.

**Run one script inside it.**

| file | what it is |
|---|---|
| `run_in_container.sh` | `exec <script>` runs any of the optimizer's scripts inside the image; `shell` gives you interactive R there. It does `cd $REPO && Rscript $*`, so paths are given relative to the optimizer root: `exec tools/watch_workers.R`. |
| `lib_apptainer.sh` | puts `apptainer` on `PATH`, loading the module if that is what it takes. Sourced by the others. |

**Submit the optimizer.**

| file | what it is |
|---|---|
| `t3opt_ceres.sh` | the `sbatch` script: one node, node-local store, `tools/watch_memory.sh` alongside the workers. |
| `lib_submit.sh` | the submission mechanism. Source it, call `submit_optimizer`. |
| `submit.local.sh.example` | copy to `submit.local.sh` (untracked) and fill in **your** account, `/project` paths and node sizing. Values only; the logic stays in the tracked files, so `git pull` never conflicts. |
| `settings.local.R.scinet` | copy to `../settings.local.R`. The Ceres values for the store and cache paths, and the guards that keep the store off network storage. |

**Prepare the indices** — run once before a long optimization.

| file | what it is |
|---|---|
| `prepare_ceres.sh` | the `sbatch` script for `tools/prepare_indices.R`. Light: it fetches accession name lists, never a dosage matrix, so it is network-bound and the wall clock is the only real constraint. |
| `submit_prepare.sh` | queues it. Arguments after `--` reach the R script (`-- --only=projects` is the one-minute version). |

Filling the `project -> accessions` and `trial -> accessions` maps removes the BrAPI
relationship wizard from every evaluation — one such call measured 440 s of a 795 s
evaluation. The job is safe to lose: every fetch is cached as it completes and synced to the
durable backup every 100 keys, so resubmitting after a wall-clock kill skips what is done.

**Submit the diagnostic.**

| file | what it is |
|---|---|
| `diagnose_ceres.sh` | the `sbatch` script for `tools/diagnose_failures.R`. |
| `submit_diagnose.sh` | queues it. |

**The coupling worth stating:** `submit_diagnose.sh` and `submit_prepare.sh` both read
`ACCOUNT` out of `submit.local.sh`. So `submit.local.sh` is a prerequisite for every submission
here, not just the optimizer's — a fresh clone that has only ever run the diagnostic still has
to create it.

**The three job names are deliberately distinct** — `t3opt`, `t3diag`, `t3prep` — because each
wrapper passes `--dependency=singleton`. That serialises repeat submissions of the *same* kind
(two diagnostics would contend for one node-local cache; two prewarms would fetch the same keys
twice and flush the durable tree against each other) without letting any of them block a
different kind.

**No container at all.**

| file | what it is |
|---|---|
| `setup_fallback_libs.R` | populates a personal R library from CRAN when the image is unavailable — 30–60 min of compiling `lme4`, `sommer`, `Matrix`. Run it from the optimizer root: `Rscript container/setup_fallback_libs.R`. |
| `FALLBACK_modules.md` | when to take that route, which modules to load, and how to verify it worked. |

The container is the better answer unless you need RStudio.
