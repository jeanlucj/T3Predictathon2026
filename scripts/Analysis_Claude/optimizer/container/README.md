# `container/` — the Apptainer image, and submitting jobs that use it

Eleven files that group into five jobs. Nothing here is optional scaffolding: every file is
referenced from somewhere else in the tree.

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

**Submit the diagnostic.**

| file | what it is |
|---|---|
| `diagnose_ceres.sh` | the `sbatch` script for `tools/diagnose_failures.R`. |
| `submit_diagnose.sh` | queues it. |

**The coupling worth stating:** `submit_diagnose.sh` reads `ACCOUNT` out of `submit.local.sh`.
So `submit.local.sh` is a prerequisite for submitting the *diagnostic* too, not just the
optimizer — a fresh clone that has only ever run the diagnostic still has to create it.

**No container at all.**

| file | what it is |
|---|---|
| `setup_fallback_libs.R` | populates a personal R library from CRAN when the image is unavailable — 30–60 min of compiling `lme4`, `sommer`, `Matrix`. Run it from the optimizer root: `Rscript container/setup_fallback_libs.R`. |
| `FALLBACK_modules.md` | when to take that route, which modules to load, and how to verify it worked. |

The container is the better answer unless you need RStudio.
