#!/bin/bash
#SBATCH --job-name=t3prep
#SBATCH --partition=ceres
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --account=CHANGEME
# Relative, a fallback only: submit_prepare.sh passes absolute paths under $OPTIMIZER_HOME/logs.
#SBATCH --output=logs/prep-%j.out
#SBATCH --error=logs/prep-%j.err
#
# tools/prepare_indices.R on SciNet Ceres, inside the Apptainer image.
#
# Submit through container/submit_prepare.sh, which supplies the account, the paths and the
# real sizing. Arguments for the R script arrive in PREP_ARGS.
#
# WHY THIS IS WORTH A JOB OF ITS OWN. Without the indices the pipeline answers "which trials
# share germplasm with this one?" by asking the BrAPI relationship wizard, and one such call
# inside .germ_overlap measured 440 s of a 795 s evaluation -- 55% of the run. Filling them
# once removes that from every evaluation of a multi-week run.
#
# Sizing is nothing like diagnose_ceres.sh: this fetches ACCESSION NAME LISTS, never a dosage
# matrix, so it is network-bound and single-threaded and 16G is already generous. The wall
# clock is the only real constraint -- ~110 project calls (about a minute) plus ~7,600 trial
# calls at the default --sleep=0.05.
#
# SAFE TO LOSE. tools/prepare_indices.R caches every fetch as it completes and syncs to the
# durable backup every 100 keys, so a job killed at the wall clock loses at most that tail;
# resubmitting skips everything already done. Raising --time is cheaper than a --limit.

set -euo pipefail

: "${OPTIMIZER_HOME:?not set -- submit through container/submit_prepare.sh}"
: "${REPO:?not set -- submit through container/submit_prepare.sh}"
export OPTIMIZER_HOME

SIF="${SIF:-$REPO/container/optimizer.sif}"
PREP_ARGS="${PREP_ARGS:-}"

# apptainer is not on PATH by default on Ceres, and a batch shell may not define `module`.
. "$REPO/container/lib_apptainer.sh"
ensure_apptainer || exit 1

# cache_dir is node-local ($TMPDIR) and starts EMPTY; tools/prepare_indices.R restores the
# durable backup into it and flushes back on the way out, which is what makes this job's work
# outlive its allocation. Both paths must be a compute node's, not a login node's.
: "${TMPDIR:?TMPDIR is unset -- not inside a SLURM job?}"
: "${SLURM_JOB_ID:?SLURM_JOB_ID is unset -- this must run as a SLURM job, not on a login node}"

mkdir -p "$OPTIMIZER_HOME/state" "$OPTIMIZER_HOME/logs"

echo "node      : $(hostname)"
echo "TMPDIR    : $TMPDIR"
echo "home      : $OPTIMIZER_HOME"
echo "image     : $SIF"
echo "mem       : ${SLURM_MEM_PER_NODE:-?} MB"
echo "args      : ${PREP_ARGS:-(none)}"
echo "build     : $(sed -n 's/^OPTIMIZER_BUILD *<- *"\(.*\)".*/\1/p' "$REPO/settings.R" 2>/dev/null || echo '?')"

# No memory sampler here, unlike the other two jobs: this one holds name lists, not genotypes,
# and a second background process would be more moving parts than the question deserves.

# The cd into $REPO is load-bearing: R reads .Renviron -- and the T3 credentials -- from the
# working directory only. Each bind maps a path onto itself so the paths settings.R computes
# stay valid inside the container.
apptainer exec \
  --bind "$REPO:$REPO" \
  --bind "$OPTIMIZER_HOME:$OPTIMIZER_HOME" \
  --bind "$TMPDIR:$TMPDIR" \
  "$SIF" \
  bash -c "cd '$REPO' && Rscript tools/prepare_indices.R $PREP_ARGS"
