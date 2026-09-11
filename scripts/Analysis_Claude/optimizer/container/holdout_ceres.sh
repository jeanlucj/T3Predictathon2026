#!/bin/bash
#SBATCH --job-name=t3hold
#SBATCH --partition=ceres
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=800G
#SBATCH --time=24:00:00
#SBATCH --account=CHANGEME
# Relative, a fallback only: submit_holdout.sh passes absolute paths under $OPTIMIZER_HOME/logs.
#SBATCH --output=logs/hold-%j.out
#SBATCH --error=logs/hold-%j.err
#
# dev/holdout_test.R on SciNet Ceres, inside the Apptainer image.
#
# Submit through container/submit_holdout.sh, which supplies the account, the paths and the
# real sizing. Arguments for the R script arrive in HOLD_ARGS.
#
# WHAT THIS IS. The held-out test: every optimized configuration AND every submission, run on
# trials the optimizer never saw. Without it an incumbent's score is a maximum over noisy means
# on the trials it was selected on, which is optimistic by construction.
#
# Sizing is a real evaluation's, not the diagnostic's: this runs run_pipeline for every cell of a
# (configuration x trial) grid, so --mem must carry N concurrent evaluations exactly as the
# optimizer job does. Budget ~13 configurations x N trials at a ~2 h mean.
#
# SAFE TO LOSE. Every cell is claimed and stored as it completes, so a wall-clock kill costs only
# the in-flight ones and resubmitting skips what is done. Results go to a store of their own --
# never the optimizer's, so a held-out evaluation can never become training data.

set -euo pipefail

: "${OPTIMIZER_HOME:?not set -- submit through container/submit_holdout.sh}"
: "${REPO:?not set -- submit through container/submit_holdout.sh}"
export OPTIMIZER_HOME

SIF="${SIF:-$REPO/container/optimizer.sif}"
N_WORKERS="${N_WORKERS:-12}"
N_THREADS="${N_THREADS:-1}"
HOLD_ARGS="${HOLD_ARGS:-}"

# apptainer is not on PATH by default on Ceres, and a batch shell may not define `module`.
. "$REPO/container/lib_apptainer.sh"
ensure_apptainer || exit 1

# cache_dir is node-local ($TMPDIR) and starts EMPTY; dev/holdout_test.R restores the
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
echo "args      : ${HOLD_ARGS:-(none)}"
echo "workers   : $N_WORKERS x $N_THREADS threads"
# The worker count lives in three places -- N_WORKERS here, --ntasks and --mem in
# submit_holdout.sh -- and nothing enforces that they agree. Each copy runs a full run_pipeline,
# whose measured worst case is 82 GB, so print the arithmetic where a mismatch is visible in the
# first ten lines of the log rather than as an OOM eight hours in.
if [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ "$N_WORKERS" -gt 0 ]; then
  MEM_PER=$(( SLURM_MEM_PER_NODE / N_WORKERS / 1024 ))
  echo "memory    : ${MEM_PER} GB per worker (measured worst case is 82 GB)"
  [ "$MEM_PER" -lt 82 ] && echo "  ** below the worst case: a heavy configuration may be OOM-killed." \
                                "Lower N_WORKERS or raise --mem. **"
fi
if [ -n "${SLURM_NTASKS:-}" ] && [ "$SLURM_NTASKS" -lt "$N_WORKERS" ]; then
  echo "  ** N_WORKERS ($N_WORKERS) exceeds --ntasks ($SLURM_NTASKS) **"
fi
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
  bash -c "cd '$REPO' && ./dev/run_holdout.sh $N_WORKERS $N_THREADS $HOLD_ARGS"
