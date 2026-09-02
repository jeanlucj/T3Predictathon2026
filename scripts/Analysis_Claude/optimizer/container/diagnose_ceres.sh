#!/bin/bash
#SBATCH --job-name=t3diag
#SBATCH --partition=ceres
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=CHANGEME
# Relative, a fallback only: submit_diagnose.sh passes absolute paths under $OPTIMIZER_HOME/logs.
#SBATCH --output=logs/diag-%j.out
#SBATCH --error=logs/diag-%j.err
#
# diagnose_failures.R on SciNet Ceres, inside the Apptainer image.
#
# Submit through container/submit_diagnose.sh, which supplies the account, the paths and the
# real sizing. Arguments for the R script arrive in DIAG_ARGS.
#
# Sizing: the panel probe replays subtask C, which holds every covering project resident and
# then builds merged copies, so it peaks at about twice --dosage-budget-gb (default 4 GB, so
# ~8 GB). diagnose_failures.R checks that against this job's --mem before doing any work and
# refuses rather than being OOM-killed. --mem=48G leaves room for the store copy and the
# observation tables on top.
#
# The wall clock is the other constraint: ten diagnose_trial calls, each doing an uncached
# wizard round trip for project discovery and then a panel load.

set -euo pipefail

: "${OPTIMIZER_HOME:?not set -- submit through container/submit_diagnose.sh}"
: "${REPO:?not set -- submit through container/submit_diagnose.sh}"
export OPTIMIZER_HOME

SIF="${SIF:-$REPO/container/optimizer.sif}"
DIAG_ARGS="${DIAG_ARGS:-}"

# apptainer is not on PATH by default on Ceres, and a batch shell may not define `module`.
. "$REPO/container/lib_apptainer.sh"
ensure_apptainer || exit 1

# The cache diagnose_failures.R restores into is node-local, and settings.R places it under
# $TMPDIR. Both must be a compute node's, not a login node's shared storage.
: "${TMPDIR:?TMPDIR is unset -- not inside a SLURM job?}"
: "${SLURM_JOB_ID:?SLURM_JOB_ID is unset -- this must run as a SLURM job, not on a login node}"

mkdir -p "$OPTIMIZER_HOME/state" "$OPTIMIZER_HOME/logs"

echo "node      : $(hostname)"
echo "TMPDIR    : $TMPDIR"
echo "home      : $OPTIMIZER_HOME"
echo "image     : $SIF"
echo "mem       : ${SLURM_MEM_PER_NODE:-?} MB"
echo "args      : ${DIAG_ARGS:-(none)}"
echo "build     : $(sed -n 's/^OPTIMIZER_BUILD *<- *"\(.*\)".*/\1/p' "$REPO/settings.R" 2>/dev/null || echo '?')"

# Per-PID sampling on the node, outside the container: store rows are written on completion, so
# a killed run leaves no other record of what it was holding. Exits with the job.
MEM_TSV="$OPTIMIZER_HOME/logs/diagmem_${SLURM_JOB_ID}.tsv"
STOP_FILE="$OPTIMIZER_HOME/state/STOP_DIAG" LOG_DIR="$OPTIMIZER_HOME/logs" \
  "$REPO/monitor_memory.sh" 30 "$MEM_TSV" > /dev/null 2>&1 &
MEM_PID=$!
echo "memory    : sampling every 30s -> $MEM_TSV (pid $MEM_PID)"
trap 'kill "$MEM_PID" 2>/dev/null' EXIT

# The cd into $REPO is load-bearing: R reads .Renviron -- and the T3 credentials -- from the
# working directory only. Each bind maps a path onto itself so the paths settings.R computes
# stay valid inside the container.
apptainer exec \
  --bind "$REPO:$REPO" \
  --bind "$OPTIMIZER_HOME:$OPTIMIZER_HOME" \
  --bind "$TMPDIR:$TMPDIR" \
  "$SIF" \
  bash -c "cd '$REPO' && Rscript diagnose_failures.R $DIAG_ARGS"
