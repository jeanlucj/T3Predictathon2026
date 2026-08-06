#!/bin/bash
#SBATCH --job-name=t3opt
#SBATCH --partition=ceres
#SBATCH --nodes=1                    # ONE node: every worker shares one SQLite store
#SBATCH --ntasks=32
#SBATCH --mem=1800G                  # ~22 workers x 82 GB measured worst case
#SBATCH --time=21-00:00:00           # 3 weeks, the ceres maximum
#SBATCH --account=CHANGEME
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#
# Optimizer on SciNet Ceres, inside the Apptainer image. See USDA_ARS_SCINET.md.
#
#   sbatch container/t3opt_ceres.sh
#
# CHAIN IT for runs longer than 3 weeks -- the store is resumable, so only the evaluations
# in flight at the cutoff are lost (median 29 min each):
#
#   jid=$(sbatch --parsable container/t3opt_ceres.sh)
#   sbatch --dependency=afterany:$jid container/t3opt_ceres.sh
#
# BEFORE THE FIRST REAL SUBMISSION, do a 30-minute shakeout:
#   sbatch --qos=debug --time=00:30:00 --ntasks=4 --mem=200G container/t3opt_ceres.sh
# Two workers for half an hour costs nothing and catches every path and credential mistake.

set -euo pipefail

# ---- paths ---------------------------------------------------------------
# OPTIMIZER_HOME must be PROJECT storage, never home: home is 30 GB and the dosage cache
# alone can exceed that.
export OPTIMIZER_HOME="${OPTIMIZER_HOME:-/project/CHANGEME/t3_optimizer}"
REPO="${REPO:-/project/CHANGEME/T3Predictathon2026/scripts/Analysis_Claude/optimizer}"
SIF="${SIF:-$OPTIMIZER_HOME/optimizer.sif}"

# The live store goes on node-local disk. SQLite's WAL coordinates through an mmap'd -shm
# file that network filesystems do not provide; settings.local.R points db_path at $TMPDIR,
# and backup_store() VACUUM INTOs to project storage from there.
#
# On Ceres SLURM sets TMPDIR to /local/bgfs/<jobid> -- 1.5 TB of node-local SSD, erased when
# the job exits. Unset means we are not inside a job, and proceeding would put the store on a
# network filesystem shared by every worker: the one configuration that corrupts it.
: "${TMPDIR:?TMPDIR is unset -- not inside a SLURM job? The store MUST be on node-local disk}"

# Workers x threads <= ntasks. Memory, not cores, is the binding constraint: raise the
# worker count only after report_memory.R confirms headroom.
N_WORKERS="${N_WORKERS:-22}"
N_THREADS="${N_THREADS:-1}"

mkdir -p "$OPTIMIZER_HOME/state" "$OPTIMIZER_HOME/logs" logs

# ---- restore the store into node-local scratch ---------------------------
# $TMPDIR is empty at job start and erased at job end, and the optimizer does NOT restore the
# store by itself -- restore_cache_from_backup() covers the cache, but there is no equivalent
# for evals.sqlite. Without this copy every chained job would start from an empty store and
# re-run work already paid for.
BACKUP="$OPTIMIZER_HOME/state/evals_backup.sqlite"
STORE="$TMPDIR/evals.sqlite"
if [ -f "$BACKUP" ]; then
  cp "$BACKUP" "$STORE"
  # VACUUM INTO emits a self-contained file with no -wal sidecar, so the copy is complete as
  # it stands. Report what came back, so a silently-empty restore cannot pass unnoticed.
  echo "restored $(sqlite3 "$STORE" 'SELECT COUNT(*) FROM evals;' 2>/dev/null || echo '?') rows from $BACKUP"
else
  echo "no backup at $BACKUP -- starting from an EMPTY store"
  echo "  (expected only on the very first job; check this if you meant to resume)"
fi

echo "node      : $(hostname)"
echo "TMPDIR    : $TMPDIR"
echo "home      : $OPTIMIZER_HOME"
echo "image     : $SIF"
echo "workers   : $N_WORKERS x $N_THREADS threads"
apptainer inspect "$SIF" 2>/dev/null | grep -iE "optimizer.build|r.version|sha" || true

# The cd into $REPO is load-bearing: R reads .Renviron -- and therefore the T3 credentials
# -- from the working directory only. Apptainer 1.1+ does not bind $PWD, hence the explicit
# binds, each mapping a path onto itself so the paths settings.R computes stay valid.
apptainer exec \
  --bind "$REPO:$REPO" \
  --bind "$OPTIMIZER_HOME:$OPTIMIZER_HOME" \
  --bind "$TMPDIR:$TMPDIR" \
  "$SIF" \
  bash -c "cd '$REPO' && ./run_workers.sh $N_WORKERS $N_THREADS"

# run_workers.sh ends in `wait`, so this returns only when every worker has stopped. On
# timeout SLURM kills the job; the next chained job resumes from the backed-up store.
