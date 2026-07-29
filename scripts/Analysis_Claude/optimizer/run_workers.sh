#!/bin/bash
# run_workers.sh
#
# Launch N optimizer workers against ONE evals.sqlite and ONE cache.
#
# Each worker runs the ordinary loop -- choose a configuration, sample a trial, evaluate,
# store -- so N workers explore N configurations at once and all of them feed the same
# archive and the same surrogate. They are independent processes, not threads: R's heap is
# per-process, so N workers cost N times the memory of one. That is the whole reason
# report_memory.R exists. Size N from its "how many workers fit" table.
#
#   ./run_workers.sh 8                  # 8 workers, 2 BLAS threads each
#   ./run_workers.sh 4 4                # 4 workers, 4 BLAS threads each
#   nohup ./run_workers.sh 8 > /dev/null 2>&1 &
#
# Stop them ALL with the usual stop-file (they share it):   touch state/STOP
# Watch them:   tail -f logs/run_w1.out    /   ./monitor_memory.sh &
#
# Two things this script exists to get right:
#
#  * BLAS THREADS. R's linear algebra will otherwise grab every core in EACH worker, so 8
#    workers x N cores oversubscribe the machine and everything slows down together. The
#    product (workers x threads) should be at most the core count.
#  * WORKER IDENTITY. OPTIMIZER_WORKER tells settings.R which worker this is. Worker 1 is
#    the "leader" and is the only one that writes the report, rsyncs the cache, and backs up
#    the store -- see settings$is_leader.
#
# PREREQUISITE: db_path must be on LOCAL disk. SQLite's WAL mode, which is what makes
# concurrent writers safe, cannot work on NFS. open_store warns if the pragma did not take;
# if you see that warning, run ONE worker or move the store (settings.local.R):
#     db_path         = "/workdir/<user>/optimizer/evals.sqlite"
#     db_backup_path  = "~/t3_optimizer/state/evals_backup.sqlite"

set -u

N_WORKERS="${1:-4}"
N_THREADS="${2:-2}"

cd "$(dirname "$0")" || exit 1
mkdir -p logs

# One thread setting for every BLAS R might be linked against.
export OMP_NUM_THREADS="$N_THREADS"
export OPENBLAS_NUM_THREADS="$N_THREADS"
export MKL_NUM_THREADS="$N_THREADS"
export VECLIB_MAXIMUM_THREADS="$N_THREADS"

STOP_FILE="${OPTIMIZER_PATH:-.}/state/STOP"
if [ -f "$STOP_FILE" ]; then
  echo "run_workers.sh: $STOP_FILE exists -- the workers would exit at once. Remove it first." >&2
  exit 1
fi

echo "run_workers.sh: starting $N_WORKERS worker(s), $N_THREADS BLAS thread(s) each"
for i in $(seq 1 "$N_WORKERS"); do
  # </dev/null so an archived-VCF download can never block on an interactive prompt.
  OPTIMIZER_WORKER="$i" nohup Rscript run_optimizer.R </dev/null > "logs/run_w${i}.out" 2>&1 &
  echo "  worker $i -> pid $! -> logs/run_w${i}.out$([ "$i" = 1 ] && echo '  (leader: writes the report, backs up the store)')"
  # Stagger the starts. The workers would otherwise hit the trial catalogue and the same
  # uncached genotyping projects simultaneously; the download lock makes that correct but
  # waiting is still wasted time, and a thundering herd on the T3 server is worth avoiding.
  sleep 20
done

echo
echo "all workers launched. stop them with:  touch $STOP_FILE"
wait
