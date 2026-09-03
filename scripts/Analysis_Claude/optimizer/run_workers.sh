#!/bin/bash
# run_workers.sh
#
# Launch N optimizer workers against ONE evals.sqlite and ONE cache.
#
# Each worker runs the ordinary loop -- choose a configuration, sample a trial, evaluate,
# store -- so N workers explore N configurations at once and all of them feed the same
# archive and the same surrogate. They are independent processes, not threads: R's heap is
# per-process, so N workers cost N times the memory of one. That is the whole reason
# tools/report_memory.R exists. Size N from its "how many workers fit" table.
#
#   ./run_workers.sh 8                  # 8 workers, 2 BLAS threads each
#   ./run_workers.sh 4 4                # 4 workers, 4 BLAS threads each
#   nohup ./run_workers.sh 8 > logs/workers.out 2>&1 &
#   OPTIMIZER_FIRST_WORKER=5 ./run_workers.sh 4   # ADD workers 5-8 to a job already running
#
# Stop them ALL with the usual stop-file (they share it):
#   touch "${OPTIMIZER_HOME:-.}/state/STOP"
# There is no need to remove it afterwards: a fresh launch clears a leftover one itself.
# Watch them:   tail -f logs/run_w1.out    /   nohup ./tools/watch_memory.sh > /dev/null 2>&1 &
#
# Two things this script exists to get right:
#
#  * BLAS THREADS. R's linear algebra will otherwise grab every core in EACH worker, so 8
#    workers x N cores oversubscribe the machine and everything slows down together. The
#    product (workers x threads) should be at most the core count.
#  * WORKER IDENTITY. OPTIMIZER_WORKER tells settings.R which worker this is. Worker 1 is
#    the "leader", and its ONLY remaining exclusive job is RESTORING the cache at startup --
#    see settings$is_leader. Both backups (store and cache) are done by every worker,
#    deliberately: leader-only meant one worker's long evaluation stalled everyone else's
#    backup. See LESSONS #25.
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

# One thread setting for every BLAS R might be linked against.
export OMP_NUM_THREADS="$N_THREADS"
export OPENBLAS_NUM_THREADS="$N_THREADS"
export MKL_NUM_THREADS="$N_THREADS"
export VECLIB_MAXIMUM_THREADS="$N_THREADS"

# Resolve the stop file by asking R. OPTIMIZER_HOME is set in .Renviron, which only R reads,
# so deriving this from the shell environment would act on ./state/STOP while the workers watch
# a different file -- clearing a stale stop below would then leave the real one in place and
# every worker would exit at once. (See tools/optimizer_paths.sh.)
. "$(dirname "$0")/tools/optimizer_paths.sh"
if [ -z "${STOP_FILE:-}" ]; then
  echo "run_workers.sh: tools/optimizer_paths.sh could not resolve the stop file from R." >&2
  echo "  Guessing ./state/STOP would act on a different file from the workers." >&2
  echo "  Check: Rscript -e 'source(\"settings.R\"); optimizer_settings()\$stop_file'" >&2
  exit 1
fi

# Worker logs go to settings$log_dir, which is <OPTIMIZER_HOME>/logs on a server and ./logs on
# a laptop (settings.R derives both from perm_dir) -- so this is unchanged locally, and on a
# cluster the logs land on durable storage and survive the node instead of vanishing with it.
LOG_DIR="${LOG_DIR:-logs}"
mkdir -p "$LOG_DIR"

# Worker ids run FIRST..FIRST+N-1. The default start of 1 is a fresh launch; set
# OPTIMIZER_FIRST_WORKER to ADD workers to a run that is already going:
#
#   ./run_workers.sh 4                             # workers 1-4 (1 is the leader)
#   OPTIMIZER_FIRST_WORKER=5 ./run_workers.sh 4    # add workers 5-8, no second leader
#
# Running `./run_workers.sh 4` twice would instead start a second worker 1 -- two leaders both
# restoring the cache into a directory the other workers are reading, ambiguous `worker` values
# in the store, and, because the redirect below truncates, the second batch wiping the first
# batch's logs.
FIRST="${OPTIMIZER_FIRST_WORKER:-1}"
LAST=$((FIRST + N_WORKERS - 1))

# The stop file is on durable storage, so one consumed by the last job outlives it and would
# halt this one before its first iteration. Clearing it HERE and not in the leader is what
# makes that work under N workers: workers 2..N test the file every iteration and would exit
# before worker 1 got as far as its own unlink. Adding workers is the opposite case -- the run
# is being stopped on purpose, and clearing it would restart what someone just halted.
if [ -f "$STOP_FILE" ]; then
  if [ "$FIRST" = 1 ]; then
    echo "run_workers.sh: clearing stale $STOP_FILE"
    rm -f "$STOP_FILE" || exit 1
  else
    echo "run_workers.sh: $STOP_FILE exists -- the run you are adding workers to is" >&2
    echo "  stopping. Remove it first if you meant to keep going." >&2
    exit 1
  fi
fi

already=$(ps -e -o args= 2>/dev/null | grep -c '[r]un_optimizer\.R')
if [ "$already" -gt 0 ]; then
  echo "run_workers.sh: note -- $already run_optimizer.R process(es) are already running"
fi

# Refuse to truncate a log that a live worker is writing to: a log touched in the last two
# minutes almost certainly belongs to a running worker with this id.
for i in $(seq "$FIRST" "$LAST"); do
  if [ -n "$(find "$LOG_DIR/run_w${i}.out" -mmin -2 2>/dev/null)" ]; then
    # Suggest the next FREE id, taken from the highest run_w<id>.out on disk -- not from the
    # process count, which says nothing about which ids are in use.
    max_id=$(ls "$LOG_DIR"/run_w*.out 2>/dev/null | sed 's/.*run_w\([0-9]*\)\.out/\1/' | sort -n | tail -1)
    echo "run_workers.sh: $LOG_DIR/run_w${i}.out was written to seconds ago -- worker $i looks" >&2
    echo "  like it is already running. To ADD workers, start above the running ids:" >&2
    echo "    OPTIMIZER_FIRST_WORKER=$(( ${max_id:-0} + 1 )) $0 $N_WORKERS $N_THREADS" >&2
    exit 1
  fi
done

echo "run_workers.sh: starting $N_WORKERS worker(s) (ids $FIRST-$LAST), $N_THREADS BLAS thread(s) each"
for i in $(seq "$FIRST" "$LAST"); do
  # </dev/null so an archived-VCF download can never block on an interactive prompt.
  OPTIMIZER_WORKER="$i" nohup Rscript run_optimizer.R </dev/null > "$LOG_DIR/run_w${i}.out" 2>&1 &
  echo "  worker $i -> pid $! -> $LOG_DIR/run_w${i}.out$([ "$i" = 1 ] && echo '  (leader: restores the cache at startup)')"
  # Stagger the starts. The workers would otherwise hit the trial catalogue and the same
  # uncached genotyping projects simultaneously; the download lock makes that correct but
  # waiting is still wasted time, and a thundering herd on the T3 server is worth avoiding.
  sleep 20
done

echo
echo "all workers launched. stop them with:  touch $STOP_FILE"
echo "(the next fresh launch clears that file itself -- no rm needed)"

# Block until every worker exits. This is for SLURM: a batch script that returns immediately
# would end the job and tear the allocation down under the workers. The workers are each
# nohup'd above, so they do NOT depend on this process staying alive -- which is why an
# interactive launch should be `nohup ./run_workers.sh N > logs/workers.out 2>&1 &`, or the
# wait below simply occupies your terminal for the length of the run.
wait
