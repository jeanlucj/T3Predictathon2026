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

# Refuse to TRUNCATE a log a live worker is writing to -- but decide "live" by asking the OS,
# not by the log's mtime.
#
# mtime alone cannot tell a worker writing now from one killed thirty seconds ago, and the
# killed case is the common one: cancel and resubmit is how a code change gets picked up.
# $LOG_DIR is on DURABLE storage and shared across nodes, so those warm logs survive the node
# change too. It is erratic as well as wrong -- workers touch their logs only when they emit a
# message, so which id looks "active" is whichever happened to log last before the kill. That
# is how a relaunch was refused on worker 5 while workers 1-4 passed. Same defect as
# docs/LESSONS.md #24: liveness judged by a timestamp. See dev/README.md thread 7.
#
# CAVEAT: liveness is per NODE while $LOG_DIR is shared, so a worker on another node is
# invisible here. --dependency=singleton makes two concurrent t3opt jobs impossible, and the
# documented add-workers route (srun --overlap --jobid=<jid>) puts you on the job's own node.
worker_alive() {                       # $1 = worker id; identified by its OPTIMIZER_WORKER
  local p
  for p in $(pgrep -u "$(id -u)" -f 'run_optimizer\.R' 2>/dev/null); do
    [ -r "/proc/$p/environ" ] || continue
    tr '\0' '\n' < "/proc/$p/environ" 2>/dev/null | grep -qx "OPTIMIZER_WORKER=$1" && return 0
  done
  return 1
}

if [ -n "${OPTIMIZER_FORCE_LAUNCH:-}" ]; then
  echo "run_workers.sh: OPTIMIZER_FORCE_LAUNCH set -- skipping the running-worker check"
elif [ "$already" -eq 0 ]; then
  # Nothing of ours is running on this node, so no log here can belong to a live worker,
  # whatever its mtime says. This is the cancel-and-resubmit case.
  warm=$(find "$LOG_DIR" -maxdepth 1 -name 'run_w*.out' -mmin -2 2>/dev/null | wc -l | tr -d ' ')
  [ "$warm" -gt 0 ] && echo "run_workers.sh: $warm recently-written log(s) belong to workers that" \
                            "are gone; truncating them"
else
  for i in $(seq "$FIRST" "$LAST"); do
    if [ -d /proc ]; then
      worker_alive "$i" || continue                       # exact: this id is not running
    else
      # No /proc (a macOS interactive run): fall back to the mtime heuristic, which is all
      # that is available. Only reached when some run_optimizer.R IS running.
      [ -n "$(find "$LOG_DIR/run_w${i}.out" -mmin -2 2>/dev/null)" ] || continue
    fi
    # Suggest the next FREE id, taken from the highest run_w<id>.out on disk -- not from the
    # process count, which says nothing about which ids are in use.
    max_id=$(ls "$LOG_DIR"/run_w*.out 2>/dev/null | sed 's/.*run_w\([0-9]*\)\.out/\1/' | sort -n | tail -1)
    echo "run_workers.sh: worker $i is already running (its log is $LOG_DIR/run_w${i}.out)." >&2
    echo "  To ADD workers to it, start above the running ids:" >&2
    echo "    OPTIMIZER_FIRST_WORKER=$(( ${max_id:-0} + 1 )) $0 $N_WORKERS $N_THREADS" >&2
    echo "  If you believe it is NOT running, re-run with OPTIMIZER_FORCE_LAUNCH=1." >&2
    exit 1
  done
fi

# Seconds between worker launches. Also the unit the "failed launch" window below is built
# from, so the two cannot drift apart.
STAGGER=20

echo "run_workers.sh: starting $N_WORKERS worker(s) (ids $FIRST-$LAST), $N_THREADS BLAS thread(s) each"
WPIDS=(); WIDS=(); T0=$SECONDS
for i in $(seq "$FIRST" "$LAST"); do
  # </dev/null so an archived-VCF download can never block on an interactive prompt.
  OPTIMIZER_WORKER="$i" nohup Rscript run_optimizer.R </dev/null > "$LOG_DIR/run_w${i}.out" 2>&1 &
  WPIDS+=("$!"); WIDS+=("$i")
  echo "  worker $i -> pid $! -> $LOG_DIR/run_w${i}.out$([ "$i" = 1 ] && echo '  (leader: restores the cache at startup)')"
  # Stagger the starts. The workers would otherwise hit the trial catalogue and the same
  # uncached genotyping projects simultaneously; the download lock makes that correct but
  # waiting is still wasted time, and a thundering herd on the T3 server is worth avoiding.
  sleep "$STAGGER"
done

echo
echo "all workers launched. stop them with:  touch $STOP_FILE"
echo "(the next fresh launch clears that file itself -- no rm needed)"

# Block until every worker exits, COLLECTING their exit statuses. This is for SLURM: a batch
# script that returns immediately would end the job and tear the allocation down under the
# workers. The workers are each nohup'd above, so they do NOT depend on this process staying
# alive -- which is why an interactive launch should be
# `nohup ./run_workers.sh N > logs/workers.out 2>&1 &`, or the wait below simply occupies your
# terminal for the length of the run.
#
# `wait` WITH NO ARGUMENTS RETURNS 0 whatever the children did. That is how two multi-day
# launches were lost: every worker died within minutes, this script exited 0, apptainer exited
# 0, and sacct recorded the job COMPLETED -- so --mail-type=FAIL could not fire and the only
# signal was t3opt quietly leaving squeue. Waiting on each pid is what makes SLURM's own
# failure reporting work. See dev/README.md, "A dead run reports success to SLURM".
n_ok=0; n_fail=0; failed_ids=""
for k in "${!WPIDS[@]}"; do
  if wait "${WPIDS[$k]}"; then
    n_ok=$((n_ok + 1))
  else
    n_fail=$((n_fail + 1)); failed_ids="$failed_ids ${WIDS[$k]}"
  fi
done
ELAPSED=$((SECONDS - T0))

if [ "$n_fail" -eq 0 ]; then
  echo "run_workers.sh: all $n_ok worker(s) exited cleanly after ${ELAPSED}s."
  exit 0
fi

echo "run_workers.sh: $n_fail of $((n_ok + n_fail)) worker(s) exited NON-ZERO (ids:$failed_ids)" >&2
for i in $failed_ids; do
  echo "  --- last 5 lines of $LOG_DIR/run_w${i}.out ---" >&2
  tail -n 5 "$LOG_DIR/run_w${i}.out" >&2 2>/dev/null || true
done

# Nothing survived: whatever else is true, this run produced no work.
if [ "$n_ok" -eq 0 ]; then
  echo "run_workers.sh: NO worker survived -- this run produced nothing." >&2
  exit 1
fi

# A FAILED LAUNCH, not attrition. Distinguishing the two is the point: a worker OOM-killed at
# hour 30 is tolerated by design (the run continues with fewer), while workers gone within
# minutes of starting means none of them ever got going. The window is the stagger -- the
# launch itself takes STAGGER * (N-1) seconds -- plus a grace period. A STOP file means the
# operator asked for this, so a quick exit is intentional rather than a fault.
EARLY_GRACE="${OPTIMIZER_EARLY_EXIT_S:-300}"
EARLY_LIMIT=$(( STAGGER * (N_WORKERS - 1) + EARLY_GRACE ))
if [ "$ELAPSED" -lt "$EARLY_LIMIT" ] && [ ! -f "$STOP_FILE" ]; then
  echo "run_workers.sh: workers were gone ${ELAPSED}s after launch (< ${EARLY_LIMIT}s) with no" >&2
  echo "  STOP file -- this is a failed LAUNCH, not attrition. Read the tails above." >&2
  exit 1
fi

echo "run_workers.sh: $n_ok worker(s) ran to completion; treating the $n_fail loss(es) as" >&2
echo "  attrition rather than a failed run. Check report_memory.R if they were OOM kills." >&2
exit 0
