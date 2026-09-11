#!/bin/bash
# run_holdout.sh -- N copies of dev/holdout_test.R against one shared holdout store.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   ./dev/run_holdout.sh 8 1 --trials=15633,15636,...
#
# Arg 1 is the worker count, arg 2 the BLAS threads each; everything after goes to the R script.
# The grid is a fixed (configuration x trial) product, so there is no search to coordinate --
# claim_eval() alone keeps two copies off the same cell, and each copy starts at a different
# offset so they do not queue behind one another on the first row.
#
# Like run_workers.sh, this blocks until every copy exits and reports a non-zero status if they
# all failed: a batch script that returns immediately would tear the allocation down under them.

set -u
cd "$(dirname "$0")/.." || exit 1     # the optimizer ROOT: R reads .Renviron from the cwd

N_WORKERS="${1:-4}"; N_THREADS="${2:-1}"; shift 2 2>/dev/null || shift $# 
export OMP_NUM_THREADS="$N_THREADS" OPENBLAS_NUM_THREADS="$N_THREADS"
export MKL_NUM_THREADS="$N_THREADS" VECLIB_MAXIMUM_THREADS="$N_THREADS"

. "./tools/optimizer_paths.sh"
LOG_DIR="${LOG_DIR:-logs}"; mkdir -p "$LOG_DIR"

# Seconds between launches, matching run_workers.sh. Not cosmetic: four processes opening one
# SQLite store within milliseconds contend on open_store()'s `PRAGMA synchronous`, and
# .with_busy_retry's five attempts can all be exhausted -- which is what makes
# tests/test_concurrency.R fail about one run in nine. Copy 1 also restores the cache while the
# others wait, so a wide stagger costs nothing here.
STAGGER="${OPTIMIZER_STAGGER:-20}"

echo "run_holdout.sh: starting $N_WORKERS copy(ies), $N_THREADS BLAS thread(s) each," \
     "${STAGGER}s apart"
PIDS=(); IDS=()
for i in $(seq 1 "$N_WORKERS"); do
  OPTIMIZER_WORKER="$i" nohup Rscript dev/holdout_test.R "$@" \
    </dev/null > "$LOG_DIR/holdout_w${i}.out" 2>&1 &
  PIDS+=("$!"); IDS+=("$i")
  echo "  copy $i -> pid $! -> $LOG_DIR/holdout_w${i}.out"
  sleep "$STAGGER"
done

n_ok=0; n_fail=0
for k in "${!PIDS[@]}"; do
  if wait "${PIDS[$k]}"; then n_ok=$((n_ok+1)); else n_fail=$((n_fail+1))
    echo "  --- tail of $LOG_DIR/holdout_w${IDS[$k]}.out ---" >&2
    tail -n 5 "$LOG_DIR/holdout_w${IDS[$k]}.out" >&2 2>/dev/null || true
  fi
done
echo "run_holdout.sh: $n_ok copy(ies) finished, $n_fail failed"
[ "$n_ok" -eq 0 ] && exit 1
exit 0
