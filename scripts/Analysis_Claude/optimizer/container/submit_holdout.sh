#!/bin/bash
# submit_holdout.sh -- queue dev/holdout_test.R as a Ceres batch job.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer/container
#   ./submit_holdout.sh -- --trials=15633,15636,...        # the held-out grid
#   ./submit_holdout.sh --time=48:00:00 -- --trials=... --k=12
#   ./submit_holdout.sh --test-only                        # sbatch's dry run; queues nothing
#
# Everything before `--` goes to sbatch; everything after it goes to dev/holdout_test.R.
#
# Run this AFTER an optimization has produced contenders. It evaluates every selected
# configuration and every submission on trials the optimizer never saw -- the measurement the
# "we beat the submissions" claim rests on. Read the result with
# `Rscript dev/holdout_test.R --report`.
#
# TRACKED, unlike submit.local.sh, because it holds no site values: the account comes from
# your submit.local.sh (or $ACCOUNT), and OPTIMIZER_HOME from .Renviron. Nothing to re-copy
# when this file changes.

set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
repo="$(dirname "$here")"

# _renviron_get, and only that: lib_submit.sh defines functions and submits nothing on source.
. "$here/lib_submit.sh"

# ---- sbatch flags before `--`, R arguments after --------------------------
sbatch_args=()
hold_args=()
seen_sep=0
for a in "$@"; do
  if [ "$seen_sep" -eq 0 ] && [ "$a" = "--" ]; then seen_sep=1; continue; fi
  if [ "$seen_sep" -eq 1 ]; then hold_args+=("$a"); else sbatch_args+=("$a"); fi
done

# ---- account: yours, wherever you keep it ---------------------------------
# submit.local.sh is the file that already holds it. Read the assignment rather than sourcing:
# that file ends by calling submit_optimizer, so sourcing it would queue the OPTIMIZER.
ACCOUNT="${ACCOUNT:-}"
if [ -z "$ACCOUNT" ] && [ -f "$here/submit.local.sh" ]; then
  ACCOUNT="$(sed -n 's/^[[:space:]]*ACCOUNT=//p' "$here/submit.local.sh" | tail -n 1 \
             | sed -e 's/[[:space:]]*$//' -e 's/^["'\'']//' -e 's/["'\'']$//')"
fi
if [ -z "$ACCOUNT" ] || [ "$ACCOUNT" = "CHANGEME" ]; then
  echo "no account: set ACCOUNT=... in the environment, or fill it in $here/submit.local.sh" >&2
  echo "  find it with: sacctmgr -Pns show user format=account,defaultaccount" >&2
  exit 1
fi

# ---- OPTIMIZER_HOME: from .Renviron, the single source --------------------
renv="$repo/.Renviron"
[ -f "$renv" ] || { echo "no .Renviron at $renv -- copy .Renviron.example first" >&2; exit 1; }
oh="$(_renviron_get "$renv" OPTIMIZER_HOME)"
case "$oh" in
  "")   echo "OPTIMIZER_HOME is not set in $renv" >&2; exit 1 ;;
  /*)   : ;;
  *'$'*|'~'*)
        echo "OPTIMIZER_HOME in $renv must be a literal absolute path, not '$oh'" >&2; exit 1 ;;
  *)    echo "OPTIMIZER_HOME in $renv must be an absolute path, got '$oh'" >&2; exit 1 ;;
esac

export REPO="$repo"
export OPTIMIZER_HOME="$oh"
export SIF="${SIF:-$here/optimizer.sif}"
# Read by holdout_ceres.sh inside the job; --export=ALL carries it there.
export HOLD_ARGS="${hold_args[*]-}"

# SLURM does not create the log directory, and a job whose output file cannot be opened dies at
# launch with the error going nowhere.
mkdir -p "$OPTIMIZER_HOME/logs"

# --dependency=singleton on a name of its own. Two prewarm jobs would fetch the same keys twice
# and flush the same durable cache tree against each other; and this must NOT queue behind, or
# hold up, the t3opt and t3diag jobs, which is what a shared name would do.
exec sbatch \
  --account="$ACCOUNT" \
  --job-name=t3hold \
  --dependency=singleton \
  --chdir="$repo" \
  --output="$OPTIMIZER_HOME/logs/hold-%j.out" \
  --error="$OPTIMIZER_HOME/logs/hold-%j.err" \
  --nodes=1 \
  --ntasks=12 \
  --mem=800G \
  --time=24:00:00 \
  --export=ALL \
  ${sbatch_args[@]+"${sbatch_args[@]}"} \
  "$here/holdout_ceres.sh"
