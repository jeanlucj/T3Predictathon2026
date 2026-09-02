#!/bin/bash
# submit_diagnose.sh -- queue diagnose_failures.R as a Ceres batch job.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer/container
#   ./submit_diagnose.sh                                  # defaults
#   ./submit_diagnose.sh -- --cached-only --trials=10253   # arguments for the R script
#   ./submit_diagnose.sh --mem=96G --time=12:00:00 -- --dosage-budget-gb=40
#   ./submit_diagnose.sh --test-only                      # what it WOULD get; queues nothing
#
# Everything before `--` goes to sbatch; everything after it goes to diagnose_failures.R.
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
diag_args=()
seen_sep=0
for a in "$@"; do
  if [ "$seen_sep" -eq 0 ] && [ "$a" = "--" ]; then seen_sep=1; continue; fi
  if [ "$seen_sep" -eq 1 ]; then diag_args+=("$a"); else sbatch_args+=("$a"); fi
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
# Read by diagnose_ceres.sh inside the job; --export=ALL carries it there.
export DIAG_ARGS="${diag_args[*]-}"

# SLURM does not create the log directory, and a job whose output file cannot be opened dies at
# launch with the error going nowhere.
mkdir -p "$OPTIMIZER_HOME/logs"

# --dependency=singleton on a name of its own: two diagnostics would contend for the same
# node-local cache, and neither should block or be blocked by the optimizer's t3opt jobs.
exec sbatch \
  --account="$ACCOUNT" \
  --job-name=t3diag \
  --dependency=singleton \
  --chdir="$repo" \
  --output="$OPTIMIZER_HOME/logs/diag-%j.out" \
  --error="$OPTIMIZER_HOME/logs/diag-%j.err" \
  --nodes=1 \
  --ntasks=4 \
  --mem=48G \
  --time=08:00:00 \
  --export=ALL \
  ${sbatch_args[@]+"${sbatch_args[@]}"} \
  "$here/diagnose_ceres.sh"
