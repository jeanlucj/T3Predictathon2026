# lib_submit_multiple.sh -- the submission mechanism for a multi-setting plan. Source it, then
# call submit_multiple.
#
#   . "$(cd "$(dirname "$0")" && pwd)/lib_submit_multiple.sh"
#   submit_multiple "$@"
#
# Same split as lib_submit.sh, and for the same reason: submit_multiple.local.sh is gitignored,
# so `git pull` never updates it and a copy made months ago would keep submitting the way it did
# then. VALUES there, LOGIC here. If you want to change HOW submission works, edit this file.
#
# Shaped like lib_submit.sh and lib_apptainer.sh: it DEFINES a function rather than doing the
# work on source, so sourcing a file never has the surprising side effect of submitting a job.
#
# It reuses _renviron_get from lib_submit.sh -- which is safe to source precisely because that
# file only defines functions.

_here_lsm="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$_here_lsm/lib_submit.sh"
unset _here_lsm

submit_multiple() {
  local here repo renv oh plan

  # container/ lives inside the optimizer directory, so the repo is this file's grandparent.
  # Derived, never configured: it cannot then point at a different checkout than the one you
  # are standing in.
  here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  repo="$(dirname "$here")"

  # .Renviron is the single source for OPTIMIZER_HOME. R overrides the inherited environment
  # from it at startup, so a value set in the shell would lose to it inside the container --
  # the shell would bind one directory while R used another, and that failure looks like a
  # permissions problem rather than a config mismatch.
  renv="$repo/.Renviron"
  [ -f "$renv" ] || { echo "no .Renviron at $renv -- copy .Renviron.example first" >&2; return 1; }
  oh="$(_renviron_get "$renv" OPTIMIZER_HOME)"
  case "$oh" in
    *'$'*|'~'*)
      echo "OPTIMIZER_HOME in $renv must be a literal path, not '$oh'" >&2
      echo "  R expands \${HOME} and ~; this shell read does not, so they would differ." >&2
      return 1 ;;
  esac
  case "$oh" in
    "")  echo "OPTIMIZER_HOME is not set in $renv" >&2; return 1 ;;
    /*)  : ;;
    *)   echo "OPTIMIZER_HOME in $renv must be an absolute path, got '$oh'" >&2; return 1 ;;
  esac

  # The plan itself. Checked HERE, on the login node, because the alternative is discovering a
  # typo'd filename after the job has queued. The deep checks -- do the settings files exist,
  # do they parse, is a holdout inside its own target domain -- need R and so run inside the
  # image as the job's first act; this is only "is there a plan at all".
  plan="${MULTIPLE_SETTINGS:-$repo/multiple_settings.R}"
  case "$plan" in /*) : ;; *) plan="$repo/$plan" ;; esac
  if [ ! -f "$plan" ]; then
    echo "no plan at $plan" >&2
    echo "  cp $repo/multiple_settings.R.example $repo/multiple_settings.R, then edit it." >&2
    echo "  Check it before submitting:  cd $repo && Rscript tools/plan_multiple.R --check" >&2
    return 1
  fi

  export REPO="$repo"
  export OPTIMIZER_HOME="$oh"
  export MULTIPLE_SETTINGS="$plan"
  # Beside the recipe that built it, matching build.sh and run_in_container.sh. Not in
  # OPTIMIZER_HOME, which is for state that cannot be regenerated.
  export SIF="${SIF:-$here/optimizer.sif}"
  export OPTIMIZER_SUBMIT_LIB=1

  # SLURM does not create the log directory, and a job whose output file cannot be opened dies
  # at launch with the error going nowhere -- because the output file is where it would have
  # gone. That failure looks exactly like "nothing happened".
  mkdir -p "$OPTIMIZER_HOME/logs"

  # A job name of its own, NOT t3opt. --dependency=singleton keys on the name, so sharing one
  # would make a plan queue behind an ordinary optimizer job and vice versa -- and the two want
  # the same store, so they must never run at once either. Keeping the names distinct while
  # both carry singleton gives: one plan at a time, one optimizer at a time, and nothing
  # stopping you from checking first that the other is not running.
  #
  # Absolute --output/--error: the #SBATCH directives in t3multi_ceres.sh are relative, so
  # SLURM would otherwise resolve them against whatever directory you submitted from.
  # "$@" comes before the script path so extra arguments reach sbatch, not the job.
  exec sbatch \
    --account="$ACCOUNT" \
    --job-name="${JOB_NAME:-t3multi}" \
    --dependency=singleton \
    --chdir="$repo" \
    --output="$OPTIMIZER_HOME/logs/multi-%j.out" \
    --error="$OPTIMIZER_HOME/logs/multi-%j.err" \
    --nodes=1 \
    --ntasks="$NTASKS" \
    --mem="$MEM" \
    --time="$TIME" \
    --export=ALL \
    "$@" \
    "$here/t3multi_ceres.sh"
}
