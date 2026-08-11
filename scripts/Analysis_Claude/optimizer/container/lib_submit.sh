# lib_submit.sh -- the submission mechanism. Source it, then call submit_optimizer.
#
#   . "$(cd "$(dirname "$0")" && pwd)/lib_submit.sh"
#   submit_optimizer "$@"
#
# WHY THIS IS SEPARATE FROM submit.local.sh. That file is gitignored, so `git pull` never
# updates it -- a copy made in August still submits the way August did. Every flag added since
# was silently absent from the copy in use, which is how a run ended up writing
# slurm-<jid>.out to the repo instead of $OPTIMIZER_HOME/logs.
#
# So the split is: submit.local.sh holds VALUES (account, sizes -- things that genuinely differ
# per site and per node), this file holds LOGIC (paths, validation, the sbatch flags). Fix a
# flag here and `git pull` delivers it; nobody has to remember to re-copy anything.
#
# Shaped like lib_apptainer.sh: it DEFINES a function rather than doing the work on source, so
# sourcing a file never has the surprising side effect of submitting a job.

# Read one KEY=VALUE out of an .Renviron, stripping surrounding whitespace and quotes.
# Deliberately does NOT expand ${VAR} or ~: R would, this does not, so a value containing
# either would mean the shell and R disagree -- see the check in submit_optimizer.
_renviron_get() {
  sed -n "s/^[[:space:]]*$2[[:space:]]*=[[:space:]]*//p" "$1" \
    | tail -n 1 | sed -e 's/[[:space:]]*$//' -e 's/^["'\'']//' -e 's/["'\'']$//'
}

submit_optimizer() {
  local here repo renv oh

  # container/ lives inside the optimizer directory, so the repo is this file's grandparent.
  # Derived, never configured: it cannot then point at a different checkout than the one you
  # are standing in.
  here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  repo="$(dirname "$here")"

  # .Renviron is the single source for OPTIMIZER_HOME. R overrides the inherited environment
  # from it at startup (see optimizer_paths.sh), so a value set in the shell would lose to it
  # inside the container -- the shell would bind and restore one directory while R used
  # another, and that failure looks like a permissions problem, not a config mismatch.
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

  export REPO="$repo"
  export OPTIMIZER_HOME="$oh"
  # Beside the recipe that built it, matching build.sh and run_in_container.sh. Not in
  # OPTIMIZER_HOME, which is for state that cannot be regenerated.
  export SIF="${SIF:-$here/optimizer.sif}"
  # Lets t3opt_ceres.sh tell a current submit.local.sh from one predating this file.
  export OPTIMIZER_SUBMIT_LIB=1

  # SLURM does not create the log directory, and a job whose output file cannot be opened dies
  # at launch with the error going nowhere -- because the output file is where it would have
  # gone. That failure looks exactly like "nothing happened".
  mkdir -p "$OPTIMIZER_HOME/logs"

  # Absolute --output/--error: the #SBATCH directives in t3opt_ceres.sh are relative, so SLURM
  # would otherwise resolve them against whatever directory you happened to submit from.
  # "$@" comes before the script path so extra arguments reach sbatch, not the job.
  exec sbatch \
    --account="$ACCOUNT" \
    --job-name="$JOB_NAME" \
    --dependency=singleton \
    --chdir="$repo" \
    --output="$OPTIMIZER_HOME/logs/slurm-%j.out" \
    --error="$OPTIMIZER_HOME/logs/slurm-%j.err" \
    --ntasks="$NTASKS" \
    --mem="$MEM" \
    --time="$TIME" \
    --export=ALL \
    "$@" \
    "$here/t3opt_ceres.sh"
}
