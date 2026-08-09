#!/bin/bash
# run_in_container.sh -- run the optimizer (or any of its scripts) inside optimizer.sif.
#
#   ./run_in_container.sh workers 8 2          # 8 workers, 2 BLAS threads each
#   ./run_in_container.sh shell                # interactive R inside the container
#   ./run_in_container.sh exec peek_failures.R # any Rscript in the optimizer directory
#
# Environment:
#   OPTIMIZER_SIF    path to the image      (default: <this dir>/optimizer.sif)
#   OPTIMIZER_HOME   durable storage root   (must already be set for remote_server mode)
#
# THREE THINGS THIS SCRIPT EXISTS TO GET RIGHT
#
#  1. THE WORKING DIRECTORY. R reads .Renviron from the working directory ONLY -- no parent
#     walk -- and T3_USERNAME/T3_PASSWORD live in optimizer/.Renviron. Apptainer 1.1+ no
#     longer bind-mounts $PWD automatically, so the repo is bound explicitly and every
#     command is run through a `cd` into it. Get this wrong and the symptom is not an error:
#     the connection silently falls back to anonymous and every catalogue fetch 401s.
#
#  2. CREDENTIALS STAY OUT OF THE IMAGE. .Renviron is bind-mounted with the repo, never
#     baked into the .sif -- credentials inside an image are credentials in every copy of it.
#
#  3. PATHS MUST MATCH INSIDE AND OUT. settings.R computes db_path, cache_dir and the rest
#     from OPTIMIZER_HOME. Binding a path to a DIFFERENT path inside the container makes
#     those computed paths point at nothing, and settings.R falls back to local mode
#     silently. So every bind below maps a path onto itself.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
# By default the optimizer directory is this folder's parent -- container/ normally lives
# inside the repo. Override when the image and the code are kept apart (e.g. the .sif on
# project storage, the clone somewhere else).
REPO="${OPTIMIZER_REPO:-$(dirname "$HERE")}"
SIF="${OPTIMIZER_SIF:-$HERE/optimizer.sif}"
MODE="${1:-workers}"; shift || true

[ -f "$SIF" ] || { echo "no image at $SIF -- run ./build.sh first" >&2; exit 1; }
# Performs Ceres's required `module load apptainer`; a no-op where it is already on PATH.
. "$HERE/lib_apptainer.sh"
ensure_apptainer || exit 1

# Check REPO really is the optimizer directory. Without this the mistake surfaces much
# later and much less helpfully -- as "cannot open file 'tests/run_all.R'" from inside the
# container, or worse, as a run that starts and behaves oddly.
for f in settings.R run_workers.sh tests/run_all.R R/store.R; do
  [ -e "$REPO/$f" ] && continue
  echo "ERROR: $REPO is not the optimizer directory (no $f)." >&2
  echo >&2
  echo "container/ normally sits INSIDE the optimizer folder, so the parent is the repo." >&2
  echo "On a fresh machine the code is usually just not cloned yet:" >&2
  echo >&2
  echo "  git clone <repo-url> <path>/T3Predictathon2026" >&2
  echo "  mv $HERE <path>/T3Predictathon2026/scripts/Analysis_Claude/optimizer/" >&2
  echo >&2
  echo "Or keep them apart and point at the clone explicitly:" >&2
  echo "  OPTIMIZER_REPO=<path>/scripts/Analysis_Claude/optimizer $0 $MODE" >&2
  exit 1
done

if [ ! -f "$REPO/.Renviron" ]; then
  echo "WARNING: no .Renviron in $REPO" >&2
  echo "  T3_USERNAME/T3_PASSWORD live there; without it every BrAPI call runs anonymous." >&2
fi

# Bind each path onto ITSELF (see note 3). Home is auto-mounted; project storage is not.
BINDS="--bind $REPO:$REPO"
[ -n "${OPTIMIZER_HOME:-}" ] && BINDS="$BINDS --bind $OPTIMIZER_HOME:$OPTIMIZER_HOME"
# The live store must be on node-local disk -- SQLite WAL cannot work over NFS.
[ -n "${TMPDIR:-}" ] && BINDS="$BINDS --bind $TMPDIR:$TMPDIR"

# --cleanenv would drop OPTIMIZER_HOME and the BLAS thread settings, so it is deliberately
# not used; APPTAINERENV_ passthrough keeps the few variables that matter explicit.
#
# TMPDIR is forced through rather than left to Apptainer's defaults because settings.local.R
# computes db_path from it. If the container saw a different TMPDIR than the host -- /tmp, say
# -- the store would land on shared disk instead of the job's local SSD, and the only visible
# sign would be a db_path that looks subtly wrong.
export APPTAINERENV_OPTIMIZER_HOME="${OPTIMIZER_HOME:-}"
[ -n "${TMPDIR:-}" ] && export APPTAINERENV_TMPDIR="$TMPDIR"

echo "image : $SIF"
echo "repo  : $REPO"
echo "binds : $BINDS"

case "$MODE" in
  workers)
    N="${1:-4}"; T="${2:-2}"
    exec apptainer exec $BINDS "$SIF" \
      bash -c "cd '$REPO' && ./run_workers.sh $N $T"
    ;;
  exec)
    [ $# -ge 1 ] || { echo "usage: $0 exec <script.R> [args]" >&2; exit 1; }
    exec apptainer exec $BINDS "$SIF" \
      bash -c "cd '$REPO' && Rscript $*"
    ;;
  shell)
    exec apptainer exec $BINDS "$SIF" \
      bash -c "cd '$REPO' && R --no-save"
    ;;
  test)
    # The checks from the plan's verification list, in one command.
    exec apptainer exec $BINDS "$SIF" bash -c "cd '$REPO' && \
      Rscript -e 'cat(R.version.string, \"\n\"); cat(\"rsync:\", Sys.which(\"rsync\"), \"\n\")' && \
      Rscript tests/run_all.R"
    ;;
  *)
    echo "usage: $0 {workers [n] [threads] | exec <script.R> | shell | test}" >&2
    exit 1 ;;
esac
