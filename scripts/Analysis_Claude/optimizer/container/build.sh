#!/bin/bash
# build.sh -- build optimizer.sif from optimizer.def.
#
# Works on any Linux host with Apptainer: the Cornell BioHPC Rocky 9 machines (Apptainer
# 1.1.2) or a SciNet COMPUTE node. There is no native macOS build, so this cannot run on
# the laptop.
#
#   ssh <biohpc-host>
#   cd <path>/optimizer/container
#   ./build.sh                       # writes ./optimizer.sif
#   ./build.sh /project/x/opt.sif    # or name the output
#
# Takes 20-45 minutes on a first build: lme4, sommer and Matrix all compile.
#
# --fakeroot is implied automatically for an unprivileged `apptainer build`, so no admin
# involvement is needed. If the host has user namespaces disabled the build fails with a
# permissions error -- that is the case Singularity Cloud Builder exists for.

set -euo pipefail
cd "$(dirname "$0")"

DEF="optimizer.def"
SIF="${1:-optimizer.sif}"

command -v apptainer >/dev/null 2>&1 || {
  echo "apptainer not on PATH." >&2
  echo "  SciNet Ceres  : should already be present, no module load needed." >&2
  echo "  BioHPC Rocky 9: present; CentOS 7 hosts have 'singularity' instead." >&2
  exit 1
}
[ -f "$DEF" ] || { echo "no $DEF in $(pwd)" >&2; exit 1; }

# Build on a compute node, not a login node. SciNet documents this explicitly, and a build
# is a heavy compile that does not belong on a shared login host either way.
case "$(hostname -s)" in
  ceres|atlas|*-login*|login*)
    echo "This looks like a login node ($(hostname -s))." >&2
    echo "Get an allocation first, e.g.:  salloc -N1 -n8 --mem=32G -t 2:00:00 -A <account>" >&2
    exit 1 ;;
esac

# Image layers land in ~/.apptainer unless redirected, and a 30 GB home fills fast.
# Ceres sets these for SLURM jobs; BioHPC does not.
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${TMPDIR:-/tmp}/apptainer_cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${TMPDIR:-/tmp}/apptainer_tmp}"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

echo "apptainer : $(apptainer --version)"
echo "cachedir  : $APPTAINER_CACHEDIR"
echo "output    : $SIF"
echo

apptainer build --force "$SIF" "$DEF"

echo
echo "== built =="
ls -lh "$SIF"
# The digest is what makes a result traceable to an exact image. Record it alongside
# OPTIMIZER_BUILD when reporting anything.
apptainer inspect "$SIF" 2>/dev/null | sed 's/^/  /' || true
echo
echo "sha256: $(sha256sum "$SIF" | cut -d' ' -f1)"
echo
echo "Now run the built-in test:   apptainer test $SIF"
