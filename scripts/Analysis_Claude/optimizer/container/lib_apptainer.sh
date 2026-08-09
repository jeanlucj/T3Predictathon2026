# lib_apptainer.sh -- put apptainer on PATH, loading the module if that is what it takes.
#
#   source "$(dirname "$0")/lib_apptainer.sh"
#   ensure_apptainer
#
# On SciNet Ceres `apptainer` is NOT on PATH until `module load apptainer`. On BioHPC's
# Rocky 9 nodes it is already there. This handles both without the caller caring.
#
# The awkward part is batch jobs. `module` is a shell FUNCTION defined by /etc/profile.d, and
# an sbatch script running under `#!/bin/bash` is neither a login nor an interactive shell, so
# that function may not exist -- `module load` then fails with "module: command not found"
# even on a host where it works perfectly from the terminal. So the module system is
# initialised first when needed, rather than assumed.

ensure_apptainer() {
  command -v apptainer >/dev/null 2>&1 && return 0

  # `command -v` reports shell functions too, which is what `module` normally is.
  if ! command -v module >/dev/null 2>&1; then
    for _init in /etc/profile.d/lmod.sh \
                 /etc/profile.d/modules.sh \
                 /usr/share/lmod/lmod/init/bash \
                 /usr/share/Modules/init/bash; do
      # shellcheck disable=SC1090
      [ -r "$_init" ] && . "$_init" && break
    done
    unset _init
  fi

  if command -v module >/dev/null 2>&1; then
    # Not fatal on its own: the re-check below is what decides.
    module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
  fi

  command -v apptainer >/dev/null 2>&1 && return 0

  echo "apptainer is not on PATH and could not be module-loaded." >&2
  echo "  tried: module load apptainer, module load singularity" >&2
  echo "  module command available: $(command -v module >/dev/null 2>&1 && echo yes || echo no)" >&2
  echo >&2
  echo "  SciNet Ceres  : module load apptainer   (needed -- it is NOT loaded by default)" >&2
  echo "  BioHPC Rocky 9: already present; CentOS 7 hosts have 'singularity' instead" >&2
  echo "  check what exists: module spider apptainer" >&2
  return 1
}
