#!/bin/bash
# check_oom.sh
#
# Did something in this allocation get killed for memory, and how close to the limit did it
# get? Reads the job's cgroup counters and the accounting record. READ-ONLY.
#
# Run it INSIDE the allocation whose processes you are asking about -- the cgroup high-water
# mark is not reset between processes, so it still holds the peak of a run that has ended, but
# it disappears with the allocation.
#
#   ./check_oom.sh              # this allocation
#   ./check_oom.sh 12345678     # sacct section for another job id as well
#
# The decisive line is `oom_kill` (cgroup v2) or `failcnt` (v1): non-zero means the kernel
# killed something here for memory. `peak` against `limit` says how much headroom there was.

set -u
JOBID="${1:-${SLURM_JOB_ID:-}}"
hr() { printf '\n%s\n%s\n' "$1" "$(printf '%.0s-' $(seq ${#1}))"; }
mb() { awk '{ if ($1+0 > 0) printf "%.1f GB\n", $1/1073741824; else print $1 }'; }

hr "cgroup (this allocation, host $(hostname -s))"
CG=$(awk -F: '$1=="0"{print $2}' /proc/self/cgroup 2>/dev/null)
if [ -n "${CG:-}" ] && [ -d "/sys/fs/cgroup${CG}" ]; then
  # Walk up: the memory limit is usually set on the job-level cgroup, not the leaf.
  d="/sys/fs/cgroup${CG}"
  while [ "$d" != "/sys/fs/cgroup" ] && [ -n "$d" ]; do
    if [ -r "$d/memory.peak" ] || [ -r "$d/memory.events" ]; then
      echo "  $d"
      [ -r "$d/memory.peak" ]   && printf '    peak      : %s' "$(mb < "$d/memory.peak")"
      [ -r "$d/memory.max" ]    && printf '    limit     : %s' "$(mb < "$d/memory.max")"
      [ -r "$d/memory.events" ] && { echo "    events    :"; sed 's/^/      /' "$d/memory.events"; }
    fi
    d=$(dirname "$d")
  done
else
  # cgroup v1 layout.
  for d in /sys/fs/cgroup/memory/slurm/uid_$(id -u)/job_${JOBID:-0} /sys/fs/cgroup/memory; do
    [ -d "$d" ] || continue
    echo "  $d"
    [ -r "$d/memory.max_usage_in_bytes" ] && printf '    peak      : %s' "$(mb < "$d/memory.max_usage_in_bytes")"
    [ -r "$d/memory.limit_in_bytes" ]     && printf '    limit     : %s' "$(mb < "$d/memory.limit_in_bytes")"
    [ -r "$d/memory.failcnt" ]            && printf '    failcnt   : %s\n' "$(cat "$d/memory.failcnt")"
  done
fi

hr "accounting"
if [ -z "${JOBID:-}" ]; then
  echo "  no job id (not in an allocation, and none given on the command line)"
else
  # MaxRSS is recorded per STEP, never on the parent job row, and only once the step ENDS --
  # so a running interactive allocation shows blanks on the row people usually look at.
  sacct -j "$JOBID" --units=G \
        --format=JobID%20,JobName%14,State%14,ExitCode,MaxRSS,AveRSS,ReqMem,Elapsed 2>&1 |
    sed 's/^/  /'
  echo
  echo "  (blank MaxRSS on the parent row is normal: it is a per-step, end-of-step figure."
  echo "   If every row is blank, this cluster's jobacct_gather does not sample RSS and the"
  echo "   cgroup section above is the only memory evidence.)"
fi

hr "kernel log"
if dmesg -T >/dev/null 2>&1; then
  dmesg -T | grep -iE "killed process|out of memory|oom-kill" | tail -20 | sed 's/^/  /' ||
    echo "  no OOM lines"
else
  echo "  dmesg not readable by this user (normal on a compute node)"
fi
