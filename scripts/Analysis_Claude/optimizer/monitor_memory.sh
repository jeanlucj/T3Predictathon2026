#!/bin/bash
# monitor_memory.sh
#
# Sample the memory footprint of the running optimizer worker(s) and append it to a TSV.
#
# READ-ONLY and EXTERNAL: it only reads `ps` and the OS memory counters, so it attaches to a
# run that is ALREADY GOING -- no restart, no code change, no interaction with the optimizer.
# (For per-EVALUATION peaks attributed to a configuration, see the peak_r_mb column that
# R/memory.R records in the store, and report_memory.R. This script answers the different
# question of what the machine as a whole is doing right now.)
#
# Usage, from the optimizer directory:
#   ./monitor_memory.sh &                    # 60 s interval, logs/memory_<host>.tsv
#   ./monitor_memory.sh 30 logs/mem.tsv      # custom interval (s) and output file
#   nohup ./monitor_memory.sh > /dev/null 2>&1 &     # survives logout
# Stop it with kill, or by creating the same STOP file the optimizer watches.
#
# Columns (tab-separated, one row per sample):
#   iso_time  n_proc  rss_total_mb  rss_max_mb  mem_avail_mb  mem_total_mb  pids_rss
# where pids_rss is a comma-separated pid:rss_mb list, so a single greedy worker can be
# told apart from N workers that are each fine.
#
# Read it back in R with:
#   readr::read_tsv("logs/memory_<host>.tsv")

set -u

INTERVAL="${1:-60}"
HOST="$(hostname -s 2>/dev/null || hostname)"
OUT="${2:-logs/memory_${HOST}.tsv}"
# The optimizer's own stop-file, so stopping the run stops the monitor too. Resolved by
# asking R: OPTIMIZER_PATH lives in .Renviron, which only R reads, so building this path from
# the shell environment would silently watch ./state/STOP on a server whose real stop file is
# under $HOME -- and the monitor would never exit. (See optimizer_paths.sh.)
. "$(dirname "$0")/optimizer_paths.sh"
STOP_FILE="${STOP_FILE:-./state/STOP}"

mkdir -p "$(dirname "$OUT")"
if [ ! -s "$OUT" ]; then
  printf 'iso_time\tn_proc\trss_total_mb\trss_max_mb\tmem_avail_mb\tmem_total_mb\tpids_rss\n' > "$OUT"
fi

# Total and available system memory in MB. MemAvailable (Linux) is the honest number -- it
# counts reclaimable page cache, which MemFree does not, and the page cache here is large
# because of the VCF/RDS traffic. macOS has no equivalent, so free pages are used instead.
sys_mem() {
  if [ -r /proc/meminfo ]; then
    awk '/^MemTotal:/ {t=$2} /^MemAvailable:/ {a=$2} END {printf "%d\t%d", a/1024, t/1024}' /proc/meminfo
  elif command -v vm_stat >/dev/null 2>&1; then
    local pgsize total
    pgsize=$(sysctl -n hw.pagesize 2>/dev/null || echo 4096)
    total=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
    vm_stat | awk -v ps="$pgsize" -v tot="$total" '
      /Pages free/          {gsub(/\./,""); f=$3}
      /Pages inactive/      {gsub(/\./,""); i=$3}
      /Pages speculative/   {gsub(/\./,""); s=$3}
      END {printf "%d\t%d", (f+i+s)*ps/1048576, tot/1048576}'
  else
    printf 'NA\tNA'
  fi
}

echo "monitor_memory.sh: sampling every ${INTERVAL}s into ${OUT} (Ctrl-C or ${STOP_FILE} to stop)" >&2

while true; do
  [ -f "$STOP_FILE" ] && { echo "monitor_memory.sh: STOP file present -- exiting" >&2; break; }

  # Every process whose command line mentions run_optimizer.R, excluding this script and the
  # grep/ps pipeline itself. RSS from ps is in KB. `ps -e -o` is the portable spelling that
  # behaves the same on Linux and macOS (`ps -C` is Linux-only).
  procs=$(ps -e -o pid=,rss=,args= 2>/dev/null \
    | grep 'run_optimizer\.R' \
    | grep -v 'monitor_memory' \
    | grep -v 'grep')

  # Two awk passes over the same text so the counts and the pid list can be placed on
  # either side of the system-memory columns (pids_rss is last, per the header).
  counts=$(printf '%s\n' "$procs" | awk 'NF { mb = $2/1024; n++; tot += mb; if (mb > mx) mx = mb }
                                         END { printf "%d\t%.0f\t%.0f", n+0, tot+0, mx+0 }')
  pids=$(printf '%s\n' "$procs" | awk 'NF { list = (list == "" ? "" : list ",") $1 ":" sprintf("%.0f", $2/1024) }
                                       END { printf "%s", (list == "" ? "-" : list) }')

  printf '%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$counts" "$(sys_mem)" "$pids" >> "$OUT"
  sleep "$INTERVAL"
done
