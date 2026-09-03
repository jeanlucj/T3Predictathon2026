#!/bin/bash
# tools/optimizer_paths.sh
#
# Resolve the optimizer's paths and export them for use in the shell.
#
#   source ./tools/optimizer_paths.sh
#   touch "$STOP_FILE"          # stops every worker
#   ls -la "$DB_PATH"
#
# WHY THIS EXISTS. `OPTIMIZER_HOME` is set in `.Renviron`, which **only R reads** -- it is not
# a shell variable unless you also export it from `.bashrc`. So in a plain shell,
# "$OPTIMIZER_HOME/state/STOP" expands to "/state/STOP": `touch` fails with permission denied,
# and you conclude you have stopped a run that is in fact still going. That is the failure
# this file exists to prevent.
#
# Asking R is also the only way to be *right*: `settings.local.R` can override `db_path`,
# `stop_file` and the rest, so reconstructing paths from `OPTIMIZER_HOME` in shell would be
# wrong on exactly the machines that matter. R applies the same layering the optimizer does.
#
# Exports: OPTIMIZER_HOME STOP_FILE DB_PATH REPORT_PATH LOG_DIR CACHE_DIR

# The optimizer ROOT, which is this script's PARENT directory -- it lives in tools/, and
# `.Renviron` and `settings.R` are at the root. Resolved from BASH_SOURCE rather than the cwd,
# so sourcing it from anywhere still finds them.
_opt_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"

# The `cd "$_opt_dir"` below is LOAD-BEARING, not tidiness. R reads `.Renviron` at the startup
# of EVERY R process -- there is no "has R run yet" state to worry about, a fresh Rscript
# reads it itself -- but it looks for `./.Renviron` in the CURRENT WORKING DIRECTORY (then
# ~/.Renviron). The optimizer's `.Renviron` lives at the ROOT, so without the cd, running this
# script from anywhere else silently misses OPTIMIZER_HOME and returns local-mode paths.
#
# `.Renviron` also OVERRIDES a variable already exported in the shell, so it is authoritative
# and cannot conflict with .bashrc.
#
# NEVER pass `--vanilla` to any R invocation in this project -- not the call below, and above
# all not `Rscript run_optimizer.R`. R's help: "--vanilla: Combine --no-save, --no-restore,
# --no-site-file, --no-init-file and --no-environ", and `--no-environ` means "Don't read the
# site and user environment files" -- i.e. skip `.Renviron`. That file carries T3_USERNAME and
# T3_PASSWORD as well as OPTIMIZER_HOME, so `--vanilla` does not merely mislocate the state
# directory: the optimizer cannot log in to T3 at all.
_opt_vals="$(cd "$_opt_dir" && Rscript -e '
  here::i_am("run_optimizer.R")
  suppressMessages(source("settings.R"))
  s <- optimizer_settings()
  cat(Sys.getenv("OPTIMIZER_HOME"), s$stop_file, s$db_path,
      s$report_path, s$log_dir, s$cache_dir, sep = "\n")
' 2>/dev/null)"

if [ -z "$_opt_vals" ]; then
  echo "tools/optimizer_paths.sh: could not read settings from R." >&2
  echo "  Run from the optimizer ROOT and check: Rscript -e 'source(\"settings.R\"); optimizer_settings()'" >&2
else
  # One path per line, read in the order cat() wrote them.
  {
    IFS= read -r OPTIMIZER_HOME
    IFS= read -r STOP_FILE
    IFS= read -r DB_PATH
    IFS= read -r REPORT_PATH
    IFS= read -r LOG_DIR
    IFS= read -r CACHE_DIR
  } <<EOF
$_opt_vals
EOF
  export OPTIMIZER_HOME STOP_FILE DB_PATH REPORT_PATH LOG_DIR CACHE_DIR
fi
unset _opt_dir _opt_vals

# Executed rather than sourced? Print what was found -- exporting would be pointless, since
# the exports die with this process.
if [ "${BASH_SOURCE[0]:-$0}" = "$0" ]; then
  echo "OPTIMIZER_HOME=$OPTIMIZER_HOME"
  echo "STOP_FILE=$STOP_FILE"
  echo "DB_PATH=$DB_PATH"
  echo "REPORT_PATH=$REPORT_PATH"
  echo "LOG_DIR=$LOG_DIR"
  echo "CACHE_DIR=$CACHE_DIR"
  echo
  echo "(these were printed, not exported -- use 'source ${0}' to get them in your shell)"
fi
