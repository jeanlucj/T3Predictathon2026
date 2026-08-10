# setup_fallback_libs.R
#
# Install the optimizer's packages into R_LIBS_USER, so a BARE `Rscript` works -- RStudio via
# Ceres OnDemand, or any session where you are not inside the container.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   module load r/4.5.3
#   Rscript setup_fallback_libs.R --dry-run     # show what would be installed
#   Rscript setup_fallback_libs.R               # 30-60 min: lme4, sommer, Matrix compile
#
# YOU PROBABLY DO NOT NEED THIS. The container already has everything, and running a script in
# it is one line:
#
#   ./container/run_in_container.sh exec peek_failures.R
#
# This exists for the cases where a container is awkward -- chiefly RStudio OnDemand, which
# hands you the module's R.
#
# VERSIONS ARE READ FROM container/optimizer.def, never restated here. Two copies of a version
# list drift, and drift means the container route and this one silently compute different
# results -- exactly what pinning exists to prevent. If the parse fails, this script stops
# rather than falling back to "whatever CRAN holds today".

# BASE R ONLY, deliberately. Every other script here starts with library(tidyverse) and
# here::i_am(); this one cannot. It runs when those packages are exactly what is missing, so
# depending on them would make it fail in the only situation it exists for. Do not "fix" this
# to match the house style.
if (!file.exists("setup_fallback_libs.R"))
  stop("run this FROM the optimizer directory:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript setup_fallback_libs.R\n",
       "  (working directory was: ", getwd(), ")")

args    <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

# ---- refuse where this makes no sense ------------------------------------
# Inside the container the packages are already present, and anything installed would land in
# the ephemeral overlay and vanish at exit.
if (dir.exists("/.singularity.d") || nzchar(Sys.getenv("APPTAINER_CONTAINER"))) {
  stop("This is running INSIDE the container, where every package is already installed.\n",
       "  Run it on the host instead -- or, more likely, you do not need it at all.")
}

lib <- Sys.getenv("R_LIBS_USER")
# R sets R_LIBS_USER to a default under ~ when it is unset, and on Ceres home is 30 GB.
if (!nzchar(lib) || startsWith(normalizePath(lib, mustWork = FALSE),
                               normalizePath(Sys.getenv("HOME"), mustWork = FALSE))) {
  stop("R_LIBS_USER is unset, or points inside your home directory (", lib, ").\n",
       "  Home is 30 GB on Ceres and this package set will not fit comfortably.\n",
       "  Add this to the .Renviron IN THIS DIRECTORY (not ~/.Renviron -- R reads only one,\n",
       "  the working directory's if it exists, and this project always runs from here):\n",
       "    R_LIBS_USER=/project/<account>/R_packages/%v\n",
       "  Then restart R -- .Renviron is read only at startup -- and check with .libPaths().")
}

# ---- read the manifest out of the container recipe -----------------------
def_path <- file.path("container", "optimizer.def")   # relative: we checked the cwd above
if (!file.exists(def_path)) stop("cannot find ", def_path)
def <- readLines(def_path, warn = FALSE)

# `pkgs <- c("a", "b", ...)` in %post, possibly wrapped over several lines.
.parse_pkgs <- function(x) {
  i <- grep("^pkgs <- c\\(", x)
  if (length(i) != 1L) return(character())
  j <- i; while (j <= length(x) && !grepl("\\)\\s*$", x[j])) j <- j + 1L
  blob <- paste(x[i:j], collapse = " ")
  unlist(regmatches(blob, gregexpr('"[^"]+"', blob))) |> gsub('"', '', x = _)
}
.parse_label <- function(x, key) {
  hit <- grep(paste0("^\\s*", key, "\\s"), x, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(paste0("^\\s*", key, "\\s+"), "", hit[1]))
}

pkgs     <- .parse_pkgs(def)
snapshot <- .parse_label(def, "CRAN\\.Snapshot")
sha_helpers <- .parse_label(def, "T3BrapiHelpers\\.Sha")
sha_brapi   <- .parse_label(def, "BrAPI\\.Sha")

# Fail loudly: a silent parse failure would install today's CRAN and quietly unpin everything.
if (!length(pkgs))            stop("could not parse the package list from ", def_path)
if (is.na(snapshot))          stop("could not parse CRAN.Snapshot from ", def_path)
if (is.na(sha_helpers) || is.na(sha_brapi))
  stop("could not parse the GitHub SHAs from ", def_path)

# The binary repo for this platform; the bare .../cran/<date> URL is source-only and much
# slower. On a cluster the codename may not match a PPM build, in which case source is fine.
codename <- if (file.exists("/etc/os-release")) {
  os <- readLines("/etc/os-release", warn = FALSE)
  hit <- grep("^VERSION_CODENAME=", os, value = TRUE)
  if (length(hit)) sub('^VERSION_CODENAME="?([^"]+)"?$', "\\1", hit[1]) else NA_character_
} else NA_character_
# Braces are load-bearing: at top level R ends the statement after the `if` branch when `else`
# starts a new line, and the file no longer parses.
repo <- if (!is.na(codename) && nzchar(codename)) {
  sprintf("https://packagemanager.posit.co/cran/__linux__/%s/%s", codename, snapshot)
} else {
  sprintf("https://packagemanager.posit.co/cran/%s", snapshot)
}

cat("\n== manifest, read from container/optimizer.def ==\n")
cat("  CRAN snapshot   :", snapshot, "\n")
cat("  repo            :", repo, "\n")
cat("  packages (", length(pkgs), "):\n    ", paste(pkgs, collapse = ", "), "\n", sep = "")
cat("  BrAPI           :", sha_brapi, "\n")
cat("  T3BrapiHelpers  :", sha_helpers, "\n")
cat("\n== target ==\n")
cat("  R_LIBS_USER     :", lib, "\n")
cat("  R               :", R.version.string, "\n")

if (dry_run) { cat("\n--dry-run: nothing installed.\n"); quit(status = 0) }

dir.create(lib, recursive = TRUE, showWarnings = FALSE)

# PUT IT ON THE SEARCH PATH, not just create it. R builds .libPaths() at STARTUP and silently
# drops entries that do not exist -- so on the first run, when R_LIBS_USER pointed at a
# directory this script had yet to create, `lib` was absent from .libPaths() for the whole
# session. install.packages(lib = lib) still wrote there correctly, and the completeness check
# below still found the files, but `remotes::install_github` a few lines on could not see the
# package it had just installed. Prepending here fixes it for this session; the next R start
# picks it up on its own, because by then the directory exists.
.libPaths(c(lib, .libPaths()))
if (!lib %in% .libPaths())
  stop("could not add ", lib, " to .libPaths() -- is it writable?")

options(repos = c(CRAN = repo),
        # PPM serves binaries only when the User-Agent identifies the platform; without this
        # it silently falls back to source and the install takes hours instead of minutes.
        HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
                                paste(getRversion(), R.version["platform"],
                                      R.version["arch"], R.version["os"])),
        Ncpus = max(1L, parallel::detectCores()))

t0 <- Sys.time()

# Install only what is absent, so re-running after an interruption costs seconds rather than
# another hour. install.packages() would otherwise re-download and rebuild everything.
todo <- setdiff(pkgs, rownames(installed.packages(lib.loc = lib)))
if (length(todo)) {
  cat("\n== installing", length(todo), "CRAN packages ==\n"); utils::flush.console()
  install.packages(todo, lib = lib)
} else {
  cat("\n== all", length(pkgs), "CRAN packages already present ==\n")
}

missing <- setdiff(pkgs, rownames(installed.packages(lib.loc = lib)))
if (length(missing))
  stop("FAILED to install: ", paste(missing, collapse = ", "))

cat("\n== installing the two GitHub packages, pinned ==\n"); utils::flush.console()
remotes::install_github(paste0("TriticeaeToolbox/BrAPI.R@", sha_brapi),
                        lib = lib, upgrade = "never")
remotes::install_github(paste0("jeanlucj/T3BrapiHelpers@", sha_helpers),
                        lib = lib, upgrade = "never")

# ---- verify, with the same assertions optimizer.def's %test makes --------
cat("\n== verifying ==\n")
for (p in c(pkgs, "BrAPI", "T3BrapiHelpers")) {
  if (p == "remotes") next
  ok <- suppressMessages(require(p, character.only = TRUE, quietly = TRUE))
  if (!ok) stop("installed but will not load: ", p)
}
# The private function R/data_access.R reaches for. A pin that drifted shows up here.
stopifnot(exists("make_row_from_trial_result", asNamespace("T3BrapiHelpers")))
for (f in c("get_all_trial_meta_data", "covariance_combiner",
            "get_lat_long_elev_from_location_vec"))
  stopifnot(exists(f, asNamespace("T3BrapiHelpers")))

cat(sprintf("\nOK -- %d packages into %s (%.0f min)\n",
            length(pkgs) + 2L, lib,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("Check it took:  Rscript -e '.libPaths(); library(here)'\n")
