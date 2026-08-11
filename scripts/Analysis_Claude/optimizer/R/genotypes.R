# genotypes.R
#
# The genotype subsystem: get a project's archived VCF from T3, parse it into a dosage matrix,
# and cache that matrix. Split out of data_access.R because it is self-contained and is the
# piece most likely to graduate into T3BrapiHelpers.
#
# The seam with data_access.R runs one way. This file USES the cache primitives (.cache_path,
# .cache_save, .cache_existing, .neg_cache_*, .with_cache_lock) and .brapi_try; nothing there
# calls in except through get_project_dosage() and .project_stat(). Keep it that way.
#
# Cannot be exercised offline beyond the download-budget tests: everything else needs a real
# archive. Run it once against a known project before trusting a change.

library(tidyverse)

# --- per-session VCF-download retry budget ---------------------------------
# Failed download attempts per project, in RAM: each retry tries less hard, and after
# vcf_max_download_attempts the project is skipped for the rest of the run -- LESSONS #10.
.vcf_download_fails <- new.env(parent = emptyenv())

# Given the prior failure count for a project, how hard to try THIS time.
.vcf_download_plan <- function(pid, settings) {
  n       <- .vcf_download_fails[[pid]] %||% 0L
  max_att <- as.integer(settings$vcf_max_download_attempts %||% 3L)
  base    <- .brapi_tries(settings)
  list(n = n,
       skip       = n >= max_att,           # give up for this session
       tries      = max(1L, base - n),      # fewer in-call retries each prior failure
       base_delay = max(0.5, 2 / (n + 1)))  # shorter backoff each prior failure
}
.note_vcf_download_fail <- function(pid)
  assign(pid, (.vcf_download_fails[[pid]] %||% 0L) + 1L, envir = .vcf_download_fails)
.clear_vcf_download_fail <- function(pid)
  if (!is.null(.vcf_download_fails[[pid]])) rm(list = pid, envir = .vcf_download_fails)

# --- dosage matrix for one genotyping project ------------------------------
# Return the accessions x markers dosage (0/1/2) for `project_id`, restricted to keep_samples.
# The whole project is extracted and cached once; later calls subset that cache at read time,
# and the raw VCF is deleted once parsed.
#
# Cache files per project:
#   * dosage_<pid>[_thin<e>]_sz<bytes>.rds -- parsed once at the densest thin the budget allows
#     (.cache_thin). Any coarser thin is a read-time column subset, never a re-download.
#   * stat_<pid>.rds -- {n_samples, n_markers}.
#   * unparseable_<pid>.rds -- negative cache; delete to retry -- LESSONS #9.
#
# A truncated download is transient and retried next run; a structural failure on a COMPLETE
# file is permanent and negative-cached -- LESSONS #7.

# The densest thin at which a project's dense matrix (n_samples x n_markers x 4 B) fits
# dosage_budget_bytes. Most projects fit at thin 1; only very large panels are thinned.
.cache_thin <- function(n_samples, n_markers, budget_bytes) {
  # as.numeric: n_samples * n_markers overflows 32-bit integer for big panels.
  fit <- ceiling(as.numeric(n_samples) * as.numeric(n_markers) * 4 / max(1, budget_bytes %||% 2e9))
  max(1, fit)
}

# The thin level encoded in a dosage cache filename (1 when there is no _thin tag).
.dosage_thin_from_name <- function(paths) {
  vapply(basename(paths), function(nm) {
    m <- regmatches(nm, regexpr("_thin[0-9]+", nm))
    if (length(m)) as.integer(sub("_thin", "", m)) else 1L
  }, integer(1), USE.NAMES = FALSE)
}

# {n_samples, n_markers} recorded when a project's VCF was parsed; NULL if it never was. The
# sanctioned way to read the stat cache from outside this file -- R/pipeline.R sizes its
# per-project thinning plan from it.
.project_stat <- function(settings, project_id) {
  p <- .cache_existing(settings, "stat", as.character(project_id))
  if (is.na(p)) return(NULL)
  st <- tryCatch(readRDS(p), error = function(e) NULL)
  if (is.list(st)) st else NULL
}

# The DENSEST cached dosage for a project (smallest thin), across the nested cache/dosage/
# folder and the legacy flat cache/. Returns list(path, thin) or NULL. Normally there is
# exactly one dosage per project; if several thin levels linger from pre-change caches, the
# densest is chosen and the rest are redundant (prune_dosage_cache.R removes them).
.find_densest_dosage <- function(settings, project_id) {
  pat  <- paste0("^dosage_", project_id, "(_thin[0-9]+)?_sz[0-9]+[.]rds$")
  dirs <- c(file.path(settings$cache_dir, "dosage"), settings$cache_dir)
  hits <- unlist(lapply(dirs, function(d) list.files(d, pattern = pat, full.names = TRUE)))
  if (!length(hits)) return(NULL)
  thins <- .dosage_thin_from_name(hits)
  i <- which.min(thins)
  list(path = hits[[i]], thin = thins[[i]])
}

get_project_dosage <- function(project_id, keep_samples, conn, settings,
                               marker_thin = 1L) {
  pid <- as.character(project_id)
  marker_thin <- max(1L, as.integer(marker_thin))

  subset_samples <- function(full) {
    # keep_samples = NULL -> the WHOLE project population (marker QC + allele freqs need it).
    if (is.null(full) || is.null(keep_samples)) return(full)
    keep <- intersect(rownames(full), as.character(keep_samples))
    if (!length(keep)) NULL else full[keep, , drop = FALSE]
  }
  # Serve thin `marker_thin` from a cache at thin `e` by keeping every k-th marker, k =
  # max(1, floor(marker_thin / e)), so the result is at least as dense as requested.
  serve <- function(full, e) {
    if (is.null(full)) return(NULL)
    k <- max(1L, as.integer(floor(marker_thin / e)))
    if (k > 1L) full <- full[, seq(1L, ncol(full), by = k), drop = FALSE]
    subset_samples(full)
  }

  # Permanent negative cache: an archive already found unparseable. Delete
  # cache/unparseable/unparseable_<pid>.rds to force a retry.
  if (.neg_cache_hit(settings, "unparseable", pid)) return(NULL)

  # The cached dosage if it is already as dense as this machine's budget wants, else NULL --
  # a cache thinned on a smaller machine would otherwise pin a bigger one to that thin. Costs
  # a one-time re-download; dosage_redensify = FALSE keeps whatever is cached.
  dense_enough <- function() {
    cd <- .find_densest_dosage(settings, pid)
    if (is.null(cd)) return(NULL)
    e_ideal <- cd$thin
    if (!isFALSE(settings$dosage_redensify) && cd$thin > 1L) {
      st <- .project_stat(settings, pid)
      if (!is.null(st))
        e_ideal <- .cache_thin(st$n_samples, st$n_markers, settings$dosage_budget_bytes)
    }
    if (e_ideal >= cd$thin) cd else NULL
  }
  serve_cached <- function(cd) if (is.null(cd)) NULL else serve(readRDS(cd$path), cd$thin)

  hit <- dense_enough()
  if (!is.null(hit)) return(serve_cached(hit))                 # cache is dense enough

  coarse <- .find_densest_dosage(settings, pid)                # present but too coarse?
  if (!is.null(coarse))
    message(sprintf(
      "project %s: cached at every %dth marker, but this budget affords denser -- re-downloading to re-densify",
      pid, coarse$thin))

  # Miss, or too coarse -> download and parse, under the lock: two workers on one project
  # would delete each other's in-flight file.
  .with_cache_lock(settings, paste0("dosage_", pid),
    ready    = function() !is.null(dense_enough()),
    on_ready = function() serve_cached(dense_enough()),
    expr     = {
      # Re-check under the lock: a worker may have finished between the fast path and here.
      hit <- dense_enough()
      if (!is.null(hit)) serve_cached(hit)
      else .fetch_and_cache_dosage(pid, project_id, conn, settings, serve)
    })
}

# Download, measure and parse a project ONCE at the densest thin the budget allows, cache it,
# and return it served at the caller's requested thin. Callers hold the project's lock.
.fetch_and_cache_dosage <- function(pid, project_id, conn, settings, serve) {
  path <- .ensure_project_vcf(project_id, conn, settings)      # download + integrity check
  rm_raw <- function() { base <- sub("[.]gz$", "", path); unlink(c(base, paste0(base, ".gz"))) }
  st <- tryCatch(.vcf_stat(path), error = function(e) e)
  if (inherits(st, "error")) {                                 # not a standard VCF
    .neg_cache_mark(settings, "unparseable", pid, conditionMessage(st))
    rm_raw(); return(NULL)
  }
  .cache_save(settings, "stat", pid, list(n_samples = st$n_samples, n_markers = st$n_markers))
  e <- .cache_thin(st$n_samples, st$n_markers, settings$dosage_budget_bytes)
  if (!is.finite(e) || e < 1) {                                # unmeasurable -> treat as unusable
    .neg_cache_mark(settings, "unparseable", pid, "could not size the VCF")
    rm_raw(); return(NULL)
  }
  e <- as.integer(e)
  if (e > 1L)
    message(sprintf(
      "project %s: %d markers x %d samples exceeds the %.1f GB dosage budget; caching at every %dth marker",
      pid, st$n_markers, st$n_samples, (settings$dosage_budget_bytes %||% 2e9) / 1e9, e))
  sz <- file.info(path)$size
  # Any coarser cache this parse supersedes. Removed only AFTER the new one is safely in
  # place (the two have different filenames), so there is never a moment with no cache --
  # which a concurrent reader would otherwise see as a cache miss and re-download.
  old <- .find_densest_dosage(settings, pid)
  d  <- .vcf_to_dosage(path, NULL, e)
  if (is.null(d)) {                                            # complete VCF, no genotypes
    .neg_cache_mark(settings, "unparseable", pid, "no usable genotypes")
    rm_raw(); return(NULL)
  }
  tag <- if (e > 1L) paste0("_thin", e) else ""
  new_path <- .cache_save(settings, "dosage", paste0(pid, tag, "_sz", sz), d)
  if (!is.null(old) && !identical(normalizePath(old$path, mustWork = FALSE),
                                  normalizePath(new_path, mustWork = FALSE)))
    unlink(old$path)
  rm_raw()                                                     # cache written; raw redundant
  serve(d, e)
}

# Ensure a complete archived VCF for a project is on disk; return its path
# (.vcf or .vcf.gz). Re-downloads a missing or incomplete file once; errors if it
# still will not validate, so a transient failure is retried rather than cached.
.ensure_project_vcf <- function(project_id, conn, settings) {
  base    <- .cache_path(settings, "raw_project", project_id, ext = "vcf")        # nested target
  legacy  <- .cache_legacy_path(settings, "raw_project", project_id, ext = "vcf") # pre-migration flat
  nested  <- c(base, paste0(base, ".gz"))
  all_loc <- c(nested, legacy, paste0(legacy, ".gz"))
  # already on disk (nested or legacy) and complete?
  p <- all_loc[file.exists(all_loc)][1]
  if (!is.na(p) && .vcf_complete(p)) return(p)

  # Per-session download budget: a project that has already failed to download this run tries
  # less hard, and past the cap is skipped outright (T3 is likely down for that archive).
  pid  <- as.character(project_id)
  plan <- .vcf_download_plan(pid, settings)
  if (plan$skip)
    stop(sprintf("VCF download %s: skipped (failed %d times this session; T3 may be down)",
                 pid, plan$n))

  # missing or incomplete -> clear any stale copy and (re)download, retrying a transient failure
  for (q in all_loc) if (file.exists(q)) unlink(q)
  dir.create(dirname(base), showWarnings = FALSE, recursive = TRUE)
  ok <- tryCatch({
    .brapi_try(function() conn$vcf_archived(output = base, genotyping_project_id = project_id),
               conn = conn, settings = settings, tries = plan$tries, base_delay = plan$base_delay,
               what = sprintf("VCF download %s (attempt %d)", pid, plan$n + 1L))
    TRUE
  }, error = function(e) FALSE)
  p <- nested[file.exists(nested)][1]
  if (!ok || is.na(p) || !.vcf_complete(p)) {
    for (q in nested) if (file.exists(q)) unlink(q)
    .note_vcf_download_fail(pid)                # remember -> less effort next time, then skip
    stop("incomplete VCF download for project ", project_id,
         " (removed; will retry next run)")
  }
  .clear_vcf_download_fail(pid)                 # success -> forget prior failures
  p
}

# A VCF is "complete enough" if it has a #CHROM header naming >=1 sample AND at
# least one variant (data) line after it. Catches header-only / truncated /
# 0-sample files that otherwise parse to an empty or tiny dosage matrix.
.vcf_complete <- function(path) {
  con <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  saw_chrom <- FALSE; n_samples <- 0L
  repeat {
    line <- tryCatch(readLines(con, 1), error = function(e) character())
    if (!length(line)) break
    if (!saw_chrom) {
      if (startsWith(line, "#CHROM")) {
        saw_chrom <- TRUE
        n_samples <- length(strsplit(line, "\t")[[1]]) - 9L
      }
      next
    }
    if (nzchar(line) && !startsWith(line, "#")) return(saw_chrom && n_samples >= 1L)  # a data line
  }
  FALSE
}

# Read up to the #CHROM line from an OPEN connection, validate the standard VCF column
# layout, and return the split header. Rejects a transposed / non-VCF export (some T3
# archives store markers as columns -- the "#CHROM" row then holds chromosome labels and
# the next row starts with POS), which no VCF parser can read.
.vcf_header <- function(con) {
  repeat {
    line <- readLines(con, 1L)
    if (!length(line)) stop("no #CHROM header line")
    if (startsWith(line, "#CHROM")) {
      h <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(h) < 10L || h[2] != "POS" ||
          !identical(h[c(3L, 4L, 5L, 9L)], c("ID", "REF", "ALT", "FORMAT")))
        stop("not a standard VCF header (transposed or malformed): first fields ",
             paste(utils::head(h, 4L), collapse = ","))
      return(h)
    }
  }
}

# Validate and measure a VCF in one cheap streaming pass (no genotype parsing): sample
# names from the header, variant count from the body. Used to pick the marker thinning
# before the heavier dosage pass, and to reject non-VCF archives early.
.vcf_stat <- function(path) {
  con <- file(path, "rt"); on.exit(close(con), add = TRUE)   # file() auto-decompresses gz/BGZF
  header <- .vcf_header(con)
  n <- 0L
  repeat {
    ls <- readLines(con, 100000L); if (!length(ls)) break
    ls <- ls[!is.na(ls)]                              # a huge file can yield undecodable lines
    n <- n + sum(!startsWith(ls, "#"))
  }
  list(samples = header[-(1:9)], n_samples = length(header) - 9L, n_markers = n)
}

# VCF -> dosage (accessions x markers, coded 0/1/2). keep_samples = NULL extracts ALL
# samples; `thin` keeps every thin-th variant. Streams the file in fixed-size chunks so
# peak memory is bounded by one chunk plus the thinned result, whatever the file size.
# gzip/BGZF decompresses transparently; a malformed variant line is skipped rather than failing
# the file, and a transposed/non-VCF header is rejected -- LESSONS #9.
.vcf_to_dosage <- function(path, keep_samples, thin = 1L) {
  con <- file(path, "rt"); on.exit(close(con), add = TRUE)
  header  <- .vcf_header(con)
  samples <- header[-(1:9)]
  keep_idx <- if (is.null(keep_samples)) seq_along(samples)
              else which(samples %in% keep_samples)
  if (!length(keep_idx)) return(NULL)
  scol   <- 9L + keep_idx
  nfield <- length(header)
  thin   <- max(1L, as.integer(thin))
  blocks <- list(); ids <- character(); seen <- 0L
  repeat {
    lines <- readLines(con, 20000L); if (!length(lines)) break
    lines <- lines[!is.na(lines)]                     # drop undecodable lines (huge files)
    lines <- lines[!startsWith(lines, "#")]; if (!length(lines)) next
    parts <- strsplit(lines, "\t", fixed = TRUE)
    good  <- lengths(parts) == nfield                 # skip malformed (ragged) variant lines
    gi    <- seen + cumsum(good); seen <- seen + sum(good)
    keeprow <- if (thin > 1L) good & ((gi - 1L) %% thin == 0L) else good
    parts <- parts[keeprow]; if (!length(parts)) next
    m  <- do.call(rbind, parts)                        # rows all have nfield fields
    ids <- c(ids, m[, 3L])
    g  <- sub(":.*", "", m[, scol, drop = FALSE])      # GT is the first ':'-subfield
    g  <- gsub("|", "/", g, fixed = TRUE)              # phased | -> unphased /
    num <- matrix(NA_integer_, nrow(g), ncol(g))       # anything but 0/0,0/1,1/0,1/1 -> NA
    num[g == "0/0"] <- 0L; num[g == "0/1" | g == "1/0"] <- 1L; num[g == "1/1"] <- 2L
    blocks[[length(blocks) + 1L]] <- num               # markers(chunk) x samples
  }
  if (!length(blocks)) return(NULL)
  mat <- t(do.call(rbind, blocks))                     # -> samples x markers
  rownames(mat) <- samples[keep_idx]; colnames(mat) <- ids
  mat
}
