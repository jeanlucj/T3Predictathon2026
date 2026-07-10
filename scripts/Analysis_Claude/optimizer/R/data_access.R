# data_access.R
#
# The real-data layer: pull from T3/Wheat over BrAPI and cache everything on disk
# so the cost of touching a given trial (its phenotypes, its genotypes) is paid
# once and reused across configurations and iterations. This is what lets the
# optimizer run for days in the background without re-downloading.
#
# It reuses the proven calls from this repo:
#   * BrAPI: conn$search("/observations", ...), conn$wizard("accessions", ...),
#     conn$wizard("genotyping_projects", ...) for downloadable project ids,
#     conn$vcf_archived(...)  (see scripts/Prediction5, build_grm_for_cv00.R,
#     expand_project_universe.R)
#   * T3BrapiHelpers: get_all_trial_meta_data, find_other_studies_evaluating_same_germplasm,
#     covariance_combiner
#
# NOTE: the exact column names in BrAPI responses can vary by server version.
# The accessors below (.obs_tibble, trial_catalog) normalize the common shapes
# and fail loudly if a response is unexpected, so the first real run surfaces any
# mismatch immediately rather than corrupting results.

library(tidyverse)

# --- tiny on-disk memorizer -------------------------------------------------
.cache_path <- function(settings, name) file.path(settings$cache_dir, paste0(name, ".rds"))

cached <- function(settings, name, expr, max_age_days = Inf) {
  p <- .cache_path(settings, name)
  if (file.exists(p)) {
    age <- as.numeric(difftime(Sys.time(), file.info(p)$mtime, units = "days"))
    if (age <= max_age_days) return(readRDS(p))
  }
  val <- force(expr)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  saveRDS(val, p)
  val
}

# --- trial catalogue -------------------------------------------------------
# All trials in the crop, with the metadata we use to sample and to select
# training trials. Refreshed weekly.
trial_catalog <- function(conn, settings) {
  cached(settings, "trial_catalog", max_age_days = 7, expr = {
    # We want only trials that measured the focal trait. get_all_trial_meta_data
    # has no trait filter (it just does conn$search("studies", commonCropNames=)),
    # so when focal_trait_db_id is set we run that same studies search ourselves
    # with an observationVariableDbIds filter -- one bulk query that returns ONLY
    # focal-trait trials (e.g. 7561 instead of ~9030 for grain yield on T3/Wheat).
    # The rows are built with the package's own make_row_from_trial_result so the
    # columns match get_all_trial_meta_data exactly.
    #
    # Columns are janitor::clean_names() snake_case: study_db_id, study_name,
    # location_name, program_name, start_date, end_date, study_type, trial_db_id,
    # common_crop_name, experimental_design, create_date. There is no year /
    # latitude / longitude column, so we derive year from start_date (POSIXct)
    # and join lat/long/elev by location (see below).
    id <- settings$focal_trait_db_id
    meta <- if (!is.null(id) && nzchar(id)) {
      make_row <- getFromNamespace("make_row_from_trial_result", "T3BrapiHelpers")
      search <- conn$search("studies", body = list(
        commonCropNames        = settings$crop_name,
        observationVariableDbIds = list(as.character(id))))
      janitor::clean_names(dplyr::bind_rows(lapply(search$combined_data, make_row)))
    } else {
      T3BrapiHelpers::get_all_trial_meta_data(conn, settings$crop_name)
    }
    meta <- tibble::as_tibble(meta)
    if ("start_date" %in% names(meta) && !("year" %in% names(meta))) {
      meta <- meta |>
        dplyr::mutate(year = suppressWarnings(as.integer(format(start_date, "%Y"))))
    }
    # Trial metadata has no coordinates; fetch lat/long/elev per location and join
    # by location_name (cached with the catalogue). One /locations search for the
    # unique locations -> columns latitude, longitude, elevation.
    if ("location_name" %in% names(meta)) {
      locs <- unique(stats::na.omit(meta$location_name))
      coords <- if (length(locs)) tryCatch(
        T3BrapiHelpers::get_lat_long_elev_from_location_vec(as.list(locs), conn, id_or_name = "name"),
        error = function(e) NULL) else NULL
      if (!is.null(coords) && nrow(coords)) {
        coords <- coords |>
          dplyr::select(location_name, latitude, longitude, elevation) |>
          dplyr::distinct(location_name, .keep_all = TRUE)
        meta <- dplyr::left_join(meta, coords, by = "location_name")
      }
    }
    meta
  })
}

# --- focal trait matching --------------------------------------------------
# A focal_trait is "name|CO_id" (e.g. "Grain yield - kg/ha|CO_321:0001218").
# Split it so trials/observations can be matched on either part (BrAPI sometimes
# returns the variable as a name, sometimes as a CO id).
.focal_trait_parts <- function(focal_trait) {
  parts <- trimws(strsplit(focal_trait, "\\|")[[1]])
  list(name = parts[1], id = if (length(parts) > 1) parts[2] else NA_character_)
}
.matches_trait <- function(x, parts) {
  hit <- grepl(parts$name, x, fixed = TRUE)
  if (!is.na(parts$id)) hit <- hit | grepl(parts$id, x, fixed = TRUE)
  hit
}

# --- target-domain filter --------------------------------------------------
# Restrict candidate focal trials to the programs / years / locations the user
# wants the pipeline optimized for (settings$target_domain). A NULL field is no
# constraint; a constraint on a column the catalogue lacks is skipped with a
# warning rather than silently dropping all trials.
.apply_target_domain <- function(cand, td) {
  if (is.null(td)) return(cand)
  if (!is.null(td$programs)) {
    if ("program_name" %in% names(cand)) cand <- dplyr::filter(cand, program_name %in% td$programs)
    else warning("target_domain$programs set but no program_name column in catalogue")
  }
  if (!is.null(td$years)) {
    if ("year" %in% names(cand))
      cand <- dplyr::filter(cand, suppressWarnings(as.integer(year)) %in% as.integer(td$years))
    else warning("target_domain$years set but no year column in catalogue")
  }
  if (!is.null(td$locations)) {
    if ("location_name" %in% names(cand)) cand <- dplyr::filter(cand, location_name %in% td$locations)
    else warning("target_domain$locations set but no location_name column in catalogue")
  }
  if (!is.null(td$trials)) {
    if ("study_name" %in% names(cand)) cand <- dplyr::filter(cand, study_name %in% td$trials)
    else warning("target_domain$trials set but no study_name column in catalogue")
  }
  cand
}

# --- sample a random focal trial -------------------------------------------
# A trial is usable if it measured the focal trait, falls within the target
# domain, and has at least min_trial_acc genotyped accessions for which we can
# assemble a training set. We sample, validate, and retry a few times.
sample_real_trial <- function(settings, conn, max_tries = 12) {
  cat <- trial_catalog(conn, settings)
  # trial_catalog() already restricts to focal-trait trials when focal_trait_db_id
  # is set; the per-trial observation check in the loop below is the safety net
  # (and the only trait filter when focal_trait_db_id is NULL).
  cand <- cat |> dplyr::filter(!is.na(study_db_id))
  cand <- .apply_target_domain(cand, settings$target_domain)
  if (!nrow(cand)) fatal("No trials match the target domain in settings", "no_trials_in_domain")

  for (i in seq_len(max_tries)) {
    row <- cand[sample.int(nrow(cand), 1), ]
    id  <- as.character(row$study_db_id)
    acc <- tryCatch(get_trial_accessions(id, conn, settings), error = function(e) character())
    if (length(acc) < settings$min_trial_acc) next
    # Confirm the trial actually measured the focal trait (cached download).
    obs <- tryCatch(get_observations(id, conn, settings), error = function(e) NULL)
    if (is.null(obs) || !nrow(obs)) next
    return(.trial_descriptor(row, acc))
  }
  sample_failed(paste0("none of ", max_tries, " sampled trials had >= ",
                       settings$min_trial_acc, " accessions with focal-trait observations"))
}

# Build the trial descriptor from a catalogue row + its accessions. Shared by
# sample_real_trial() (random) and build_trial_descriptor() (a specific id, for
# canary checks and diagnose_trial).
.trial_descriptor <- function(row, acc) {
  list(
    id          = as.character(row$study_db_id),
    study_name  = as.character(row$study_name %||% row$study_db_id),
    accessions  = acc,
    n_acc       = length(acc),
    program     = as.character(row$program_name %||% NA),
    location    = as.character(row$location_name %||% NA),
    year        = suppressWarnings(as.integer(row$year %||% NA)),
    # lat/long/elev are joined into the catalogue by location_name in
    # trial_catalog(); they enable .trial_similarity()'s environmental branch.
    # Still NA for locations with no coordinates on record.
    lat         = suppressWarnings(as.numeric(row$latitude %||% NA)),
    long        = suppressWarnings(as.numeric(row$longitude %||% NA)),
    elev        = suppressWarnings(as.numeric(row$elevation %||% NA)),
    start_date  = as.character(row$start_date %||% NA))
}

# Descriptor for ONE specified trial id (not sampled). Used by the diagnostics.
# Looks the trial up in the catalogue; if it is not there (e.g. a focal trial
# whose trait filter excluded it) it still fetches accessions so the trial can be
# diagnosed. Errors only if the trial has no accessions at all.
build_trial_descriptor <- function(study_id, conn, settings) {
  id  <- as.character(study_id)
  cat <- trial_catalog(conn, settings)
  row <- cat |> dplyr::filter(as.character(study_db_id) == id)
  if (!nrow(row)) row <- tibble::tibble(study_db_id = id)   # minimal stub
  acc <- get_trial_accessions(id, conn, settings)
  if (!length(acc)) stop("trial ", id, " has no accessions (check the id / connection)")
  .trial_descriptor(row[1, , drop = FALSE], acc)
}

# --- phenotypes ------------------------------------------------------------
# Focal-trait observations for a set of study ids, as a tidy tibble
# (study_id, germplasm_name, value, unit_id, rep, block, col, row). One cache per
# study, obs_<sid>, holds the JOINED all-trait table (the /observations records
# joined to the /observationunits records on observationUnitDbId) -- so finding
# obs_<sid> avoids BOTH network searches. The join is here because /observations
# carries value + observationUnitDbId while rep/block/coordinates live only on
# /observationunits. The cache is trait-independent; the focal trait is filtered
# at read time, so changing focal_trait reuses it.
get_observations <- function(study_ids, conn, settings) {
  parts <- .focal_trait_parts(settings$focal_trait)
  study_ids <- unique(as.character(study_ids))
  purrr::map_dfr(study_ids, function(sid) {
    all_obs <- cached(settings, paste0("obs_", sid), max_age_days = 30, expr = {
      o_resp <- tryCatch(
        conn$search("/observations", body = list(studyDbIds = list(sid), pageSize = 100000L)),
        error = function(e) NULL)
      u_resp <- tryCatch(
        conn$search("/observationunits", body = list(studyDbIds = list(sid), pageSize = 100000L)),
        error = function(e) NULL)
      dplyr::left_join(.obs_tibble(o_resp, sid), .obsunits_tibble(u_resp), by = "unit_id")
    })
    all_obs |>
      dplyr::filter(.matches_trait(trait, parts)) |>
      dplyr::select(-trait)
  },
  .progress = "Observations from study ids")
}

# Normalize a BrAPI /observations search response to our tibble, keeping ALL
# numeric traits (with a `trait` column); get_observations() selects the focal
# trait. The records are a list under $combined_data (the auto-paginated result
# of conn$search), each a nested list with fields observationVariableName,
# germplasmName, value, observationUnitDbId, studyDbId, ... (NOT a data frame --
# $data holds search metadata, not the rows). rep/block/coordinates are NOT on the
# /observations records; we keep observationUnitDbId here and get_observations()
# joins them from the /observationunits query (see .obsunits_tibble).
.obs_tibble <- function(resp, sid) {
  empty <- tibble::tibble(study_id = character(), germplasm_name = character(),
                          trait = character(), value = numeric(), unit_id = character())
  if (is.null(resp)) return(empty)
  recs <- resp$combined_data %||% resp$data
  if (is.null(recs) || !length(recs)) return(empty)
  g1 <- function(o, ...) { for (k in c(...)) if (!is.null(o[[k]])) return((o[[k]])[1]); NULL }
  ch <- function(x) as.character(x %||% NA_character_)
  tibble::tibble(
    study_id       = sid,
    germplasm_name = vapply(recs, function(o) ch(g1(o, "germplasmName", "germplasmDbId")), character(1)),
    trait          = vapply(recs, function(o) ch(g1(o, "observationVariableName", "observationVariableDbId")), character(1)),
    value          = vapply(recs, function(o) suppressWarnings(as.numeric(g1(o, "value", "observationValue") %||% NA)), numeric(1)),
    unit_id        = vapply(recs, function(o) ch(g1(o, "observationUnitDbId")), character(1))
  ) |> dplyr::filter(is.finite(value))
}

# Normalize a BrAPI /observationunits search response to a per-plot table keyed by
# observationUnitDbId (matches the observations' observationUnitDbId). Captures the
# field position -- col = observationUnitPosition$positionCoordinateX,
# row = observationUnitPosition$positionCoordinateY (kept for possible spatial
# analysis, not yet used) -- and rep / block, which live in
# observationUnitPosition$observationLevelRelationships as {levelName, levelCode}
# entries (levelName "rep" -> the rep number in levelCode; "block" -> block number).
.obsunits_tibble <- function(resp) {
  empty <- tibble::tibble(unit_id = character(), col = character(), row = character(),
                          rep = character(), block = character())
  if (is.null(resp)) return(empty)
  recs <- resp$combined_data %||% resp$data
  if (is.null(recs) || !length(recs)) return(empty)
  ch  <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x[[1]])
  lvl <- function(rels, name) {          # levelCode where levelName == name
    if (is.null(rels) || !length(rels)) return(NA_character_)
    for (r in rels) if (identical(tolower(ch(r$levelName)), name)) return(ch(r$levelCode))
    NA_character_
  }
  out <- purrr::map_dfr(recs, function(o) {
    pos <- o$observationUnitPosition
    tibble::tibble(
      unit_id = ch(o$observationUnitDbId),
      col     = ch(pos$positionCoordinateX),
      row     = ch(pos$positionCoordinateY),
      rep     = lvl(pos$observationLevelRelationships, "rep"),
      block   = lvl(pos$observationLevelRelationships, "block"))
  },
  .progress = "Extract Observation Unit Position")
  dplyr::distinct(out, unit_id, .keep_all = TRUE)
}

# --- accessions for a trial ------------------------------------------------
get_trial_accessions <- function(study_id, conn, settings) {
  cached(settings, paste0("acc_", study_id), max_age_days = 30, expr = {
    w <- conn$wizard("accessions", list(trials = list(as.character(study_id))))
    unique(as.character(w$data$names))
  })
}

# --- genotyping projects covering a set of accessions ----------------------
# Returns downloadable genotyping_project_ids -- the id space conn$vcf_archived()
# accepts. NOTE: a genotyping PROJECT id is NOT a genotyping PROTOCOL id; the two
# id spaces are distinct and vcf_archived() needs the project id. We use the
# breeder-search wizard's `genotyping_projects` category (the proven path in
# expand_project_universe.R), batching accessions so large sets don't overload a
# single request. (get_geno_protocol_from_germ_vec returns protocol ids and must
# NOT be used here.)
projects_for_accessions <- function(accessions, conn, settings) {
  key <- paste0("proj_", substr(rlang::hash(sort(unique(accessions))), 1, 12))
  cached(settings, key, max_age_days = 30, expr = {
    acc <- unique(as.character(accessions))
    batches <- split(acc, ceiling(seq_along(acc) / 500L))
    ids <- unlist(purrr::map(batches, function(b) {
      w <- tryCatch(conn$wizard("genotyping_projects", list(accessions = b)),
                    error = function(e) NULL)
      if (is.null(w)) character() else as.character(w$data$ids)
    },
    .progress = "Genotyping projects for accessions"), use.names = FALSE)
    unique(stats::na.omit(ids))
  })
}

# --- dosage matrix for one genotyping project ------------------------------
# Return the accessions x markers dosage (0/1/2) for `project_id`, restricted to
# keep_samples. Extracts the WHOLE project (all samples) once and caches that,
# keyed by project (+ marker_thin + the raw VCF's byte size); every later call --
# for any trial / any keep_samples -- reads that one cache and subsets at read
# time. So a project's VCF is downloaded and parsed exactly once; afterwards the
# large raw VCF is redundant and is DELETED to reclaim space.
#
# Robustness (a flaky/partial download used to poison the cache permanently):
#   * .ensure_project_vcf validates the VCF is COMPLETE before extraction; a
#     truncated file is re-downloaded once, else an error (retried next run);
#   * the full-dosage cache name carries the VCF byte size, so a later re-download
#     at a different size (were the VCF still present) would not be shadowed by a
#     stale dosage. The cache is looked up by pattern (ignoring size) so it works
#     WITHOUT the VCF -- which is what lets the VCF be deleted.
get_project_dosage <- function(project_id, keep_samples, conn, settings,
                               marker_thin = 1L) {
  thin_tag <- if (marker_thin > 1) paste0("_thin", marker_thin) else ""
  # Look up the cached FULL-project dosage without needing the raw VCF.
  pat  <- paste0("^dosage_", project_id, thin_tag, "_sz[0-9]+[.]rds$")
  hits <- list.files(settings$cache_dir, pattern = pat, full.names = TRUE)
  full <- if (length(hits)) {
    readRDS(hits[[which.max(file.info(hits)$mtime)]])          # newest, if several
  } else {
    path <- .ensure_project_vcf(project_id, conn, settings)    # download + validate
    sz   <- file.info(path)$size
    d    <- .vcf_to_dosage(path, NULL, marker_thin)            # NULL -> ALL samples
    if (!is.null(d)) {
      saveRDS(d, .cache_path(settings, paste0("dosage_", project_id, thin_tag, "_sz", sz)))
      # full-project dosage now cached -> the (large) raw VCF is redundant; delete it.
      base <- sub("[.]gz$", "", path)
      unlink(c(base, paste0(base, ".gz")))
    }
    d
  }
  if (is.null(full)) return(NULL)
  keep <- intersect(rownames(full), as.character(keep_samples))
  if (!length(keep)) return(NULL)
  full[keep, , drop = FALSE]
}

# Ensure a complete archived VCF for a project is on disk; return its path
# (.vcf or .vcf.gz). Re-downloads a missing or incomplete file once; errors if it
# still will not validate, so a transient failure is retried rather than cached.
.ensure_project_vcf <- function(project_id, conn, settings) {
  base  <- file.path(settings$cache_dir, paste0("raw_project_", project_id, ".vcf"))
  found <- function() { p <- c(base, paste0(base, ".gz")); p[file.exists(p)][1] }
  p <- found()
  if (!is.na(p) && .vcf_complete(p)) return(p)
  # missing or incomplete -> (re)download once
  for (q in c(base, paste0(base, ".gz"))) if (file.exists(q)) unlink(q)
  conn$vcf_archived(output = base, genotyping_project_id = project_id)
  p <- found()
  if (is.na(p) || !.vcf_complete(p)) {
    for (q in c(base, paste0(base, ".gz"))) if (file.exists(q)) unlink(q)
    stop("incomplete VCF download for project ", project_id,
         " (removed; will retry next run)")
  }
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

# VCF -> dosage (accessions x markers). keep_samples = NULL extracts ALL samples
# (used to cache the whole project once); a character vector subsets to those
# samples. Optional genome-wide marker thinning for very large projects.
# (Condensed from build_grm_for_cv00.R; same encoding and the same transposition fix.)
.vcf_to_dosage <- function(path, keep_samples, thin = 1L) {
  read_samples <- function(p) {
    con <- if (grepl("[.]gz$", p)) gzfile(p, "rt") else file(p, "rt"); on.exit(close(con))
    repeat { line <- readLines(con, 1); if (!length(line)) stop("no #CHROM in ", p)
      if (startsWith(line, "#CHROM")) return(strsplit(line, "\t")[[1]][-(1:9)]) }
  }
  samples  <- read_samples(path)
  keep_cols <- if (is.null(keep_samples)) 9 + seq_along(samples)
               else 9 + which(samples %in% keep_samples)
  if (!length(keep_cols)) return(NULL)
  read_path <- path
  if (thin > 1) {
    tmp <- tempfile(fileext = ".vcf"); on.exit(unlink(tmp), add = TRUE)
    ci <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt"); co <- file(tmp, "wt")
    seen <- 0L
    repeat { ls <- readLines(ci, 50000L); if (!length(ls)) break
      hdr <- startsWith(ls, "#"); gi <- seen + cumsum(!hdr)
      writeLines(ls[hdr | (!hdr & ((gi - 1L) %% thin == 0L))], co); seen <- seen + sum(!hdr) }
    close(ci); close(co); read_path <- tmp
  }
  vcf <- vcfR::read.vcfR(read_path, cols = c(1:9, keep_cols), verbose = FALSE)
  if (is.null(vcf@gt) || ncol(vcf@gt) < 2) return(NULL)
  gt <- vcf@gt[, -1, drop = FALSE]          # drop=FALSE: keep a matrix even for 1 sample
  gt_dim <- dim(gt)
  if (any(grepl(":", gt))) gt <- sub(":.*", "", gt)
  gt <- gsub("\\|", "/", gt)
  dim(gt) <- gt_dim                          # restore dims that the string ops may strip
  num <- matrix(NA_integer_, nrow = nrow(gt), ncol = ncol(gt))
  num[gt %in% "0/0"] <- 0L; num[gt %in% c("0/1", "1/0")] <- 1L; num[gt %in% "1/1"] <- 2L
  mat <- t(num)                                    # transpose AFTER filling (the bugfix)
  colnames(mat) <- vcf@fix[, "ID"]; rownames(mat) <- colnames(vcf@gt)[-1]
  mat
}
