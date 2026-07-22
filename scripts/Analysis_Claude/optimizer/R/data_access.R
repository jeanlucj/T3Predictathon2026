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
#   * T3BrapiHelpers: get_all_trial_meta_data, covariance_combiner
#
# NOTE: the exact column names in BrAPI responses can vary by server version.
# The accessors below (.obs_tibble, trial_catalog) normalize the common shapes
# and fail loudly if a response is unexpected, so the first real run surfaces any
# mismatch immediately rather than corrupting results.

library(tidyverse)

# --- tiny on-disk memorizer -------------------------------------------------
# --- category-partitioned cache paths --------------------------------------
# A cache entry lives at   cache/<category>/<category>[_<identifier>].<ext>
# so the (large, flat) cache is navigable by kind. `identifier = NULL` is a singleton
# category (e.g. "trial_catalog") whose file is just cache/<category>/<category>.<ext>.
# READS fall back to the pre-migration FLAT path (cache/<category>[_<identifier>].<ext>) so
# an un-migrated or half-migrated cache still hits; WRITES always use the nested path.
# (Run migrate_cache_layout.R to relocate the existing flat files once.)
.cache_stem <- function(category, identifier = NULL)
  if (is.null(identifier)) category else paste0(category, "_", identifier)

.cache_path <- function(settings, category, identifier = NULL, ext = "rds")   # nested (write)
  file.path(settings$cache_dir, category, paste0(.cache_stem(category, identifier), ".", ext))

.cache_legacy_path <- function(settings, category, identifier = NULL, ext = "rds")  # flat
  file.path(settings$cache_dir, paste0(.cache_stem(category, identifier), ".", ext))

# Path of an existing cache file for this key: nested if present, else legacy flat, else NA.
.cache_existing <- function(settings, category, identifier = NULL, ext = "rds") {
  np <- .cache_path(settings, category, identifier, ext)
  if (file.exists(np)) return(np)
  lp <- .cache_legacy_path(settings, category, identifier, ext)
  if (file.exists(lp)) return(lp)
  NA_character_
}

# Write `value` to the nested cache path, creating the category subfolder. Returns the path.
.cache_save <- function(settings, category, identifier = NULL, value, ext = "rds") {
  p <- .cache_path(settings, category, identifier, ext)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  saveRDS(value, p); invisible(p)
}

# --- BrAPI auth + retry ----------------------------------------------------
# How many attempts a flaky BrAPI call gets. Single source of the default (was duplicated
# as `settings$brapi_tries %||% 4L` at every call site).
.brapi_tries <- function(settings = NULL) as.integer((settings$brapi_tries %||% 4L))

# Log a connection in from environment credentials (T3_USERNAME / T3_PASSWORD, e.g. from
# .Renviron -- see .Renviron.example). Never prompts interactively (which would hang a
# background run); errors clearly if the credentials are absent. conn$login() mutates the
# connection in place, storing the bearer token on it.
t3_login <- function(conn, settings = NULL) {
  user <- Sys.getenv("T3_USERNAME"); pass <- Sys.getenv("T3_PASSWORD")
  if (!nzchar(user) || !nzchar(pass))
    stop("T3 login needs T3_USERNAME and T3_PASSWORD in the environment ",
         "(copy .Renviron.example to .Renviron, fill it in, restart R).")
  conn$login(username = user, password = pass)
  invisible(conn)
}

# Construct a connection AND log it in -- the single place a connection is made.
t3_connect <- function(settings) {
  conn <- BrAPI::createBrAPIConnection(settings$brapi_host, is_breedbase = TRUE)
  t3_login(conn, settings)
  conn
}

# The T3 server surfaces an unauthenticated call NOT as an R error but as a WARNING
# ("Unauthorized (HTTP 401)." / "You must login ...") plus a returned response whose
# $status reason is "Unauthorized" and whose data is empty. So auth failure is detected two
# ways (either suffices): the warning text, and the response status.
.is_auth_warning <- function(msg)
  grepl("unauthor|must login|permission to access|\\b401\\b", msg, ignore.case = TRUE)

.response_auth_failed <- function(r) {
  if (!is.list(r) || is.null(r$status)) return(FALSE)
  # $status varies by call type: a flat httr::http_status list (category/reason/message,
  # e.g. from conn$wizard) or a per-page LIST of those (from a paginated conn$search).
  # Flatten either shape to its character values and scan -- no structural assumptions.
  msgs <- tryCatch(as.character(unlist(r$status, use.names = FALSE)),
                   error = function(e) character())
  any(grepl("unauthor|\\b401\\b", msgs, ignore.case = TRUE))
}

# Retry a BrAPI network call through a flaky server, and re-authenticate on a 401.
# `thunk` is zero-arg (re-evaluated each attempt). Transport errors (timeout/refused) are
# caught and retried with exponential backoff + jitter, up to `tries`, then re-raised.
# An UNAUTHORIZED result triggers one re-login (t3_login) via `conn`/`settings` and an
# immediate free retry (it does not consume a transient attempt); if login does not resolve
# it, it is treated like a transient failure and finally raised. Pass conn = settings = NULL
# (the default) to skip re-auth -- then behaviour is the pure retry loop as before.
.brapi_try <- function(thunk, conn = NULL, settings = NULL, tries = NULL,
                       base_delay = 2, what = "BrAPI call") {
  if (is.null(tries)) tries <- .brapi_tries(settings)
  tries <- max(1L, as.integer(tries))
  last <- NULL; relogins <- 0L; attempt <- 0L
  while (attempt < tries) {
    attempt <- attempt + 1L
    auth_fail <- FALSE
    r <- withCallingHandlers(
      tryCatch(thunk(), error = function(e) { last <<- e; e }),
      warning = function(w) {
        if (.is_auth_warning(conditionMessage(w))) {          # swallow only auth warnings
          auth_fail <<- TRUE; invokeRestart("muffleWarning")  # (others propagate normally)
        }
      })
    is_err   <- inherits(r, "error")
    auth_fail <- auth_fail || (!is_err && .response_auth_failed(r))
    if (!is_err && !auth_fail) return(r)                      # success

    # 401 -> one free re-login, then retry the same thunk with the refreshed token.
    if (auth_fail && !is.null(conn) && !is.null(settings) && relogins < 1L) {
      relogins <- relogins + 1L
      if (tryCatch({ t3_login(conn, settings); TRUE }, error = function(e) { last <<- e; FALSE })) {
        message(sprintf("%s: was unauthorized -- logged in, retrying", what))
        attempt <- attempt - 1L; next
      }
    }
    if (attempt < tries) {
      msg <- if (is_err) conditionMessage(r) else "unauthorized"
      message(sprintf("%s failed (%d/%d): %s -- retrying", what, attempt, tries, msg))
      Sys.sleep(base_delay * 2^(attempt - 1L) + stats::runif(1))
    }
  }
  if (!is.null(last)) stop(last)
  stop(sprintf("%s failed: unauthorized and login did not resolve it", what))
}

# --- per-session VCF-download retry budget ---------------------------------
# A VCF download that keeps timing out (T3 choking on a big archive) must not be re-stormed
# on every trial that covers it. This tracks failed download attempts per project IN RAM
# (like .overlap_memo): each subsequent attempt tries less hard, and after
# settings$vcf_max_download_attempts the project is skipped for the rest of the run. The
# counter resets on a new run, so a transiently-unavailable project is retried fresh.
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

# On-disk memoizer. `valid(val)` gates the WRITE: only a result that passes it is
# persisted, so a soft failure (a 200-with-empty response, or a degraded result missing
# data a sub-fetch failed to supply) is returned for the current call but NOT cached --
# the next call retries instead of serving the bad answer for `max_age_days`. A hard error
# in `expr` propagates before the write is reached, so it is never cached either. (The
# default accepts everything, preserving old behaviour.)
cached <- function(settings, category, identifier = NULL, expr, max_age_days = Inf,
                   valid = function(v) TRUE) {
  hit <- .cache_existing(settings, category, identifier)
  if (!is.na(hit)) {
    age <- as.numeric(difftime(Sys.time(), file.info(hit)$mtime, units = "days"))
    if (age <= max_age_days) return(readRDS(hit))
  }
  val <- force(expr)
  if (isTRUE(valid(val))) .cache_save(settings, category, identifier, val)
  val
}

# --- trial catalogue -------------------------------------------------------
# All trials in the crop, with the metadata we use to sample and to select
# training trials. Refreshed weekly.
trial_catalog <- function(conn, settings) {
  # Cache only a NON-EMPTY catalogue that, if it has locations, also has coordinates --
  # so a truncated/empty studies search or a failed lat/long fetch is retried, not stored
  # coordinate-less for 7 days (which would silently zero out environmental similarity).
  cat_valid <- function(m) is.data.frame(m) && nrow(m) > 0 &&
    (!("location_name" %in% names(m)) || "latitude" %in% names(m))
  cached(settings, "trial_catalog", max_age_days = 7, valid = cat_valid, expr = {   # singleton (no identifier)
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
      search <- .brapi_try(function() conn$search("studies", body = list(
        commonCropNames        = settings$crop_name,
        observationVariableDbIds = list(as.character(id)))),
        conn = conn, settings = settings, what = "studies search")
      janitor::clean_names(dplyr::bind_rows(lapply(search$combined_data, make_row)))
    } else {
      .brapi_try(function() T3BrapiHelpers::get_all_trial_meta_data(conn, settings$crop_name),
                 conn = conn, settings = settings, what = "trial metadata")
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
        .brapi_try(function() T3BrapiHelpers::get_lat_long_elev_from_location_vec(
          as.list(locs), conn, id_or_name = "name"),
          conn = conn, settings = settings, what = "location coords"),
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
    hit <- .cache_existing(settings, "obs", sid)
    all_obs <- if (!is.na(hit) &&
                   as.numeric(difftime(Sys.time(), file.info(hit)$mtime, units = "days")) <= 30) {
      readRDS(hit)
    } else {
      o_resp <- tryCatch(.brapi_try(function() conn$search(
        "/observations", body = list(studyDbIds = list(sid), pageSize = 100000L)),
        conn = conn, settings = settings, what = "observations search"), error = function(e) NULL)
      u_resp <- tryCatch(.brapi_try(function() conn$search(
        "/observationunits", body = list(studyDbIds = list(sid), pageSize = 100000L)),
        conn = conn, settings = settings, what = "observationunits search"), error = function(e) NULL)
      tb <- dplyr::left_join(.obs_tibble(o_resp, sid), .obsunits_tibble(u_resp), by = "unit_id")
      # Cache ONLY when BOTH fetches SUCCEEDED. A transient error returns NULL, which
      # .obs_tibble/.obsunits_tibble turn into an empty tibble -- caching that would store
      # "this trial has no phenotypes" for 30 days (the no_focal_obs / all-NA rep/block
      # signatures). A genuinely empty but SUCCESSFUL response is a real answer and is cached.
      if (!is.null(o_resp) && !is.null(u_resp)) .cache_save(settings, "obs", sid, tb)
      tb
    }
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
  # A real trial always has accessions; an empty result means a soft failure (200 with
  # empty data). Don't cache it -- retry next call. (A hard wizard error propagates and is
  # not cached either.)
  cached(settings, "acc", as.character(study_id), max_age_days = 30,
         valid = function(a) length(a) > 0, expr = {
    w <- .brapi_try(function() conn$wizard("accessions", list(trials = list(as.character(study_id)))),
                    conn = conn, settings = settings, what = "accessions wizard")
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
  acc <- unique(as.character(accessions))
  key <- substr(rlang::hash(sort(acc)), 1, 12)       # identifier: hash of the accession set
  hit <- .cache_existing(settings, "proj", key)
  if (!is.na(hit) &&
      as.numeric(difftime(Sys.time(), file.info(hit)$mtime, units = "days")) <= 30)
    return(readRDS(hit))

  # Cache ONLY a result assembled from wizard calls that all SUCCEEDED. A transient
  # server error (HTTP 500 / timeout) used to be swallowed to an empty result and then
  # cached for 30 days -- permanently hiding a trial's genotype data (its projects read
  # as zero). Now a failed batch sets `failed`, we still return what we have so the
  # current run proceeds, but we do NOT persist it, so the next run retries.
  batches <- split(acc, ceiling(seq_along(acc) / 500L))
  failed  <- FALSE
  ids <- unlist(purrr::map(batches, function(b) {
    w <- tryCatch(.brapi_try(function() conn$wizard("genotyping_projects", list(accessions = b)),
                             conn = conn, settings = settings, what = "genotyping_projects wizard"),
                  error = function(e) { failed <<- TRUE; NULL })
    if (is.null(w)) character() else as.character(w$data$ids)
  }, .progress = "Genotyping projects for accessions"), use.names = FALSE)
  ids <- unique(stats::na.omit(ids))
  if (!failed) .cache_save(settings, "proj", key, ids)
  ids
}

# --- dosage matrix for one genotyping project ------------------------------
# Return the accessions x markers dosage (0/1/2) for `project_id`, restricted to
# keep_samples. Extracts the WHOLE project (all samples) once and caches that; every
# later call -- for any trial / any keep_samples -- reads that one cache and subsets at
# read time. So a project's VCF is downloaded and parsed exactly once; afterwards the
# large raw VCF is redundant and is DELETED to reclaim space.
#
# Cache files per project (all under settings$cache_dir):
#   * dosage_<pid>[_thin<e>]_sz<bytes>.rds -- the project's dosage, parsed ONCE at the
#     densest thin `e` the memory budget allows (see .cache_thin). A config's requested
#     marker_thin is NOT baked in here; it is applied at read time by column-subsetting
#     this cache, so a project is never re-downloaded for a different thin. `_thin<e>` is
#     absent when e == 1 (full markers -- the case for every project that fits the budget).
#   * stat_<pid>.rds -- {n_samples, n_markers}, written on first extraction.
#   * unparseable_<pid>.rds -- a negative cache: an archive we found is not a usable VCF
#     (transposed/non-VCF layout, or no genotypes). Future calls skip it instead of
#     re-downloading and re-failing every run. Delete it to force a retry.
#
# Robustness (a flaky/partial download used to poison the cache permanently):
#   * .ensure_project_vcf validates the VCF is COMPLETE (not truncated) before use; a
#     truncated file is re-downloaded once, else an error (retried next run). That is a
#     TRANSIENT failure -- not negative-cached. A STRUCTURAL failure on a complete file
#     (.vcf_stat rejects it) is permanent -> negative cache.

# The densest thin a project can be cached at while its dense integer matrix
# (n_samples x n_markers x 4 bytes) fits settings$dosage_budget_bytes. A project is parsed
# ONCE at this thin; any coarser requested thin is derived by column-subsetting the cache
# (see get_project_dosage), so the budget is the SOLE control on cached marker density --
# set it deliberately. A normal project fits at thin 1 (full markers); only a very large
# panel (e.g. 7.5M markers x 683 samples = >20 GB) is thinned to fit.
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
  # Serve the requested thin from a cache at thin `e` by keeping every k-th marker, where
  # k = max(1, floor(marker_thin / e)) -- the largest multiple of e not exceeding the
  # request (so the served matrix is as dense as, or denser than, requested). When
  # marker_thin < e (a denser matrix than the budget allowed for the cache), k = 1 and the
  # cache is served as-is: it is already the densest available.
  serve <- function(full, e) {
    if (is.null(full)) return(NULL)
    k <- max(1L, as.integer(floor(marker_thin / e)))
    if (k > 1L) full <- full[, seq(1L, ncol(full), by = k), drop = FALSE]
    subset_samples(full)
  }

  # Permanent negative cache: an archive already found unparseable. Delete
  # cache/unparseable/unparseable_<pid>.rds to force a retry.
  if (!is.na(.cache_existing(settings, "unparseable", pid))) return(NULL)

  # A single cached dosage (the densest one) serves every requested thin by column-subset --
  # no re-download when a different config asks for a different thin.
  cd <- .find_densest_dosage(settings, pid)
  if (!is.null(cd)) return(serve(readRDS(cd$path), cd$thin))

  # Miss -> fetch, measure, and parse ONCE at the densest thin the budget allows.
  path <- .ensure_project_vcf(project_id, conn, settings)      # download + integrity check
  rm_raw <- function() { base <- sub("[.]gz$", "", path); unlink(c(base, paste0(base, ".gz"))) }
  st <- tryCatch(.vcf_stat(path), error = function(e) e)
  if (inherits(st, "error")) {                                 # not a standard VCF
    .cache_save(settings, "unparseable", pid, list(reason = conditionMessage(st), when = Sys.time()))
    rm_raw(); return(NULL)
  }
  .cache_save(settings, "stat", pid, list(n_samples = st$n_samples, n_markers = st$n_markers))
  e <- .cache_thin(st$n_samples, st$n_markers, settings$dosage_budget_bytes)
  if (!is.finite(e) || e < 1) {                                # unmeasurable -> treat as unusable
    .cache_save(settings, "unparseable", pid, list(reason = "could not size the VCF", when = Sys.time()))
    rm_raw(); return(NULL)
  }
  e <- as.integer(e)
  if (e > 1L)
    message(sprintf(
      "project %s: %d markers x %d samples exceeds the %.1f GB dosage budget; caching at every %dth marker",
      pid, st$n_markers, st$n_samples, (settings$dosage_budget_bytes %||% 2e9) / 1e9, e))
  sz <- file.info(path)$size
  d  <- .vcf_to_dosage(path, NULL, e)
  if (is.null(d)) {                                            # complete VCF, no genotypes
    .cache_save(settings, "unparseable", pid, list(reason = "no usable genotypes", when = Sys.time()))
    rm_raw(); return(NULL)
  }
  tag <- if (e > 1L) paste0("_thin", e) else ""
  .cache_save(settings, "dosage", paste0(pid, tag, "_sz", sz), d)
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
# peak memory is bounded by one chunk plus the (thinned) result, regardless of file size.
# Robust to the real T3 archives seen: gzip/BGZF is decompressed transparently by file();
# a malformed variant line (wrong field count) is skipped rather than failing the whole
# file; a transposed/non-VCF header is rejected. Encoding is byte-identical to the former
# vcfR-based reader (tests/test_subtasks.R checks this against a hand-built VCF).
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
