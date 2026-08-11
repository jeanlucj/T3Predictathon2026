# data_access.R
#
# The real-data layer: pull phenotypes and genotypes from T3/Wheat over BrAPI and cache them on
# disk, so a trial costs one download however many configurations touch it.
#
# BrAPI response column names vary by server version; .obs_tibble and trial_catalog normalise
# the common shapes and fail loudly on anything else.

library(tidyverse)

# --- category-partitioned cache paths --------------------------------------
# An entry lives at cache/<category>/<category>[_<identifier>].<ext>; identifier = NULL is a
# singleton category. Reads also try the pre-migration flat path, so an un-migrated cache still
# hits; writes always use the nested one. migrate_cache_layout.R relocates old files.
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
# Atomic: temp file then rename, so a reader never sees a half-written RDS -- LESSONS #24.
.cache_save <- function(settings, category, identifier = NULL, value, ext = "rds") {
  p <- .cache_path(settings, category, identifier, ext)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  # The pid keeps two writers apart; the leading dot keeps the temp out of the dosage_*.rds glob.
  tmp <- file.path(dirname(p), sprintf(".tmp%d_%s", Sys.getpid(), basename(p)))
  ok <- tryCatch({ saveRDS(value, tmp); file.rename(tmp, p) }, error = function(e) e)
  if (inherits(ok, "error") || !isTRUE(ok)) {
    unlink(tmp)
    if (inherits(ok, "error")) stop(ok)
    stop("could not write cache file ", p)
  }
  invisible(p)
}

# --- cross-worker locking --------------------------------------------------
# Workers share cache_dir, and a VCF download is the one operation they must not overlap on:
# each deletes any existing copy before downloading -- LESSONS #24.
#
# dir.create() is the primitive: atomic on POSIX, and FALSE rather than success when the
# directory exists, so exactly one caller wins.
.lock_dir <- function(settings, key)
  file.path(settings$cache_dir, "locks", paste0(key, ".lock"))

.acquire_lock <- function(settings, key) {
  p <- .lock_dir(settings, key)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  # Break a lock whose holder died. Judged on the lock's mtime, so lock_stale_minutes must
  # exceed a realistic download+parse or a live holder gets robbed.
  stale <- as.numeric(settings$lock_stale_minutes %||% 90)
  if (dir.exists(p)) {
    age <- as.numeric(difftime(Sys.time(), file.info(p)$mtime, units = "mins"))
    if (is.finite(age) && age > stale) {
      message(sprintf("lock %s is %.0f min old (> %g) -- assuming its holder died, breaking it",
                      key, age, stale))
      unlink(p, recursive = TRUE)
    }
  }
  if (!suppressWarnings(dir.create(p, showWarnings = FALSE))) return(FALSE)
  writeLines(sprintf("host=%s pid=%d worker=%s time=%s", Sys.info()[["nodename"]],
                     Sys.getpid(), settings$worker_id %||% "1", format(Sys.time())),
             file.path(p, "owner"))
  TRUE
}

.release_lock <- function(settings, key) invisible(unlink(.lock_dir(settings, key), recursive = TRUE))

# Run `expr` under the lock for `key`. If another worker holds it, wait for `ready()` to
# become TRUE (its result landing in the cache) rather than duplicating the work, and return
# `on_ready()`. Falls through to doing the work itself if the wait times out -- a slow peer
# must not mean this worker never makes progress.
.with_cache_lock <- function(settings, key, expr, ready = function() FALSE,
                             on_ready = function() NULL) {
  # Check before taking the lock: a free lock says nothing about whether the work is still
  # needed.
  if (isTRUE(ready())) return(on_ready())
  held <- .acquire_lock(settings, key)
  if (!held) {
    wait_s <- 60 * as.numeric(settings$lock_wait_minutes %||% 60)
    message(sprintf("waiting for another worker to finish %s (up to %.0f min)", key, wait_s / 60))
    waited <- 0
    while (waited < wait_s && !held) {
      Sys.sleep(15); waited <- waited + 15
      # The good case: the holder finished and its result is in the cache.
      if (isTRUE(ready())) {
        message(sprintf("%s: another worker produced it after %.0f min -- using that",
                        key, waited / 60))
        return(on_ready())
      }
      held <- .acquire_lock(settings, key)
    }
    if (!held) {
      # Proceed without it rather than fail: duplicated work is wasteful, never wrong.
      message(sprintf("%s: gave up waiting -- proceeding without the lock", key))
      return(force(expr))
    }
    message(sprintf("%s: took the lock after %.0f min of waiting", key, waited / 60))
  }
  on.exit(.release_lock(settings, key), add = TRUE)
  force(expr)
}

# --- BrAPI auth + retry ----------------------------------------------------
# How many attempts a flaky BrAPI call gets. Single source of the default.
.brapi_tries <- function(settings = NULL) as.integer((settings$brapi_tries %||% 4L))

# Log in from T3_USERNAME / T3_PASSWORD in the environment. Never prompts -- a prompt hangs a
# background run. conn$login() mutates the connection, storing the bearer token on it.
t3_login <- function(conn, settings = NULL) {
  user <- Sys.getenv("T3_USERNAME"); pass <- Sys.getenv("T3_PASSWORD")
  if (!nzchar(user) || !nzchar(pass))
    # Its own class: a setup error, so .brapi_try fails fast instead of retrying.
    stop(structure(
      class = c("t3_missing_credentials", "error", "condition"),
      list(message = paste("T3 login needs T3_USERNAME and T3_PASSWORD in the environment.",
             "Copy .Renviron.example to .Renviron, fill it in, and RESTART R --",
             ".Renviron is read only at startup, so editing it mid-session does nothing.",
             "Check with Sys.getenv(\"T3_USERNAME\")."),
           call = NULL)))
  conn$login(username = user, password = pass)
  # login() returns normally on a rejected password, leaving auth_token NULL -- LESSONS #27.
  if (!nzchar(conn$auth_token %||% ""))
    stop(structure(
      class = c("t3_bad_credentials", "error", "condition"),
      list(message = paste0(
             "T3 login was REJECTED for user '", user, "' -- the server issued no token ",
             "(look for an 'Incorrect Password' warning above). Fix T3_USERNAME / ",
             "T3_PASSWORD in .Renviron ON THIS MACHINE and RESTART R -- .Renviron is read ",
             "only at startup."),
           call = NULL)))
  invisible(conn)
}

# Construct a connection AND log it in -- the single place a connection is made.
t3_connect <- function(settings) {
  conn <- BrAPI::createBrAPIConnection(settings$brapi_host, is_breedbase = TRUE)
  t3_login(conn, settings)
  conn
}

# T3 surfaces an unauthenticated call as a WARNING plus an empty response, never an error, so
# auth failure is detected two ways: the warning text and the response status.
.is_auth_warning <- function(msg)
  grepl("unauthor|must login|permission to access|\\b401\\b", msg, ignore.case = TRUE)

.response_auth_failed <- function(r) {
  # [["status"]], not $status: $ on a tibble lacking the column warns.
  if (!is.list(r) || is.null(r[["status"]])) return(FALSE)
  # status is a flat httr::http_status list (wizard) or a per-page list of them (search), so
  # flatten either shape and scan.
  msgs <- tryCatch(as.character(unlist(r[["status"]], use.names = FALSE)),
                   error = function(e) character())
  any(grepl("unauthor|\\b401\\b", msgs, ignore.case = TRUE))
}

# Retry a BrAPI call through a flaky server, and re-authenticate on a 401. `thunk` is zero-arg,
# re-evaluated each attempt; transport errors retry with exponential backoff + jitter up to
# `tries`. A 401 triggers one re-login and a free retry that does not consume an attempt.
# conn = settings = NULL skips re-auth and leaves the pure retry loop.
.brapi_try <- function(thunk, conn = NULL, settings = NULL, tries = NULL,
                       base_delay = 2, what = "BrAPI call") {
  if (is.null(tries)) tries <- .brapi_tries(settings)
  tries <- max(1L, as.integer(tries))
  last <- NULL; relogins <- 0L; max_relogin <- 2L; attempt <- 0L
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

    # 401 -> re-login and retry the same thunk with the refreshed token. Bounded by
    # max_relogin so a token that stays unauthorized after login cannot loop forever.
    if (auth_fail && !is.null(conn) && !is.null(settings) && relogins < max_relogin) {
      relogins <- relogins + 1L
      lr <- tryCatch({ t3_login(conn, settings); TRUE }, error = function(e) { last <<- e; e })
      if (isTRUE(lr)) {
        message(sprintf("%s: was unauthorized -- logged in, retrying", what))
        attempt <- attempt - 1L; next                         # free retry (fresh token)
      }
      message(sprintf("%s: re-login FAILED: %s", what, conditionMessage(lr)))
      # A setup error, not a transient one: fail fast rather than burn the retry budget.
      if (inherits(lr, c("t3_missing_credentials", "t3_bad_credentials"))) stop(lr)
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

# --- cache backup / restore ------------------------------------------------
# Back up the cache to durable storage. Additive (cache files are write-once), so --delete is
# omitted and the transient raw_project/ VCFs are excluded. No-op without cache_backup_dir or
# rsync.
#
# Self-throttling, and any worker may do it: `min_age_minutes` is the interval, keyed on a stamp
# file every worker can see. One at a time, because this is an rsync over thousands of files
# rather than the store's millisecond copy -- LESSONS #25.
sync_cache_to_backup <- function(settings, quiet = TRUE,
                                 min_age_minutes = settings$cache_sync_minutes %||% 0) {
  dst <- settings$cache_backup_dir
  if (is.null(dst) || !nzchar(dst)) return(invisible(FALSE))
  src <- settings$cache_dir
  if (is.null(src) || !dir.exists(src)) return(invisible(FALSE))
  rsync <- Sys.which("rsync")
  if (!nzchar(rsync)) { message("cache backup: rsync not on PATH -- skipping"); return(invisible(FALSE)) }
  dir.create(dst, showWarnings = FALSE, recursive = TRUE)

  # backup_age_minutes() (R/store.R) is just "age of a file in minutes"; a missing stamp reads
  # as infinitely old, so the first call always syncs.
  stamp <- file.path(dst, ".last_sync")
  if (min_age_minutes > 0 && backup_age_minutes(stamp) < min_age_minutes)
    return(invisible(FALSE))
  file.create(stamp)                      # claim BEFORE the work, not after

  code <- suppressWarnings(system2(rsync,
    c("-a", "--exclude", "raw_project/", paste0(src, "/"), paste0(dst, "/")),
    stdout = if (quiet) FALSE else "", stderr = if (quiet) FALSE else ""))
  if (code != 0) message(sprintf("cache backup -> %s failed (rsync exit %d)", dst, code))
  invisible(code == 0)
}

# Fill the work cache from its durable backup, additively: rsync copies only what is missing,
# so this works whether the cache is empty or partly populated. Safe to call any time.
restore_cache_from_backup <- function(settings) {
  bak <- settings$cache_backup_dir
  if (is.null(bak) || !nzchar(bak) || !dir.exists(bak)) return(invisible(FALSE))
  src <- settings$cache_dir
  rsync <- Sys.which("rsync")
  if (!nzchar(rsync)) return(invisible(FALSE))
  dir.create(src, showWarnings = FALSE, recursive = TRUE)
  before <- length(list.files(src, recursive = TRUE))
  suppressWarnings(system2(rsync, c("-a", paste0(bak, "/"), paste0(src, "/")),
                           stdout = FALSE, stderr = FALSE))
  added <- length(list.files(src, recursive = TRUE)) - before
  if (added > 0L) message(sprintf("cache restore: +%d file(s) from backup %s", added, bak))
  invisible(TRUE)
}

# On-disk memoizer. `valid(val)` gates the WRITE, so a soft failure -- a 200 with an empty
# body, a degraded result -- is returned to this caller but not cached, and the next call
# retries. Default accepts everything.
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
    # get_all_trial_meta_data has no trait filter, so with focal_trait_db_id set we run the
    # same studies search ourselves with observationVariableDbIds -- LESSONS #1. Rows are built
    # with the package's own make_row_from_trial_result, so the columns match it exactly.
    #
    # Columns are clean_names() snake_case and carry no year or coordinates: year is derived
    # from start_date, lat/long/elev joined by location -- LESSONS #2.
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
# Restrict candidate focal trials to settings$target_domain. A NULL field is no constraint; a
# constraint on a column the catalogue lacks warns rather than dropping every trial.
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

# Descriptor for one specified trial id, for the diagnostics. A trial absent from the catalogue
# (the trait filter excluded it) still gets its accessions fetched, so it can be diagnosed.
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
# Focal-trait observations for a set of study ids: study_id, germplasm_name, value, unit_id,
# rep, block, col, row. One cache per study holds the JOINED all-trait table -- /observations
# carries the value, /observationunits the rep/block/coordinates -- so a hit avoids both
# searches. Trait-independent: the focal trait is filtered at read time.
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

# Normalise a /observations response, keeping ALL numeric traits; get_observations() selects
# the focal one. Records are the nested list under $combined_data, not $data -- LESSONS #5.
# observationUnitDbId is kept so get_observations() can join rep/block/coordinates.
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

# Normalise a /observationunits response to a per-plot table keyed by observationUnitDbId.
# col/row come from observationUnitPosition (kept for spatial analysis, not yet used); rep and
# block from its observationLevelRelationships, as {levelName, levelCode} entries.
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

# --- local wizards: shared primary maps, locally derived inversions --------
#
# Three BrAPI wizard calls answer every "who relates to what" question here:
#   trial   -> accessions   get_trial_accessions       (cache/acc/acc_<id>.rds)
#   project -> accessions   get_project_accessions     (cache/proj_acc/proj_acc_<id>.rds)
#   accessions -> trials / projects
#
# The last two are INVERSIONS of the first two, so only the first two are fetched. Both are
# stored one file PER KEY, so N workers writing distinct files never contend and any worker's
# fetch is immediately available to the others.
#
# The inversions are memoised in RAM, not persisted: they cost ~1 s to rebuild, and a persisted
# copy would need invalidating.
.index_memo <- new.env(parent = emptyenv())

# Source files for a primary map, honouring both cache layouts as .cache_existing does.
.index_files <- function(settings, category) {
  pat <- paste0("^", category, "_.*[.]rds$")
  fs <- c(list.files(file.path(settings$cache_dir, category), pattern = pat, full.names = TRUE),
          list.files(settings$cache_dir, pattern = pat, full.names = TRUE))
  fs[!duplicated(basename(fs))]                   # nested wins, as in .cache_existing
}

# Invert a primary map to `accession -> keys`, memoised per session. attr(, "keys") records
# which keys were seen, so a key absent from a later tabulation reads as "zero overlap" rather
# than "unknown", which would send every lookup back to the wizard.
.inverted_index <- function(settings, category) {
  fs <- .index_files(settings, category)
  if (!length(fs)) return(NULL)
  # Signature: count AND newest mtime. Count alone misses a REPLACED file. Formatted to a
  # string and compared with identical(), never all.equal().
  # The memo key includes cache_dir: two settings objects pointing at different caches must
  # not share an entry, or a coincidental count+mtime match would serve the wrong index.
  mkey <- paste0(category, "@", settings$cache_dir)
  # Signature over full PATHS plus size and mtime. Count and newest-mtime are both invariant
  # when a flat cache file is replaced by its nested twin, which would serve a stale index for
  # the session; the path set is not. Compared with identical(), never all.equal().
  inf  <- file.info(fs)
  sig  <- rlang::hash(list(fs, inf$size,
                           format(inf$mtime, "%Y-%m-%d %H:%M:%OS6")))
  memo <- .index_memo[[mkey]]
  if (!is.null(memo) && identical(memo$sig, sig)) return(memo$index)

  ids  <- sub(paste0("^", category, "_"), "", sub("[.]rds$", "", basename(fs)))
  vals <- lapply(fs, function(f) tryCatch(as.character(readRDS(f)), error = function(e) character()))
  keep <- lengths(vals) > 0
  if (!any(keep)) return(NULL)
  index <- split(rep(ids[keep], lengths(vals[keep])), unlist(vals[keep], use.names = FALSE))
  index <- lapply(index, unique)      # one key listed twice must not count twice
  attr(index, "keys") <- ids[keep]
  assign(mkey, list(sig = sig, index = index), envir = .index_memo)
  index
}

.trial_index   <- function(settings) .inverted_index(settings, "acc")
.project_index <- function(settings) .inverted_index(settings, "proj_acc")

# Keys ATTEMPTED but yielding nothing. get_trial_accessions / get_project_accessions
# deliberately do not cache an empty result (a soft failure must be retried), so a genuinely
# empty trial or project would otherwise hold coverage below 100% for ever.
.attempted <- function(settings, category) {
  p <- .cache_existing(settings, "attempted", category)
  if (is.na(p)) return(character())
  tryCatch(as.character(readRDS(p)), error = function(e) character())
}

.note_attempted <- function(settings, category, ids) {
  cur <- union(.attempted(settings, category), as.character(ids))
  .cache_save(settings, "attempted", category, cur)
  invisible(cur)
}

# Can this question be answered locally? Only if every key the UNIVERSE can offer is either
# indexed or known-attempted. Self-healing: a refreshed universe that gains keys drops
# coverage and the wizard resumes until the pre-warm is re-run -- there is no window in which
# a local answer is silently incomplete.
.index_covers <- function(idx, universe, settings, category) {
  if (is.null(idx) || !length(universe)) return(FALSE)
  all(universe %in% union(attr(idx, "keys") %||% character(), .attempted(settings, category)))
}

# Every genotyping project in the crop -- the coverage denominator for the projects index,
# as trial_catalog is for trials. One unfiltered wizard call; 110 ids on T3/Wheat.
.all_project_ids <- function(conn, settings) {
  cached(settings, "all_projects", NULL, max_age_days = 1,
         valid = function(v) length(v) > 0, expr = {
    w <- .brapi_try(function() conn$wizard("genotyping_projects", list()),
                    conn = conn, settings = settings, what = "all genotyping_projects wizard")
    unique(as.character(w$data$ids))
  })
}

# Accessions genotyped in ONE project. The mirror of get_trial_accessions, and the primary
# map the accession -> projects inversion is built from. The wizard answers this direction
# directly, which is what makes the projects index cheap: 110 calls rather than one per
# accession.
get_project_accessions <- function(project_id, conn, settings) {
  cached(settings, "proj_acc", as.character(project_id), max_age_days = 30,
         valid = function(a) length(a) > 0, expr = {
    w <- .brapi_try(function() conn$wizard("accessions",
                      list(genotyping_projects = list(as.character(project_id)))),
                    conn = conn, settings = settings, what = "accessions-by-project wizard")
    unique(as.character(w$data$names))
  })
}

# --- genotyping projects covering a set of accessions ----------------------
# Returns downloadable genotyping_project_ids -- the id space conn$vcf_archived() accepts, and
# NOT protocol ids -- LESSONS #4. Uses the wizard's `genotyping_projects` category, batching
# accessions so a large set does not overload one request.
projects_for_accessions <- function(accessions, conn, settings) {
  acc <- unique(as.character(accessions))
  if (!length(acc)) return(character())

  # Local answer when the index covers every project in the crop -- equivalent to the wizard,
  # and keyed on projects (~110) rather than on the accession set, so it hits once built.
  idx <- tryCatch(.project_index(settings), error = function(e) NULL)
  universe <- tryCatch(.all_project_ids(conn, settings), error = function(e) character())
  if (.index_covers(idx, universe, settings, "proj_acc")) {
    .note_geno_once("project_discovery_mode", sprintf(
      "project discovery: LOCAL (index covers all %d genotyping projects) -- no wizard call",
      length(universe)))
    return(unique(unlist(idx[intersect(acc, names(idx))], use.names = FALSE)))
  }

  # Fallback: ask the wizard. Reached until prewarm_indices.R has filled the projects map.
  .note_geno_once("project_discovery_mode",
                  "project discovery: WIZARD -- run prewarm_indices.R to answer this locally")
  # Cache only when every batch SUCCEEDED: an empty result from a swallowed HTTP 500, once
  # cached, hides a trial's genotypes entirely -- LESSONS #7. On failure return what we have
  # so this run proceeds, but do not persist it.
  batches <- split(acc, ceiling(seq_along(acc) / 500L))
  failed  <- FALSE
  ids <- unlist(purrr::map(batches, function(b) {
    w <- tryCatch(.brapi_try(function() conn$wizard("genotyping_projects", list(accessions = b)),
                             conn = conn, settings = settings, what = "genotyping_projects wizard"),
                  error = function(e) { failed <<- TRUE; NULL })
    if (is.null(w)) character() else as.character(w$data$ids)
  }, .progress = "Genotyping projects for accessions"), use.names = FALSE)
  unique(stats::na.omit(ids))
}

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
  if (!is.na(.cache_existing(settings, "unparseable", pid))) return(NULL)

  # The cached dosage if it is already as dense as this machine's budget wants, else NULL --
  # a cache thinned on a smaller machine would otherwise pin a bigger one to that thin. Costs
  # a one-time re-download; dosage_redensify = FALSE keeps whatever is cached.
  dense_enough <- function() {
    cd <- .find_densest_dosage(settings, pid)
    if (is.null(cd)) return(NULL)
    e_ideal <- cd$thin
    if (!isFALSE(settings$dosage_redensify) && cd$thin > 1L) {
      st <- tryCatch(readRDS(.cache_existing(settings, "stat", pid)), error = function(e) NULL)
      if (is.list(st))
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
  # Any coarser cache this parse supersedes. Removed only AFTER the new one is safely in
  # place (the two have different filenames), so there is never a moment with no cache --
  # which a concurrent reader would otherwise see as a cache miss and re-download.
  old <- .find_densest_dosage(settings, pid)
  d  <- .vcf_to_dosage(path, NULL, e)
  if (is.null(d)) {                                            # complete VCF, no genotypes
    .cache_save(settings, "unparseable", pid, list(reason = "no usable genotypes", when = Sys.time()))
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
