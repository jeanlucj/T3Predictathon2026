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

# On-disk memoizer, and the ONE place a soft failure is kept out of the cache. Two ways to
# refuse the write, for the two things the decision can depend on:
#   valid(val)      the VALUE is not worth keeping (empty, degraded, missing a column).
#   no_cache(val)   something the value cannot show -- a failed fetch whose empty result is
#                   indistinguishable from a real one -- says do not keep it.
# Either way the value is returned to this caller and the next call retries.
cached <- function(settings, category, identifier = NULL, expr, max_age_days = Inf,
                   valid = function(v) TRUE) {
  hit <- .cache_existing(settings, category, identifier)
  if (!is.na(hit)) {
    age <- as.numeric(difftime(Sys.time(), file.info(hit)$mtime, units = "days"))
    if (age <= max_age_days) return(readRDS(hit))
  }
  val <- force(expr)
  vetoed <- isTRUE(attr(val, ".no_cache"))
  if (vetoed) attr(val, ".no_cache") <- NULL       # the caller must not see the marker
  if (!vetoed && isTRUE(valid(val))) .cache_save(settings, category, identifier, val)
  val
}

# Return this value but do not keep it. For use inside a cached() expr.
no_cache <- function(value) {
  if (is.null(value))
    stop("no_cache() cannot mark NULL; gate on the value instead: valid = function(v) !is.null(v)")
  structure(value, .no_cache = TRUE)
}

# Permanent negative cache: a key whose failure is STRUCTURAL rather than transient, so it must
# not be retried every run. Delete the file to force a retry -- LESSONS #9.
.neg_cache_hit  <- function(settings, category, key)
  !is.na(.cache_existing(settings, category, as.character(key)))
.neg_cache_mark <- function(settings, category, key, reason)
  .cache_save(settings, category, as.character(key), list(reason = reason, when = Sys.time()))

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
# Catalogue rows a focal trial may be drawn from: the crop's focal-trait trials restricted to
# the target domain. The universe both sample_real_trial() and eligible_trial_ids() work from.
.eligible_trials <- function(conn, settings) {
  cand <- trial_catalog(conn, settings) |> dplyr::filter(!is.na(study_db_id))
  .apply_target_domain(cand, settings$target_domain)
}

# The eligible trials' study ids. A config evaluated on all of them has nothing left to learn
# from; choose_config() drops it from the replication backlog.
eligible_trial_ids <- function(conn, settings) {
  unique(as.character(.eligible_trials(conn, settings)$study_db_id))
}

# When the catalogue this universe was resolved from was last refreshed.
.trial_catalog_ts <- function(settings) {
  f <- tryCatch(.cache_existing(settings, "trial_catalog", NULL), error = function(e) NA)
  if (is.na(f)) NA_character_ else format(file.info(f)$mtime, "%Y-%m-%d %H:%M:%S")
}

# The trial universe for this run, resolved ONCE and shared through the store. The leader
# resolves it from the catalogue and records the run; every other worker reads it back. The
# universe is a set of trial ids, so a catalogue that later reports a trial differently cannot
# change what this run searches. Returns list(run_id, universe, resumed); a NULL universe means
# unconstrained.
pin_trial_universe <- function(con, conn, settings, leader = TRUE, wait_s = 600) {
  build  <- settings$build %||% OPTIMIZER_BUILD
  run_id <- run_id_for(settings, build)

  row <- read_run(con, run_id)
  if (!is.null(row)) return(list(run_id = run_id, universe = run_universe(row), resumed = TRUE))

  if (isTRUE(settings$simulate)) {
    record_run(con, run_id, settings, build)
    return(list(run_id = run_id, universe = NULL, resumed = FALSE))
  }

  if (!leader) {
    waited <- 0
    while (is.null(row) && waited < wait_s) {
      Sys.sleep(5); waited <- waited + 5
      row <- read_run(con, run_id)
    }
    if (is.null(row))
      fatal(paste0("no run row after ", round(wait_s / 60), " min: the leader never recorded ",
                   "run ", run_id), "no_run_row")
    return(list(run_id = run_id, universe = run_universe(row), resumed = TRUE))
  }

  cand <- .eligible_trials(conn, settings)
  ids  <- unique(as.character(cand$study_db_id))
  if (!length(ids))
    fatal("no trials match the target domain in settings", "no_trials_in_domain")
  # The catalogue must account for every named trial before a multi-day search is pinned to it:
  # a trial that is absent here is one the run can never reach.
  want <- settings$target_domain$trials
  if (!is.null(want)) {
    missing <- setdiff(as.character(want), as.character(cand$study_name))
    if (length(missing))
      fatal(paste0(length(missing), " of ", length(want), " target_domain$trials are not in ",
                   "the resolved universe: ", paste(missing, collapse = ", ")),
            "domain_not_covered")
  }
  record_run(con, run_id, settings, build, ids,
             cand$study_name[match(ids, as.character(cand$study_db_id))],
             catalog_ts = .trial_catalog_ts(settings))
  list(run_id = run_id, universe = ids, resumed = FALSE)
}

sample_real_trial <- function(settings, conn, max_tries = 12) {
  # trial_catalog() already restricts to focal-trait trials when focal_trait_db_id is set; the
  # per-trial observation check in the loop below is the safety net (and the only trait filter
  # when focal_trait_db_id is NULL).
  cand <- .eligible_trials(conn, settings)
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
    study_name  = as.character(row$study_name %||% NA_character_),
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
    all_obs <- cached(settings, "obs", sid, max_age_days = 30, expr = {
      o_resp <- tryCatch(.brapi_try(function() conn$search(
        "/observations", body = list(studyDbIds = list(sid), pageSize = 100000L)),
        conn = conn, settings = settings, what = "observations search"), error = function(e) NULL)
      u_resp <- tryCatch(.brapi_try(function() conn$search(
        "/observationunits", body = list(studyDbIds = list(sid), pageSize = 100000L)),
        conn = conn, settings = settings, what = "observationunits search"), error = function(e) NULL)
      tb <- dplyr::left_join(.obs_tibble(o_resp, sid), .obsunits_tibble(u_resp), by = "unit_id")
      # Keep it only when BOTH fetches SUCCEEDED. A transient error returns NULL, which the
      # tibble builders turn into an empty table -- storing that means "this trial has no
      # phenotypes" for 30 days (LESSONS #7). A successful but genuinely empty response is a
      # real answer and is kept.
      if (is.null(o_resp) || is.null(u_resp)) no_cache(tb) else tb
    })
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

  # Fallback: ask the wizard. Reached until tools/prepare_indices.R has filled the projects map.
  .note_geno_once("project_discovery_mode",
                  "project discovery: WIZARD -- run tools/prepare_indices.R to answer this locally")
  # Deliberately NOT cached: the key would be an arbitrary accession set, which almost never
  # recurs, so a cache of it would hold stale answers without ever hitting. The projects index
  # above is what makes this question cheap. A failed batch yields nothing and the run proceeds
  # on what the others returned.
  batches <- split(acc, ceiling(seq_along(acc) / 500L))
  ids <- unlist(purrr::map(batches, function(b) {
    w <- tryCatch(.brapi_try(function() conn$wizard("genotyping_projects", list(accessions = b)),
                             conn = conn, settings = settings, what = "genotyping_projects wizard"),
                  error = function(e) NULL)
    if (is.null(w)) character() else as.character(w$data$ids)
  }, .progress = "Genotyping projects for accessions"), use.names = FALSE)
  unique(stats::na.omit(ids))
}

