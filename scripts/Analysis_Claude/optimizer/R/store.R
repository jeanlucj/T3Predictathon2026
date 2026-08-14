# store.R
#
# Persistent, resumable results store (SQLite). Every evaluation -- one
# configuration run on one trial under one CV scheme -- is a row. The optimizer's
# entire state is reconstructable from this one table, so the background process
# can be killed and resumed with no loss of progress. Configurations are stored
# as JSON (NA-preserving) plus a stable hash for de-duplication.

library(tidyverse)

# --- NA-preserving config <-> JSON ----------------------------------------
# jsonlite maps NA to null, which fromJSON then drops, silently corrupting a
# config on reload. We instead encode NA as a sentinel and decode each key back
# to its schema type (numeric params -> numeric, everything else -> character).
.NA_SENTINEL <- "__NA__"

config_to_json <- function(cfg) {
  cfg <- cfg[canonical_keys()]                       # stable key order
  enc <- lapply(cfg, function(v) if (length(v) == 0 || is.na(v)) .NA_SENTINEL else v)
  # digits = NA keeps full numeric precision so round-trips are lossless; the
  # canonical order makes the JSON string itself a stable config identity.
  jsonlite::toJSON(enc, auto_unbox = TRUE, digits = NA)
}

config_from_json <- function(json) {
  raw <- jsonlite::fromJSON(json, simplifyVector = TRUE)
  schema <- feature_schema()
  out <- lapply(names(raw), function(key) {
    v <- raw[[key]]
    if (length(v) == 0 || identical(as.character(v), .NA_SENTINEL)) return(NA)
    if (!is.null(schema[[key]]) && schema[[key]]$kind == "numeric") as.numeric(v) else as.character(v)
  })
  names(out) <- names(raw)
  out[canonical_keys()]
}

# --- store lifecycle -------------------------------------------------------
# SEVERAL WORKERS MAY SHARE ONE STORE. Three pragmas make that safe, issued before anything
# else touches the database -- LESSONS #24:
#   journal_mode = WAL   readers never block the writer. A property of the FILE, and it cannot
#                        work on a network filesystem, so the store belongs on local disk with
#                        a backup to durable storage. Verified, because failure is silent.
#   busy_timeout         without it SQLite returns SQLITE_BUSY immediately rather than waiting.
#   synchronous = NORMAL the backup provides durability; the unit of loss is one evaluation.
open_store <- function(path, busy_timeout_ms = 60000) {
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  DBI::dbGetQuery(con, sprintf("PRAGMA busy_timeout = %d", as.integer(busy_timeout_ms)))
  jm <- tryCatch(DBI::dbGetQuery(con, "PRAGMA journal_mode = WAL")$journal_mode[1],
                 error = function(e) NA_character_)
  if (!identical(tolower(jm %||% ""), "wal"))
    warning("could not put the store in WAL mode (got '", jm %||% "?", "') at ", path,
            ".\nThis is what a network filesystem looks like. ONE worker only, or move ",
            "db_path to local disk (/workdir) and back it up with db_backup_path.",
            call. = FALSE)
  # Retried: another process VACUUMing or writing holds a lock, and a bare dbExecute reports
  # "database is locked" rather than waiting.
  .with_busy_retry(function() DBI::dbExecute(con, "PRAGMA synchronous = NORMAL"))
  .with_busy_retry(function() DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS evals (
      id            INTEGER PRIMARY KEY AUTOINCREMENT,
      config_hash   TEXT,
      config_json   TEXT,
      trial_id      TEXT,
      study_name    TEXT,
      program_name  TEXT,
      location_name TEXT,
      year          INTEGER,
      scheme        TEXT,
      score         REAL,
      n_test        INTEGER,
      status        TEXT,
      reason        TEXT,
      detail        TEXT,
      seconds       REAL,
      ts            TEXT
    )"))
  .with_busy_retry(function()
    DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_hash ON evals(config_hash)"))
  # In-flight work, so a worker can see what the others are RUNNING -- `evals` holds only what
  # finished, and selection reading it alone let N workers start the same (config, trial).
  # The PRIMARY KEY is the interlock: two workers issue the same INSERT OR IGNORE, SQLite
  # permits one row, and the winner is the one that inserted it. See claim_eval().
  .with_busy_retry(function() DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS claims (
      config_hash TEXT,
      trial_id    TEXT,
      scheme      TEXT,
      worker      TEXT,
      host        TEXT,
      pid         INTEGER,
      ts          TEXT,
      PRIMARY KEY (config_hash, trial_id, scheme)
    )"))
  # Claims whose owner died, kept after the claim itself is reclaimed. Without this the record
  # of a killed worker survives only until the next choose_trial purges it -- seconds. No
  # primary key: the same pair may legitimately be lost more than once.
  .with_busy_retry(function() DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS claims_reaped (
      config_hash TEXT,
      trial_id    TEXT,
      scheme      TEXT,
      worker      TEXT,
      host        TEXT,
      pid         INTEGER,
      ts          TEXT,
      reaped_ts   TEXT
    )"))
  # Migrations: columns that older stores lack. `detail` holds the failure funnel. The four
  # domain columns reuse the catalogue's names, so .apply_target_domain filters sampling and
  # evals identically. peak_rss_mb is the figure to size a machine from; peak_r_mb is kept for
  # the ratio (R/memory.R). dosage_budget and em_df_method record HIDDEN axes -- neither is a
  # config parameter, so without them rows made under different settings would be averaged
  # together.
  have <- DBI::dbListFields(con, "evals")
  add <- c(detail = "TEXT", study_name = "TEXT", program_name = "TEXT",
           location_name = "TEXT", year = "INTEGER",
           peak_r_mb = "REAL", rss_mb = "REAL", worker = "TEXT", dosage_budget = "REAL",
           peak_rss_mb = "REAL", em_df_method = "TEXT", build = "TEXT")
  for (col in names(add)) {
    if (col %in% have) next
    # Two workers starting together both see the column missing and both try to add it.
    # The loser gets "duplicate column name" and must shrug, not die -- so swallow the
    # error and confirm the column exists afterwards either way.
    tryCatch(
      DBI::dbExecute(con, sprintf("ALTER TABLE evals ADD COLUMN %s %s", col, add[[col]])),
      error = function(e) {
        if (!(col %in% DBI::dbListFields(con, "evals"))) stop(e)   # a real failure
      })
  }
  con
}

close_store <- function(con) invisible(DBI::dbDisconnect(con))

# Copy the live store to durable storage: WAL forces it onto local disk, which on a cluster is
# the disk that gets wiped.
#
# VACUUM INTO runs against a live database and emits a self-contained file with no WAL sidecar,
# so the backup is directly usable. It refuses to overwrite, hence temp-then-rename -- which
# also means an interrupted backup never replaces a good one. `dest` is a FILE path.
backup_store <- function(con, dest) {
  if (is.null(dest) || !nzchar(dest)) return(invisible(FALSE))
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(dest, ".tmp", Sys.getpid())
  # PITFALL: file.rename does NOT raise an error when it fails -- it returns FALSE and signals
  # a warning. Handling only `error` here (as this once did) meant a backup that never happened
  # reported nothing at all, for the length of a run, and the whole reason db_path may sit on a
  # wipeable /workdir was quietly void. So the FALSE is checked explicitly.
  err <- NULL
  ok <- tryCatch({
    unlink(tmp)
    DBI::dbExecute(con, sprintf("VACUUM INTO '%s'", gsub("'", "''", tmp)))
    isTRUE(file.rename(tmp, dest))
  }, error = function(e) { err <<- conditionMessage(e); FALSE })
  if (!isTRUE(ok)) {
    unlink(tmp)
    message("store backup -> ", dest, " FAILED",
            if (!is.null(err)) paste0(": ", err) else "",
            # The likeliest cause of a bare rename failure, and worth naming outright.
            if (dir.exists(dest))
              " -- it is a DIRECTORY; db_backup_path must name a FILE, e.g. <dir>/evals_backup.sqlite"
            else "")
  }
  invisible(isTRUE(ok))
}

# Copy a store to `dest` WITH its -wal/-shm sidecars, and return dest. WAL keeps recent writes
# in the sidecar, so the main file alone can be 0 bytes -- copying it without them loses
# everything. Read the copy, never the live file.
.copy_store_with_sidecars <- function(path, dest) {
  invisible(file.copy(path, dest, overwrite = TRUE))
  for (ext in c("-wal", "-shm"))
    if (file.exists(paste0(path, ext)))
      invisible(file.copy(paste0(path, ext), paste0(dest, ext), overwrite = TRUE))
  dest
}

# Read-only shape of any store file: rows, distinct configs and trials, and the breakdowns that
# say whether a run will resume or restart. NULL when the file does not exist. Reads a sidecar
# copy, so it is safe against a live store.
store_summary <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  tmp <- .copy_store_with_sidecars(path, file.path(tempdir(), "summary.sqlite"))
  con <- DBI::dbConnect(RSQLite::SQLite(), tmp)
  on.exit({ DBI::dbDisconnect(con); unlink(paste0(tmp, c("", "-wal", "-shm"))) }, add = TRUE)
  e <- tryCatch(tibble::as_tibble(DBI::dbReadTable(con, "evals")), error = function(err) NULL)
  if (is.null(e)) return(NULL)
  col <- function(nm) if (nm %in% names(e)) e[[nm]] else rep(NA_character_, nrow(e))
  list(
    path = path, bytes = file.info(path)$size, mtime = file.mtime(path),
    rows = nrow(e), evals = e,
    n_config = dplyr::n_distinct(e$config_hash), n_trial = dplyr::n_distinct(e$trial_id),
    by_scheme = table(col("scheme"), useNA = "ifany"),
    by_status = table(col("status"), useNA = "ifany"),
    by_build  = table(col("build"),  useNA = "ifany"),
    ts_range  = if (nrow(e)) range(col("ts"), na.rm = TRUE) else c(NA_character_, NA_character_))
}

# Age of the backup in minutes; Inf when there is none, or when no backup is configured. What
# the report reads to say whether backups are actually happening.
backup_age_minutes <- function(dest) {
  if (is.null(dest) || !nzchar(dest) || !file.exists(dest)) return(Inf)
  as.numeric(difftime(Sys.time(), file.mtime(dest), units = "mins"))
}

# Rows in a store FILE; NA when it is absent or unreadable. Opens it directly rather than a
# copy: a VACUUM INTO backup has no -wal sidecar, and rename is atomic, so a reader sees either
# the old file or the new one.
.stored_rows <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NA_integer_)
  con <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), path), error = function(e) NULL)
  if (is.null(con)) return(NA_integer_)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  tryCatch(as.integer(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM evals")$n),
           error = function(e) NA_integer_)
}

# --- which store a READ-ONLY script should open ----------------------------
# On a cluster db_path is node-local scratch belonging to the job that wrote it, so a
# diagnostic run from any other allocation finds nothing there -- the durable copy is
# db_backup_path. Every read-only script used to hand-roll `if (!file.exists(db_path)) stop()`,
# which on Ceres reported "no store" while a perfectly good backup sat on /project.
#
# Order: an explicit path, then --store=, then the live store IF IT HAS ROWS, then the backup.
resolve_read_store <- function(explicit = NULL, what = "this script") {
  if (.on_login_node()) stop(.login_node_message(what), call. = FALSE)

  want <- function(p, why) {
    if (!file.exists(p)) stop("no store at ", p, " (", why, ")", call. = FALSE)
    message("store: ", p, " (", why, ")")
    p
  }
  # A flag is not a path. report_*.R pass argv[1] straight in for their documented positional
  # form, so `--store=<path>` would otherwise be taken as a filename by whichever check ran
  # first, and both spellings have to work in the same script.
  if (!is.null(explicit) && !is.na(explicit) && nzchar(explicit) &&
      !grepl("^--", explicit))
    return(want(explicit, "given on the command line"))
  flag <- grep("^--store=", commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(flag)) return(want(sub("^--store=", "", flag[1]), "from --store="))

  # No tryCatch: on a login node cluster_scratch_paths() raises the message that says how to
  # fix it, and swallowing that is what used to turn a clear diagnosis into "no store at ...".
  s <- optimizer_settings()

  # Rows, not existence. open_store() creates a full empty schema at whatever path it is given,
  # so a node that merely STARTED a job would otherwise shadow a good backup with an empty file.
  live <- s$db_path
  if (!is.null(live) && nzchar(live) && isTRUE(.stored_rows(live) > 0)) {
    message("store: ", live, " (live)")
    return(live)
  }
  bak <- s$db_backup_path
  if (!is.null(bak) && nzchar(bak) && file.exists(bak)) {
    age <- backup_age_minutes(bak)
    message("store: ", bak, "\n  Reading the BACKUP -- the live store (", live,
            ") is node-local to the job that wrote it and is not on this node.\n  ",
            .stored_rows(bak), " rows, written ",
            if (!is.finite(age)) "at an unknown time"
            else if (age < 90) sprintf("%.0f min ago", age) else sprintf("%.1f h ago", age / 60),
            ". Anything evaluated since that backup is NOT included.")
    return(bak)
  }
  stop("no store to read.\n  live  : ", live %||% "(not configured)",
       if (!is.null(live) && file.exists(live)) "  (exists but holds no rows)" else "  (absent)",
       "\n  backup: ", if (is.null(bak) || !nzchar(bak)) "(not configured)" else paste0(bak, "  (absent)"),
       "\n  Pass one explicitly with --store=<path> if it is somewhere else.", call. = FALSE)
}

# Merge the durable backup INTO the working store, the way restore_cache_from_backup() merges
# the cache: additive, so nothing already on the work disk is lost, and safe to call repeatedly.
# The leader calls it at startup, since db_path is node-local scratch that starts empty.
#
# Redundancy is keyed on (config_hash, trial_id, scheme) -- the pipeline is deterministic in
# those three, so a cell already present is the same answer. `scheme` is in the key because CV0
# and CV00 are different tasks scored separately. A collision keeps the LOCAL row: on a fresh
# node there are none, and mid-job the local row is the one this build computed.
#
# NOT keyed on build / dosage_budget / em_df_method, which do change the result (see the
# migration note in open_store). Skipped rows differing on one are counted and reported, so the
# mixing stays visible rather than being silently resolved either way.
restore_store_from_backup <- function(settings) {
  src <- settings$db_backup_path
  if (is.null(src) || !nzchar(src) || !file.exists(src)) return(invisible(NULL))
  dir.create(dirname(settings$db_path), showWarnings = FALSE, recursive = TRUE)
  con <- open_store(settings$db_path)               # also creates the schema ATTACH needs
  on.exit(close_store(con), add = TRUE)

  ok <- tryCatch({
    DBI::dbExecute(con, sprintf("ATTACH DATABASE '%s' AS bak", gsub("'", "''", src)))
    TRUE
  }, error = function(e) { message("store restore: cannot read ", src, " (", conditionMessage(e), ")"); FALSE })
  if (!ok) return(invisible(NULL))
  on.exit(try(DBI::dbExecute(con, "DETACH DATABASE bak"), silent = TRUE), add = TRUE)

  # Only the columns both schemas have: a backup from an older build lacks the migrated ones,
  # and those are exactly the backups worth rescuing. `id` is AUTOINCREMENT and is reassigned.
  cols <- setdiff(intersect(DBI::dbListFields(con, "evals"),
                            DBI::dbGetQuery(con, "SELECT * FROM bak.evals LIMIT 0") |> names()),
                  "id")
  q    <- paste(sprintf('"%s"', cols), collapse = ", ")
  n_bak <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM bak.evals")$n
  before <- n_evals(con)

  # MAX(id) collapses duplicate cells already in the backup: nothing ever prevented them, and a
  # NOT EXISTS against main cannot see rows the same statement is inserting.
  ins <- .with_busy_retry(function() DBI::dbExecute(con, sprintf(
    "INSERT INTO evals (%s)
     SELECT %s FROM bak.evals b
      WHERE b.id IN (SELECT MAX(id) FROM bak.evals GROUP BY config_hash, trial_id, scheme)
        AND NOT EXISTS (SELECT 1 FROM evals m
                         WHERE m.config_hash = b.config_hash
                           AND m.trial_id    = b.trial_id
                           AND m.scheme IS   b.scheme)",
    q, paste(sprintf('b."%s"', cols), collapse = ", "))))

  # Cells present on both sides but computed under different hidden axes -- worth seeing,
  # because both rows describe the same (config, trial, scheme) yet need not agree.
  hidden <- if (all(c("build", "dosage_budget") %in% cols)) {
    DBI::dbGetQuery(con,
      "SELECT COUNT(*) AS n FROM bak.evals b JOIN evals m
         ON m.config_hash = b.config_hash AND m.trial_id = b.trial_id AND m.scheme IS b.scheme
        WHERE IFNULL(m.build,'') != IFNULL(b.build,'')
           OR IFNULL(m.dosage_budget,-1) != IFNULL(b.dosage_budget,-1)")$n
  } else 0L

  message(sprintf("store restore: +%d row(s) from %s (%d of %d already present%s)",
                  ins, src, n_bak - ins, n_bak,
                  if (hidden > 0) sprintf("; %d differ on build/dosage_budget", hidden) else ""))
  invisible(list(backup_rows = n_bak, inserted = ins, skipped = n_bak - ins,
                 skipped_hidden_axis = hidden, before = before, after = n_evals(con)))
}

# Append one evaluation. `score` may be NA when the pipeline failed (status
# records why); failures are recorded too, so the optimizer learns not to revisit
# configurations that reliably break.
store_eval <- function(con, cfg, trial_id, scheme, score, n_test, status,
                       reason = NA_character_, detail = NA_character_,
                       seconds = NA_real_,
                       study_name = NA_character_, program_name = NA_character_,
                       location_name = NA_character_, year = NA_integer_,
                       peak_r_mb = NA_real_, rss_mb = NA_real_, peak_rss_mb = NA_real_,
                       worker = NA_character_, dosage_budget = NA_real_,
                       em_df_method = NA_character_, build = NA_character_) {
  invisible(.with_busy_retry(function() DBI::dbExecute(con,
    "INSERT INTO evals
       (config_hash, config_json, trial_id, study_name, program_name, location_name, year,
        scheme, score, n_test, status, reason, detail, seconds, ts,
        peak_r_mb, rss_mb, worker, dosage_budget, peak_rss_mb, em_df_method, build)
     VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
    params = list(config_hash(cfg), config_to_json(cfg), trial_id,
                  study_name %||% NA_character_, program_name %||% NA_character_,
                  location_name %||% NA_character_,
                  suppressWarnings(as.integer(year %||% NA_integer_)),
                  scheme,
                  ifelse(is.finite(score), score, NA_real_),
                  as.integer(n_test), status, reason,
                  detail %||% NA_character_,
                  as.numeric(seconds), format(Sys.time(), tz = "UTC"),
                  as.numeric(peak_r_mb %||% NA_real_), as.numeric(rss_mb %||% NA_real_),
                  as.character(worker %||% NA_character_),
                  as.numeric(dosage_budget %||% NA_real_),
                  as.numeric(peak_rss_mb %||% NA_real_),
                  as.character(em_df_method %||% NA_character_),
                  as.character(build %||% NA_character_)))))
}

# --- in-flight claims ------------------------------------------------------
# The pipeline is deterministic in (config, trial, scheme), so two workers running the same
# triple burn a node-hour to recompute a known number and then count it twice in that config's
# mean. Selection cannot prevent that on its own: it reads committed rows and is blind to what
# is running. These three functions are the missing piece.
#
# A config being picked by several workers at once is FINE and wanted -- a config owed more
# than one repeat should be spread across workers. Only the (config, trial) pair is exclusive.

# Is `pid` still running on this host? Used to decide whether a claim's owner died, which is
# the ONLY thing that makes a claim stale -- never elapsed time. A claim held for three days is
# a three-day evaluation, and expiring it would start a duplicate of work that is still going.
.pid_alive <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (!is.finite(pid) || pid <= 0) return(FALSE)
  if (dir.exists("/proc")) return(dir.exists(file.path("/proc", pid)))   # Linux: exact, cheap
  identical(suppressWarnings(system2("ps", c("-p", pid), stdout = FALSE, stderr = FALSE)), 0L)
}

# Drop claims whose owning process is gone. Only rows from THIS host are judged: a pid from
# another machine says nothing here, and all workers share one node by design, so that path is
# rare. PID reuse across a reboot is covered by clear_claims() at leader startup.
.purge_dead_claims <- function(con) {
  here <- DBI::dbGetQuery(con, "SELECT rowid, pid FROM claims WHERE host = ?",
                          params = list(Sys.info()[["nodename"]]))
  if (!nrow(here)) return(invisible(0L))
  dead <- here$rowid[!vapply(here$pid, .pid_alive, logical(1))]
  if (!length(dead)) return(invisible(0L))
  ids <- paste(as.integer(dead), collapse = ",")
  # COPY BEFORE DELETE. A dead claim is the only record that a worker was killed
  # mid-evaluation, and this purge runs on every choose_trial -- so the evidence is gone within
  # seconds of the death unless it is kept. `claims_reaped` is that history: what peek_workers.R
  # reports as "workers that died", and the answer to "has anything broken since I launched it".
  .with_busy_retry(function() DBI::dbExecute(con, sprintf(
    "INSERT INTO claims_reaped (config_hash, trial_id, scheme, worker, host, pid, ts, reaped_ts)
       SELECT config_hash, trial_id, scheme, worker, host, pid, ts, ?
         FROM claims WHERE rowid IN (%s)", ids),
    params = list(format(Sys.time(), tz = "UTC"))))
  invisible(.with_busy_retry(function() DBI::dbExecute(con,
    sprintf("DELETE FROM claims WHERE rowid IN (%s)", ids))))
}

# Take (config_hash, trial_id, scheme) unless a live worker holds it OR it is already evaluated
# under this build. TRUE means this worker owns it and must release_claim() when done.
#
# BOTH conditions are tested in the one statement, and that is the point. Testing only "is
# anyone running it" leaves a hole: the holder commits its row and releases, and a worker whose
# selection read the store just before that commit then claims a pair that is already done. The
# window is milliseconds and it fired repeatedly in a 4-worker simulate run. A caller's read is
# always stale by the time it acts on it, so the check has to happen where the write happens.
#
# Scoped to `build`: a build that invalidated a row (BUILD_CHANGES) must be able to recompute
# that pair. Which rows a build retires is R logic, not SQL, so this blocks only same-build
# duplicates -- the build-aware exclusion is choose_trial's `seen`, which does apply it.
claim_eval <- function(con, config_hash, trial_id, scheme, worker = NA_character_,
                       build = NA_character_) {
  take <- function() .with_busy_retry(function() DBI::dbExecute(con,
    "INSERT OR IGNORE INTO claims (config_hash, trial_id, scheme, worker, host, pid, ts)
     SELECT ?,?,?,?,?,?,?
      WHERE NOT EXISTS (SELECT 1 FROM evals
                         WHERE config_hash = ? AND trial_id = ? AND scheme IS ? AND build IS ?)",
    params = list(config_hash, as.character(trial_id), scheme,
                  as.character(worker %||% NA_character_), Sys.info()[["nodename"]],
                  Sys.getpid(), format(Sys.time(), tz = "UTC"),
                  config_hash, as.character(trial_id), scheme,
                  as.character(build %||% NA_character_))))
  if (isTRUE(take() == 1L)) return(TRUE)
  # Refused. Usually a live worker holds it or it is already done -- neither is recoverable.
  # A dead holder is, but it is the rare case, so the liveness sweep is paid for here rather
  # than on every claim: it costs a `ps` per claimed row where /proc is unavailable.
  .purge_dead_claims(con)
  isTRUE(take() == 1L)
}

release_claim <- function(con, config_hash, trial_id, scheme) {
  invisible(.with_busy_retry(function() DBI::dbExecute(con,
    "DELETE FROM claims WHERE config_hash = ? AND trial_id = ? AND scheme IS ?",
    params = list(config_hash, as.character(trial_id), scheme))))
}

# Every claim, as (config_hash, trial_id) pairs -- what selection must avoid.
active_claims <- function(con, scheme = NULL) {
  .purge_dead_claims(con)
  if (is.null(scheme))
    tibble::as_tibble(DBI::dbGetQuery(con, "SELECT config_hash, trial_id FROM claims"))
  else
    tibble::as_tibble(DBI::dbGetQuery(
      con, "SELECT config_hash, trial_id FROM claims WHERE scheme IS ?", params = list(scheme)))
}

# Every claim with all its columns, plus whether its owner is still running. For LOOKING at the
# table -- what each worker is doing, and which claims belong to a process that died.
#
# Deliberately does NOT purge, which is the whole reason it is separate from active_claims():
# a claim held by a dead pid IS the evidence that a worker was killed mid-evaluation, and
# purging while displaying would delete it as it was read. Adds hours_held so a long-running
# evaluation is visible as such -- length is not a fault (LESSONS #28).
claims_snapshot <- function(con) {
  cl <- tibble::as_tibble(DBI::dbGetQuery(con, "SELECT * FROM claims ORDER BY worker, ts"))
  if (!nrow(cl)) return(cl)
  here <- Sys.info()[["nodename"]]
  cl$alive <- vapply(seq_len(nrow(cl)), function(i)
    # A pid from another host says nothing here, so it is reported as unknown rather than dead.
    if (!identical(cl$host[i], here)) NA else .pid_alive(cl$pid[i]), logical(1))
  cl$hours_held <- as.numeric(difftime(Sys.time(), as.POSIXct(cl$ts, tz = "UTC"), units = "hours"))
  cl
}

# Deaths recorded by .purge_dead_claims and clear_claims, newest first -- the run's history of
# workers lost mid-evaluation.
reaped_claims <- function(con) {
  tibble::as_tibble(DBI::dbGetQuery(
    con, "SELECT * FROM claims_reaped ORDER BY reaped_ts DESC"))
}

# Wipe the table. The leader calls this at startup, when no worker can be mid-evaluation, which
# is what clears claims left by a job the scheduler killed on a durable db_path. (On a cluster
# db_path is node-local scratch and is wiped with the job, so there is nothing to clear.)
#
# Anything still here was held when the last job died, so it is recorded as a death like any
# other -- that is precisely the "the scheduler killed us mid-evaluation" case worth keeping.
clear_claims <- function(con) {
  .with_busy_retry(function() DBI::dbExecute(con,
    "INSERT INTO claims_reaped (config_hash, trial_id, scheme, worker, host, pid, ts, reaped_ts)
       SELECT config_hash, trial_id, scheme, worker, host, pid, ts, ? FROM claims",
    params = list(format(Sys.time(), tz = "UTC"))))
  invisible(.with_busy_retry(function() DBI::dbExecute(con, "DELETE FROM claims")))
}

# Belt-and-braces around the busy_timeout set in open_store: retry a write that still comes
# back locked. An evaluation costs minutes, so losing one to a lock that a few seconds of
# waiting would clear is a bad trade. Anything that is NOT a lock error is re-raised at once.
.with_busy_retry <- function(thunk, tries = 5L, base_delay = 1) {
  for (attempt in seq_len(tries)) {
    r <- tryCatch(thunk(), error = function(e) e)
    if (!inherits(r, "error")) return(r)
    if (!grepl("database is locked|database table is locked|SQLITE_BUSY",
               conditionMessage(r), ignore.case = TRUE)) stop(r)
    if (attempt == tries) stop(r)
    Sys.sleep(base_delay * 2^(attempt - 1L) + stats::runif(1))
  }
}

# All evaluations as a tibble (config_json kept as text; parse on demand).
read_evals <- function(con) {
  tibble::as_tibble(DBI::dbReadTable(con, "evals"))
}

n_evals <- function(con) {
  DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM evals")$n
}

# Configs already tried (by hash), to avoid re-running identical pipelines.
tried_hashes <- function(con) {
  DBI::dbGetQuery(con, "SELECT DISTINCT config_hash FROM evals")$config_hash
}

# (config, trial, scheme) cells holding more than one row. The pipeline is deterministic in
# that triple, so every one of these is a duplicated evaluation -- the thing claim_eval()
# prevents. Reported so the count can be watched: it should stop growing.
duplicate_cells <- function(con) {
  DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM (SELECT 1 FROM evals
                          GROUP BY config_hash, trial_id, scheme HAVING COUNT(*) > 1)")$n
}
