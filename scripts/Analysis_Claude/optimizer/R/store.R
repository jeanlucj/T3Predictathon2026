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
  DBI::dbExecute(con, "PRAGMA synchronous = NORMAL")
  DBI::dbExecute(con, "
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
    )")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_hash ON evals(config_hash)")
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

# Age of the backup in minutes; Inf when there is none, or when no backup is configured.
# The throttle keys on this rather than on a per-process timestamp so that every worker shares
# one interval -- see should_backup_now(). Also what the report reads to say whether backups
# are actually happening.
backup_age_minutes <- function(dest) {
  if (is.null(dest) || !nzchar(dest) || !file.exists(dest)) return(Inf)
  as.numeric(difftime(Sys.time(), file.mtime(dest), units = "mins"))
}

# Is it time to back the store up? A predicate rather than an inline condition so the rule is
# testable on its own.
#
# Two things it deliberately does NOT do: consult is_leader, and track its own last-backup
# time. The caller runs this between evaluations that can last many hours, so gating on one
# worker makes the interval a floor of that worker's slowest evaluation; and a per-process
# timestamp would let N workers each honour the interval separately. See LESSONS #25.
should_backup_now <- function(settings) {
  dest   <- settings$db_backup_path
  db_min <- settings$db_backup_minutes %||% 0
  if (is.null(dest) || !nzchar(dest) || db_min <= 0) return(FALSE)
  backup_age_minutes(dest) >= db_min
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
