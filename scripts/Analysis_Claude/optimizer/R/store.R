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
# SEVERAL WORKERS MAY SHARE ONE STORE. Three settings make that safe, and they must be
# issued before anything else touches the database:
#   journal_mode = WAL  readers never block the writer and vice versa. WAL is a property of
#                       the FILE (it persists), and it CANNOT work on a network filesystem --
#                       it coordinates through an mmap'd -shm index, which NFS does not
#                       provide. That is why the store belongs on local disk (/workdir) with
#                       a periodic backup to durable storage; see settings$db_backup_path.
#                       We verify the pragma actually took and warn loudly if it did not,
#                       because the failure is silent and the consequence is corruption.
#   busy_timeout        without it SQLite returns SQLITE_BUSY *immediately* on contention
#                       rather than waiting. 60 s is enormous relative to the workload here
#                       (one small INSERT per multi-minute evaluation).
#   synchronous = NORMAL the durable copy is made by the backup, and a fsync per INSERT buys
#                       nothing against a workload whose unit of loss is one evaluation.
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
  # Migrations: add columns that stores created before them will lack. `detail`
  # holds the failure funnel; the four domain columns (study_name/program_name/
  # location_name/year) record the trial's target-domain attributes so the
  # surrogate can be trained on only the in-domain slice of this shared store.
  # They deliberately reuse the catalogue's column names so the exact same
  # target-domain predicate (.apply_target_domain) filters both sampling and evals.
  # peak_rss_mb is this evaluation's true peak RSS and the figure to size a machine from;
  # peak_r_mb is R's heap peak, an UNDER-estimate kept for the ratio (R/memory.R). `worker` names
  # which concurrent worker produced the row; dosage_budget records the marker-density
  # budget in force, WITHOUT which rows made at different densities are silently
  # incomparable (density is not a config parameter -- see settings$dosage_budget_bytes).
  # em_df_method is the same kind of hidden axis for em_combine: it names how each partial
  # covariance's EM weight was derived. NULL on rows written before 2026-07-31 means the old
  # accession-count weighting, which scores differently -- see EM_COMBINE_COMPARISON.md.
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

# Copy the live store to durable storage. Needed because WAL forces the working store onto
# LOCAL disk (see open_store), which on a cluster is the disk that gets wiped -- the store is
# the one file whose loss costs real work.
#
# `VACUUM INTO` is the right tool: it runs against a live, concurrently-written database and
# emits a single self-contained file with no WAL sidecar, so the backup is directly usable.
# It refuses to overwrite, hence the write-to-temp-then-rename (which also means a backup
# interrupted half-way never replaces a good one).
# `dest` is a FILE path, filename included (e.g. .../state/evals_backup.sqlite), not a
# directory -- see the paths block in settings.R.
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
