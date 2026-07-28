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
open_store <- function(path) {
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
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
  have <- DBI::dbListFields(con, "evals")
  add <- c(detail = "TEXT", study_name = "TEXT", program_name = "TEXT",
           location_name = "TEXT", year = "INTEGER")
  for (col in names(add)) {
    if (!(col %in% have))
      DBI::dbExecute(con, sprintf("ALTER TABLE evals ADD COLUMN %s %s", col, add[[col]]))
  }
  con
}

close_store <- function(con) invisible(DBI::dbDisconnect(con))

# Append one evaluation. `score` may be NA when the pipeline failed (status
# records why); failures are recorded too, so the optimizer learns not to revisit
# configurations that reliably break.
store_eval <- function(con, cfg, trial_id, scheme, score, n_test, status,
                       reason = NA_character_, detail = NA_character_,
                       seconds = NA_real_,
                       study_name = NA_character_, program_name = NA_character_,
                       location_name = NA_character_, year = NA_integer_) {
  invisible(DBI::dbExecute(con,
    "INSERT INTO evals
       (config_hash, config_json, trial_id, study_name, program_name, location_name, year,
        scheme, score, n_test, status, reason, detail, seconds, ts)
     VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
    params = list(config_hash(cfg), config_to_json(cfg), trial_id,
                  study_name %||% NA_character_, program_name %||% NA_character_,
                  location_name %||% NA_character_,
                  suppressWarnings(as.integer(year %||% NA_integer_)),
                  scheme,
                  ifelse(is.finite(score), score, NA_real_),
                  as.integer(n_test), status, reason,
                  detail %||% NA_character_,
                  as.numeric(seconds), format(Sys.time(), tz = "UTC"))))
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
