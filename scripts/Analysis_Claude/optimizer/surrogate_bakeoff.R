# surrogate_bakeoff.R
#
# Which surrogate design predicts a held-out configuration best, ON YOUR DATA? Scores the
# competing designs by REPEATED cross-validation over configurations and reports each with an
# error bar and a paired test.
#
# WHY THIS EXISTS. On 2026-07-31 a single 5-fold split suggested trial_id blocking improved
# the surrogate (0.108 -> 0.135); twelve splits showed the difference was -0.001 (p 0.84).
# One split is noise. Simulations that day flipped sign three times for the same reason --
# under-replication. So: never decide this from one split, and always print how large a
# difference the current store can actually resolve.
#
# THE ARMS
#   A  pooled    per-config mean scores, forest on config features only (pre-0.7.5 behaviour)
#   B  blocked   one row per (config, trial), trial_id a factor, marginalised when scoring
#                (0.7.5; settings$surrogate_block_trial)
#   C  merf      mixed-effects random forest -- trial as a SHRUNKEN random effect OUTSIDE the
#                tree: fit the forest to y - Zb, refit b on y - f(X), iterate. This is the
#                actual analogue of fitting a many-level factor as random in an LMM, where the
#                group costs one variance parameter instead of competing for tree structure
#                (Hajjem, Bellavance & Larocque 2014). Needs lme4; skipped if absent.
#
# Read-only, no network. Copies the store AND its -wal/-shm sidecars and reads the copy.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript surrogate_bakeoff.R
#
# Options:
#   --reps=12        how many independent CV repetitions (more = tighter error bars)
#   --folds=5        folds per repetition
#   --ntree=200      trees per forest
#   --arms=A,B,C     which arms to score
#   --scheme=CV00    restrict to one CV scheme; default: the most common in the store

if (!file.exists("surrogate_bakeoff.R"))
  stop("run this FROM the optimizer directory:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript surrogate_bakeoff.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("surrogate_bakeoff.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", n, "="), "", h[1]) }
o_reps  <- as.integer(opt("reps", "12"))
o_folds <- as.integer(opt("folds", "5"))
o_ntree <- as.integer(opt("ntree", "200"))
o_arms  <- trimws(strsplit(opt("arms", "A,B,C"), ",")[[1]])
o_scheme <- opt("scheme")

s <- optimizer_settings()
if (!file.exists(s$db_path)) stop("no store at ", s$db_path)
tmp <- file.path(tempdir(), "bakeoff.sqlite")
invisible(file.copy(s$db_path, tmp, overwrite = TRUE))
for (ext in c("-wal", "-shm"))
  if (file.exists(paste0(s$db_path, ext)))
    invisible(file.copy(paste0(s$db_path, ext), paste0(tmp, ext), overwrite = TRUE))
con <- DBI::dbConnect(RSQLite::SQLite(), tmp)
e   <- tibble::as_tibble(DBI::dbReadTable(con, "evals")); DBI::dbDisconnect(con)

if (is.null(o_scheme)) o_scheme <- names(sort(table(e$scheme), decreasing = TRUE))[1]
d <- e |>
  dplyr::filter(scheme == o_scheme, status == "ok", is.finite(score)) |>
  dplyr::mutate(z = atanh(pmin(pmax(score, -0.999), 0.999)), w = pmax(n_test - 3, 1))
if (nrow(d) < 30) stop("only ", nrow(d), " scored rows for scheme ", o_scheme,
                       " -- too few to cross-validate")

F      <- configs_to_features(lapply(d$config_json, config_from_json))
trials <- sort(unique(d$trial_id))
hashes <- unique(d$config_hash)
# The target: each config's own weighted mean, i.e. what we would like to predict for a
# configuration we have not tried.
truth <- d |> dplyr::group_by(config_hash) |>
  dplyr::summarise(y = sum(w * z) / sum(w), .groups = "drop")

# --- the three fits ---------------------------------------------------------
fit_A <- function(tr) {
  ag <- d[tr, ] |> dplyr::group_by(config_hash) |>
    dplyr::summarise(y = sum(w * z) / sum(w), .groups = "drop")
  m <- fit_surrogate(F[match(ag$config_hash, d$config_hash), , drop = FALSE], ag$y,
                     ntree = o_ntree, min_obs = 8)
  if (is.null(m)) return(NULL)
  function(i) predict_surrogate(m, F[i, , drop = FALSE])$mean
}
fit_B <- function(tr) {
  FB <- F; FB$trial_id <- factor(d$trial_id, levels = trials)
  m <- fit_surrogate(FB[tr, , drop = FALSE], d$z[tr], ntree = o_ntree, min_obs = 8)
  if (is.null(m)) return(NULL)
  function(i) predict_surrogate_marginal(m, F[i, , drop = FALSE], trials)$mean
}
fit_C <- function(tr, iters = 3L) {
  if (!requireNamespace("lme4", quietly = TRUE)) return(NULL)
  dt <- d[tr, ]; Xt <- F[tr, , drop = FALSE]
  if (dplyr::n_distinct(dt$trial_id) < 2) return(NULL)
  y <- dt$z; m <- NULL
  for (it in seq_len(iters)) {
    # Random effect for trial, on what the forest has NOT explained.
    lm_d <- data.frame(y = y, trial_id = factor(dt$trial_id))
    re <- tryCatch(suppressMessages(lme4::lmer(y ~ 1 + (1 | trial_id), data = lm_d)),
                   error = function(err) NULL)
    if (is.null(re)) return(NULL)
    b <- lme4::ranef(re)$trial_id[as.character(dt$trial_id), 1]
    m <- fit_surrogate(Xt, dt$z - b, ntree = o_ntree, min_obs = 8)
    if (is.null(m)) return(NULL)
    fitted_rf <- rowMeans(vapply(m$trees, function(tr2)
      as.numeric(stats::predict(tr2, newdata = Xt)), numeric(nrow(Xt))))
    y <- dt$z - fitted_rf                     # residual for the next random-effect fit
  }
  function(i) predict_surrogate(m, F[i, , drop = FALSE])$mean
}
FITS <- list(A = fit_A, B = fit_B, C = fit_C)
LABEL <- c(A = "pooled  (config means, no trial)",
           B = "blocked (trial_id in forest, marginalised)",
           C = "merf    (trial as random effect outside)")

# --- repeated CV ------------------------------------------------------------
arms <- intersect(o_arms, names(FITS))
res  <- matrix(NA_real_, nrow = length(arms), ncol = o_reps, dimnames = list(arms, NULL))
cat("scheme ", o_scheme, " | rows ", nrow(d), " | configs ", length(hashes),
    " | trials ", length(trials), "\n", sep = "")
cat("running ", o_reps, " x ", o_folds, "-fold CV over configurations...\n\n", sep = "")
for (rep_i in seq_len(o_reps)) {
  set.seed(rep_i * 97)
  folds <- sample(rep_len(seq_len(o_folds), length(hashes)))
  preds <- lapply(arms, function(a) numeric(0)); names(preds) <- arms
  obs <- numeric(0)
  for (k in seq_len(o_folds)) {
    te <- hashes[folds == k]; tr <- !(d$config_hash %in% te)
    fns <- lapply(arms, function(a) FITS[[a]](tr)); names(fns) <- arms
    for (h in te) {
      i1 <- which(d$config_hash == h)[1]
      for (a in arms)
        preds[[a]] <- c(preds[[a]], if (is.null(fns[[a]])) NA_real_ else fns[[a]](i1))
      obs <- c(obs, truth$y[truth$config_hash == h])
    }
  }
  for (a in arms) res[a, rep_i] <- suppressWarnings(cor(preds[[a]], obs, use = "complete.obs"))
}

# --- report -----------------------------------------------------------------
cat(sprintf("%-45s %8s %7s %10s %8s\n", "arm", "mean cor", "sd", "vs A", "p"))
base <- res["A", ]
for (a in arms) {
  m <- mean(res[a, ], na.rm = TRUE); sdv <- stats::sd(res[a, ], na.rm = TRUE)
  if (a == "A") {
    cat(sprintf("%-45s %8.3f %7.3f %10s %8s\n", LABEL[a], m, sdv, "--", "--"))
  } else {
    dif <- res[a, ] - base
    p <- tryCatch(stats::t.test(dif)$p.value, error = function(err) NA_real_)
    cat(sprintf("%-45s %8.3f %7.3f %+10.3f %8.3f\n", LABEL[a], m, sdv, mean(dif, na.rm = TRUE), p))
  }
}
# The line that stops a null result being read as "no effect" when it means "not enough data".
sd_dif <- if (length(arms) > 1) {
  max(apply(res[setdiff(arms, "A"), , drop = FALSE], 1,
            function(r) stats::sd(r - base, na.rm = TRUE)))
} else {
  NA_real_
}
cat(sprintf("\nsmallest difference this store can resolve (2 x paired SE): %s\n",
            if (is.finite(sd_dif)) sprintf("~%.3f", 2 * sd_dif / sqrt(o_reps)) else "n/a"))
cat("a difference smaller than that is NOT evidence of no effect -- it is too little data.\n")
