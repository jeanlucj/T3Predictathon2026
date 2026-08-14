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
#   A  pooled    per-config mean scores, forest on config features only
#   B  blocked   one row per (config, trial), trial_id a factor, marginalised when scoring.
#                What choose_config() uses.
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
#   --n-grid=5       points on the learning curve (see below)
#   --curve-reps=8   CV repetitions per curve point
#   --no-curve       full-store comparison only
#   --out=logs       where the TSV and PNG go
#
# THE LEARNING CURVE. Each point re-runs the whole comparison on the FIRST n evaluations by
# timestamp, so the x-axis is the history of the run: "had I stopped at n, how good would each
# surrogate have been?" Writes <out>/surrogate_bakeoff.tsv (tidy, so the figure can be redrawn
# without re-running) and <out>/surrogate_bakeoff.png.

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
o_curve  <- !("--no-curve" %in% args)
o_ngrid  <- as.integer(opt("n-grid", "5"))
o_creps  <- as.integer(opt("curve-reps", "8"))
o_out    <- opt("out", "logs")

store_path <- resolve_read_store(what = "surrogate_bakeoff.R")
tmp <- .copy_store_with_sidecars(store_path, file.path(tempdir(), "bakeoff.sqlite"))
con <- DBI::dbConnect(RSQLite::SQLite(), tmp)
e   <- tibble::as_tibble(DBI::dbReadTable(con, "evals")); DBI::dbDisconnect(con)

if (is.null(o_scheme)) o_scheme <- names(sort(table(e$scheme), decreasing = TRUE))[1]
d <- e |>
  dplyr::filter(scheme == o_scheme, status == "ok", is.finite(score)) |>
  dplyr::mutate(z = atanh(pmin(pmax(score, -0.999), 0.999)), w = pmax(n_test - 3, 1))
if (nrow(d) < 30) stop("only ", nrow(d), " scored rows for scheme ", o_scheme,
                       " -- too few to cross-validate")

# Ordered by time: the learning curve subsamples the FIRST n rows, so the x-axis reads as the
# history of the run rather than as an abstract sample size.
d <- d[order(d$ts), ]
F_all <- configs_to_features(lapply(d$config_json, config_from_json))

# --- one cross-validated comparison over a subset of rows -------------------
# `rows` indexes d. Everything the three arms need is derived from the subset, and the fits
# are defined INSIDE so they close over it -- the full-store report and every point on the
# learning curve then run the identical code path.
cv_at <- function(rows, reps, folds, ntree) {
  d0 <- d[rows, , drop = FALSE]
  F0 <- F_all[rows, , drop = FALSE]
  trials <- sort(unique(d0$trial_id))
  hashes <- unique(d0$config_hash)
  # Too few configurations to cross-validate: a fit on a handful of rows still returns a
  # number, and that number is noise dressed as a measurement. Refuse instead.
  if (length(hashes) < 20 || nrow(d0) < 30)
    return(matrix(NA_real_, nrow = length(arms), ncol = reps, dimnames = list(arms, NULL)))
  truth <- d0 |> dplyr::group_by(config_hash) |>
    dplyr::summarise(y = sum(w * z) / sum(w), .groups = "drop")

  fit_A <- function(tr) {
    ag <- d0[tr, ] |> dplyr::group_by(config_hash) |>
      dplyr::summarise(y = sum(w * z) / sum(w), .groups = "drop")
    m <- fit_surrogate(F0[match(ag$config_hash, d0$config_hash), , drop = FALSE], ag$y,
                       ntree = ntree, min_obs = 8)
    if (is.null(m)) return(NULL)
    function(i) predict_surrogate(m, F0[i, , drop = FALSE])$mean
  }
  fit_B <- function(tr) {
    FB <- F0; FB$trial_id <- factor(d0$trial_id, levels = trials)
    m <- fit_surrogate(FB[tr, , drop = FALSE], d0$z[tr], ntree = ntree, min_obs = 8)
    if (is.null(m)) return(NULL)
    function(i) predict_surrogate_marginal(m, F0[i, , drop = FALSE], trials)$mean
  }
  fit_C <- function(tr, iters = 3L) {
    if (!requireNamespace("lme4", quietly = TRUE)) return(NULL)
    dt <- d0[tr, ]; Xt <- F0[tr, , drop = FALSE]
    if (dplyr::n_distinct(dt$trial_id) < 2) return(NULL)
    y <- dt$z; m <- NULL
    for (it in seq_len(iters)) {
      lm_d <- data.frame(y = y, trial_id = factor(dt$trial_id))
      re <- tryCatch(suppressMessages(lme4::lmer(y ~ 1 + (1 | trial_id), data = lm_d)),
                     error = function(err) NULL)
      if (is.null(re)) return(NULL)
      b <- lme4::ranef(re)$trial_id[as.character(dt$trial_id), 1]
      m <- fit_surrogate(Xt, dt$z - b, ntree = ntree, min_obs = 8)
      if (is.null(m)) return(NULL)
      fitted_rf <- rowMeans(vapply(m$trees, function(tr2)
        as.numeric(stats::predict(tr2, newdata = Xt)), numeric(nrow(Xt))))
      y <- dt$z - fitted_rf
    }
    function(i) predict_surrogate(m, F0[i, , drop = FALSE])$mean
  }
  FITS <- list(A = fit_A, B = fit_B, C = fit_C)

  out <- matrix(NA_real_, nrow = length(arms), ncol = reps, dimnames = list(arms, NULL))
  for (rep_i in seq_len(reps)) {
    set.seed(rep_i * 97)
    fold <- sample(rep_len(seq_len(folds), length(hashes)))
    preds <- lapply(arms, function(a) numeric(0)); names(preds) <- arms
    obs <- numeric(0)
    for (k in seq_len(folds)) {
      te <- hashes[fold == k]; tr <- !(d0$config_hash %in% te)
      fns <- lapply(arms, function(a) FITS[[a]](tr)); names(fns) <- arms
      for (h in te) {
        i1 <- which(d0$config_hash == h)[1]
        for (a in arms)
          preds[[a]] <- c(preds[[a]], if (is.null(fns[[a]])) NA_real_ else fns[[a]](i1))
        obs <- c(obs, truth$y[truth$config_hash == h])
      }
    }
    for (a in arms) out[a, rep_i] <- suppressWarnings(cor(preds[[a]], obs, use = "complete.obs"))
  }
  out
}

LABEL <- c(A = "pooled  (config means, no trial)",
           B = "blocked (trial_id in forest, marginalised)",
           C = "merf    (trial as random effect outside)")
arms <- intersect(o_arms, c("A", "B", "C"))

cat("scheme ", o_scheme, " | rows ", nrow(d), " | configs ",
    dplyr::n_distinct(d$config_hash), " | trials ", dplyr::n_distinct(d$trial_id), "\n", sep = "")
cat("running ", o_reps, " x ", o_folds, "-fold CV over configurations...\n\n", sep = "")
res <- cv_at(seq_len(nrow(d)), o_reps, o_folds, o_ntree)

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

# --- learning curve ---------------------------------------------------------
if (o_curve) {
  grid <- unique(round(nrow(d) * c(0.2, 0.35, 0.5, 0.7, 1)))
  if (is.finite(o_ngrid) && o_ngrid >= 2)
    grid <- unique(round(seq(0.2, 1, length.out = o_ngrid) * nrow(d)))
  cat("\n\n=== learning curve: ", o_creps, " x ", o_folds, "-fold CV at ",
      length(grid), " sizes ===\n", sep = "")
  # The confound, printed where it cannot be skipped. Early rows are seeds and random-init
  # draws; later rows are acquisition picks concentrated in one region. So a flat or falling
  # segment can mean the SEARCH NARROWED, not that the surrogate stopped learning.
  cat("NOTE: rows are in time order, so this mixes 'more data' with 'differently distributed\n",
      "      data' -- early evaluations are seeds/random, later ones are acquisition picks.\n\n", sep = "")

  long <- NULL
  SHORT <- c(A = "pooled", B = "blocked", C = "merf")
  cat(sprintf("%8s %8s %s\n", "n", "configs",
              paste(sprintf("%22s", SHORT[arms]), collapse = "")))
  for (n in grid) {
    r  <- cv_at(seq_len(n), o_creps, o_folds, o_ntree)
    nc <- dplyr::n_distinct(d$config_hash[seq_len(n)])
    cat(sprintf("%8d %8d %s\n", n, nc,
        paste(vapply(arms, function(a) {
          m <- mean(r[a, ], na.rm = TRUE)
          if (!is.finite(m)) sprintf("%22s", "-- too few configs") else
            sprintf("%16.3f +-%4.3f", m, stats::sd(r[a, ], na.rm = TRUE))
        }, character(1)), collapse = "")))
    utils::flush.console()
    long <- rbind(long, data.frame(n = n, configs = nc,
                                   arm = rep(arms, each = ncol(r)),
                                   rep = rep(seq_len(ncol(r)), times = length(arms)),
                                   cor = as.vector(t(r))))
  }

  dir.create(o_out, showWarnings = FALSE, recursive = TRUE)
  tsv <- file.path(o_out, "surrogate_bakeoff.tsv")
  readr::write_tsv(long, tsv)
  cat("\nwrote ", tsv, "\n", sep = "")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    sm <- long |>
      dplyr::filter(is.finite(cor)) |>
      dplyr::group_by(n, arm) |>
      dplyr::summarise(m = mean(cor), se = stats::sd(cor) / sqrt(dplyr::n()), .groups = "drop") |>
      dplyr::mutate(arm = factor(arm, levels = arms,
                                 labels = trimws(sub("\\s+\\(.*", "", LABEL[arms]))))
    p <- ggplot2::ggplot(sm, ggplot2::aes(n, m, colour = arm, fill = arm)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = m - se, ymax = m + se),
                           alpha = 0.15, colour = NA) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 1.6) +
      ggplot2::labs(x = "evaluations", y = "cor(predicted, observed)",
                    colour = NULL, fill = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "bottom")
    png <- file.path(o_out, "surrogate_bakeoff.png")
    ggplot2::ggsave(png, p, width = 6.5, height = 4.2, dpi = 150)
    cat("wrote ", png, "\n", sep = "")
  } else {
    cat("(ggplot2 not installed -- TSV written, no figure)\n")
  }

  # The ceiling belongs in the text, not on the figure: it moves with the store's variance
  # components, and a mis-drawn reference line would be worse than none.
  cat(sprintf("\nFor scale: with one trial per configuration the attainable cor is ~0.21\n"))
  cat("(sd_config 0.036 against sd_trial 0.078 + sd_resid 0.151); it rises with replication.\n")
}
