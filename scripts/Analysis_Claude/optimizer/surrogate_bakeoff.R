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
#   --curve-min=100  smallest learning-curve point, in ROWS (see below)
#   --curve-ratios=10,15,25,40,60
#                    the ladder's rungs within each decade; must all be >= 10 and < 100
#   --curve-reps=8   CV repetitions per curve point
#   --no-curve       full-slice comparison only
#   --all-trials     score every trial in the store, not just the run's pinned universe
#   --all-builds     score rows from every build, not just the current one
#   --trials=a,b     name the universe explicitly, for a store written before runs were
#                    recorded (no run row means no universe to filter to)
#   --exclude-trials=c  drop named trials from whatever universe applies
#   --out=logs       where the TSVs and PNG go
#
# THE SLICE. The store is a shared archive across target domains and builds, so by default only
# the rows the surrogate will actually be fitted to are scored: the run's pinned trial universe,
# at the current build. Rows outside it can be the majority, and since the arms differ only in
# how they model the TRIAL effect, an unbalanced trial factor decides the answer before the arms
# do. The header prints the largest trial's share of rows for that reason.
#
# THE LEARNING CURVE. Each point re-runs the whole comparison on the FIRST n evaluations by
# timestamp, so the x-axis is the history of the run: "had I stopped at n, how good would each
# surrogate have been?" Points are ABSOLUTE row counts on a log ladder -- the rungs of
# --curve-ratios through the powers of ten, from --curve-min up to the size of the slice, which
# is always kept as the last point. Absolute rather than fractional so two runs' curves can be
# laid side by side. Writes <out>/surrogate_bakeoff_main.tsv (per-rep values for the
# full-slice comparison, always) and, with the curve, <out>/surrogate_bakeoff.tsv (tidy, so the
# figure can be redrawn without re-running) plus <out>/surrogate_bakeoff.png.

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
o_cmin   <- as.integer(opt("curve-min", "100"))
o_crat   <- as.numeric(trimws(strsplit(opt("curve-ratios", "10,15,25,40,60"), ",")[[1]]))
# --n-grid meant "how many points"; --curve-min means "the smallest point, in rows". Reading
# one as the other would silently produce a curve nobody asked for.
if (any(grepl("^--n-grid=", args)))
  stop("--n-grid is gone: the curve is now a ladder from a minimum row count.\n",
       "  use --curve-min=<rows> (default 100), and --curve-ratios to reshape the ladder.")
if (!length(o_crat) || any(!is.finite(o_crat)) || any(o_crat < 10) || any(o_crat >= 100))
  stop("--curve-ratios must span ONE decade -- every value >= 10 and < 100, e.g. 10,15,25,40,60.\n",
       "  the ladder is built by multiplying them through the powers of ten.")
o_creps  <- as.integer(opt("curve-reps", "8"))
o_out    <- opt("out", "logs")
o_alltr  <- "--all-trials" %in% args
o_allbld <- "--all-builds" %in% args
o_utr    <- { v <- opt("trials");         if (is.null(v)) NULL else trimws(strsplit(v, ",")[[1]]) }
o_xtr    <- { v <- opt("exclude-trials"); if (is.null(v)) NULL else trimws(strsplit(v, ",")[[1]]) }

s_set <- optimizer_settings()
BUILD <- s_set$build %||% OPTIMIZER_BUILD
store_path <- resolve_read_store(what = "surrogate_bakeoff.R")
tmp <- .copy_store_with_sidecars(store_path, file.path(tempdir(), "bakeoff.sqlite"))
con <- open_store(tmp)
e   <- read_evals(con)

if (is.null(o_scheme)) o_scheme <- names(sort(table(e$scheme), decreasing = TRUE))[1]

# The store is a shared archive across domains and builds; a surrogate comparison is only
# meaningful on the slice the surrogate will actually be fitted to. Rows outside it can be the
# majority, and then they, not the arms, decide the answer.
d <- e |> dplyr::filter(scheme == o_scheme, status == "ok", is.finite(score))
n_scheme <- nrow(d)

# The universe comes from the run these ROWS belong to. Re-deriving a run id from the current
# settings would silently return "unpinned" whenever the bakeoff is run under settings that
# differ from the run's in any of run_id_for()'s keys, and score every trial without saying so.
universe <- local({
  rid <- unique(stats::na.omit(d$run_id %||% NA_character_))
  if (!length(rid)) rid <- run_id_for(s_set, BUILD)
  u <- unique(unlist(lapply(rid, function(r) run_universe(read_run(con, r)))))
  if (length(u)) return(as.character(u))
  rows <- tibble::as_tibble(DBI::dbReadTable(con, "runs"))
  if (nrow(rows) == 1L) run_universe(rows) else NULL     # one run in the store is unambiguous
})
if (!is.null(o_utr)) universe <- o_utr
if (!is.null(o_xtr))
  universe <- setdiff(universe %||% unique(as.character(d$trial_id)), o_xtr)
close_store(con)
if (!o_alltr)  d <- filter_evals_to_universe(d, universe)
if (!o_allbld) d <- filter_evals_to_build(d, BUILD)
cat(sprintf("slice: %d scored %s row(s) -> %d after %s\n", n_scheme, o_scheme, nrow(d),
            paste(c(if (!o_alltr) sprintf("universe (%s)",
                                          if (length(universe))
                                            sprintf("%d trials%s", length(universe),
                                                    if (!is.null(o_utr) || !is.null(o_xtr)) ", given on the command line" else "")
                                          else "unpinned -- no constraint"),
                    if (!o_allbld) paste("build", BUILD)) %||% "no filtering",
                  collapse = " + ")))
if (!o_alltr && !length(universe))
  cat("  no run row in this store, so every trial is scored: name the universe with --trials=,\n",
      "  or drop what does not belong with --exclude-trials=.\n", sep = "")
d <- d |> dplyr::mutate(z = atanh(pmin(pmax(score, -0.999), 0.999)), w = pmax(n_test - 3, 1))
if (nrow(d) < 30) stop("only ", nrow(d), " scored rows for scheme ", o_scheme,
                       " in this slice -- too few to cross-validate (--all-trials widens it)")

# Ordered by time: the learning curve subsamples the FIRST n rows, so the x-axis reads as the
# history of the run rather than as an abstract sample size.
d <- d[order(d$ts), ]
F_all <- configs_to_features(lapply(d$config_json, config_from_json))

# --- one cross-validated comparison over a subset of rows -------------------
# `rows` indexes d. Everything the three arms need is derived from the subset, and the fits
# are defined INSIDE so they close over it -- the full-store report and every point on the
# learning curve then run the identical code path.
# `label` names this pass in the progress lines; NULL silences them. Progress is line-oriented
# cat + flush.console(), never \r, because this is normally run under nohup into a log file --
# the same idiom as prewarm_indices.R.
cv_at <- function(rows, reps, folds, ntree, label = NULL) {
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
  # NOTE on all three closures: `i` is a VECTOR of row indices and the result is a vector.
  # Both predict functions already take multiple rows -- predict_surrogate vapply()s over the
  # trees once for the whole newdata frame, and predict_surrogate_marginal is written for n_c
  # configs. Calling them per config paid the 200-tree loop once per config instead of once per
  # fold: measured at 5,210 rows / 1,259 configs, 250 configs took 25.1 s one at a time against
  # 0.15 s batched, for bit-identical values.
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
  n_fold <- reps * folds                      # the unit of progress
  done <- 0L
  t_start <- Sys.time()
  for (rep_i in seq_len(reps)) {
    set.seed(rep_i * 97)
    fold <- sample(rep_len(seq_len(folds), length(hashes)))
    preds <- lapply(arms, function(a) numeric(0)); names(preds) <- arms
    obs <- numeric(0)
    for (k in seq_len(folds)) {
      te <- hashes[fold == k]; tr <- !(d0$config_hash %in% te)
      fns <- lapply(arms, function(a) FITS[[a]](tr)); names(fns) <- arms
      # One call per arm for the whole fold. match() also replaces a which(==)[1] scan per
      # config, which was O(rows x configs) on its own.
      i1 <- match(te, d0$config_hash)
      for (a in arms)
        preds[[a]] <- c(preds[[a]],
                        if (is.null(fns[[a]])) rep(NA_real_, length(i1)) else fns[[a]](i1))
      obs <- c(obs, truth$y[match(te, truth$config_hash)])

      # After every fold: an ETA from the folds actually done, so the first line lands within
      # one fold of starting rather than after the whole pass. The first fold is also the
      # estimate -- there is no separate calibration run.
      done <- done + 1L
      if (!is.null(label)) {
        el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
        cat(sprintf("  [%s] fold %3d/%-3d %3.0f%%  %.1f s/fold  elapsed %.1f min  ETA %.1f min\n",
                    label, done, n_fold, 100 * done / n_fold, el / done, el / 60,
                    el / done * (n_fold - done) / 60))
        utils::flush.console()
      }
    }
    for (a in arms) out[a, rep_i] <- suppressWarnings(cor(preds[[a]], obs, use = "complete.obs"))
  }
  out
}

LABEL <- c(A = "pooled  (config means, no trial)",
           B = "blocked (trial_id in forest, marginalised)",
           C = "merf    (trial as random effect outside)")
SHORT_ARM <- c(A = "pooled", B = "blocked", C = "merf")
arms <- intersect(o_arms, c("A", "B", "C"))

cat("scheme ", o_scheme, " | rows ", nrow(d), " | configs ",
    dplyr::n_distinct(d$config_hash), " | trials ", dplyr::n_distinct(d$trial_id), "\n", sep = "")
# The comparison is between ways of modelling the TRIAL effect, so the design of the trial
# factor decides what it can measure. One trial holding most of the rows makes the answer a
# statement about that trial.
local({
  per <- d |> dplyr::group_by(trial_id) |>
    dplyr::summarise(n = dplyr::n(), n_cfg = dplyr::n_distinct(config_hash), .groups = "drop")
  top <- per[which.max(per$n), ]
  cat(sprintf("design: largest trial %s holds %.0f%% of rows | configs per trial median %d (range %d-%d)%s\n",
              top$trial_id, 100 * top$n / nrow(d), stats::median(per$n_cfg),
              min(per$n_cfg), max(per$n_cfg),
              if (top$n / nrow(d) > 0.25) "   <-- one trial dominates; arms are not comparable on this slice" else ""))
})
cat("running ", o_reps, " x ", o_folds, "-fold CV over configurations...\n\n", sep = "")
res <- cv_at(seq_len(nrow(d)), o_reps, o_folds, o_ntree, label = "main")
cat("\n")

# --- report -----------------------------------------------------------------
cat(sprintf("%-45s %8s %7s\n", "arm", "mean cor", "sd"))
for (a in arms)
  cat(sprintf("%-45s %8.3f %7.3f\n", LABEL[a], mean(res[a, ], na.rm = TRUE),
              stats::sd(res[a, ], na.rm = TRUE)))

# Every pair, because which comparison matters depends on what is deployed: A vs B says
# whether modelling the trial pays at all, B vs C which way of modelling it is better. Each
# gets its OWN resolution figure -- reps are paired, so the SE of a difference has nothing to
# do with the spread of either arm on its own.
if (length(arms) > 1) {
  cat(sprintf("\n%-14s %10s %8s %8s %12s\n", "contrast", "diff", "SE", "p", "resolvable"))
  for (pair in utils::combn(arms, 2, simplify = FALSE)) {
    a1 <- pair[1]; a2 <- pair[2]
    dif <- res[a2, ] - res[a1, ]
    ok  <- sum(is.finite(dif))
    se  <- stats::sd(dif, na.rm = TRUE) / sqrt(max(1L, ok))
    p   <- tryCatch(stats::t.test(dif)$p.value, error = function(err) NA_real_)
    cat(sprintf("%-14s %+10.3f %8.3f %8.3f %12s\n",
                paste0(SHORT_ARM[a2], "-", SHORT_ARM[a1]),
                mean(dif, na.rm = TRUE), se, p,
                if (is.finite(se)) sprintf("~%.3f", 2 * se) else "n/a"))
  }
  cat("\na difference smaller than its `resolvable` column is NOT evidence of no effect --\n")
  cat("it is too little data. Reps are paired, so read the contrast rows, not the sds above.\n")
}

# The per-rep values for the MAIN comparison, so a contrast can be recomputed, or a further
# test run, without paying for the fits again.
dir.create(o_out, showWarnings = FALSE, recursive = TRUE)
main_tsv <- file.path(o_out, "surrogate_bakeoff_main.tsv")
readr::write_tsv(data.frame(arm = rep(arms, each = ncol(res)),
                            rep = rep(seq_len(ncol(res)), times = length(arms)),
                            cor = as.vector(t(res[arms, , drop = FALSE])),
                            n = nrow(d), scheme = o_scheme,
                            trials = dplyr::n_distinct(d$trial_id)), main_tsv)
cat("\nwrote ", main_tsv, "\n", sep = "")

# --- learning curve ---------------------------------------------------------
# Sizes to score at: the ratios through the powers of ten, from n_min up to n_max, which is
# always the last point. The 1.2 margins drop a rung that would land within a fifth of either
# end -- next to n_min it is a duplicate, next to n_max it is a duplicate of the most expensive
# fit in the run.
.curve_grid <- function(n_max, n_min = 100L, ratios = c(10, 15, 25, 40, 60)) {
  n_max <- as.integer(n_max); n_min <- max(1L, as.integer(n_min))
  if (n_max <= n_min) return(n_max)
  rungs <- sort(unique(as.vector(outer(ratios, 10^(0:ceiling(log10(n_max)))))))
  rungs <- rungs[rungs > n_min * 1.2 & rungs < n_max / 1.2]
  as.integer(sort(unique(c(n_min, rungs, n_max))))
}

if (o_curve) {
  grid <- .curve_grid(nrow(d), o_cmin, o_crat)
  cat("\n\n=== learning curve: ", o_creps, " x ", o_folds, "-fold CV at ",
      length(grid), " sizes ===\n", sep = "")
  # The confound, printed where it cannot be skipped. Early rows are seeds and random-init
  # draws; later rows are acquisition picks concentrated in one region. So a flat or falling
  # segment can mean the SEARCH NARROWED, not that the surrogate stopped learning.
  cat("NOTE: rows are in time order, so this mixes 'more data' with 'differently distributed\n",
      "      data' -- early evaluations are seeds/random, later ones are acquisition picks.\n",
      "      The +- is the spread across reps, which re-split the SAME rows: it is split noise\n",
      "      at fixed n and carries none of the sampling error, so the small-n points are less\n",
      "      certain than their bars suggest.\n\n", sep = "")

  long <- NULL
  cat(sprintf("%8s %8s %s\n", "n", "configs",
              paste(sprintf("%22s", SHORT_ARM[arms]), collapse = "")))
  for (n in grid) {
    r  <- cv_at(seq_len(n), o_creps, o_folds, o_ntree, label = sprintf("curve n=%d", n))
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
