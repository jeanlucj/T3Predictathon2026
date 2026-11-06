# dev/measure_churn.R
#
# Is the search converging, or is it churning? Replays the ranking machinery over prefixes of
# each run's evaluation stream and prints the three series that separate the two, plus the
# diagnostics that say WHY. Everything here is arithmetic on rows the store already holds --
# no network, no genotypes, nothing re-evaluated -- so a 70 MB store answers in minutes.
#
# This is the measurement behind docs/SEARCH_ALGORITHM.md sections 9-11. Read that first for
# what the numbers mean; this script exists so they can be re-derived rather than trusted.
#
# CONVERGENCE AND CHURN LOOK IDENTICAL IN ANY ONE OF THE THREE SERIES, which is why all three
# are printed:
#
#   TURNOVER    does the leader keep changing? A converged search has a stable top set. On its
#               own this proves nothing: an early search SHOULD churn.
#   DIRECTION   is the leader's score still rising? Turnover while the score climbs is the
#               search working. Turnover after the score goes flat is churn.
#   RESOLUTION  can the estimator tell the top configurations apart? `in_play` counts the
#               configurations whose optimistic bound still reaches the leader -- the
#               UNCAPPED contender set, before the top-8 cap. When that stays large while the
#               score is flat, the ranking is reshuffling noise.
#
# Churn is all three at once: turnover high, direction flat, resolution poor. The summary block
# after each run's tables states which pattern that run shows.
#
# Runs through the optimizer's OWN aggregate_scores(), .contenders(), fit_surrogate_merf() and
# candidate generators, so a change to the search changes these numbers too. Nothing here
# re-implements a rule.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript dev/measure_churn.R
#
# Options:
#   --store=<path>   store to read      (default: resolve_read_store(), as the tools do)
#   --runs=a,b       run_ids to analyse (default: every run with >= --min-evals rows)
#   --min-evals=500  skip runs shorter than this in the default run selection
#   --checkpoints=25 prefix points per run in the replay
#   --z=1            contender z; defaults to the live contender_z setting
#   --no-ei          skip section 4 (the only section that fits a surrogate -- the slow one)
#   --no-cause       skip section 3 (refits the variance components six times per run)
#   --seed=1         seed for section 4's candidate draws
#
# COST: one mixed-model fit per checkpoint per run, so it scales with checkpoints x rows, not
# with the whole store. Seconds for a 1,500-row run at 25 checkpoints, a couple of minutes for a
# 25,000-row one. --no-cause and --no-ei remove the two heaviest sections.
#
# SAFETY: reads a COPY of the store with its -wal/-shm sidecars, never the live file, so it is
# safe to run beside a live search. Read-only throughout.

# The optimizer ROOT, not this script's directory: `.Renviron` -- and so the T3 credentials --
# is read from the WORKING DIRECTORY only, with no parent walk. `settings.R` marks that root.
if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript dev/measure_churn.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
options(width = 150)                    # the resolution table is wide and unreadable wrapped
here::i_am("dev/measure_churn.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", n, "="), "", h[1]) }

s        <- optimizer_settings()
o_chk    <- max(3L, suppressWarnings(as.integer(opt("checkpoints", "25"))))
o_min    <- suppressWarnings(as.integer(opt("min-evals", "500")))
o_z      <- suppressWarnings(as.numeric(opt("z", as.character(s$contender_z %||% 1))))
o_seed   <- suppressWarnings(as.integer(opt("seed", "1")))
o_no_ei  <- "--no-ei" %in% args
o_no_cz  <- "--no-cause" %in% args

hr  <- function(t) cat("\n\n", t, "\n", strrep("=", nchar(t)), "\n", sep = "")
sec <- function(t) cat("\n-- ", t, "\n", sep = "")
pr  <- function(d, k = 4) print(as.data.frame(
  dplyr::mutate(d, dplyr::across(where(is.numeric), ~round(.x, k)))), row.names = FALSE)

# ---- load, from a copy ------------------------------------------------------
store_path <- resolve_read_store(opt("store"), "dev/measure_churn.R")
con  <- open_store(.copy_store_with_sidecars(
          store_path, file.path(tempdir(), "measure_churn.sqlite")))
evals_all <- read_evals(con)
runs_tbl  <- if ("runs" %in% DBI::dbListTables(con))
               tryCatch(tibble::as_tibble(DBI::dbReadTable(con, "runs")),
                        error = function(e) NULL) else NULL
close_store(con)

if (!("run_id" %in% names(evals_all)))
  stop("this store predates the run_id column, so its rows cannot be split by run")

cat("store : ", store_path, "  (", nrow(evals_all), " rows)\n", sep = "")
cat("z     : ", o_z, "   checkpoints: ", o_chk, "\n", sep = "")

# ---- which runs -------------------------------------------------------------
# Several runs can share one store -- it is a global archive, and a target-domain switch starts
# a new run_id against the same file. Mixing them would blend two different searches.
by_run <- evals_all |>
  dplyr::mutate(run_id = dplyr::coalesce(run_id, "<none>")) |>
  dplyr::count(run_id, name = "n_eval")
o_runs <- opt("runs")
want <- if (!is.null(o_runs)) {
  trimws(strsplit(o_runs, ",")[[1]])
} else {
  by_run$run_id[by_run$n_eval >= o_min]
}
want <- want[want %in% by_run$run_id]
if (!length(want)) stop("no run in this store matches --runs/--min-evals")

sec("runs in this store")
pr(by_run |>
   dplyr::mutate(analysed = run_id %in% want) |>
   dplyr::left_join(
     if (is.null(runs_tbl)) tibble::tibble(run_id = character(), scheme = character())
     else dplyr::select(runs_tbl, run_id, scheme), by = "run_id") |>
   dplyr::arrange(dplyr::desc(n_eval)))
if (!is.null(runs_tbl) && !all(want %in% runs_tbl$run_id))
  cat("note: ", sum(!want %in% runs_tbl$run_id), " of these runs have no `runs` row",
      " (they predate the table), so the declared universe is unknown --\n",
      "      every trial figure below is the REALISED set, the trials actually evaluated.\n",
      sep = "")

# ---- the pieces the sections share -----------------------------------------
# Variance components and the per-evaluation residual scale, off whatever aggregate_scores()
# actually fitted. NA throughout when it fell back to `pooled`: no fit, no components, and
# hence no `se` and no contenders at all.
vcomp <- function(agg) {
  v  <- attr(agg, "var_comps")
  wb <- attr(agg, "median_weight")
  g  <- function(nm) if (!is.null(v) && nm %in% names(v)) as.numeric(v[[nm]]) else NA_real_
  sdr <- g("sd_resid"); sdc <- g("sd_config")
  # lmer reports sd_resid at unit weight; the weights are n_test - 3, so this is the noise of
  # ONE evaluation at the median accession count -- the scale a config difference competes with.
  per_eval <- if (is.finite(sdr) && is.finite(wb) && wb > 0) sdr / sqrt(wb) else NA_real_
  list(sd_trial = g("sd_trial"), sd_config = sdc, sd_resid = sdr,
       resid_per_eval = per_eval,
       # Evaluations one configuration needs before the data outweigh the shrinkage prior --
       # i.e. before the BLUP is more measurement than assumption. Compare it to the number of
       # TRIALS available: n* above the universe size means the design cannot resolve configs
       # however long it runs.
       n_star = if (is.finite(per_eval) && is.finite(sdc) && sdc > 0)
                  (per_eval / sdc)^2 else NA_real_,
       estimator = attr(agg, "estimator") %||% NA_character_,
       note = attr(agg, "estimator_note") %||% NA_character_)
}

# Rows in the order the search saw them. `ts` is only second-resolution, so `id` breaks ties --
# without it a checkpoint can split a batch of simultaneous inserts arbitrarily.
run_rows <- function(run) {
  e <- evals_all |>
    dplyr::mutate(run_id = dplyr::coalesce(run_id, "<none>")) |>
    dplyr::filter(run_id == run)
  if ("id" %in% names(e)) dplyr::arrange(e, ts, id) else dplyr::arrange(e, ts)
}

# ============================================================================
# 1. the three churn series
# ============================================================================
series <- function(ev, run) {
  cuts <- unique(round(seq(min(200L, nrow(ev)), nrow(ev), length.out = o_chk)))
  out <- list(); prev_top <- NULL; prev_lead <- NA_character_
  for (cz in cuts) {
    agg <- aggregate_scores(ev[seq_len(cz), , drop = FALSE], s)
    sc  <- agg |> dplyr::filter(is.finite(mean_score)) |>
                  dplyr::arrange(dplyr::desc(mean_score))
    if (!nrow(sc)) next
    v    <- vcomp(agg)
    top  <- utils::head(sc$config_hash, 8L)
    # Uncapped vs capped: the gap IS the cap. `in_play` configurations have an optimistic bound
    # that still reaches the leader; only the best 8 of them can ever be served, and a
    # configuration pushed out is owed nothing (SEARCH_ALGORITHM.md section 7).
    in_play <- .contenders(agg, z = o_z, k = nrow(agg))
    capped  <- .contenders(agg, z = o_z, k = 8L)
    out[[length(out) + 1L]] <- tibble::tibble(
      n_eval = cz, n_cfg = nrow(sc),
      lead = substr(sc$config_hash[1], 1, 8), lead_score = sc$mean_score[1],
      lead_nok = sc$n_ok[1], lead_changed = !identical(sc$config_hash[1], prev_lead),
      median_score = stats::median(sc$mean_score),
      prev_lead_rank = if (is.na(prev_lead)) NA_integer_
                       else match(prev_lead, sc$config_hash),
      top8_jaccard = if (is.null(prev_top)) NA_real_
                     else length(intersect(top, prev_top)) / length(union(top, prev_top)),
      se_med = stats::median(agg$se, na.rm = TRUE),
      sd_trial = v$sd_trial, sd_config = v$sd_config,
      resid_per_eval = v$resid_per_eval, n_star = v$n_star,
      in_play = length(in_play), capped = length(capped),
      estimator = v$estimator)
    prev_top <- top; prev_lead <- sc$config_hash[1]
  }
  dplyr::bind_rows(out)
}

# ============================================================================
# 2. is one evaluation worth anything? se against replication
# ============================================================================
se_by_reps <- function(ev) {
  agg <- aggregate_scores(ev, s)
  v   <- vcomp(agg)
  # Two different ways this table cannot exist, and they mean opposite things: no fit at all
  # (no se, so .contenders() returns nothing and contender replication is off), versus a fit
  # whose sd_config came back at zero (every se is zero, so EVERY config is a contender).
  if (!identical(v$estimator, "blup"))
    return(paste0("no random-effects fit (estimator = ", v$estimator %||% "?",
                  "): no se, so .contenders() returns nothing at all"))
  if (!is.finite(v$sd_config) || v$sd_config <= 0)
    return(paste0("sd_config fitted to ", signif(v$sd_config, 3),
                  ": every se is zero, so se/sd_config is undefined and every ",
                  "configuration is a contender"))
  agg |>
    dplyr::filter(is.finite(se)) |>
    dplyr::mutate(bin = cut(n_ok, c(0, 1, 2, 3, 5, 8, 12, Inf),
                            labels = c("1", "2", "3", "4-5", "6-8", "9-12", "13+"))) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(configs = dplyr::n(), se_median = stats::median(se),
                     # se in units of sd_config: 1.0 means the estimate carries no more
                     # information than the prior, so the config is indistinguishable from
                     # the population it was drawn from.
                     se_over_sd_config = stats::median(se) / v$sd_config,
                     .groups = "drop") |>
    dplyr::filter(configs > 0)
}

# ============================================================================
# 3. why sd_config falls: dilution, noise, or a homogeneous population
# ============================================================================
# sd_config is the between-configuration sd -- the signal the search is trying to resolve, and
# simultaneously the shrinkage prior and the unit `se` is judged in. Three things make it fall,
# and they have different remedies, so they are separated by measurement rather than argued:
#
#   DILUTION      singly-evaluated configs shrink to the mean and flatten the estimate.
#                 Signature: sd_config RISES when the fit is restricted to replicated configs.
#   NOISE         the early estimate is unstable on few rows. Signature: only the first
#                 checkpoints move; section 1's series settles.
#   HOMOGENEITY   the configurations evaluated late genuinely do not differ, because the
#                 search stopped drawing anything new. Signature: fitting each HALF of the
#                 run alone gives a lower sd_config for the second half.
#
# Only the third is a defect, and only when the leader's score is not rising with it.
sd_cause <- function(ev) {
  agg0 <- aggregate_scores(ev, s)
  nok  <- dplyr::select(agg0, config_hash, n_ok)
  thr_tab <- purrr::map_dfr(c(1L, 2L, 3L, 4L, 6L), function(thr) {
    keep <- nok$config_hash[nok$n_ok >= thr]
    d <- dplyr::filter(ev, config_hash %in% keep)
    v <- if (nrow(d)) vcomp(aggregate_scores(d, s))
         else list(sd_config = NA_real_, resid_per_eval = NA_real_)
    tibble::tibble(min_n_ok = thr, configs = length(keep), rows = nrow(d),
                   sd_config = v$sd_config, resid_per_eval = v$resid_per_eval)
  })
  # Each half refitted alone, restricted to configurations with two rows IN that half, so the
  # two fits are comparably replicated and the only difference is WHICH configs the search was
  # drawing at the time.
  half <- floor(nrow(ev) / 2)
  half_tab <- purrr::map_dfr(c("first", "second"), function(lab) {
    idx <- if (lab == "first") seq_len(half) else (half + 1L):nrow(ev)
    d <- ev[idx, , drop = FALSE] |>
      dplyr::group_by(config_hash) |> dplyr::filter(dplyr::n() >= 2L) |> dplyr::ungroup()
    v <- if (nrow(d)) vcomp(aggregate_scores(d, s))
         else list(sd_config = NA_real_, resid_per_eval = NA_real_)
    tibble::tibble(half = lab, configs = dplyr::n_distinct(d$config_hash), rows = nrow(d),
                   sd_config = v$sd_config, resid_per_eval = v$resid_per_eval)
  })
  list(thr = thr_tab, half = half_tab)
}

# ============================================================================
# 4. does Expected Improvement still explore?
# ============================================================================
# Replays phase 4 of choose_config() on the run's final state. Two numbers decide it:
#
#   candidates at the sd floor   predict_surrogate() floors the between-tree sd at
#                                0.25 * sd(y_obs). When every candidate is at the floor they
#                                all carry the SAME uncertainty, and EI -- strictly increasing
#                                in the mean at fixed sd -- collapses to argmax(mean): pure
#                                greedy exploitation, with the exploration term inert.
#                                (predict_surrogate_marginal, the `blocked` path, applies no
#                                floor, so this failure is specific to merf and pooled.)
#   Spearman(EI, mean)           1.000 confirms the collapse directly.
#
# The candidate groups are built here rather than via propose_candidates() because that
# function returns one undifferentiated pool, and the question is WHICH group the pick came
# from. The generators, the dedup and the >= 5 novelty rule are the same, in the same order.
ei_probe <- function(ev, run) {
  set.seed(o_seed)
  agg    <- aggregate_scores(ev, s)
  scored <- dplyr::filter(agg, is.finite(mean_score))
  if (!nrow(scored)) return(NULL)
  rows <- dplyr::filter(ev, is.finite(score))
  if (nrow(rows)) {
    n_cfg <- tapply(rows$config_hash, rows$trial_id, function(x) length(unique(x)))
    rows  <- rows[n_cfg[as.character(rows$trial_id)] >= 2L, , drop = FALSE]
  }
  if (!nrow(rows) || dplyr::n_distinct(rows$trial_id) < 2L) return(NULL)
  feats <- configs_to_features(lapply(rows$config_json, config_from_json))
  surro <- fit_surrogate_merf(feats, rows$score, rows$trial_id,
                              ntree = s$ntree, min_obs = .surrogate_min_obs(s))
  if (is.null(surro)) return(NULL)

  # get_parents(), not get_elites(): the parent pool is the in-play set since 2026-11-04
  # (LESSONS #42), and this section's whole point is what the LIVE generator produces.
  elites <- get_parents(agg, s)
  if (length(elites) < 2L) return(NULL)
  grp <- c(rep("crossover", s$n_cross), rep("mutation", s$n_mut), rep("random", s$n_rand))
  cands <- c(
    lapply(seq_len(s$n_cross), function(i) { p <- sample(length(elites), 2)
      repair_config(crossover(elites[[p[1]]], elites[[p[2]]])) }),
    lapply(seq_len(s$n_mut), function(i)
      mutate_config(elites[[sample(length(elites), 1)]], n_subtasks = sample(1:2, 1))),
    lapply(seq_len(s$n_rand), function(i) sample_config()))
  h    <- vapply(cands, config_hash, character(1))
  dedup <- tibble::tibble(grp = grp, hash = h) |>
    dplyr::group_by(grp) |>
    dplyr::summarise(drawn = dplyr::n(), distinct = dplyr::n_distinct(hash), .groups = "drop")
  keep <- !duplicated(h); cands <- cands[keep]; grp <- grp[keep]; h <- h[keep]
  nov  <- !(h %in% unique(ev$config_hash))
  if (sum(nov) >= 5L) { cands <- cands[nov]; grp <- grp[nov] }

  cf   <- configs_to_features(cands)
  p    <- predict_surrogate(surro, cf)
  flr  <- 0.25 * stats::sd(surro$y_obs)
  best <- max(scored$mean_score)
  ei   <- expected_improvement(p$mean, p$sd, best, s$ei_xi)
  pick <- which.max(ei)
  list(dedup = dedup,
       by_grp = tibble::tibble(grp = grp, pred = p$mean, ei = ei) |>
         dplyr::group_by(grp) |>
         dplyr::summarise(candidates = dplyr::n(), pred_median = stats::median(pred),
                          pred_max = max(pred), ei_max = max(ei), .groups = "drop"),
       floor = flr, at_floor = sum(p$sd <= flr + 1e-12), n_cand = length(cands),
       greedy = identical(pick, which.max(p$mean)),
       spearman = suppressWarnings(stats::cor(ei, p$mean, method = "spearman")),
       best = best, pred_max = max(p$mean), pick_grp = grp[pick])
}

# ============================================================================
# 5. has the population collapsed to one lineage?
# ============================================================================
# A frozen parameter is one the search has stopped varying. Crossover swaps whole subtask
# BLOCKS, so once every elite shares a block it can no longer produce a new value for it, and
# only mutation or a fresh random draw can -- neither of which wins if EI has gone greedy.
# Counting distinct values per key over the run's distinct configurations shows how far that
# has gone.
lineage <- function(ev) {
  cfgs <- ev |> dplyr::distinct(config_hash, .keep_all = TRUE) |> dplyr::pull(config_json) |>
    lapply(function(j) tryCatch(config_from_json(j), error = function(e) NULL))
  cfgs <- cfgs[!vapply(cfgs, is.null, logical(1))]
  if (!length(cfgs)) return(NULL)
  keys <- canonical_keys()
  tibble::tibble(key = keys,
    distinct = vapply(keys, function(k) dplyr::n_distinct(vapply(cfgs, function(cf)
paste(cf[[k]] %||% NA), character(1))), integer(1)),
    modal_frac = vapply(keys, function(k) { v <- vapply(cfgs, function(cf)
paste(cf[[k]] %||% NA), character(1)); max(table(v)) / length(v) }, numeric(1))) |>
    dplyr::arrange(dplyr::desc(modal_frac), distinct) |>
    dplyr::mutate(n_configs = length(cfgs))
}

# ============================================================================
# per run
# ============================================================================
for (run in want) {
  ev <- run_rows(run)
  hr(paste0("RUN ", run, "  --  ", nrow(ev), " evaluations, ",
            sum(is.finite(ev$score)), " scored, ",
            dplyr::n_distinct(ev$trial_id), " trials realised"))

  d <- series(ev, run)
  if (!nrow(d)) { cat("no checkpoint produced a usable estimate -- nothing to measure\n"); next }

  sec("1a. TURNOVER -- does the top set keep changing?")
  pr(dplyr::select(d, n_eval, n_cfg, lead, lead_changed, top8_jaccard, prev_lead_rank), 3)
  sec("1b. DIRECTION -- is the leader's score still rising?")
  pr(dplyr::select(d, n_eval, lead_score, lead_nok, median_score) |>
     dplyr::mutate(gain_from_first = lead_score - d$lead_score[1]), 4)
  sec("1c. RESOLUTION -- can the estimator tell them apart?")
  pr(dplyr::select(d, n_eval, sd_trial, sd_config, resid_per_eval, se_med, n_star,
                   in_play, capped, estimator), 4)

  last <- d[nrow(d), ]
  cat("\nsummary\n")
  cat("  leader changed at ", sum(d$lead_changed, na.rm = TRUE), " of ", nrow(d),
      " checkpoints; ", dplyr::n_distinct(d$lead), " distinct leaders\n", sep = "")
  cat("  median top-8 Jaccard between consecutive checkpoints: ",
      round(stats::median(d$top8_jaccard, na.rm = TRUE), 3), "\n", sep = "")
  cat("  the previous checkpoint's leader, ranked at the next: median ",
      stats::median(d$prev_lead_rank, na.rm = TRUE), " (IQR ",
      paste(round(stats::quantile(d$prev_lead_rank, c(.25, .75), na.rm = TRUE)),
            collapse = "-"), ")\n", sep = "")
  cat("  leader's score: ", round(d$lead_score[1], 4), " -> ", round(last$lead_score, 4),
      ";  range over the last half of the checkpoints: ", round(
        max(utils::tail(d$lead_score, max(2, floor(nrow(d) / 2)))) -
        min(utils::tail(d$lead_score, max(2, floor(nrow(d) / 2)))), 4),
      ", against a median se of ", round(last$se_med, 4), "\n", sep = "")
  cat("  in play at the end: ", last$in_play, " of ", last$n_cfg,
      " configurations within z=", o_z, " of the leader, of which ", last$capped,
      " can be served\n", sep = "")
  cat("  n* = ", round(last$n_star, 1), " evaluations to outweigh the prior, against ",
      dplyr::n_distinct(ev$trial_id), " trials realised",
      if (is.finite(last$n_star) && last$n_star > dplyr::n_distinct(ev$trial_id))
        "  <-- EXCEEDS the trials available: se stays near sd_config at any budget"
      else "", "\n", sep = "")
  # The verdict, stated rather than left to the reader: all three series at once.
  flat <- {
    tailv <- utils::tail(d$lead_score, max(2, floor(nrow(d) / 2)))
    is.finite(last$se_med) && (max(tailv) - min(tailv)) < last$se_med }
  degenerate <- !is.finite(last$sd_config) || last$sd_config <= 0
  cat("  VERDICT: ",
      if (degenerate)
        paste0("DEGENERATE -- sd_config fitted to zero, so every se is zero, every one of the ",
               last$n_cfg, " configurations is in play, and the contender rule has no ",
               "yardstick left. Nothing this run evaluates can discriminate.")
      else if (stats::median(d$top8_jaccard, na.rm = TRUE) < 0.6 && flat)
        "CHURN -- the top set turns over while the score sits inside one se"
      else if (stats::median(d$top8_jaccard, na.rm = TRUE) < 0.6)
        "turnover, but the score is still moving: the search is working"
      else "stable top set", "\n", sep = "")

  sec("2. se against replication (final state)")
  t2 <- se_by_reps(ev)
  if (is.character(t2)) cat(t2, "\n", sep = "") else pr(t2, 5)

  if (!o_no_cz) {
    cz <- sd_cause(ev)
    sec("3a. sd_config restricted to replicated configs -- rises => dilution")
    pr(cz$thr, 5)
    sec("3b. sd_config in each half of the run alone -- falls => homogeneity")
    pr(cz$half, 5)
  }

  if (!o_no_ei) {
    sec("4. Expected Improvement at the final state")
    e4 <- ei_probe(ev, run)
    if (is.null(e4)) cat("no merf fit at the final state -- phase 4 cannot be replayed\n")
    else {
      cat("distinct candidates each generator produced (crossover collapses when the elites",
          " share a block):\n", sep = "")
      pr(e4$dedup)
      cat("\nsd floor = 0.25*sd(y_obs) = ", round(e4$floor, 4), "; at the floor: ",
          e4$at_floor, " of ", e4$n_cand, " candidates (",
          round(100 * e4$at_floor / e4$n_cand), "%)\n", sep = "")
      cat("argmax(EI) == argmax(mean): ", e4$greedy,
          "    Spearman(EI, mean) = ", round(e4$spearman, 4), "\n", sep = "")
      cat("EI reference `best` = ", round(e4$best, 4),
          " (a shrunken BLUP) vs max predicted mean = ", round(e4$pred_max, 4),
          " (a raw score) -- the two are on different scales\n", sep = "")
      cat("the pick came from: ", e4$pick_grp, "\n", sep = "")
      pr(e4$by_grp, 4)
    }
  }

  sec("5. frozen parameters -- modal_frac is the share of configurations carrying the commonest value")
  t5 <- lineage(ev)
  if (is.null(t5)) cat("no parseable configs\n") else {
    cat(t5$n_configs[1], " distinct configurations; ", sum(t5$modal_frac >= 0.9), " of ",
        nrow(t5), " parameters effectively frozen (one value in >= 90% of them), ",
        sum(t5$distinct <= 1), " never varied at all\n", sep = "")
    pr(utils::head(dplyr::select(t5, key, distinct, modal_frac), 14), 3)
  }
}

cat("\n\ndone. See docs/SEARCH_ALGORITHM.md sections 9-11 for what these mean.\n")
