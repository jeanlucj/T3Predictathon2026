# dev/sweep_contender_z.R
#
# How large is the "could still be best" set at each contender_z, AT EVERY STAGE of a run?
#
# `contender_z` decides two things today and would decide a third under the unification in
# SEARCH_ALGORITHM.md section 7:
#
#   IN PLAY     ucb_i >= max_j(mean_j) under the POINT bar: the set that cannot be ruled out.
#               Read off .in_play(rule = "point"), so the interval arithmetic is the production
#               one (formed on the z scale, endpoints transformed -- SEARCH_ALGORITHM.md s7).
#   SERVED      the best 8 of that set by UCB -- the only ones that earn extra evaluations.
#   PARENTS     the first n_parents of it are also the parents of crossover and mutation. This
#               script measured the asymmetry that motivated that unification, when parents were
#               the top 8 by mean_score with no uncertainty consulted; `overlap` below is how
#               many parents the change moves. Adopted 2026-11-04 (LESSONS #42).
#
# It measures the POINT bar throughout, which is no longer the production rule. For the racing
# bar, the resolution floor and the reference-population sd_config, see sweep_contender_rule.R.
#
# The question a single final-state number cannot answer: is a FIXED z usable, or does the set
# size drift so much over a run that z would have to adapt? `se` shrinks as the population
# homogenizes (LESSONS #33, section 11), and z is denominated in `se`, so the yardstick moves
# under the search. This replays the real estimator over prefixes of each run's evaluation
# stream and reports the set at every z at every stage.
#
# Three things per run per z:
#   in_play       how many configurations are in play
#   overlap       how many of the UCB top 8 are ALSO in the mean-ranked top 8 -- i.e. how much
#                 the proposed parent pool would differ from today's
#   med_n_ok      median replication of the UCB top 8: does UCB ranking promote singletons?
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript dev/sweep_contender_z.R
#
# Options:
#   --store=<path>    store to read      (default: resolve_read_store(), as the tools do)
#   --runs=a,b        run_ids to analyse (default: every run with >= --min-evals rows)
#   --min-evals=500   skip runs shorter than this in the default run selection
#   --checkpoints=25  prefix points per run in the replay
#   --zs=0,0.25,...   the z grid
#   --band=0.10       "usable" upper bound for in_play, as a fraction of the scored population
#
# COST: one mixed-model fit per checkpoint per run -- the z grid is free, since every z reads
# the same fit. Same order as dev/measure_churn.R section 1.
#
# SAFETY: reads a COPY of the store with its -wal/-shm sidecars, never the live file. Read-only.

if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript dev/sweep_contender_z.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
options(width = 170)
here::i_am("dev/sweep_contender_z.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(n, d = NULL) {
  h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", n, "="), "", h[1])
}

s      <- optimizer_settings()
o_min  <- suppressWarnings(as.integer(opt("min-evals", "500")))
o_chk  <- max(3L, suppressWarnings(as.integer(opt("checkpoints", "25"))))
o_band <- suppressWarnings(as.numeric(opt("band", "0.10")))
o_zs   <- suppressWarnings(as.numeric(strsplit(opt("zs", "0,0.1,0.25,0.5,0.75,1,1.5,2"),
                                               ",")[[1]]))
o_zs   <- sort(unique(o_zs[is.finite(o_zs) & o_zs >= 0]))
CAP    <- 8L                      # the hardcoded cap in .contenders()' call sites

hr <- function(t) cat("\n\n", t, "\n", strrep("=", nchar(t)), "\n", sep = "")
pr <- function(d, k = 4) print(as.data.frame(
  dplyr::mutate(d, dplyr::across(where(is.numeric), ~round(.x, k)))), row.names = FALSE)

store_path <- resolve_read_store(opt("store"), "dev/sweep_contender_z.R")
con  <- open_store(.copy_store_with_sidecars(
          store_path, file.path(tempdir(), "sweep_contender_z.sqlite")))
evals_all <- read_evals(con)
close_store(con)
if (!("run_id" %in% names(evals_all)))
  stop("this store predates the run_id column, so its rows cannot be split by run")

cat("store      : ", store_path, "  (", nrow(evals_all), " rows)\n", sep = "")
cat("live z     : ", s$contender_z %||% 1, "    cap: ", CAP, "\n", sep = "")
cat("z grid     : ", paste(o_zs, collapse = ", "), "\n", sep = "")
cat("checkpoints: ", o_chk, "    usable band: ", CAP, " <= in_play <= ",
    o_band, " * scored\n", sep = "")

by_run <- evals_all |>
  dplyr::mutate(run_id = dplyr::coalesce(run_id, "<none>")) |>
  dplyr::count(run_id, name = "n_eval")
want <- {
  if (!is.null(opt("runs"))) strsplit(opt("runs"), ",")[[1]]
  else by_run$run_id[by_run$n_eval >= o_min]
}
want <- want[want %in% by_run$run_id]
if (!length(want)) stop("no run has >= ", o_min, " rows; pass --runs= or lower --min-evals=")

# Rows in the order the search saw them. `ts` is second-resolution, so `id` breaks ties.
run_rows <- function(run) {
  e <- evals_all |>
    dplyr::mutate(run_id = dplyr::coalesce(run_id, "<none>")) |>
    dplyr::filter(run_id == run)
  if ("id" %in% names(e)) dplyr::arrange(e, ts, id) else dplyr::arrange(e, ts)
}

# One fit per checkpoint; every z reads that same fit.
sweep_run <- function(ev, run) {
  cuts <- unique(round(seq(min(200L, nrow(ev)), nrow(ev), length.out = o_chk)))
  out <- list()
  for (cz in cuts) {
    agg <- aggregate_scores(ev[seq_len(cz), , drop = FALSE], s)
    sc  <- agg |> dplyr::filter(is.finite(mean_score)) |>
                  dplyr::arrange(dplyr::desc(mean_score))
    if (!nrow(sc)) next
    vc <- attr(agg, "var_comps")
    # today's parent pool, for the overlap column
    el <- utils::head(sc$config_hash, CAP)
    # The production rule compares an OPTIMISTIC bound against a POINT estimate of the maximum
    # (SEARCH_ALGORITHM.md section 7). max() over noisy estimates is the most upward-biased
    # order statistic in the set, so a lucky under-replicated config raises the bar for
    # everyone. `ip_lcb` is the racing alternative -- optimistic vs PESSIMISTIC,
    # ucb_i >= max_j(mean_j - z*se_j) -- whose bar is discounted by the leader's own
    # uncertainty. Measured here, not adopted; the question is which one is stable.
    a2 <- agg[is.finite(agg$mean_score) & is.finite(agg$se %||% NA_real_), , drop = FALSE]
    for (z in o_zs) {
      # rule = "point" explicitly: .contenders() defaults to the racing bar since 2026-11-04,
      # and this script's whole point is the comparison, so neither column may inherit a default.
      ip   <- .contenders(agg, z = z, k = nrow(agg), rule = "point")
      top  <- .contenders(agg, z = z, k = CAP, rule = "point")
      nok  <- sc$n_ok[match(top, sc$config_hash)]
      n_lcb <- if (!nrow(a2)) NA_integer_ else
                 length(.in_play(agg, z = z, rule = "racing"))
      out[[length(out) + 1L]] <- tibble::tibble(
        run = run, n_eval = cz, n_cfg = nrow(sc),
        sd_config = if (is.null(vc)) NA_real_ else vc[["sd_config"]],
        lead_nok = sc$n_ok[1], lead_gap = sc$mean_score[1] - stats::median(sc$mean_score),
        z = z, in_play = length(ip), pct = 100 * length(ip) / nrow(sc),
        in_play_lcb = n_lcb,
        overlap = length(intersect(top, el)),
        med_n_ok = suppressWarnings(stats::median(nok, na.rm = TRUE)),
        n_single = sum(nok == 1L, na.rm = TRUE))
    }
  }
  dplyr::bind_rows(out)
}

all <- list()
for (run in want) {
  ev <- run_rows(run)
  hr(paste0("RUN ", run, "   ", nrow(ev), " evaluations"))
  d <- sweep_run(ev, run)
  if (!nrow(d)) { cat("no checkpoint produced a usable estimate\n"); next }
  all[[run]] <- d
  deg <- d |> dplyr::filter(z == min(o_zs)) |>
              dplyr::summarise(n = dplyr::n(), z0 = sum(in_play == n_cfg))
  if (deg$z0 > 0)
    cat("NOTE: sd_config fits to 0 at ", deg$z0, " of ", deg$n, " checkpoints. There every BLUP ",
        "ties,\n      every se is 0, and the whole population is in play AT EVERY z -- the ",
        "z column\n      is inoperative on those rows (SEARCH_ALGORITHM.md section 4.6).\n",
        sep = "")

  cat("\n-- in_play, by checkpoint (rows) and z (columns)\n")
  pr(d |> dplyr::select(n_eval, n_cfg, sd_config, z, in_play) |>
       tidyr::pivot_wider(names_from = z, values_from = in_play, names_prefix = "z="))
  cat("\n-- overlap of the UCB top 8 with the mean-ranked top 8 (out of 8)\n")
  pr(d |> dplyr::select(n_eval, z, overlap) |>
       tidyr::pivot_wider(names_from = z, values_from = overlap, names_prefix = "z="))
  cat("\n-- median n_ok of the UCB top 8\n")
  pr(d |> dplyr::select(n_eval, z, med_n_ok) |>
       tidyr::pivot_wider(names_from = z, values_from = med_n_ok, names_prefix = "z="))
  cat("\n-- in_play under the RACING rule (ucb_i >= max of the LOWER bounds), and the leader\n")
  pr(d |> dplyr::select(n_eval, lead_nok, lead_gap, z, in_play_lcb) |>
       tidyr::pivot_wider(names_from = z, values_from = in_play_lcb, names_prefix = "z="))
}

if (!length(all)) stop("nothing measured")
A <- dplyr::bind_rows(all)

hr("Is a FIXED z usable? in_play across the stages of each run")
cat("usable = the set is bigger than the cap (so the cap has something to ration) and still\n",
    "means something (<= ", o_band, " of the scored population). drift = max/min over the\n",
    "checkpoints: 1 would be a set size that does not move at all.\n", sep = "")
pr(A |> dplyr::group_by(run, z) |>
     dplyr::summarise(chk = dplyr::n(),
       ip_min = min(in_play), ip_med = stats::median(in_play), ip_max = max(in_play),
       drift = max(in_play) / pmax(min(in_play), 1),
       pct_med = stats::median(pct),
       usable = 100 * mean(in_play >= CAP & in_play <= o_band * n_cfg),
       overlap_med = stats::median(overlap),
       .groups = "drop") |>
     dplyr::arrange(run, z))

hr("The same, pooled over the runs where z does anything (sd_config > 0)")
pr(A |> dplyr::filter(is.finite(sd_config), sd_config > 0) |>
     dplyr::group_by(z) |>
     dplyr::summarise(runs = dplyr::n_distinct(run), chk = dplyr::n(),
       ip_med = stats::median(in_play), pct_med = stats::median(pct),
       usable = 100 * mean(in_play >= CAP & in_play <= o_band * n_cfg),
       overlap_med = stats::median(overlap),
       med_n_ok = stats::median(med_n_ok, na.rm = TRUE),
       singles = stats::median(n_single, na.rm = TRUE),
       .groups = "drop"))

hr("Which RULE is stable? step-to-step instability of the in-play set, same fits")
cat("jump = median |log2(in_play[t+1] / in_play[t])| over consecutive checkpoints: 0 means a\n",
    "set size that does not move, 1 means it typically doubles or halves between checkpoints.\n",
    "point = the production rule (vs max(mean_score)); racing = vs max(mean - z*se).\n", sep = "")
jump <- function(x) { x <- x[is.finite(x) & x > 0]
  if (length(x) < 3) NA_real_ else stats::median(abs(diff(log2(x)))) }
pr(A |> dplyr::filter(is.finite(sd_config), sd_config > 0) |>
     dplyr::group_by(run, z) |> dplyr::arrange(n_eval, .by_group = TRUE) |>
     dplyr::summarise(point = jump(in_play), racing = jump(in_play_lcb),
                      ip_med = stats::median(in_play),
                      lcb_med = stats::median(in_play_lcb), .groups = "drop") |>
     dplyr::filter(z %in% c(0.25, 0.5, 1)) |> dplyr::arrange(z, run))
cat("\noverlap_med is out of 8: it is how many parents would CHANGE if get_elites() were\n",
    "replaced by the contender set at that z. 8 = no change, 0 = a completely new pool.\n",
    sep = "")
