# dev/sweep_contender_rule.R
#
# What the 2026-11-04 contender settings actually do, replayed over the real runs.
#
# Four changes landed together and three of them can only be judged against a store:
#
#   contender_rule = "racing"     the bar is max_j(mean_j - z*se_j), not max_j(mean_j)
#   contender_z    = 0.25         halved, because z does double duty under racing
#   resolution_delta = ?          the externally chosen difference worth resolving
#   sd_config_source = "reference"  shrink on the random draws' spread, not the search's
#
# This script replays the estimator over prefixes of each run's evaluation stream and, at every
# checkpoint, reports:
#
#   SD_CONFIG   fitted (whole evaluated population) beside reference (the surrogate-independent
#               draws). The reference value is the one that should NOT fall as the search
#               concentrates; `ratio` is how far apart the two priors are.
#   IN PLAY     set size under each (rule, z, delta) combination, as a count and as a fraction
#               of the scored population.
#   STABILITY   `jump`, the median |log2(in_play[t+1]/in_play[t])| over consecutive checkpoints:
#               0 means a set size that never moves, 1 that it typically doubles or halves.
#               This is the quantity the rule change is FOR.
#
# Every (rule, z, delta) combination is read off the same two fits per checkpoint, so the grid
# is nearly free; the cost is two mixed-model fits per checkpoint per run.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript dev/sweep_contender_rule.R
#
# Options:
#   --store=<path>     store to read       (default: resolve_read_store(), as the tools do)
#   --runs=a,b         run_ids to analyse  (default: every run with >= --min-evals rows)
#   --min-evals=500    skip runs shorter than this in the default run selection
#   --checkpoints=12   prefix points per run in the replay
#   --zs=0.25,0.5,1    the z grid
#   --deltas=0,0.002,0.005,0.01,0.02   the resolution-floor grid
#   --band=0.10        "usable" upper bound for in_play, as a fraction of the scored population
#
# SAFETY: reads a COPY of the store with its -wal/-shm sidecars, never the live file. Read-only.

if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript dev/sweep_contender_rule.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
options(width = 180)
here::i_am("dev/sweep_contender_rule.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(n, d = NULL) {
  h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", n, "="), "", h[1])
}
nums <- function(x) { v <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
                      sort(unique(v[is.finite(v) & v >= 0])) }

o_min  <- suppressWarnings(as.integer(opt("min-evals", "500")))
o_chk  <- max(3L, suppressWarnings(as.integer(opt("checkpoints", "12"))))
o_band <- suppressWarnings(as.numeric(opt("band", "0.10")))
o_zs   <- nums(opt("zs", "0.25,0.5,1"))
o_dl   <- nums(opt("deltas", "0,0.002,0.005,0.01,0.02"))
CAP    <- 8L

# Two settings objects differing in ONE field, so the comparison is exactly that field.
s_fit <- optimizer_settings(extra = list(sd_config_source = "fitted"))
s_ref <- optimizer_settings(extra = list(sd_config_source = "reference"))

hr <- function(t) cat("\n\n", t, "\n", strrep("=", nchar(t)), "\n", sep = "")
pr <- function(d, k = 4) print(as.data.frame(
  dplyr::mutate(d, dplyr::across(where(is.numeric), ~round(.x, k)))), row.names = FALSE)

store_path <- resolve_read_store(opt("store"), "dev/sweep_contender_rule.R")
con  <- open_store(.copy_store_with_sidecars(
          store_path, file.path(tempdir(), "sweep_contender_rule.sqlite")))
evals_all <- read_evals(con)
close_store(con)
if (!("run_id" %in% names(evals_all)))
  stop("this store predates the run_id column, so its rows cannot be split by run")

cat("store      : ", store_path, "  (", nrow(evals_all), " rows)\n", sep = "")
cat("live spec  : rule ", s_ref$contender_rule, "  z ", s_ref$contender_z,
    "  delta ", s_ref$resolution_delta, "  sd_config_source ", s_ref$sd_config_source,
    "\n", sep = "")
cat("grids      : z ", paste(o_zs, collapse = ", "), "   delta ",
    paste(o_dl, collapse = ", "), "\n", sep = "")
cat("checkpoints: ", o_chk, "    usable band: ", CAP, " <= in_play <= ", o_band,
    " * scored\n", sep = "")
# How much of the store can even be read as a reference population. `source` has only been
# stored since 2026-11-03, so on these runs .reference_hashes() takes its arrival-order path --
# which is what the `ref_src` column below records.
cat("source col : ", sum(!is.na(evals_all$source %||% NA)), " of ", nrow(evals_all),
    " rows carry a phase\n", sep = "")

by_run <- evals_all |>
  dplyr::mutate(run_id = dplyr::coalesce(run_id, "<none>")) |>
  dplyr::count(run_id, name = "n_eval")
want <- {
  if (!is.null(opt("runs"))) strsplit(opt("runs"), ",")[[1]]
  else by_run$run_id[by_run$n_eval >= o_min]
}
want <- want[want %in% by_run$run_id]
if (!length(want)) stop("no run has >= ", o_min, " rows; pass --runs= or lower --min-evals=")

run_rows <- function(run) {
  e <- evals_all |>
    dplyr::mutate(run_id = dplyr::coalesce(run_id, "<none>")) |>
    dplyr::filter(run_id == run)
  if ("id" %in% names(e)) dplyr::arrange(e, ts, id) else dplyr::arrange(e, ts)
}

sd_of <- function(agg, nm) {
  v <- attr(agg, "var_comps")
  if (is.null(v) || !(nm %in% names(v))) NA_real_ else as.numeric(v[[nm]])
}

sweep_run <- function(ev, run) {
  cuts <- unique(round(seq(min(200L, nrow(ev)), nrow(ev), length.out = o_chk)))
  comp <- list(); grid <- list()
  for (cz in cuts) {
    pre <- ev[seq_len(cz), , drop = FALSE]
    af  <- aggregate_scores(pre, s_fit)
    ar  <- aggregate_scores(pre, s_ref)
    nsc <- sum(is.finite(ar$mean_score))
    if (!nsc) next
    # Which path .reference_hashes() took, and on how many configurations.
    w   <- .fisher_weights(pre, .min_test(s_ref), .weight_cap(s_ref))
    rh  <- .reference_hashes(w, s_ref)
    tagged <- any(!is.na(w$source %||% NA) & (w$source %in% .REFERENCE_SOURCES))
    comp[[length(comp) + 1L]] <- tibble::tibble(
      run = run, n_eval = cz, n_cfg = nsc,
      sd_fit = sd_of(af, "sd_config"), sd_ref = sd_of(ar, "sd_config_ref"),
      ratio  = sd_of(ar, "sd_config_ref") / sd_of(af, "sd_config"),
      n_ref  = length(rh), ref_src = if (tagged) "source" else "arrival",
      se_fit = stats::median(af$se, na.rm = TRUE),
      se_ref = stats::median(ar$se, na.rm = TRUE))
    for (est in c("fitted", "reference")) {
      a <- if (identical(est, "fitted")) af else ar
      for (z in o_zs) for (dl in o_dl) for (rl in c("point", "racing")) {
        n <- length(.in_play(a, z = z, delta = dl, rule = rl))
        grid[[length(grid) + 1L]] <- tibble::tibble(
          run = run, n_eval = cz, n_cfg = nsc, est = est, rule = rl, z = z, delta = dl,
          in_play = n, pct = 100 * n / nsc)
      }
    }
  }
  list(comp = dplyr::bind_rows(comp), grid = dplyr::bind_rows(grid))
}

# |log2| of consecutive ratios: scale-free, so a set of 10 and a set of 400 are comparable.
jump <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 3L) return(NA_real_)
  stats::median(abs(diff(log2(x))), na.rm = TRUE)
}

comps <- list(); grids <- list()
for (run in want) {
  ev <- run_rows(run)
  hr(paste0("RUN ", run, "   ", nrow(ev), " evaluations"))
  r <- sweep_run(ev, run)
  if (!nrow(r$comp)) { cat("no checkpoint produced a usable estimate\n"); next }
  comps[[run]] <- r$comp; grids[[run]] <- r$grid
  cat("\n-- sd_config: the whole evaluated population vs the reference draws\n")
  pr(r$comp, 4)
}
if (!length(grids)) stop("nothing to summarise")
comp <- dplyr::bind_rows(comps); grid <- dplyr::bind_rows(grids)

# ---------------------------------------------------------------------------------------------
hr("1. Does the reference estimate hold up where the fitted one collapses?")
cat("sd_fit is the REML estimate on everything evaluated; sd_ref is the same estimator on the\n",
    "surrogate-independent draws alone. A ratio above 1 means the fitted prior is tighter than\n",
    "the configuration space warrants, which is the collapse SEARCH_ALGORITHM.md section 11\n",
    "describes. NA sd_ref means no reference estimate was available and the run keeps the\n",
    "fitted value.\n", sep = "")
pr(comp |>
     dplyr::group_by(run) |>
     dplyr::summarise(chk = dplyr::n(),
                      ref_src   = dplyr::first(ref_src),
                      n_ref_med = stats::median(n_ref),
                      sd_fit_first = dplyr::first(sd_fit), sd_fit_last = dplyr::last(sd_fit),
                      sd_ref_med   = stats::median(sd_ref, na.rm = TRUE),
                      ratio_med    = stats::median(ratio, na.rm = TRUE),
                      n_fit_zero   = sum(sd_fit <= 1e-9, na.rm = TRUE),
                      n_ref_na     = sum(is.na(sd_ref)),
                      se_ratio_med = stats::median(se_ref / se_fit, na.rm = TRUE),
                      .groups = "drop"), 4)

# ---------------------------------------------------------------------------------------------
hr("2. Set size and stability, by rule x z x delta")
cat("`jump` is the median |log2| change in in_play between consecutive checkpoints -- the\n",
    "quantity the rule change is for. `usable` is the share of checkpoints with\n", CAP,
    " <= in_play <= ", o_band, " * scored. Pooled over runs whose sd_config is not pinned to 0.\n",
    sep = "")
live <- comp |> dplyr::group_by(run) |>
  dplyr::summarise(dead = all(sd_fit <= 1e-9, na.rm = TRUE), .groups = "drop")
for (e in unique(grid$est)) {
  cat("\n-- sd_config_source = ", e, ", runs with a non-zero fitted sd_config\n", sep = "")
  pr(grid |>
       dplyr::filter(est == e, run %in% live$run[!live$dead]) |>
       dplyr::group_by(rule, z, delta) |>
       dplyr::summarise(ip_med = stats::median(in_play), pct_med = stats::median(pct),
                        jump   = stats::median(vapply(split(in_play, run), jump, numeric(1)),
                                               na.rm = TRUE),
                        usable = mean(in_play >= CAP & pct <= 100 * o_band),
                        .groups = "drop") |>
       dplyr::arrange(rule, z, delta), 3)
}

# ---------------------------------------------------------------------------------------------
hr("3. Per-run jump at the live spec vs the old one")
old_spec <- list(rule = "point", z = 1, delta = 0, est = "fitted")
new_spec <- list(rule = s_ref$contender_rule, z = s_ref$contender_z,
                 delta = s_ref$resolution_delta, est = s_ref$sd_config_source)
pick <- function(sp, tag) grid |>
  dplyr::filter(est == sp$est, rule == sp$rule, z == sp$z, delta == sp$delta) |>
  dplyr::group_by(run) |>
  dplyr::summarise(!!paste0("ip_", tag)   := stats::median(in_play),
                   !!paste0("pct_", tag)  := stats::median(pct),
                   !!paste0("jump_", tag) := jump(in_play), .groups = "drop")
if (all(c(old_spec$z, new_spec$z) %in% o_zs) &&
    all(c(old_spec$delta, new_spec$delta) %in% o_dl)) {
  pr(dplyr::left_join(pick(old_spec, "old"), pick(new_spec, "new"), by = "run"), 3)
  cat("\nold = point rule, z 1, delta 0, fitted sd_config (the settings through 2026-11-03)\n")
  cat("new = ", new_spec$rule, " rule, z ", new_spec$z, ", delta ", new_spec$delta, ", ",
      new_spec$est, " sd_config\n", sep = "")
} else {
  cat("pass --zs and --deltas covering both specs to get this table:\n  old z ",
      old_spec$z, " delta ", old_spec$delta, "   new z ", new_spec$z, " delta ",
      new_spec$delta, "\n", sep = "")
}

cat("\nTwo fits per checkpoint; every rule x z x delta cell is read off those, so the grid is\n",
    "nearly free. Read `jump` beside ip_med: a larger set is mechanically steadier, so a rule\n",
    "is only better when it wins at COMPARABLE set size.\n", sep = "")
