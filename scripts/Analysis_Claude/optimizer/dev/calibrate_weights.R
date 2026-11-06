# dev/calibrate_weights.R
#
# What weight should one evaluation carry? Solves for the saturation cap in
#
#     1/w_i = 1/(n_test_i - 3)  +  1/C                     C = 1 / sigma_int^2
#
# the weight implied by TWO independent error sources in a single score:
#
#   SAMPLING ERROR   variance 1/(n_test - 3) on the Fisher-z scale. Shrinks as the trial grows.
#                    This is the only term the current weighting models (R/optimizer.R:112).
#   CONFIG x TRIAL   the configuration genuinely performing differently on THIS trial than on
#   INTERACTION      an average one; variance sigma_int^2. Does NOT shrink with trial size -- a
#                    bigger trial measures that trial's idiosyncrasy more precisely, it does not
#                    make the idiosyncrasy smaller.
#
# C is the weight a trial of unbounded size would receive, so it reads directly as a trial size:
# C = 50 means "no trial may count for more than one with 53 usable accessions".
#
# HOW C IS IDENTIFIED. lmer reports sd_resid at UNIT weight, so the fitted per-row residual
# variance is sd_resid^2 / w_i. The weight above is correct exactly when that equals the true
# row variance sigma_int^2 + 1/(n-3) = 1/w_i -- that is, when
#
#     sd_resid == 1.
#
# So sd_resid is the calibration readout, not a nuisance parameter. It rises monotonically with
# C (bigger C, bigger weights, bigger unit-weight residual), so the crossing is found by
# fitting a grid of C and interpolating on log-log. The uncapped fit (C = Inf, today's
# behaviour) is included as the baseline: if ITS sd_resid is already ~1, no cap is warranted and
# n_test - 3 is the whole story.
#
# Runs the optimizer's own .fisher_weights() and .blup_scores(), overriding only the `.w`
# column, so the fit being calibrated is the fit in production.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript dev/calibrate_weights.R
#
# Options:
#   --store=<path>   store to read      (default: resolve_read_store(), as the tools do)
#   --runs=a,b       run_ids to analyse (default: every run with >= --min-evals rows)
#   --min-evals=500  skip runs shorter than this in the default run selection
#   --caps=5,10,...  the C grid (default 5,10,20,50,100,300,1000); Inf is always added
#
# COST: one mixed-model fit per cap per run. Minutes for a 50,000-row store.
#
# SAFETY: reads a COPY of the store with its -wal/-shm sidecars, never the live file, so it is
# safe to run beside a live search. Read-only throughout.

if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript dev/calibrate_weights.R\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
options(width = 150)
here::i_am("dev/calibrate_weights.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", n, "="), "", h[1]) }

s     <- optimizer_settings()
o_min <- suppressWarnings(as.integer(opt("min-evals", "500")))
o_cap <- suppressWarnings(as.numeric(strsplit(opt("caps", "5,10,20,50,100,300,1000"),
                                              ",")[[1]]))
o_cap <- sort(unique(c(o_cap[is.finite(o_cap) & o_cap > 0], Inf)))
min_t <- .min_test(s)

hr <- function(t) cat("\n\n", t, "\n", strrep("=", nchar(t)), "\n", sep = "")
pr <- function(d, k = 4) print(as.data.frame(
  dplyr::mutate(d, dplyr::across(where(is.numeric), ~round(.x, k)))), row.names = FALSE)

# ---- load, from a copy ------------------------------------------------------
store_path <- resolve_read_store(opt("store"), "dev/calibrate_weights.R")
con  <- open_store(.copy_store_with_sidecars(
          store_path, file.path(tempdir(), "calibrate_weights.sqlite")))
evals_all <- read_evals(con)
runs_tbl  <- if ("runs" %in% DBI::dbListTables(con))
               tryCatch(tibble::as_tibble(DBI::dbReadTable(con, "runs")),
                        error = function(e) NULL) else NULL
close_store(con)

if (!("run_id" %in% names(evals_all)))
  stop("this store predates the run_id column, so its rows cannot be split by run")

cat("store   : ", store_path, "  (", nrow(evals_all), " rows)\n", sep = "")
cat("min_test: ", min_t, "   (settings$min_test_acc -- the estimator's usability floor)\n",
    sep = "")
cat("caps    : ", paste(o_cap, collapse = ", "), "\n", sep = "")

by_run <- evals_all |>
  dplyr::mutate(run_id = dplyr::coalesce(run_id, "<none>")) |>
  dplyr::count(run_id, name = "n_eval")
want <- {
  if (!is.null(opt("runs"))) strsplit(opt("runs"), ",")[[1]]
  else by_run$run_id[by_run$n_eval >= o_min]
}
want <- want[want %in% by_run$run_id]
if (!length(want)) stop("no run has >= ", o_min, " rows; pass --runs= or lower --min-evals=")

scheme_of <- function(run) {
  if (is.null(runs_tbl)) return(NA_character_)
  x <- runs_tbl$scheme[runs_tbl$run_id == run]
  if (length(x)) x[1] else NA_character_
}

# ---- one fit at one cap ----------------------------------------------------
# `.blup_scores()` filters to .usable itself and reads `.w`, so capping the weight here is the
# whole intervention: everything else is the production fit.
fit_at <- function(d, cap) {
  d$.w <- if (is.infinite(cap)) d$.w_raw else 1 / (1 / d$.w_raw + 1 / cap)
  b <- tryCatch(.blup_scores(d), error = function(e) NULL)
  if (is.null(b)) return(NULL)
  vc <- attr(b, "var_comps")
  w  <- d$.w[d$.usable]
  tibble::tibble(
    cap       = cap,
    sd_trial  = vc[["sd_trial"]],
    sd_config = vc[["sd_config"]],
    sd_resid  = vc[["sd_resid"]],
    w_med     = stats::median(w),
    w_ratio   = max(w) / min(w),
    # the per-evaluation residual scale and the replication n* that SEARCH_ALGORITHM 4.5-4.6
    # report, recomputed at this cap so the consequence of the choice is visible
    resid_per_eval = vc[["sd_resid"]] / sqrt(stats::median(w)),
    n_star    = if (is.finite(vc[["sd_config"]]) && vc[["sd_config"]] > 0)
                  (vc[["sd_resid"]] / sqrt(stats::median(w)))^2 / vc[["sd_config"]]^2
                else Inf)
}

# Where does sd_resid cross 1? Monotone increasing in cap, so interpolate on log-log between
# the bracketing grid points. Returns NA with a reason when the grid does not bracket it.
solve_cap <- function(tab) {
  t2 <- dplyr::filter(tab, is.finite(cap), is.finite(sd_resid), sd_resid > 0) |>
        dplyr::arrange(cap)
  if (nrow(t2) < 2) return(list(cap = NA_real_, why = "fewer than 2 usable grid fits"))
  if (all(t2$sd_resid > 1))
    return(list(cap = NA_real_, why = paste0("sd_resid > 1 at every cap (min ",
                round(min(t2$sd_resid), 3), " at C=", t2$cap[1],
                "): the true cap is BELOW the grid")))
  if (all(t2$sd_resid < 1)) {
    # The uncapped fit is the C = Inf end of the same curve: if IT is above 1 the crossing is
    # simply above the grid, which is a different conclusion from "no cap helps".
    u <- tab$sd_resid[is.infinite(tab$cap)]
    if (length(u) && is.finite(u) && u > 1)
      return(list(cap = NA_real_, why = paste0(
        "crossing is ABOVE the top finite cap (sd_resid ", round(max(t2$sd_resid), 3),
        " at C=", t2$cap[nrow(t2)], ", ", round(u, 3),
        " uncapped) -- extend --caps= to bracket it")))
    return(list(cap = NA_real_, why = paste0("sd_resid < 1 at every cap (max ",
                round(max(t2$sd_resid), 3), " at C=", t2$cap[nrow(t2)],
                ") and uncapped -- the weight shape cannot be calibrated here; check whether ",
                "n_test is near-constant in this run, which makes the fit invariant to C")))
  }
  i <- max(which(t2$sd_resid <= 1))
  if (i == nrow(t2)) return(list(cap = t2$cap[i], why = "at the top grid point"))
  lo <- t2[i, ]; hi <- t2[i + 1L, ]
  lc <- log(lo$cap) + (log(1) - log(lo$sd_resid)) *
        (log(hi$cap) - log(lo$cap)) / (log(hi$sd_resid) - log(lo$sd_resid))
  list(cap = exp(lc), why = paste0("interpolated between C=", lo$cap, " and C=", hi$cap))
}

# ---- per run ---------------------------------------------------------------
summary_rows <- list()

for (run in want) {
  d <- evals_all |>
    dplyr::mutate(run_id = dplyr::coalesce(run_id, "<none>")) |>
    dplyr::filter(run_id == run) |>
    .fisher_weights(min_t, Inf)   # uncapped here; this script supplies the cap itself
  d$.w_raw <- d$.w                      # n_test - 3, floored at 1 -- the current weight

  hr(paste0("RUN ", run, "   scheme ", scheme_of(run)))
  nu <- sum(d$.usable)
  cat("rows ", nrow(d), "   usable ", nu, "   trials ",
      dplyr::n_distinct(d$trial_id[d$.usable]), "   configs ",
      dplyr::n_distinct(d$config_hash[d$.usable]), "\n", sep = "")
  if (nu) {
    q <- stats::quantile(d$.n_test[d$.usable], c(0, .5, 1), na.rm = TRUE)
    cat("n_test   min ", q[1], "  median ", q[2], "  max ", q[3], "\n", sep = "")
  }
  skip <- .blup_skip_reason(d)
  if (!is.null(skip)) { cat("SKIPPED: ", skip, "\n", sep = ""); next }

  tab <- purrr::map(o_cap, ~fit_at(d, .x)) |> purrr::compact() |> dplyr::bind_rows()
  if (!nrow(tab)) { cat("SKIPPED: no cap produced a fit\n"); next }

  cat("\n-- sd_resid against the cap (calibrated where sd_resid = 1)\n")
  pr(tab)

  sol <- solve_cap(tab)
  cat("\nsolved C : ", if (is.na(sol$cap)) "--" else signif(sol$cap, 4),
      "   (", sol$why, ")\n", sep = "")
  if (!is.na(sol$cap)) {
    f <- fit_at(d, sol$cap)
    cat("saturates at n_test = ", round(sol$cap + 3, 1),
        "   sigma_int = ", signif(sqrt(1 / sol$cap), 3), " (z scale)\n", sep = "")
    if (!is.null(f)) {
      cat("at that cap:\n"); pr(f)
      summary_rows[[run]] <- dplyr::mutate(f, run_id = run, scheme = scheme_of(run),
                                           .before = 1)
    }
  }
  # the uncapped baseline, stated plainly next to it
  u <- dplyr::filter(tab, is.infinite(cap))
  if (nrow(u))
    cat("uncapped (today): sd_resid ", round(u$sd_resid, 3),
        "  -- ", if (u$sd_resid > 1.15) "too HIGH: rows carry variance the weight does not model, so a cap is warranted"
                 else if (u$sd_resid < 0.85) "too LOW: n_test-3 already over-states the noise of large trials"
                 else "within 15% of 1: n_test-3 is close to right, a cap would change little",
        "\n", sep = "")
}

hr("SUMMARY: each run at its own calibrated cap")
if (!length(summary_rows)) cat("no run calibrated\n") else {
  pr(dplyr::bind_rows(summary_rows))
  cat("\nw_ratio is the largest-to-smallest weight within the run AT that cap --",
      "\nthe quantity the cap exists to limit. Compare it against the uncapped ratio above.\n")
}
cat("\n")
