# dev/holdout_test.R
#
# Do the optimized configurations beat the five submissions on trials the optimizer NEVER SAW?
#
# The design this serves: choose a target domain of trials resembling the one you care about,
# optimize on it, then test the optimized configurations AND the submissions on held-out trials.
# Without that last step an incumbent's score is a maximum over noisy means on the trials it was
# selected on -- optimistic by construction (docs/BACKGROUND.md sec. 4). This is the measurement
# the project's headline claim rests on.
#
#   cd <repo>/scripts/Analysis_Claude/optimizer
#   Rscript dev/holdout_test.R --trials=10673,10675,Big6_2026_URB       # evaluate the grid
#   Rscript dev/holdout_test.R --trials=... --dry-run                   # what it WOULD run
#   Rscript dev/holdout_test.R --report                                 # analyse what is stored
#
# Options:
#   --trials=<ids or names>  the HELD-OUT trials. Names are resolved to ids (T3 renames them).
#   --select=mean|ucb        how the tested configurations are chosen; default `mean`.
#   --k=8                    how many of them.
#   --store=<path>           the optimizer store to read configurations from.
#   --out=<path>             where results go; default $OPTIMIZER_HOME/state/holdout.sqlite.
#   --dry-run                resolve and print the grid, evaluate nothing.
#   --report                 analysis only, no evaluation.
#
# WHY `mean` AND NOT `ucb` BY DEFAULT. .contenders() is an ELIMINATION rule -- `ucb >= max(mean)`
# asks "which configurations can I not yet rule out?", which is the right question for deciding
# where to spend more evaluations (its only functional use is .replication_backlog()). It is not
# a ranking: incumbent_config() and get_elites() both sort on mean_score. At equal replication
# the two orders coincide, but a UCB set can admit an under-replicated configuration that is
# merely uncertain rather than good -- which would bias this comparison downward. `mean` is what
# you would actually deploy. `--select=ucb` reproduces the report's set.
#
# SAFETY: results go to a SEPARATE store, never the optimizer's. Beyond keeping the run's slice
# counts clean, it makes it impossible for a held-out evaluation to become training data if the
# target domain is ever widened.

# The optimizer ROOT, not this script's directory: `.Renviron` -- and so the T3
# credentials -- is read from the WORKING DIRECTORY only, with no parent walk. `settings.R`
# is the marker for that root. here::i_am() below also halts from the wrong place, but only
# when it cannot find the project at all, and a cwd inside the project is not one of those
# cases: here() still resolves while .Renviron silently does not.
if (!file.exists("settings.R"))
  stop("run this from the optimizer ROOT, so R reads ./.Renviron:\n",
       "  cd <repo>/scripts/Analysis_Claude/optimizer\n",
       "  Rscript dev/holdout_test.R --trials=<ids>\n",
       "  (working directory was: ", getwd(), ")")

suppressMessages(library(tidyverse))
here::i_am("dev/holdout_test.R")
source(here::here("settings.R"))
for (f in list.files(here::here("R"), pattern = "[.]R$", full.names = TRUE)) source(f)

args <- commandArgs(trailingOnly = TRUE)
opt  <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^--", n, "="), "", h[1]) }
o_trials <- opt("trials"); o_sel <- opt("select", "mean")
o_k      <- suppressWarnings(as.integer(opt("k", "8")))
o_dry    <- "--dry-run" %in% args
o_report <- "--report"  %in% args
if (!o_sel %in% c("mean", "ucb")) stop("--select must be `mean` or `ucb`")

s      <- optimizer_settings()
scheme <- s$optimize_scheme
hr <- function(t) cat("\n\n", t, "\n", strrep("=", nchar(t)), "\n", sep = "")

# ---- where the results live ------------------------------------------------
# The working store goes on NODE-LOCAL disk and is backed up to durable storage, exactly as the
# optimizer's does. SQLite's WAL needs an mmap'd -shm file that a network filesystem does not
# provide, so a store under $OPTIMIZER_HOME cannot coordinate N writers -- open_store() says so
# ("ONE worker only"). This script runs N copies, so the distinction is load-bearing.
#
# `hs` is a settings COPY with just those two paths repointed: restore_store_from_backup() reads
# only db_path and db_backup_path, so it drives the holdout pair unchanged.
o_out <- opt("out")
if (!is.null(o_out)) {                      # explicit path: single file, caller's problem
  work_db <- o_out; back_db <- NULL
} else if (!is.null(s$db_backup_path) && nzchar(s$db_backup_path %||% "")) {
  work_db <- file.path(dirname(s$db_path), "holdout.sqlite")          # node-local, per cluster_scratch_paths
  back_db <- file.path(dirname(s$db_backup_path), "holdout_backup.sqlite")
} else {                                    # laptop: one disk, no NFS problem
  work_db <- file.path(dirname(s$db_path %||% "state/x"), "holdout.sqlite"); back_db <- NULL
}
hs <- modifyList(s, list(db_path = work_db, db_backup_path = back_db))
out_db <- work_db

# ---- the configurations under test -----------------------------------------
# Only needed to BUILD the grid. --report reads the group off each stored row, so it must not
# require the optimizer's store to be present -- the analysis is often run somewhere else, long
# after the evaluations.
grid_cfg <- list(); grid_grp <- character(); tested <- list(); seeds <- list()
universe <- NULL
if (!o_report) {
# Read the optimizer's store through a copy, as the tools do: it may be live.
src <- resolve_read_store(opt("store"), "dev/holdout_test.R")
con_src <- open_store(.copy_store_with_sidecars(src, file.path(tempdir(), "holdout_src.sqlite")))
src_evals <- read_evals(con_src)
rid  <- run_id_for(s, s$build %||% OPTIMIZER_BUILD)
rrow <- tryCatch(read_run(con_src, rid), error = function(e) NULL)
universe <- run_universe(rrow)
close_store(con_src)

slice <- src_evals |>
  filter_evals_to_scheme(scheme) |>
  filter_evals_to_build(s$build %||% OPTIMIZER_BUILD)
if (length(universe)) slice <- filter_evals_to_universe(slice, universe)
agg <- aggregate_scores(slice)
if (!nrow(agg)) stop("no configurations in the optimizer store's own domain/scheme/build slice")

picked <- if (identical(o_sel, "ucb")) {
  .contenders(agg, s$contender_z %||% 1, k = o_k)
} else {
  agg |> dplyr::filter(is.finite(mean_score)) |>
    dplyr::arrange(dplyr::desc(mean_score)) |> head(o_k) |> dplyr::pull(config_hash)
}
tested <- lapply(picked, function(h) config_from_json(agg$config_json[match(h, agg$config_hash)]))
names(tested) <- paste0("opt", seq_along(tested))
seeds <- seed_configs(scheme)
# A seed that the optimizer also ranks highly would otherwise be counted in both groups and
# make the comparison partly self-referential.
seed_h <- vapply(seeds, config_hash, character(1))
dupe   <- intersect(seed_h, picked)
if (length(dupe)) {
  keep <- !(picked %in% dupe)
  message("note: ", sum(!keep), " selected configuration(s) ARE seeds; kept in the seed group only")
  tested <- tested[keep]; picked <- picked[keep]
}
grid_cfg <- c(tested, seeds)
grid_grp <- c(rep("contender", length(tested)), rep("seed", length(seeds)))
}

# ---- the held-out trials ---------------------------------------------------
conn <- NULL
resolve_trials <- function(spec) {
  want <- trimws(strsplit(spec, ",")[[1]])
  conn <<- t3_connect(s)
  cat0 <- trial_catalog(conn, s)
  ids  <- as.character(cat0$study_db_id); nms <- as.character(cat0$study_name)
  out <- vapply(want, function(w) {
    if (w %in% ids) return(w)
    i <- match(w, nms)
    if (!is.na(i)) return(ids[i])
    NA_character_
  }, character(1))
  out <- unname(out)          # vapply names its result; identical() below would see the name
  bad <- want[is.na(out)]
  if (length(bad))
    stop("not in the focal-trait catalogue (no grain-yield data, or the name changed): ",
         paste(bad, collapse = ", "),
         "\n  dev/check_domain_coverage.R explains which of the two it is.", call. = FALSE)
  for (i in seq_along(want)) if (!identical(want[i], out[i]))
    cat("  resolved ", want[i], " -> ", out[i], "\n", sep = "")
  out
}

if (!o_report) {
  if (is.null(o_trials)) stop("--trials= is required unless --report")
  hr("held-out trials")
  new_ids <- resolve_trials(o_trials)
  # THE guard. A trial the optimizer trained on is not held out, and nothing downstream would
  # notice -- the result would simply be wrong in the optimistic direction.
  overlap <- intersect(new_ids, as.character(universe %||% character()))
  if (length(overlap))
    stop("these trials are IN the optimizer's pinned universe, so they are not held out:\n  ",
         paste(overlap, collapse = ", "),
         "\n  The whole point of this test is trials the optimizer never saw.", call. = FALSE)
  cat("  ", length(new_ids), " trial(s), none in the run's ",
      length(universe %||% character()), "-trial universe\n", sep = "")
  cat("  grid: ", length(grid_cfg), " configurations x ", length(new_ids), " trials = ",
      length(grid_cfg) * length(new_ids), " evaluations\n", sep = "")
  cat("  selected by ", o_sel, ": ", length(tested), " optimized + ", length(seeds), " seeds\n", sep = "")
  if (o_dry) {
    hr("configurations (dry run -- nothing evaluated)")
    for (i in seq_along(grid_cfg))
      cat("\n--- ", names(grid_cfg)[i], " (", grid_grp[i], ") ---\n",
          paste(format_config(grid_cfg[[i]]), collapse = "\n"), "\n", sep = "")
    quit(save = "no")
  }
}

# ---- evaluate the grid -----------------------------------------------------
# Claims make the grid safe for N concurrent copies and resumable after a kill, exactly as
# run_optimizer.R:60-95 uses them -- the only differences are that the (config, trial) pairs are
# fixed rather than chosen, and that everything lands in `out_db`.
dir.create(dirname(out_db), showWarnings = FALSE, recursive = TRUE)
# Merge whatever a previous job backed up, so a new node resumes rather than repeats. No-op when
# there is no durable path (laptop) or nothing backed up yet.
if (!is.null(back_db)) restore_store_from_backup(hs)
con <- open_store(out_db)
on.exit(close_store(con), add = TRUE)

if (!o_report) {
  cat("  working store -> ", out_db, "\n", sep = "")
  if (!is.null(back_db)) cat("  backed up to  -> ", back_db, "\n", sep = "")

  # The node-local cache starts EMPTY on a fresh node, so without this every VCF the durable
  # cache already holds is downloaded again. Copy 1 restores; the others wait on a ready file of
  # OUR OWN -- settings$cache_ready_file belongs to the optimizer, and a holdout job unlinking it
  # while an optimizer job is starting would be a real cross-job fault.
  ready <- if (!is.null(s$cache_backup_dir) && nzchar(s$cache_backup_dir %||% ""))
             file.path(dirname(s$db_backup_path %||% "state/x"), ".holdout_cache_ready") else NULL
  w <- suppressWarnings(as.integer(s$worker_id %||% "1")); if (!is.finite(w)) w <- 1L
  if (w == 1L) {
    if (!is.null(ready)) unlink(ready)
    restore_cache_from_backup(s)
    if (!is.null(ready)) file.create(ready)
  } else if (!is.null(ready)) {
    waited <- 0
    while (!file.exists(ready) && waited < 1800) { Sys.sleep(10); waited <- waited + 10 }
    if (waited > 0) message(sprintf("copy %d: waited %.0f s for the cache restore", w, waited))
  }
  # EVERY copy flushes on the way out, not just copy 1: whichever is killed last should still
  # have contributed what it downloaded (run_optimizer.R:156 does the same, for the same reason).
  on.exit(sync_cache_to_backup(s, min_age_minutes = 2), add = TRUE)

  if (is.null(conn)) conn <- t3_connect(s)
  done <- read_evals(con)
  done_key <- paste(done$config_hash, done$trial_id)
  # Offset by worker id so N copies start on different rows rather than contending for the
  # first cell -- the same reason run_optimizer.R offsets its seed pick.
  grid <- expand_grid(ci = seq_along(grid_cfg), ti = seq_along(new_ids))
  grid <- grid[c(seq(w, nrow(grid)), seq_len(max(0, w - 1))), , drop = FALSE]

  n_run <- 0L; n_skip <- 0L
  for (r in seq_len(nrow(grid))) {
    cfg <- grid_cfg[[grid$ci[r]]]; tid <- new_ids[grid$ti[r]]
    h   <- config_hash(cfg)
    if (paste(h, tid) %in% done_key) { n_skip <- n_skip + 1L; next }
    if (!claim_eval(con, h, tid, scheme, s$worker_id %||% NA_character_,
                    s$build %||% OPTIMIZER_BUILD)) { n_skip <- n_skip + 1L; next }
    trial <- tryCatch(build_trial_descriptor(tid, conn, s), error = function(e) NULL)
    if (is.null(trial)) {
      release_claim(con, h, tid, scheme)
      message("  ", tid, ": could not build a descriptor -- skipping"); next
    }
    cat(sprintf("[%s] %s (%s) on %s ...\n", format(Sys.time(), "%H:%M"),
                names(grid_cfg)[grid$ci[r]], grid_grp[grid$ci[r]], tid))
    ev <- evaluate_config_on_trial(cfg, trial, scheme, s, conn)
    store_eval(con, cfg, tid, scheme, ev$score, ev$n_test, ev$status, ev$reason,
               ev$detail %||% NA_character_, ev$seconds,
               study_name    = trial$study_name %||% NA_character_,
               program_name  = trial$program    %||% NA_character_,
               location_name = trial$location   %||% NA_character_,
               year          = trial$year       %||% NA_integer_,
               peak_rss_mb   = ev$peak_rss_mb %||% NA_real_,
               peak_r_mb     = ev$peak_r_mb   %||% NA_real_,
               rss_mb        = ev$rss_mb      %||% NA_real_,
               worker        = s$worker_id %||% NA_character_,
               dosage_budget = s$dosage_budget_bytes %||% NA_real_,
               em_df_method  = "effective_n",
               build         = s$build %||% OPTIMIZER_BUILD,
               # The GROUP is the one thing `evals` has no column for, and it is what the whole
               # analysis turns on. run_id is free text and unused here, so it carries it.
               run_id        = paste0("holdout:", grid_grp[grid$ci[r]]))
    release_claim(con, h, tid, scheme)
    # After every evaluation, as run_optimizer.R:270 does: the working store is node-local and
    # dies with the allocation, so the backup is the only durable record.
    if (!is.null(back_db)) backup_store(con, back_db)
    sync_cache_to_backup(s)
    cat(sprintf("     -> %s  score %s  n_test %s\n", ev$status,
                if (is.finite(ev$score)) sprintf("%+.3f", ev$score) else "--",
                ev$n_test %||% "--"))
    n_run <- n_run + 1L
  }
  cat(sprintf("\n%d evaluated, %d already done or claimed elsewhere\n", n_run, n_skip))
}

# ---- analysis --------------------------------------------------------------
e <- read_evals(con)
if (!nrow(e)) { cat("\nnothing stored yet.\n"); quit(save = "no") }
e$group <- sub("^holdout:", "", e$run_id %||% NA_character_)
e <- dplyr::filter(e, group %in% c("contender", "seed"))
lab <- e |> dplyr::distinct(config_hash, group) |>
  dplyr::arrange(group != "contender", config_hash) |>
  dplyr::mutate(label = paste0(ifelse(group == "contender", "opt", "seed"),
                               ave(seq_along(group), group, FUN = seq_along)))
e$label <- lab$label[match(e$config_hash, lab$config_hash)]

hr("1. completion: which cells scored")
grid_tbl <- e |>
  dplyr::mutate(ok = status == "ok" & is.finite(score)) |>
  dplyr::group_by(label, group) |>
  dplyr::summarise(trials = dplyr::n(), scored = sum(ok),
                   rate = scored / trials, .groups = "drop") |>
  dplyr::arrange(group != "contender", dplyr::desc(rate))
print(grid_tbl, n = Inf)
byt <- e |> dplyr::mutate(ok = status == "ok" & is.finite(score)) |>
  dplyr::group_by(study_name, trial_id) |>
  dplyr::summarise(configs = dplyr::n(), scored = sum(ok), .groups = "drop") |>
  dplyr::arrange(scored)
cat("\nby trial (a trial no configuration can score says nothing about either group):\n")
print(byt, n = Inf)
if (nrow(e[e$status != "ok", ]))
  print(dplyr::count(dplyr::filter(e, status != "ok"), status, reason, sort = TRUE), n = Inf)

ok <- dplyr::filter(e, is.finite(score))
if (dplyr::n_distinct(ok$group) < 2) { cat("\nboth groups need scores before the comparison.\n"); quit(save = "no") }

hr("2. the two distributions (mean over held-out trials, per configuration)")
per_cfg <- ok |> dplyr::group_by(label, group) |>
  dplyr::summarise(n = dplyr::n(), mean = mean(score), sd = stats::sd(score), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean))
print(per_cfg, n = Inf)
cs <- per_cfg$mean[per_cfg$group == "contender"]; ss <- per_cfg$mean[per_cfg$group == "seed"]
cat(sprintf("\n  optimized: n=%d mean %+.4f   seeds: n=%d mean %+.4f   difference %+.4f\n",
            length(cs), mean(cs), length(ss), mean(ss), mean(cs) - mean(ss)))
inv <- sum(outer(cs, ss, "<="))
cat(sprintf("  %d of %d optimized/seed pairs inverted\n", inv, length(cs) * length(ss)))
wt <- suppressWarnings(stats::wilcox.test(cs, ss, alternative = "greater"))
tt <- stats::t.test(cs, ss, alternative = "greater")
cat(sprintf("  Wilcoxon one-sided p = %.4f;  Welch t one-sided p = %.4f\n", wt$p.value, tt$p.value))

hr("3. primary: mixed model, configuration nested in group")
# group is a FIXED effect; configuration is RANDOM and supplies the error term for testing it.
# Without (1|config_hash), group would be tested against evaluation-level noise -- pseudo-
# replication across ~20 evaluations per configuration. A FIXED config term would alias group;
# a random one does not. Failures are simply absent rows: lmer uses what exists and the trial
# effect adjusts for WHICH trials each configuration completed.
fw <- .fisher_weights(ok)
fit <- tryCatch(suppressMessages(lme4::lmer(
        .z ~ group + (1 | trial_id) + (1 | config_hash), data = fw, weights = .w)),
        error = function(err) { cat("  fit failed: ", conditionMessage(err), "\n"); NULL })
if (!is.null(fit)) {
  co <- summary(fit)$coefficients
  i  <- grep("^group", rownames(co))
  est <- -co[i, 1]; se <- co[i, 2]        # negate: `contender` is the reference level
  cat(sprintf("  optimized - seeds = %+.4f  (SE %.4f, 95%% CI %+.4f to %+.4f) on the z scale\n",
              est, se, est - 1.96 * se, est + 1.96 * se))
  cat(sprintf("  on the correlation scale: %+.4f\n", tanh(est)))
  v <- as.data.frame(lme4::VarCorr(fit))
  cat(sprintf("  sd_trial %.3f  sd_config %.3f  sd_resid %.3f (unit weight)\n",
              v$sdcor[v$grp == "trial_id"][1], v$sdcor[v$grp == "config_hash"][1],
              v$sdcor[v$grp == "Residual"][1]))
}

hr("4. how much to trust it")
# Contenders are built by crossover and mutation from shared elites, so they are NOT independent
# draws. The number of distinct METHOD SKELETONS among them is the honest handle on how much
# less than n they are worth.
# A malformed config_json must not take down a report whose other sections succeeded, so an
# unreadable one becomes "?" rather than an error.
skel <- function(h) {
  j <- e$config_json[match(h, e$config_hash)]
  cfg <- tryCatch(config_from_json(j), error = function(err) list())
  one <- function(x) { x <- as.character(x); if (length(x) == 1L) x else "?" }
  paste(vapply(names(SUBTASKS), function(st)
    one(cfg[[paste0(st, ".method")]] %||% "?"), character(1)), collapse = "|")
}
for (g in c("contender", "seed")) {
  hs <- unique(e$config_hash[e$group == g])
  cat(sprintf("  %-10s %d configuration(s), %d distinct method skeleton(s)\n",
              g, length(hs), length(unique(vapply(hs, skel, character(1))))))
}
cat("\n  Contenders sharing a skeleton are not independent evidence: SMAC builds them from the\n")
cat("  same elites. Read the tests above against the SKELETON count, not the configuration count.\n")

inc <- per_cfg |> dplyr::filter(group == "contender") |> dplyr::slice_max(mean, n = 1)
bs  <- per_cfg |> dplyr::filter(group == "seed")      |> dplyr::slice_max(mean, n = 1)
if (nrow(inc) && nrow(bs))
  cat(sprintf("\n  DESCRIPTIVE, not the claim -- best optimized %s %+.4f vs best seed %s %+.4f.\n",
              inc$label, inc$mean, bs$label, bs$mean),
      "  A maximum over configurations is biased upward; the group comparison above is the result.\n",
      sep = "")
cat("\n\ndone.\n")
