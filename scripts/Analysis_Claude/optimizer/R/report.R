# report.R
#
# Make the run interpretable, not a black box: the learning curve, the current
# best pipeline, and which subtask METHODS most raise the score. Written after
# every checkpoint so you can watch progress while it runs in the background.

library(tidyverse)

# A human-readable line for one configuration.
format_config <- function(cfg) {
  parts <- vapply(names(SUBTASKS), function(st) {
    m <- cfg[[paste0(st, ".method")]]
    ps <- names(SUBTASKS[[st]]$params)
    pv <- vapply(ps, function(p) {
      v <- cfg[[paste0(st, ".", p)]]
      if (length(v) == 0 || is.na(v)) NA_character_
      else if (is.numeric(v)) paste0(p, "=", signif(v, 3))
      else paste0(p, "=", v)
    }, character(1))
    pv <- pv[!is.na(pv)]
    paste0("  ", st, ": ", m, if (length(pv)) paste0(" (", paste(pv, collapse = ", "), ")") else "")
  }, character(1))
  paste(parts, collapse = "\n")
}

# Marginal effect of each subtask METHOD: mean score of all configs using it,
# minus the overall mean. Positive = that method tends to help. This is the
# "which submitted ideas actually win" view.
method_importance <- function(evals) {
  if (!nrow(evals)) return(tibble::tibble())
  cfgs <- lapply(evals$config_json, config_from_json)
  feats <- configs_to_features(cfgs)
  overall <- mean(evals$score[is.finite(evals$score)], na.rm = TRUE)
  purrr::map_dfr(names(SUBTASKS), function(st) {
    key <- paste0(st, ".method")
    tibble::tibble(subtask = st,
                   method  = as.character(feats[[key]]),
                   score   = evals$score) |>
      dplyr::filter(is.finite(score)) |>
      dplyr::group_by(subtask, method) |>
      dplyr::summarise(n = dplyr::n(), mean_score = mean(score),
                       delta = mean(score) - overall, .groups = "drop")
  }) |>
    dplyr::arrange(subtask, dplyr::desc(mean_score))
}

# Failure-log analysis. Every (config, trial, scheme) attempt is one eval row;
# status is "ok" | "infeasible" | "error". This summarises (a) how attempts break
# down by status, (b) which infeasibility reasons dominate, and (c) the failure
# rate of each subtask METHOD -- i.e. which (often demanding) method choices fail
# most often. Returns a list of tibbles; write_report() renders a compact view.
failure_summary <- function(evals) {
  empty <- list(by_status = tibble::tibble(), by_reason = tibble::tibble(),
                by_method = tibble::tibble(), suspect = tibble::tibble())
  if (!nrow(evals)) return(empty)

  by_status <- evals |>
    dplyr::count(status, name = "n") |>
    dplyr::arrange(dplyr::desc(n))

  by_reason <- evals |>
    dplyr::filter(status %in% c("infeasible", "suspect", "error")) |>
    dplyr::count(status, reason, name = "n") |>
    dplyr::arrange(dplyr::desc(n))

  # The rows that most likely indicate a BUG (data that should be visible isn't):
  # the actual trial + reason + funnel, so they can be drilled into directly.
  suspect <- if ("detail" %in% names(evals)) {
    evals |>
      dplyr::filter(status == "suspect") |>
      dplyr::distinct(trial_id, reason, detail) |>
      utils::head(50)
  } else tibble::tibble()

  # Per-subtask-method failure rate (failed = anything that did not score "ok").
  cfgs  <- lapply(evals$config_json, config_from_json)
  feats <- configs_to_features(cfgs)
  failed <- evals$status != "ok"
  by_method <- purrr::map_dfr(names(SUBTASKS), function(st) {
    key <- paste0(st, ".method")
    tibble::tibble(subtask = st,
                   method  = as.character(feats[[key]]),
                   failed  = failed) |>
      dplyr::group_by(subtask, method) |>
      dplyr::summarise(n = dplyr::n(), n_fail = sum(failed),
                       fail_rate = mean(failed), .groups = "drop")
  }) |>
    dplyr::arrange(subtask, dplyr::desc(fail_rate))

  list(by_status = by_status, by_reason = by_reason, by_method = by_method,
       suspect = suspect)
}

# Write a Markdown snapshot to disk.
write_report <- function(con, settings) {
  all_evals <- read_evals(con)
  # Report on this optimization's own domain + scheme (what the surrogate actually
  # learns from); keep the global count for context.
  td      <- if (isTRUE(settings$simulate)) NULL else settings$target_domain
  evals   <- filter_evals_to_domain(all_evals, td) |>
               filter_evals_to_scheme(settings$optimize_scheme) |>
               filter_evals_to_build(settings$build %||% OPTIMIZER_BUILD)
  n_other <- nrow(all_evals) - nrow(evals)
  agg   <- aggregate_scores(evals)
  inc   <- incumbent_config(agg, settings$incumbent_min_reps)

  ok <- evals |> dplyr::filter(is.finite(score))
  # Best-so-far learning curve over evaluation order.
  curve <- ok |>
    dplyr::arrange(id) |>
    dplyr::mutate(running_best = cummax(score))

  seed_hashes <- vapply(seed_configs(settings$optimize_scheme), config_hash, character(1))
  seed_scores <- agg |>
    dplyr::filter(config_hash %in% seed_hashes, is.finite(mean_score)) |>
    dplyr::summarise(best_seed = suppressWarnings(max(mean_score))) |>
    dplyr::pull(best_seed)
  best_seed <- if (length(seed_scores) && is.finite(seed_scores)) seed_scores else NA_real_

  lines <- c(
    "# Optimizer report",
    paste0("_", format(Sys.time(), tz = "UTC", usetz = TRUE), "_"),
    "",
    paste0("- optimizer build: ", settings$build %||% OPTIMIZER_BUILD),
    # Whether the store is actually being backed up. This belongs in the report because the
    # report is NOT leader-gated: it keeps updating while a long evaluation runs, so it is the
    # one artefact that can report a stalled backup while the stall is happening.
    local({
      bp <- settings$db_backup_path
      if (is.null(bp) || !nzchar(bp)) return(NULL)   # local mode: no backup configured
      db_min <- settings$db_backup_minutes %||% 0
      age <- backup_age_minutes(bp)
      paste0("- store backup: ",
             if (!is.finite(age)) "_never_"
             else if (age < 90) sprintf("%.0f min ago", age)
             else sprintf("%.1f h ago", age / 60),
             # Only a FINITE age can be stale: a run's first report legitimately finds no
             # backup yet, and flagging that would be a false alarm on every fresh start.
             if (is.finite(age) && db_min > 0 && age > 2 * db_min) "  ** STALE **" else "")
    }),
    # Which estimator ranked the configs, and -- when the random-effects fit ran -- the
    # variance components. sd_trial vs sd_config is the number that says whether adjusting for
    # trials is worth anything on this data; it used to have to be computed by hand.
    local({
      est <- attr(agg, "estimator") %||% "pooled"
      vc  <- attr(agg, "var_comps")
      paste0("- config score estimator: ", est,
             if (!is.null(vc) && all(is.finite(vc)))
               sprintf("  (sd_trial %.3f, sd_config %.3f, sd_resid %.3f)",
                       vc[["sd_trial"]], vc[["sd_config"]], vc[["sd_resid"]])
             else "")
    }),
    paste0("- optimized scheme: ", settings$optimize_scheme),
    paste0("- evaluations: ", nrow(evals),
           " (", sum(is.finite(evals$score)), " scored, ",
           sum(!is.finite(evals$score)), " failed)",
           if (n_other > 0) paste0("; ", n_other, " more in the store are out of this scheme/domain") else ""),
    paste0("- distinct configurations: ", nrow(agg)),
    paste0("- best single-trial score: ",
           ifelse(nrow(ok), signif(max(ok$score), 3), "NA")),
    paste0("- best seed (submission) mean score: ",
           ifelse(is.finite(best_seed), signif(best_seed, 3), "NA")),
    "",
    "## Incumbent (best mean over its trials)",
    if (is.null(inc)) "_none yet_" else
      c(paste0("mean score ", signif(inc$mean_score, 3),
               " over ", inc$n_ok, " trials",
               if (is.finite(best_seed)) paste0("  (vs best seed ", signif(best_seed, 3),
                 if (inc$mean_score > best_seed) " -- BEATS submissions)" else ")") else ""),
        "```", format_config(inc$config), "```"),
    "",
    "## Subtask-method importance (mean score, delta vs overall)"
  )

  imp <- method_importance(evals)
  if (nrow(imp)) {
    imp_tbl <- imp |>
      dplyr::mutate(line = sprintf("  %-13s %-20s n=%-4d  mean=%+.3f  delta=%+.3f",
                                   subtask, method, n, mean_score, delta)) |>
      dplyr::pull(line)
    lines <- c(lines, "```", imp_tbl, "```")
  }

  # Failure log: status breakdown, dominant infeasibility reasons, and the
  # failure rate of each subtask method (which method choices break most often).
  fs <- failure_summary(evals)
  if (nrow(fs$by_status)) {
    lines <- c(lines, "", "## Failure log",
               "```",
               "by status:",
               sprintf("  %-12s %d", fs$by_status$status, fs$by_status$n))
    if (nrow(fs$by_reason)) {
      lines <- c(lines, "", "by reason (infeasible / error):",
                 sprintf("  %-11s %-28s %d",
                         fs$by_reason$status, fs$by_reason$reason, fs$by_reason$n))
    }
    if (nrow(fs$by_method)) {
      lines <- c(lines, "", "failure rate by subtask method:",
                 sprintf("  %-13s %-20s n=%-4d  fail=%-4d  rate=%.2f",
                         fs$by_method$subtask, fs$by_method$method,
                         fs$by_method$n, fs$by_method$n_fail, fs$by_method$fail_rate))
    }
    lines <- c(lines, "```")

    # Suspected bugs get their own section: these failed with a funnel cliff
    # (data that should be visible was not), so they are the ones to investigate
    # with diagnose_trial() before trusting the "infeasible" verdict.
    n_suspect <- fs$by_status$n[fs$by_status$status == "suspect"]
    n_suspect <- if (length(n_suspect)) n_suspect else 0L
    if (n_suspect > 0) {
      lines <- c(lines, "",
                 paste0("### ⚠ Suspected bugs (", n_suspect,
                        " evals -- funnel cliff, NOT necessarily real infeasibility)"))
      if (nrow(fs$suspect)) {
        lines <- c(lines, "```",
                   sprintf("  trial=%-12s %-26s %s",
                           fs$suspect$trial_id, fs$suspect$reason,
                           ifelse(is.na(fs$suspect$detail), "", fs$suspect$detail)),
                   "```",
                   "_Investigate with_ `diagnose_trial(<trial_id>, settings)`.")
      }
    }
  }

  if (nrow(curve) > 1) {
    # A coarse text sparkline of the running best.
    rb <- curve$running_best
    idx <- unique(round(seq(1, length(rb), length.out = min(20, length(rb)))))
    lines <- c(lines, "", "## Running best (over evaluation order)",
               "```",
               paste(sprintf("%.3f", rb[idx]), collapse = "  "),
               "```")
  }

  # ATOMIC: every worker writes this file (see run_optimizer.R -- keying it to the leader's
  # iteration count meant the report only refreshed when worker 1 finished an evaluation,
  # which with 8 workers is an eighth of the activity and can be hours apart). Concurrent
  # writeLines() to one path would interleave, and a `cat report.md` could catch a half-written
  # file, so render to a sibling temp and rename into place.
  tmp <- paste0(settings$report_path, ".tmp", Sys.getpid())
  writeLines(lines, tmp)
  if (!isTRUE(file.rename(tmp, settings$report_path))) {
    unlink(tmp)
    message("report write -> ", settings$report_path, " FAILED")
  }
  invisible(settings$report_path)
}
