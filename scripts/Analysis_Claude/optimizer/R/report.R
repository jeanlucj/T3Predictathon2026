# report.R
#
# The Markdown snapshot of a run: score precision, current best pipeline, and which subtask
# METHODS most raise the score. Rewritten at every checkpoint, so a background run can be
# watched from the file.

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

# NA prints as "--" rather than "NA": a blank-looking cell in a precision table reads as zero.
fmt_num <- function(x, digits) if (!is.finite(x)) "--" else formatC(x, format = "f", digits = digits)

# Marginal effect of each subtask METHOD: mean score of the configs using it, minus the overall
# mean. Positive = the method tends to help.
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

# Pull one key out of the funnel string evaluate.R stores in `evals$detail`, which
# funnel_string() (R/conditions.R) writes as "k=v; k=v". Vectorised over a column; NA for a row
# with no funnel (scoring statuses carry a reason instead) or a key the funnel does not hold.
#
# Shared so the funnel is parsed in ONE place: tools/inspect_failures.R and
# dev/analyze_failures.R both read the same six fields, and two parsers of one format drift.
.funnel_get <- function(detail, key) {
  unname(vapply(detail, function(x) {
    if (is.na(x)) return(NA_real_)
    p <- strsplit(trimws(strsplit(x, ";", fixed = TRUE)[[1]]), "=", fixed = TRUE)
    v <- stats::setNames(vapply(p, function(z) z[2] %||% NA_character_, character(1)),
                         vapply(p, `[`, character(1), 1))
    suppressWarnings(as.numeric(v[key]))
  }, numeric(1)))
}

# Failure-log analysis over the eval rows: the breakdown by status, the dominant infeasibility
# reasons, and each subtask METHOD's failure rate. Returns a list of tibbles; write_report()
# renders them.
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

  # The rows that most likely indicate a BUG: trial + reason + funnel, ready to drill into.
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
  # This run's own universe + scheme -- what the surrogate learns from. n_other keeps the
  # global count for context.
  evals   <- filter_evals_to_universe(all_evals, settings$trial_universe) |>
               filter_evals_to_scheme(settings$optimize_scheme) |>
               filter_evals_to_build(settings$build %||% OPTIMIZER_BUILD)
  n_other <- nrow(all_evals) - nrow(evals)
  agg   <- aggregate_scores(evals)
  inc   <- incumbent_config(agg, settings$incumbent_min_reps)

  ok <- evals |> dplyr::filter(is.finite(score))

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
    # Whether the store is actually being backed up. It belongs here because the report is not
    # leader-gated, so it keeps updating while a long evaluation runs -- the one artefact that
    # can report a stalled backup as it happens. "Rows behind" is the direct measure; a backup
    # lands within a second of a row, so a lag that persists means backups are failing.
    local({
      bp <- settings$db_backup_path
      if (is.null(bp) || !nzchar(bp)) return(NULL)
      age    <- backup_age_minutes(bp)
      behind <- nrow(all_evals) - (.stored_rows(bp) %||% NA_integer_)
      paste0("- store backup: ",
             if (!is.finite(age)) "_never_"
             else if (age < 90) sprintf("%.0f min ago", age)
             else sprintf("%.1f h ago", age / 60),
             if (is.na(behind)) ""
             else if (behind <= 0) ", up to date"
             else sprintf(", %d row(s) behind", behind),
             if (isTRUE(is.finite(age) && !is.na(behind) && behind > 0 && age > 5))
               "  ** STALE **" else "")
    }),
    # The backup runs on every loop iteration, storing rows or not, so its age says the workers
    # are alive and nothing about whether they are producing. This is the one that says that.
    local({
      if (!nrow(all_evals)) return("- newest row: _none_")
      last <- suppressWarnings(max(all_evals$ts, na.rm = TRUE))
      mins <- suppressWarnings(as.numeric(difftime(Sys.time(),
                                 as.POSIXct(last, tz = "UTC"), units = "mins")))
      paste0("- newest row: ",
             if (!is.finite(mins)) last
             else if (mins < 90) sprintf("%.0f min ago", mins)
             else sprintf("%.1f h ago", mins / 60),
             if (isTRUE(is.finite(mins) && mins > 120)) "  ** NOTHING STORED **" else "")
    }),
    # Which estimator ranked the configs and, when the random-effects fit ran, the variance
    # components. sd_resid stays where lmer puts it -- at unit weight, i.e. n_test = 4, since
    # the fit uses n_test - 3 as inverse-variance weights. Rescaling it to a typical evaluation
    # is the Score precision block's job, which does it for two populations rather than one.
    #
    # A missing decomposition is stated, never left blank: no fit also means no `se`, and
    # .contenders() drops every config without one -- so a silent "pooled" here is the visible
    # end of contender replication having stopped.
    local({
      est  <- attr(agg, "estimator") %||% "pooled"
      vc   <- attr(agg, "var_comps")
      note <- attr(agg, "estimator_note")
      paste0("- config score estimator: ", est,
             if (!is.null(vc) && all(is.finite(vc)))
               paste0(sprintf("  (sd_trial %.3f, sd_config %.3f, sd_resid %.3f at unit weight)",
                              vc[["sd_trial"]], vc[["sd_config"]], vc[["sd_resid"]]),
                      if (!is.null(note)) paste0("  _fit warned: ", note, "_") else "")
             else paste0("  ** no variance components, and no `se` so contender replication is",
                         " idle: ", note %||% "reason not recorded", " **"))
    }),
    # How much is still undecided. 1 means the field is settled at this contender_z; every
    # contender being domain-covered means no further evidence about the leaders is obtainable.
    local({
      z    <- settings$contender_z %||% 1
      cand <- .contenders(agg, z, k = 8L)
      if (!length(cand)) return(NULL)
      # The listed set is capped at 8; the COUNT must not be, or a field of 280 reads as 8.
      n_all <- length(.contenders(agg, z, k = .Machine$integer.max))
      seen <- evals |> dplyr::filter(config_hash %in% cand) |>
        dplyr::group_by(config_hash) |>
        dplyr::summarise(n_trial = dplyr::n_distinct(trial_id), .groups = "drop")
      nu <- length(settings$trial_universe %||% character())
      paste0("- contenders: ", n_all, " within z = ", z,
             if (n_all > length(cand)) sprintf("  (listing the top %d)", length(cand)) else "",
             if (n_all == 1L) "  ** the field is settled at contender_z; raise it to continue **"
             else if (nu > 0 && all(seen$n_trial >= nu))
               sprintf(paste0("  ** each has been run on all %d trials in the domain",
                              " (some runs may fail); new configurations can",
                              " change the answer **"), nu)
             else "")
    }),
    # What replication still owes, and whether workers are colliding. `base` is the unrationed
    # config_replication floor; `contender` is the tier replicate_every rations, and it is empty
    # whenever the estimator above supplied no `se`. `duplicated cells` counts (config, trial,
    # scheme) triples with more than one row -- work done twice, which claim_eval() now
    # prevents, so the number should stop growing.
    local({
      bl <- tryCatch(.replication_backlog(evals, agg, settings, settings$trial_universe),
                     error = function(e) NULL)
      if (is.null(bl)) return(NULL)
      dup <- tryCatch(duplicate_cells(con), error = function(e) NA_integer_)
      n_claim <- tryCatch(nrow(active_claims(con, settings$optimize_scheme)),
                          error = function(e) NA_integer_)
      paste0("- replication backlog: ", length(bl$base), " base, ", length(bl$extra),
             " contender; ", n_claim, " in flight",
             if (is.finite(dup) && dup > 0) sprintf("; %d duplicated cell(s) in the store", dup)
             else "")
    }),
    paste0("- optimized scheme: ", settings$optimize_scheme),
    paste0("- evaluations: ", nrow(evals),
           " (", sum(is.finite(evals$score)), " scored, ",
           sum(!is.finite(evals$score)), " failed)",
           if (n_other > 0) paste0("; ", n_other, " more in the store are out of this scheme/domain") else ""),
    paste0("- distinct configurations: ", nrow(agg)),
    paste0("- best seed (submission) mean score: ",
           ifelse(is.finite(best_seed), signif(best_seed, 3), "NA")),
    "",
    # How much evidence the leaders rest on, against the field. Four columns that read
    # left-to-right as a chain: how many trials, how many accessions each, what that buys per
    # evaluation, and the configuration's own SE. The contender row should beat the all-configs
    # row on every one; if it does not, the leaders are lucky rather than well measured, and a
    # contender tier no better replicated than the base means contender replication has stalled.
    #
    # `se` here is on the Fisher-z scale, as .blup_scores() produces it, while mean_score is
    # back-transformed to the correlation scale. At the scores this project sees (r ~ 0.15) the
    # two differ by ~2%, which is why the report does not convert.
    local({
      vc <- attr(agg, "var_comps")
      if (is.null(vc) || !all(is.finite(vc)) || !nrow(agg)) return(NULL)
      cand <- .contenders(agg, settings$contender_z %||% 1, k = 8L)
      # The same weighting the fit applied, from the same function -- see .fisher_weights().
      w <- .fisher_weights(evals)
      row <- function(label, hashes) {
        a  <- if (is.null(hashes)) agg else dplyr::filter(agg, config_hash %in% hashes)
        wi <- if (is.null(hashes)) w else dplyr::filter(w, config_hash %in% hashes)
        wi <- dplyr::filter(wi, .usable)
        mw <- suppressWarnings(stats::median(wi$.w, na.rm = TRUE))
        if (!is.finite(mw) || mw < 1) mw <- 1
        nev  <- stats::median(a$n_ok, na.rm = TRUE)
        nois <- vc[["sd_resid"]] / sqrt(mw)          # per-evaluation residual sd
        sprintf("  %-22s %18s %17s %17s %14s %14s", label,
                fmt_num(nev, 0),
                fmt_num(stats::median(wi$.n_test, na.rm = TRUE), 0),
                fmt_num(nois, 3),
                fmt_num(nois / sqrt(max(1, nev)), 3),
                fmt_num(stats::median(a$se, na.rm = TRUE), 3))
      }
      c("## Score precision", "```",
        sprintf("  %-22s %18s %17s %17s %14s %14s", "",
                "med. evals/conf.", "med. acc./trial",
                "noise per eval.", "SE of the mean",
                "posterior SD"),
        row("all configurations", NULL),
        if (length(cand)) row(sprintf("top %d contenders", length(cand)), cand),
        "```", "")
    }),
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

    # The same rows, one column per contender, so the two tables read as one grid: which of
    # these methods do the leaders actually use? Built from `imp` itself rather than from
    # SUBTASKS, so the rows and their order cannot drift from the table above.
    z    <- settings$contender_z %||% 1
    cand <- .contenders(agg, z, k = 8L)
    if (length(cand)) {
      cc <- agg[match(cand, agg$config_hash), , drop = FALSE]
      # One column per contender: its method for each subtask, looked up per row of `imp`.
      # Every contender's methods are present in `imp` -- method_importance keeps finite-score
      # rows and a contender has at least one -- but mark rather than assume, since that
      # invariant lives in another function.
      cols <- lapply(seq_len(nrow(cc)), function(i) {
        cfg <- tryCatch(config_from_json(cc$config_json[i]), error = function(e) list())
        vapply(seq_len(nrow(imp)), function(r)
          if (identical(as.character(cfg[[paste0(imp$subtask[r], ".method")]] %||% ""),
                        imp$method[r])) "X" else ".", character(1))
      })
      # trimws on the right: %-3s pads the final column, and trailing spaces in a fenced
      # block are invisible noise in a diff of two reports.
      pad <- function(x) trimws(paste(sprintf("%-3s", x), collapse = " "), which = "right")
      hdr <- sprintf("  %-13s %-20s %s", "subtask", "method",
                     pad(paste0("C", seq_along(cand))))
      body <- vapply(seq_len(nrow(imp)), function(r)
        sprintf("  %-13s %-20s %s", imp$subtask[r], imp$method[r],
                pad(vapply(cols, `[`, character(1), r))),
        character(1))
      # Hash prefixes, not full hashes: eight 32-character ids cannot share a readable line.
      # `n a/b` is usable evaluations over trials ATTEMPTED: the domain bullet above counts
      # attempts (a trial run that failed is still a trial used up), while n_ok counts the rows
      # the estimate rests on. They differ whenever a configuration failed on some trial, and
      # printing only one of them made the two look contradictory.
      att <- evals |> dplyr::filter(config_hash %in% cand) |>
        dplyr::group_by(config_hash) |>
        dplyr::summarise(n_try = dplyr::n_distinct(trial_id), .groups = "drop")
      ntry <- att$n_try[match(cc$config_hash, att$config_hash)]
      vc_l <- attr(agg, "var_comps")
      # se_fixed = sd_resid / sqrt(n): the standard error of the mean, which is what "how well
      # do we know this" means. The shrunk se beside it is the ranking device .contenders() uses.
      wl   <- attr(agg, "median_weight") %||% NA_real_
      if (!is.finite(wl)) wl <- 1
      se_fix <- if (!is.null(vc_l) && all(is.finite(vc_l)))
        vc_l[["sd_resid"]] / sqrt(wl) / sqrt(pmax(1, cc$n_ok)) else rep(NA_real_, nrow(cc))
      leg <- sprintf("  C%-2d %s  mean %+.3f  se %s (shrunk %.3f)  n %d/%s",
                     seq_len(nrow(cc)), substr(cc$config_hash, 1, 8), cc$mean_score,
                     vapply(se_fix, fmt_num, character(1), 3), cc$se,
                     cc$n_ok, ifelse(is.na(ntry), "?", ntry))
      lines <- c(lines, "",
                 sprintf("## Contender methods (top %d by UCB at contender_z = %g)",
                         length(cand), z),
                 "```", hdr, body, "",
                 "  n = usable evaluations / trials attempted",
                 leg, "```")

      # --- within-configuration spread, non-pooled -------------------------
      # What the pooled variance components cannot say: is one contender more trial-sensitive
      # than another? Only for contenders with enough evaluations for an sd to mean anything.
      wcfg <- attr(agg, "within_config")
      if (!is.null(wcfg)) {
        wr <- wcfg[match(cc$config_hash, wcfg$config_hash), , drop = FALSE]
        wr$label <- paste0("C", seq_len(nrow(cc)))
        keep <- !is.na(wr$n_eval) & wr$n_eval >= 6L
        if (any(keep)) {
          wk <- wr[keep, , drop = FALSE]
          wtbl <- sprintf("  %-6s %6s %18s %24s %16s", wk$label, wk$n_eval,
                          vapply(wk$sd_raw, fmt_num, character(1), 3),
                          vapply(wk$sd_adj, fmt_num, character(1), 3),
                          vapply(wk$sd_adj / sqrt(wk$n_eval), fmt_num, character(1), 3))
          lines <- c(lines, "",
                     "## Within-configuration spread among evaluations",
                     "```",
                     sprintf("  %-6s %6s %18s %24s %16s", "", "evals", "sd across trials",
                             "sd net of trial effect", "SE of the mean"),
                     wtbl,
                     if (sum(!keep)) sprintf("  (%d contender(s) with fewer than 6 evaluations omitted)",
                                             sum(!keep)) else NULL,
                     "```")
        } else {
          lines <- c(lines, "",
                     "## Within-configuration spread among evaluations",
                     "```",
                     "  no contender has 6 or more evaluations yet",
                     "```")
        }
      }
    }
  }

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

    # Suspected bugs get their own section: these failed with a funnel cliff, so investigate
    # them with diagnose_trial() before trusting the "infeasible" verdict.
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

  # ATOMIC, because every worker writes this file: concurrent writeLines() to one path would
  # interleave and a reader could catch a half-written file. Render to a sibling temp and rename.
  tmp <- paste0(settings$report_path, ".tmp", Sys.getpid())
  writeLines(lines, tmp)
  if (!isTRUE(file.rename(tmp, settings$report_path))) {
    unlink(tmp)
    message("report write -> ", settings$report_path, " FAILED")
  }
  invisible(settings$report_path)
}
