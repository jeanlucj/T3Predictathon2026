# conditions.R
#
# Two kinds of "this didn't work", deliberately distinguished so the optimizer
# knows whether to keep going or stop:
#
#   * infeasible(code)  -- this ONE (trial, config, scheme) cannot be evaluated:
#       too little genotyped overlap, no training trials, too few markers, etc.
#       This is EXPECTED. Many random focal trials are simply not predictable
#       with a given (often demanding) method configuration. evaluate.R catches
#       it, records a failed eval (status="infeasible", reason=code) in the
#       store -- the trial x method-config failure log -- and the loop carries on
#       with a freshly sampled trial AND a freshly chosen configuration.
#
#   * fatal(message)    -- the whole RUN cannot proceed: a settings problem (no
#       trials match the target domain), or a programming inconsistency (a method
#       in the config space that the pipeline does not implement). Re-raised past
#       evaluate.R so run_optimizer.R halts cleanly instead of spinning forever
#       re-hitting the same wall.
#
#   * sample_failed()   -- could not sample ANY usable focal trial after the
#       retries. Not attributable to a configuration, so it is NOT written to the
#       per-config failure log; the loop tolerates a few in a row (transient
#       network) and escalates to a clean stop only if they persist.
#
# The `code` on an infeasible condition is a short machine-readable token (e.g.
# "insufficient_geno_overlap") so failure_summary() can group on it and reveal
# which method choices break most often.

# Raise: this (trial, config, scheme) is not evaluable. `detail` adds context to
# the logged message but the groupable identity is `code`. `funnel` is a named
# numeric vector of the data cardinalities at each pipeline stage (n accessions,
# n genotyped, n overlap, ...), stored with the failure so a cliff -- many in,
# ~zero out -- is visible after the fact. `suspect = TRUE` marks a failure whose
# funnel looks like a BUG (data that should be visible isn't) rather than genuine
# infeasibility; evaluate.R records it under status "suspect" for review.
infeasible <- function(code, detail = NULL, funnel = NULL, suspect = FALSE) {
  msg <- if (is.null(detail)) code else paste0(code, ": ", detail)
  stop(structure(
    class = c("optimizer_infeasible", "error", "condition"),
    list(message = msg, code = code, detail = detail,
         funnel = funnel, suspect = isTRUE(suspect), call = sys.call(-1))))
}

# Encode a funnel (named numeric) as a compact "k=v; k=v" string for the store.
funnel_string <- function(funnel) {
  if (is.null(funnel) || !length(funnel)) return(NA_character_)
  paste(sprintf("%s=%s", names(funnel), funnel), collapse = "; ")
}

# Raise: the run cannot proceed. Halts the optimizer.
fatal <- function(message, code = "fatal") {
  stop(structure(
    class = c("optimizer_fatal", "error", "condition"),
    list(message = message, call = sys.call(-1), code = code)))
}

# Raise: no usable focal trial could be sampled (trial-sampling failure).
sample_failed <- function(detail = NULL) {
  msg <- paste0("no usable trial sampled",
                if (!is.null(detail)) paste0(": ", detail) else "")
  stop(structure(
    class = c("optimizer_sample_failed", "error", "condition"),
    list(message = msg, call = sys.call(-1), code = "no_usable_trial_sampled")))
}
