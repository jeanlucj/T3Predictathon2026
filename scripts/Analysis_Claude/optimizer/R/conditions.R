# conditions.R
#
# Three kinds of "this didn't work", distinguished so the loop knows whether to continue:
#
#   * infeasible(code)  this ONE (trial, config, scheme) cannot be evaluated -- too little
#       genotyped overlap, no training trials, too few markers. EXPECTED and common.
#       evaluate.R records it and the loop moves on with a fresh trial and configuration.
#   * fatal(message)    the RUN cannot proceed: a settings problem, or a method in the config
#       space the pipeline does not implement. Re-raised past evaluate.R so the loop halts
#       instead of re-hitting the same wall.
#   * sample_failed()   no usable focal trial after the retries. Not attributable to a
#       configuration, so it is not written to the per-config failure log.
#
# `code` is a short machine-readable token, so failure_summary() can group on it.

# Raise: this (trial, config, scheme) is not evaluable. `detail` adds context; `code` is the
# groupable identity. `funnel` is the data cardinality at each pipeline stage, stored with the
# failure so a cliff -- many in, ~zero out -- is visible afterwards. `suspect = TRUE` marks a
# funnel that looks like a BUG rather than genuine infeasibility; evaluate.R records it under
# its own status.
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
