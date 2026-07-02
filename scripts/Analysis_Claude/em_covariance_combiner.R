# em_covariance_combiner.R
#
# Wishart-EM covariance combiner (Akdemir), trimmed to the functions exercised
# when compute_sampling_cov = FALSE (the sampling-covariance routines are
# omitted). Sourced by both build_grm_for_cv00.R (global combine) and
# build_per_trial_grms.R (per-focal-trial combines) so the implementation -
# including the incremental, low-memory M-step - lives in one place.
#
# Index convention (matches the original genomic_analysis.Rmd): var_indices are
# the 1-based positions of each partial matrix's accessions within the combined
# accession ordering (i.e. match(colnames, combined_names)). n_vars then becomes
# max(index) + 1, which leaves a phantom variable at psi position 1 that is
# always "missing"; callers drop it with result$psi[-1, -1].

initialize_psi <- function(partial_covs, var_indices, n_vars) {
  psi <- diag(n_vars)
  for (i in 0:(n_vars - 1)) {
    variances <- numeric(0)
    for (j in seq_along(partial_covs)) {
      cov <- partial_covs[[j]]
      idx <- var_indices[[j]]
      if (i %in% idx) {
        i_local <- match(i, idx) - 1
        variances <- c(variances, cov[i_local + 1, i_local + 1])
      }
    }
    if (length(variances) > 0) psi[i + 1, i + 1] <- mean(variances)
  }
  psi
}

compute_conditional_expectation <- function(ya, psi, obs_idx, missing_idx) {
  obs_idx_r <- obs_idx + 1
  missing_idx_r <- missing_idx + 1
  n_total <- nrow(psi)
  result <- matrix(0, nrow = n_total, ncol = n_total)
  if (length(missing_idx) == 0) {
    result[obs_idx_r, obs_idx_r] <- ya
    return(result)
  }
  psi_aa <- psi[obs_idx_r, obs_idx_r, drop = FALSE]
  psi_aa <- psi_aa + diag(1e-5, nrow(psi_aa))
  psi_ab <- psi[obs_idx_r, missing_idx_r, drop = FALSE]
  psi_bb <- psi[missing_idx_r, missing_idx_r, drop = FALSE]
  psi_aa_inv <- solve(psi_aa)
  b_matrix <- t(psi_ab) %*% psi_aa_inv
  result[obs_idx_r, obs_idx_r] <- ya
  exp_yab <- b_matrix %*% ya
  result[missing_idx_r, obs_idx_r] <- exp_yab
  result[obs_idx_r, missing_idx_r] <- t(exp_yab)
  result[missing_idx_r, missing_idx_r] <-
    psi_bb - b_matrix %*% psi_ab + b_matrix %*% ya %*% t(b_matrix)
  result
}

compute_log_likelihood <- function(psi, partial_covs, var_indices, degrees_freedom) {
  log_lik <- 0.0
  for (i in seq_along(partial_covs)) {
    ya <- partial_covs[[i]]
    idx_r <- var_indices[[i]] + 1
    nu <- degrees_freedom[i]
    psi_subset <- psi[idx_r, idx_r, drop = FALSE] + diag(1e-6, length(idx_r))
    eig <- try(eigen(psi_subset, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
    if (inherits(eig, "try-error") || any(eig <= 0)) return(-Inf)
    logdet_psi <- sum(log(eig))
    trace_term <- sum(diag(solve(psi_subset) %*% ya))
    log_lik <- log_lik - 0.5 * nu * (logdet_psi + trace_term)
  }
  log_lik
}

EMCovarianceCombiner <- function(partial_covs, var_indices, degrees_freedom = NULL,
                                 max_iter = 100, tol = 1e-6, track_loglik = TRUE) {
  if (length(partial_covs) == 0 || length(var_indices) == 0) stop("Empty input provided")
  if (length(partial_covs) != length(var_indices)) stop("Mismatched lengths")
  n_vars <- max(sapply(var_indices, function(idx) max(idx))) + 1
  if (is.null(degrees_freedom)) degrees_freedom <- rep(100, length(partial_covs))
  loglik_path <- if (track_loglik) numeric(0) else NULL
  psi <- initialize_psi(partial_covs, var_indices, n_vars)
  if (track_loglik) {
    loglik_path <- c(loglik_path,
                     compute_log_likelihood(psi, partial_covs, var_indices, degrees_freedom))
  }
  pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
  for (iter in 1:max_iter) {
    setTxtProgressBar(pb, iter)
    psi_old <- psi
    # Accumulate the E-step expectations directly into psi_sum instead of first
    # storing all length(partial_covs) full n_vars x n_vars matrices in a list.
    # That list was the peak-memory hot spot (~n_projects x n_vars^2); summing
    # incrementally keeps at most one expectation matrix live at a time and is
    # mathematically identical.
    psi_sum <- matrix(0, nrow = n_vars, ncol = n_vars)
    total_nu <- sum(degrees_freedom)
    for (i in seq_along(partial_covs)) {
      idx <- var_indices[[i]]
      missing_idx <- setdiff(0:(n_vars - 1), idx)
      exp_y <- compute_conditional_expectation(partial_covs[[i]], psi, idx, missing_idx)
      psi_sum <- psi_sum + degrees_freedom[i] * exp_y
    }
    psi <- psi_sum / total_nu
    if (track_loglik) {
      loglik_path <- c(loglik_path,
                       compute_log_likelihood(psi, partial_covs, var_indices, degrees_freedom))
    }
    if (max(abs(psi - psi_old)) < tol) {
      message("\nConverged at iteration ", iter)
      break
    }
  }
  close(pb)
  list(psi = psi, sampling_cov = NULL, loglik_path = loglik_path)
}
