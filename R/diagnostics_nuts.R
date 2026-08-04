# ==============================================================================
# NUTS diagnostics
# ==============================================================================
#
# Post-sampling diagnostic checks for the No-U-Turn Sampler (NUTS).
# Covers energy-based warmup assessment and per-chain summary
# statistics (divergences, tree depth, E-BFMI).
# ==============================================================================


# ------------------------------------------------------------------------------
# integrated_act
# ------------------------------------------------------------------------------
# Integrated autocorrelation time, tau = n / ESS, via the package's own ESS
# (AR spectral density, coda compatible).
#
# The OLS standard error of the energy trend assumes independent draws. An MCMC
# energy trace is autocorrelated, so the naive |t| is inflated by sqrt(tau).
# This is what makes the trend criterion fire far more often than its nominal
# level suggests; see check_warmup_complete() for how it is reported.
#
# @param x  Numeric vector.
#
# Returns: Numeric scalar >= 1, or NA_real_ if ESS is not computable.
# ------------------------------------------------------------------------------
integrated_act = function(x) {
  n = length(x)
  if(n < 20) {
    return(1)
  }
  ess = .compute_ess_cpp(array(x, dim = c(n, 1L, 1L)))[1]
  if(!is.finite(ess) || ess <= 0) {
    return(NA_real_)
  }
  max(1, n / ess)
}


# ------------------------------------------------------------------------------
# check_warmup_complete
# ------------------------------------------------------------------------------
# Assess whether warmup is complete using energy stationarity.
#
# Splits each chain's energy trace at the midpoint and checks three
# criteria: (1) significant linear trend in energy, (2) first-half
# E-BFMI below 0.3, and (3) first-half / second-half variance ratio
# above 2.0. Any triggered criterion flags the chain.
#
# @param energy_mat  Numeric matrix (chains x iterations) of energy
#   values from the NUTS sampler.
#
# Returns: A named list with one element per diagnostic, each a
#   vector of length nchains:
#   - warmup_incomplete:  Logical, TRUE if any criterion triggered.
#   - energy_slope:       Numeric, OLS slope of energy on iteration.
#   - slope_significant:  Logical, |t| > 2.58 for the slope.
#   - ebfmi_first_half:   Numeric, E-BFMI for the first half.
#   - ebfmi_second_half:  Numeric, E-BFMI for the second half.
#   - var_ratio:          Numeric, var(first_half) / var(second_half).
#   - energy_tau:         Numeric, integrated autocorrelation time of the trend
#       residuals, estimated on the SECOND half. A residual transient lives in
#       the first half and would inflate tau there, so the second half gives the
#       stationary value.
#   - slope_t_corrected:  Numeric, the slope t-statistic with its standard error
#       widened by sqrt(energy_tau).
#
# energy_tau and slope_t_corrected are reported only; they do not enter
# warmup_incomplete. The flag stays deliberately sensitive because it exists to
# send the user to R-hat and ESS, and nothing else in the package does. They let
# a user who investigates separate a real drift from an autocorrelation artifact:
# a corrected |t| well under 2.58 alongside tau around 3 means the trend is an
# artifact of autocorrelation, not residual warmup.
# ------------------------------------------------------------------------------
check_warmup_complete = function(energy_mat) {
  nchains = nrow(energy_mat)
  n = ncol(energy_mat)

  if(n < 20) {
    return(list(
      warmup_incomplete = rep(FALSE, nchains),
      energy_slope = rep(NA_real_, nchains),
      slope_significant = rep(FALSE, nchains),
      ebfmi_first_half = rep(NA_real_, nchains),
      ebfmi_second_half = rep(NA_real_, nchains),
      var_ratio = rep(NA_real_, nchains),
      energy_tau = rep(NA_real_, nchains),
      slope_t_corrected = rep(NA_real_, nchains)
    ))
  }

  # Per-chain NA result for chains with too little usable energy to assess.
  na_result = list(
    warmup_incomplete = FALSE,
    energy_slope = NA_real_,
    slope_significant = FALSE,
    ebfmi_first_half = NA_real_,
    ebfmi_second_half = NA_real_,
    var_ratio = NA_real_,
    energy_tau = NA_real_,
    slope_t_corrected = NA_real_
  )

  results = lapply(seq_len(nchains), function(chain) {
    # An interrupted or degenerate run leaves the energy trace partly or
    # fully NA; assess the finite draws only.
    energy = energy_mat[chain, ]
    energy = energy[is.finite(energy)]
    n_chain = length(energy)

    # Too few usable values to split and regress: return NA diagnostics
    # rather than fitting lm/var on empty data.
    if(n_chain < 20) {
      return(na_result)
    }

    mid = floor(n_chain / 2)
    first_half = energy[1:mid]
    second_half = energy[(mid + 1):n_chain]

    # Linear trend in energy
    time_idx = seq_len(n_chain)
    trend_lm = stats::lm(energy ~ time_idx)
    slope = stats::coef(trend_lm)[2]
    slope_se = summary(trend_lm)$coefficients[2, 2]
    slope_significant = abs(slope / slope_se) > 2.58

    # Reported alongside the flag, not part of it: the same statistic with its
    # standard error corrected for autocorrelation. tau comes from the second
    # half of the residuals, which is stationary even when the first half is not.
    trend_resid = stats::residuals(trend_lm)
    energy_tau = integrated_act(trend_resid[(mid + 1):n_chain])
    slope_t_corrected = if(is.finite(energy_tau)) {
      unname(slope / (slope_se * sqrt(energy_tau)))
    } else {
      NA_real_
    }

    # E-BFMI per half
    ebfmi_first = mean(diff(first_half)^2) / stats::var(first_half)
    ebfmi_second = mean(diff(second_half)^2) / stats::var(second_half)

    # Variance ratio
    var_ratio = stats::var(first_half) / stats::var(second_half)

    # Flag if any criterion triggered
    warmup_incomplete = slope_significant || ebfmi_first < 0.3 || var_ratio > 2.0

    list(
      warmup_incomplete = warmup_incomplete,
      energy_slope = slope,
      slope_significant = slope_significant,
      ebfmi_first_half = ebfmi_first,
      ebfmi_second_half = ebfmi_second,
      var_ratio = var_ratio,
      energy_tau = energy_tau,
      slope_t_corrected = slope_t_corrected
    )
  })

  list(
    warmup_incomplete = sapply(results, `[[`, "warmup_incomplete"),
    energy_slope = sapply(results, `[[`, "energy_slope"),
    slope_significant = sapply(results, `[[`, "slope_significant"),
    ebfmi_first_half = sapply(results, `[[`, "ebfmi_first_half"),
    ebfmi_second_half = sapply(results, `[[`, "ebfmi_second_half"),
    var_ratio = sapply(results, `[[`, "var_ratio"),
    energy_tau = sapply(results, `[[`, "energy_tau"),
    slope_t_corrected = sapply(results, `[[`, "slope_t_corrected")
  )
}

# ------------------------------------------------------------------------------
# summarize_nuts_diagnostics
# ------------------------------------------------------------------------------
# Combine and summarize NUTS diagnostics across chains.
#
# Extracts treedepth, divergence, and energy traces from a list of
# chain outputs. Computes per-chain E-BFMI, runs the warmup check,
# and optionally prints a human-readable issues summary.
#
# @param out  List of chain outputs. Each element is a named list
#   that must contain "treedepth__", "divergent__", and "energy__".
#   Chains without these fields are silently dropped.
# @param nuts_max_depth  Integer scalar: the maximum tree depth used
#   during sampling. Iterations that reached this depth are counted
#   as tree-depth hits (default: 10).
# @param verbose  Logical scalar: if TRUE (the default), print a
#   summary of any detected issues to the console.
#
# Returns: An invisible named list with:
#   - treedepth:  Integer matrix (chains x iterations).
#   - divergent:  Integer matrix (chains x iterations), 0/1.
#   - energy:     Numeric matrix (chains x iterations).
#   - accept_prob: Numeric matrix (chains x iterations) of mean
#       per-trajectory Metropolis acceptance (Stan's accept_stat__).
#   - ebfmi:      Numeric vector of per-chain E-BFMI values.
#   - warmup_check: Output of check_warmup_complete().
#   - summary:    List with total_divergences, max_tree_depth_hits,
#       min_ebfmi, mean_accept_prob, and warmup_incomplete (logical).
# ------------------------------------------------------------------------------
summarize_nuts_diagnostics = function(out, nuts_max_depth = 10, verbose = TRUE) {
  nuts_chains = Filter(function(chain) {
    all(c("treedepth__", "divergent__", "energy__") %in% names(chain))
  }, out)

  if(length(nuts_chains) == 0) {
    stop("No NUTS diagnostics found in output.")
  }

  # Combine fields into matrices (chains x iterations). Count fields
  # (treedepth, divergence flags) are integer per the return contract;
  # energy and acceptance probabilities are real-valued.
  combine_diag = function(field, integer = FALSE) {
    coerce = if(integer) as.integer else as.numeric
    do.call(rbind, lapply(nuts_chains, function(chain) coerce(chain[[field]])))
  }

  treedepth_mat = combine_diag("treedepth__", integer = TRUE)
  divergent_mat = combine_diag("divergent__", integer = TRUE)
  energy_mat = combine_diag("energy__")

  accept_prob_mat = if("accept_prob__" %in% names(nuts_chains[[1]])) {
    combine_diag("accept_prob__")
  } else {
    matrix(NA_real_, nrow = nrow(divergent_mat), ncol = ncol(divergent_mat))
  }

  # E-BFMI per chain
  compute_ebfmi = function(energy) {
    mean(diff(energy)^2) / stats::var(energy)
  }
  ebfmi_per_chain = apply(energy_mat, 1, compute_ebfmi)

  warmup_check = check_warmup_complete(energy_mat)

  # Summaries
  n_total = nrow(divergent_mat) * ncol(divergent_mat)
  total_divergences = sum(divergent_mat)
  max_tree_depth_hits = sum(treedepth_mat == nuts_max_depth)
  min_ebfmi = min(ebfmi_per_chain)
  low_ebfmi_chains = which(ebfmi_per_chain < 0.2)

  divergence_rate = total_divergences / n_total
  depth_hit_rate = max_tree_depth_hits / n_total

  # Build the issue list regardless of verbose so has_issues is always known
  # (verbose only governs whether the block is printed). The vignette pointer
  # is emitted once by the output builder as a shared footer, not per issue.
  issues = character(0)

  if(total_divergences > 0) {
    if(divergence_rate > 0.001) {
      issues = c(issues, sprintf(
        "Divergences: %d (%.2f%%) - increase target acceptance or use adaptive-metropolis",
        total_divergences, 100 * divergence_rate
      ))
    } else {
      issues = c(issues, sprintf(
        "Divergences: %d (%.3f%%) - check R-hat and ESS",
        total_divergences, 100 * divergence_rate
      ))
    }
  }

  if(max_tree_depth_hits > 0) {
    if(depth_hit_rate > 0.01) {
      issues = c(issues, sprintf(
        "Tree depth: %d hits (%.1f%%) - consider max_depth > %d",
        max_tree_depth_hits, 100 * depth_hit_rate, nuts_max_depth
      ))
    } else {
      issues = c(issues, sprintf(
        "Tree depth: %d hits (%.2f%%) - check ESS",
        max_tree_depth_hits, 100 * depth_hit_rate
      ))
    }
  }

  if(length(low_ebfmi_chains) > 0) {
    issues = c(issues, sprintf(
      "E-BFMI: %.3f in chain%s %s",
      min_ebfmi,
      if(length(low_ebfmi_chains) > 1) "s" else "",
      paste(low_ebfmi_chains, collapse = ", ")
    ))
  }

  incomplete_chains = which(warmup_check$warmup_incomplete)
  if(length(incomplete_chains) > 0) {
    # Name the criterion that actually fired. A low first-half E-BFMI is a
    # statement about how the sampler moves through the energy landscape, not
    # about residual warmup, so reporting it as "energy not stationary" would
    # misdescribe it.
    trigger = if(any(warmup_check$ebfmi_first_half[incomplete_chains] < 0.3,
      na.rm = TRUE
    )) {
      "low first-half E-BFMI"
    } else {
      "energy not stationary"
    }
    issues = c(issues, sprintf(
      "Warmup may be incomplete: %s in chain%s %s - check R-hat and ESS",
      trigger,
      if(length(incomplete_chains) > 1) "s" else "",
      paste(incomplete_chains, collapse = ", ")
    ))
  }

  if(verbose && length(issues) > 0 && isTRUE(getOption("bgms.verbose", TRUE))) {
    cat("NUTS issues:\n")
    for(issue in issues) {
      cat("  -", issue, "\n")
    }
  }

  mean_accept_prob = mean(accept_prob_mat, na.rm = TRUE)

  invisible(list(
    treedepth = treedepth_mat,
    divergent = divergent_mat,
    energy = energy_mat,
    accept_prob = accept_prob_mat,
    ebfmi = ebfmi_per_chain,
    warmup_check = warmup_check,
    has_issues = length(issues) > 0,
    summary = list(
      total_divergences = total_divergences,
      max_tree_depth_hits = max_tree_depth_hits,
      min_ebfmi = min_ebfmi,
      mean_accept_prob = mean_accept_prob,
      warmup_incomplete = any(warmup_check$warmup_incomplete)
    )
  ))
}
