# ==============================================================================
# NUTS diagnostics
# ==============================================================================
#
# Post-sampling diagnostic checks for the No-U-Turn Sampler (NUTS).
# Covers energy-based warmup assessment and per-chain summary
# statistics (divergences, tree depth, E-BFMI).
# ==============================================================================


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
      var_ratio = rep(NA_real_, nchains)
    ))
  }

  mid = floor(n / 2)

  results = lapply(seq_len(nchains), function(chain) {
    energy = energy_mat[chain, ]
    first_half = energy[1:mid]
    second_half = energy[(mid + 1):n]

    # Linear trend in energy
    time_idx = seq_len(n)
    trend_lm = stats::lm(energy ~ time_idx)
    slope = stats::coef(trend_lm)[2]
    slope_se = summary(trend_lm)$coefficients[2, 2]
    slope_significant = abs(slope / slope_se) > 2.58

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
      var_ratio = var_ratio
    )
  })

  list(
    warmup_incomplete = sapply(results, `[[`, "warmup_incomplete"),
    energy_slope = sapply(results, `[[`, "energy_slope"),
    slope_significant = sapply(results, `[[`, "slope_significant"),
    ebfmi_first_half = sapply(results, `[[`, "ebfmi_first_half"),
    ebfmi_second_half = sapply(results, `[[`, "ebfmi_second_half"),
    var_ratio = sapply(results, `[[`, "var_ratio")
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

  if(verbose) {
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
        "E-BFMI: %.3f in chain%s %s - see vignette('diagnostics') for guidance",
        min_ebfmi,
        if(length(low_ebfmi_chains) > 1) "s" else "",
        paste(low_ebfmi_chains, collapse = ", ")
      ))
    }

    incomplete_chains = which(warmup_check$warmup_incomplete)
    if(length(incomplete_chains) > 0) {
      issues = c(issues, sprintf(
        "Warmup incomplete: energy not stationary in chain%s %s - increase warmup",
        if(length(incomplete_chains) > 1) "s" else "",
        paste(incomplete_chains, collapse = ", ")
      ))
    }

    if(length(issues) > 0 && isTRUE(getOption("bgms.verbose", TRUE))) {
      cat("NUTS issues:\n")
      for(issue in issues) {
        cat("  -", issue, "\n")
      }
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
    summary = list(
      total_divergences = total_divergences,
      max_tree_depth_hits = max_tree_depth_hits,
      min_ebfmi = min_ebfmi,
      mean_accept_prob = mean_accept_prob,
      warmup_incomplete = any(warmup_check$warmup_incomplete)
    )
  ))
}
