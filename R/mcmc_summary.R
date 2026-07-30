# Summary utilities for spike-and-slab MCMC output


# ------------------------------------------------------------------
# ensure_summaries
# ------------------------------------------------------------------
# Lazily computes MCMC diagnostics (ESS, Rhat, MCSE) and stores
# results in fit$cache (an environment with reference semantics).
# On first call the summaries are computed from raw chain samples;
# subsequent calls return immediately.
#
# @param fit  A bgms or bgmCompare object with a $cache environment.
#
# Returns: invisible(NULL). Results are stored in fit$cache.
# ------------------------------------------------------------------
ensure_summaries = function(fit) {
  cache = get_fit_cache(fit)
  if(is.null(cache)) {
    return(invisible(NULL))
  }
  if(isTRUE(cache$summaries_computed)) {
    return(invisible(NULL))
  }

  raw = cache$raw
  edge_selection = cache$edge_selection
  names_main = cache$names_main
  edge_names = cache$edge_names
  is_continuous = cache$is_continuous
  model_type = cache$model_type

  if(identical(model_type, "compare")) {
    names_all = cache$names_all
    summary_list = summarize_fit_compare(
      fit = raw,
      main_effect_indices = cache$main_effect_indices,
      pairwise_effect_indices = cache$pairwise_effect_indices,
      num_variables = cache$num_variables,
      num_groups = cache$num_groups,
      difference_selection = cache$difference_selection,
      param_names_main = names_all$main_baseline,
      param_names_pairwise = names_all$pairwise_baseline,
      param_names_main_diff = names_all$main_diff,
      param_names_pairwise_diff = names_all$pairwise_diff,
      param_names_indicators = names_all$indicators
    )

    cache$posterior_summary_main_baseline = summary_list$main_baseline
    cache$posterior_summary_pairwise_baseline = summary_list$pairwise_baseline
    cache$posterior_summary_main_differences = summary_list$main_differences
    cache$posterior_summary_pairwise_differences = summary_list$pairwise_differences
    if(!is.null(raw[[1]][["rb_inclusion_samples"]])) {
      # Report the Rao-Blackwellized inclusion probability with continuous
      # ESS/Rhat (the default estimator for the printed summary).
      cache$posterior_summary_indicator = summarize_rb_inclusion(
        raw, names_all$indicators,
        keep_parameter_col = TRUE
      )
    } else {
      cache$posterior_summary_indicator = summary_list$indicators
    }
  } else {
    if(isTRUE(is_continuous)) {
      # GGM raw pairwise samples are on the precision (off-diagonal) scale;
      # the association scale is precision * -0.5 (as in coef() /
      # extract_pairwise_interactions()). Transform the draws before
      # summarizing so the whole pipeline -- including the selection-aware
      # mixture summary -- reports mean, sd, and mcse on the association
      # scale. Zeros stay zeros, so the spike/slab split is unaffected.
      raw = lapply(raw, function(chain) {
        chain$pairwise_samples = -0.5 * chain$pairwise_samples
        chain
      })
    }
    summary_list = summarize_fit(raw, edge_selection = edge_selection)
    main_summary = summary_list$main[, -1]
    pairwise_summary = summary_list$pairwise[, -1]

    rownames(main_summary) = names_main
    rownames(pairwise_summary) = edge_names

    if(identical(model_type, "mixed_mrf")) {
      n_main = cache$n_main
      n_quad = cache$n_quad
      main_rows = seq_len(n_main)
      quad_rows = n_main + seq_len(n_quad)
      cache$posterior_summary_main = main_summary[main_rows, , drop = FALSE]
      # Recompute quadratic summary on the residual variance scale:
      # raw samples store negative association diagonal; transform to
      # residual variance = -1 / (2 * diag).
      array3d_main = combine_chains(raw, "main_samples")
      array3d_rv = -1 / (2 * array3d_main[, , quad_rows, drop = FALSE])
      rv_summary = summarize_manual(raw, array3d = array3d_rv)[, -1]
      rownames(rv_summary) = sub(
        " \\(precision diag\\)$", " (residual variance)",
        names_main[quad_rows]
      )
      cache$posterior_summary_quadratic = rv_summary
    } else if(isTRUE(is_continuous)) {
      cache$posterior_summary_main = NULL
      # Recompute quadratic summary on the residual variance scale:
      # raw samples store precision diagonal; transform to
      # residual variance = 1 / precision.
      array3d_main = combine_chains(raw, "main_samples")
      array3d_rv = 1 / array3d_main
      rv_summary = summarize_manual(raw, array3d = array3d_rv)[, -1]
      rownames(rv_summary) = sub(
        " \\(precision\\)$", " (residual variance)",
        names_main
      )
      cache$posterior_summary_quadratic = rv_summary
    } else {
      cache$posterior_summary_main = main_summary
    }
    cache$posterior_summary_pairwise = pairwise_summary

    if(edge_selection) {
      if(!is.null(raw[[1]][["rb_inclusion_samples"]])) {
        # Report the Rao-Blackwellized inclusion probability with continuous
        # ESS/Rhat (the default estimator for the printed summary).
        cache$posterior_summary_indicator = summarize_rb_inclusion(raw, edge_names)
      } else {
        indicator_summary = summary_list$indicator[, -1]
        rownames(indicator_summary) = edge_names
        cache$posterior_summary_indicator = indicator_summary
      }
    }
  }

  cache$summaries_computed = TRUE
  invisible(NULL)
}


# Combine MCMC chains into a 3D array [niter x nchains x nparam]
combine_chains = function(fit, component) {
  nchains = length(fit)
  samples_list = lapply(fit, function(x) x[[component]])
  niter = nrow(samples_list[[1]])
  nparam = ncol(samples_list[[1]])
  array3d = array(NA_real_, dim = c(niter, nchains, nparam))
  for(i in seq_len(nchains)) {
    array3d[, i, ] = samples_list[[i]]
  }
  array3d
}

# Split each chain into its first and second half, turning m chains into 2m
# sub-chains. Rhat computed on the split array is the split-R-hat of Gelman et
# al. / Vehtari et al.: it detects within-chain drift that whole-chain
# Gelman-Rubin cannot, and the vignette's Rhat < 1.01 guideline refers to this
# split form. An odd iteration count drops the middle draw so both halves match;
# chains shorter than two iterations are returned unchanged.
split_chains = function(array3d) {
  niter = dim(array3d)[1]
  nchains = dim(array3d)[2]
  nparam = dim(array3d)[3]
  half = niter %/% 2L
  if(half < 1L) {
    return(array3d)
  }
  first_idx = seq_len(half)
  second_idx = (niter - half + 1L):niter
  out = array(NA_real_, dim = c(half, 2L * nchains, nparam))
  for(i in seq_len(nchains)) {
    out[, 2L * i - 1L, ] = array3d[first_idx, i, ]
    out[, 2L * i, ] = array3d[second_idx, i, ]
  }
  out
}

# Compute ESS and Rhat for a single [niter x nchains] draws matrix.
# Used only by summarize_slab() where the draws are variable-length.
compute_rhat_ess = function(draws) {
  if(!is.matrix(draws)) draws = matrix(draws, ncol = 1)
  arr = array(draws, dim = c(nrow(draws), ncol(draws), 1L))
  ess = .compute_ess_cpp(arr)[1]
  rhat = .compute_rhat_cpp(arr)[1]
  list(ess = ess, rhat = rhat)
}

# Basic summarizer for continuous parameters
summarize_manual = function(fit, component = c("main_samples", "pairwise_samples"), param_names = NULL, array3d = NULL) {
  component = match.arg(component) # Add options later
  if(is.null(array3d)) array3d = combine_chains(fit, component)
  nparam = dim(array3d)[3]

  # Batch computation via C++
  ess = .compute_ess_cpp(array3d)
  rhat = .compute_rhat_cpp(split_chains(array3d))

  # Vectorized mean and sd across all iterations and chains
  pooled = matrix(array3d, nrow = dim(array3d)[1] * dim(array3d)[2], ncol = nparam)
  means = colMeans(pooled)
  sds = apply(pooled, 2, sd)
  mcse = sds / sqrt(ess)

  result = cbind(mean = means, mcse = mcse, sd = sds, n_eff = ess, Rhat = rhat)

  if(is.null(param_names)) {
    data.frame(parameter = paste0("parameter [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(parameter = param_names, result, check.names = FALSE)
  }
}

# Summarize binary indicator variables
summarize_indicator = function(fit, component = c("indicator_samples"), param_names = NULL, array3d = NULL) {
  component = match.arg(component) # Add options later
  if(is.null(array3d)) array3d = combine_chains(fit, component)
  nparam = dim(array3d)[3]

  # Batch indicator ESS + transition counts via C++
  ind_stats = .compute_indicator_ess_cpp(array3d)
  batch_rhat = .compute_rhat_cpp(split_chains(array3d))

  result = cbind(ind_stats[, c("mean", "mcse", "sd", "n00", "n01", "n10", "n11", "n_eff_mixt"), drop = FALSE], Rhat = batch_rhat)
  colnames(result)[4:7] = c("n0->0", "n0->1", "n1->0", "n1->1")

  # Where n_eff_mixt is NA (constant chain), Rhat should also be NA
  result[is.na(result[, "n_eff_mixt"]), "Rhat"] = NA_real_

  if(is.null(param_names)) {
    data.frame(parameter = paste0("indicator [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(
      parameter = paste0(param_names, "- indicator"),
      result, check.names = FALSE
    )
  }
}

# Summarize the Rao-Blackwellized inclusion draws J (continuous, in [0, 1]) with
# the standard continuous machinery: a lower-variance inclusion-probability mean
# plus MCSE/ESS/split-Rhat on the RB draws, which the binary indicator draws
# cannot give. Columns that were never updated (unselected indicators, e.g. main
# differences when main_difference_selection = FALSE) are all NA and get an NA
# row rather than being fed to the ESS/Rhat kernels. Columns whose RB draws are
# constant to double precision keep their mean but report NA for the RB
# mcse/n_eff/Rhat (see the masking block below).
summarize_rb_inclusion = function(raw, param_names, keep_parameter_col = FALSE) {
  array3d = combine_chains(raw, "rb_inclusion_samples")
  nparam = dim(array3d)[3]
  cols = c("mean", "mcse", "sd", "n_eff", "Rhat")
  mat = matrix(NA_real_,
    nrow = nparam, ncol = length(cols),
    dimnames = list(NULL, cols)
  )

  keep = vapply(
    seq_len(nparam),
    function(k) any(is.finite(array3d[, , k])),
    logical(1)
  )
  if(any(keep)) {
    sub = summarize_manual(raw, array3d = array3d[, , keep, drop = FALSE])
    mat[keep, ] = as.matrix(sub[, cols, drop = FALSE])
  }

  ind_stats = .compute_indicator_ess_cpp(combine_chains(raw, "indicator_samples"))

  # Mask the RB precision/convergence columns only where the RB draws are
  # constant to double precision: there is no variability from which an MCSE, an
  # ESS, or a split-R-hat could be formed, and the mean is at its numerical
  # bound. Everywhere else the RB draws vary and the continuous machinery
  # applies as it does to any smooth quantity -- including on edges whose
  # indicator never flipped, where the RB chain is still informative and the
  # MCSE is what a boundary-stable fragility check needs. Draws whose variance
  # falls below the autocovariance kernel's numerical floor pick up an NA ESS
  # there rather than here, with the same meaning. Chains stuck constant at
  # different values are not masked: their ESS is NA and their split-R-hat is
  # +Inf, which is the alarm. The directional flip counts stay in the table.
  constant = is.finite(mat[, "sd"]) & mat[, "sd"] == 0
  mat[constant, c("mcse", "n_eff", "Rhat")] = NA_real_

  full = cbind(
    mat[, c("mean", "mcse", "sd", "n_eff", "Rhat"), drop = FALSE],
    ind_stats[, c("n01", "n10"), drop = FALSE]
  )
  # Keep the directional transition counts whole; their asymmetry is
  # decision-relevant and no symmetric summary recovers it.
  colnames(full)[colnames(full) == "n01"] = "n0->1"
  colnames(full)[colnames(full) == "n10"] = "n1->0"

  if(keep_parameter_col) {
    data.frame(parameter = param_names, full, check.names = FALSE, row.names = NULL)
  } else {
    out = as.data.frame(full, check.names = FALSE)
    rownames(out) = param_names
    out
  }
}

# Summarize slab values where indicators are 1
summarize_slab = function(fit, component = c("pairwise_samples"), param_names = NULL, array3d = NULL, array3d_ind = NULL) {
  component = match.arg(component) # Add options later
  if(is.null(array3d)) array3d = combine_chains(fit, component)
  nparam = dim(array3d)[3]
  result = matrix(NA, nparam, 5)
  colnames(result) = c("mean", "mcse", "sd", "n_eff", "Rhat")

  for(j in seq_len(nparam)) {
    draws = array3d[, , j]
    vec = as.vector(draws)
    if(!is.null(array3d_ind)) {
      selected = as.vector(array3d_ind[, , j]) == 1
    } else {
      selected = vec != 0
    }
    vec = vec[selected]
    n_total = length(vec)

    if(n_total >= 1) {
      result[j, "mean"] = mean(vec)
    }
    if(n_total > 10) {
      sdev = sd(vec)
      est = compute_rhat_ess(vec) ## draws
      mcse = sdev / sqrt(est$ess)
      result[j, c("sd", "mcse", "n_eff", "Rhat")] = c(sdev, mcse, est$ess, est$rhat)
    }
  }

  if(is.null(param_names)) {
    data.frame(parameter = paste0("weight [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(
      parameter = paste0(param_names, "- weight"),
      result, check.names = FALSE
    )
  }
}

# Derived composite ESS for a model-averaged (spike-and-slab) weight.
# n_eff = Var(weight) / composite_MCSE^2, the only calibrated ESS for the
# mixture-output row (a raw-chain ESS on the effect inflates in
# inclusion-dominated cells). The composite MCSE^2 splits into an inclusion part
# (mu_cond^2 * var_p, with var_p the RB J-chain MCSE^2 of the inclusion
# probability) and a slab part (p_hat^2 * var_mu, from the included-only draws).
# share_incl is the inclusion part's share of the composite MCSE^2 -- the
# bottleneck: whether the Monte Carlo error is dominated by inclusion or slab
# uncertainty. Inputs may be vectors (one entry per edge).
derived_weight_ess = function(post_var, p_hat, var_p, mu_cond, var_mu) {
  # A constant / never-updated inclusion chain contributes no inclusion error.
  var_p[!is.finite(var_p)] = 0
  part_incl = mu_cond^2 * var_p
  part_slab = p_hat^2 * var_mu
  mcse2 = part_incl + part_slab
  ok = is.finite(mcse2) & mcse2 > 0
  list(
    mcse = ifelse(ok, sqrt(mcse2), NA_real_),
    n_eff = ifelse(ok & is.finite(post_var), post_var / mcse2, NA_real_),
    share_incl = ifelse(ok, part_incl / mcse2, NA_real_)
  )
}

# Combined summary for pairwise parameters with selection
summarize_pair = function(fit,
                          indicator_component = c("indicator_samples"),
                          slab_component = c("pairwise_samples"),
                          param_names = NULL,
                          summ_ind = NULL,
                          summ_slab = NULL,
                          array3d_id = NULL,
                          array3d_pw = NULL) {
  indicator_component = match.arg(indicator_component) # Add options later
  slab_component = match.arg(slab_component) # Add options later

  if(is.null(array3d_id)) array3d_id = combine_chains(fit, indicator_component)
  if(is.null(array3d_pw)) array3d_pw = combine_chains(fit, slab_component)
  if(is.null(summ_slab)) summ_slab = summarize_slab(fit, component = slab_component, array3d = array3d_pw, array3d_ind = array3d_id)
  nparam = dim(array3d_pw)[3]

  # Posterior mean and variance of the model-averaged weight, straight from the
  # raw effect chain (exact; the spike contributes exact zeros).
  pooled_pw = matrix(array3d_pw, nrow = dim(array3d_pw)[1] * dim(array3d_pw)[2], ncol = nparam)
  eap = colMeans(pooled_pw)
  post_var = apply(pooled_pw, 2, stats::var)

  # Inclusion probability and its Monte Carlo error from the RB J-chain; fall
  # back to the binary indicator's transition-based MCSE when RB draws are
  # unavailable (pre-0.2.0.0 fits).
  if(!is.null(fit[[1]][["rb_inclusion_samples"]])) {
    rb = summarize_manual(fit, array3d = combine_chains(fit, "rb_inclusion_samples"))
    p_hat = rb$mean
    p_mcse = rb$mcse
  } else {
    if(is.null(summ_ind)) summ_ind = summarize_indicator(fit, component = indicator_component, array3d = array3d_id)
    p_hat = summ_ind$mean
    p_mcse = summ_ind$mcse
  }

  comp = derived_weight_ess(
    post_var = post_var, p_hat = p_hat, var_p = p_mcse^2,
    mu_cond = summ_slab$mean, var_mu = summ_slab$mcse^2
  )

  rhat = .compute_rhat_cpp(split_chains(array3d_pw))
  names_out = if(is.null(param_names)) {
    paste0("weight [", seq_len(nparam), "]")
  } else {
    paste0(param_names, "- weight")
  }

  data.frame(
    parameter = names_out,
    mean = eap, mcse = comp$mcse, sd = sqrt(post_var),
    n_eff = comp$n_eff, share_incl = comp$share_incl, Rhat = rhat,
    check.names = FALSE
  )
}

# Unified summary dispatcher for either model type
summarize_fit = function(fit, edge_selection = FALSE) {
  main_summary = summarize_manual(fit, component = "main_samples")

  if(!edge_selection) {
    pair_summary = summarize_manual(fit, component = "pairwise_samples")
    return(list(main = main_summary, pairwise = pair_summary))
  }

  # Build 3D arrays once; reused by all summary functions below
  array3d_ind = combine_chains(fit, "indicator_samples")
  array3d_pw = combine_chains(fit, "pairwise_samples")

  # Compute indicator and slab summaries once (ind_summary backs the RB-absent
  # fallback and is returned for the raw inclusion table).
  ind_summary = summarize_indicator(fit, component = "indicator_samples", array3d = array3d_ind)
  slab_summary = summarize_slab(fit, component = "pairwise_samples", array3d = array3d_pw, array3d_ind = array3d_ind)

  # The derived composite ESS handles always-included edges naturally: with no
  # inclusion uncertainty its inclusion part vanishes and n_eff reduces to the
  # slab ESS, so no special-casing of fully selected edges is needed.
  pair_summary = summarize_pair(fit,
    indicator_component = "indicator_samples",
    slab_component = "pairwise_samples",
    summ_ind = ind_summary,
    summ_slab = slab_summary,
    array3d_id = array3d_ind,
    array3d_pw = array3d_pw
  )

  list(main = main_summary, pairwise = pair_summary, indicator = ind_summary)
}


# NOTE: SBM posterior-summary helpers (summarize_alloc_pairs,
# find_representative_clustering, compute_p_k_given_t, posterior_summary_SBM)
# moved to mcmc_summary_sbm.R (cleanup S4).


summarize_manual_compare = function(fit_or_array,
                                    component = c("main_samples", "pairwise_samples"),
                                    param_names = NULL) {
  component = match.arg(component)

  # allow either a fit list or a pre-combined 3D array
  if(is.array(fit_or_array)) {
    array3d = fit_or_array
  } else {
    array3d = combine_chains(fit_or_array, component)
  }

  nparam = dim(array3d)[3]

  # Batch computation via C++
  ess = .compute_ess_cpp(array3d)
  rhat = .compute_rhat_cpp(split_chains(array3d))

  # Vectorized mean and sd across all iterations and chains
  pooled = matrix(array3d, nrow = dim(array3d)[1] * dim(array3d)[2], ncol = nparam)
  means = colMeans(pooled)
  sds = apply(pooled, 2, sd)
  mcse = sds / sqrt(ess)

  result = cbind(mean = means, mcse = mcse, sd = sds, n_eff = ess, Rhat = rhat)

  if(is.null(param_names)) {
    data.frame(parameter = paste0("param [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(parameter = param_names, result, check.names = FALSE)
  }
}


summarize_indicator_compare = function(fit, component = "indicator_samples", param_names = NULL) {
  array3d = combine_chains(fit, component)
  nparam = dim(array3d)[3]

  # Batch indicator ESS + transition counts via C++
  ind_stats = .compute_indicator_ess_cpp(array3d)
  batch_rhat = .compute_rhat_cpp(split_chains(array3d))

  result = cbind(ind_stats[, c("mean", "mcse", "sd", "n00", "n01", "n10", "n11", "n_eff_mixt"), drop = FALSE], Rhat = batch_rhat)
  colnames(result)[4:7] = c("n0->0", "n0->1", "n1->0", "n1->1")

  result[is.na(result[, "n_eff_mixt"]), "Rhat"] = NA_real_

  if(is.null(param_names)) {
    data.frame(parameter = paste0("indicator [", seq_len(nparam), "]"), result, check.names = FALSE)
  } else {
    data.frame(parameter = param_names, result, check.names = FALSE)
  }
}


# Summarize one effect with spike-and-slab draws. n_eff is the derived composite
# ESS for the model-averaged weight (see derived_weight_ess); share_incl is the
# inclusion part's share of the composite MCSE^2. The inclusion Monte Carlo error
# comes from the RB J-chain (draws_rb) when available, else from the binary
# indicator's transition-based MCSE.
summarize_mixture_effect = function(draws_pw, draws_id, name, draws_rb = NULL) {
  # Handle case where single-chain extraction returns a vector
  # (dimension gets dropped when extracting [, , idx] from array with nchains=1)
  if(is.null(dim(draws_pw))) {
    draws_pw = matrix(draws_pw, ncol = 1L)
  }
  if(is.null(dim(draws_id))) {
    draws_id = matrix(draws_id, ncol = 1L)
  }

  nchains = ncol(draws_pw)
  niter = nrow(draws_pw)

  ## --- slab part (included-only draws) ---
  vec = as.vector(draws_pw)
  vec = vec[vec != 0]
  if(length(vec) > 10) {
    eap_slab = mean(vec)
    est_slab = compute_rhat_ess(vec) # treat as single chain
    mcse_slab = sqrt(var(vec)) / sqrt(est_slab$ess)
  } else {
    eap_slab = NA_real_
    mcse_slab = NA_real_
  }

  ## --- inclusion part: RB J-chain when available, else the binary indicator ---
  if(!is.null(draws_rb) && any(is.finite(draws_rb))) {
    rb_array = array(draws_rb, dim = c(niter, nchains, 1L))
    pooled_rb = as.vector(rb_array)
    p_hat = mean(pooled_rb, na.rm = TRUE)
    ess_rb = .compute_ess_cpp(rb_array)[1]
    sd_rb = stats::sd(pooled_rb)
    p_mcse = if(is.finite(ess_rb) && ess_rb > 0 && is.finite(sd_rb)) sd_rb / sqrt(ess_rb) else 0
  } else {
    id_stats = .compute_indicator_ess_cpp(array(draws_id, dim = c(niter, nchains, 1L)))
    p_hat = id_stats[1, "mean"]
    p_mcse = id_stats[1, "mcse"]
  }

  ## --- combined summaries: exact posterior mean/variance from the raw effect
  ## chain (includes the spike zeros), composite MCSE and derived ESS ---
  pooled_pw = as.vector(draws_pw)
  posterior_mean = mean(pooled_pw)
  post_var = var(pooled_pw)

  comp = derived_weight_ess(
    post_var = post_var, p_hat = p_hat, var_p = p_mcse^2,
    mu_cond = eap_slab, var_mu = mcse_slab^2
  )

  pw_array = array(draws_pw, dim = c(niter, nchains, 1L))
  rhat = if(nchains > 1) .compute_rhat_cpp(split_chains(pw_array))[1] else NA_real_

  data.frame(
    parameter = name,
    mean = posterior_mean,
    mcse = comp$mcse,
    sd = sqrt(post_var),
    n_eff = comp$n_eff,
    share_incl = comp$share_incl,
    Rhat = rhat,
    check.names = FALSE
  )
}


# --- indicator index helpers (1-based) ---
indicator_row_starts = function(V) {
  # positions where each "row i" (i..V) starts in the flattened (i,j) list
  starts = integer(V)
  starts[1L] = 1L
  if(V > 1L) {
    for(i in 2L:V) {
      # previous row length = V - (i-1) + 1
      starts[i] = starts[i - 1L] + (V - (i - 1L) + 1L)
    }
  }
  starts
}


summarize_main_diff_compare = function(
  fit,
  main_effect_indices,
  num_groups,
  param_names = NULL
) {
  main_effect_samples = combine_chains(fit, "main_samples")
  indicator_samples = combine_chains(fit, "indicator_samples")
  rb_samples = if(!is.null(fit[[1]][["rb_inclusion_samples"]])) {
    combine_chains(fit, "rb_inclusion_samples")
  } else {
    NULL
  }

  V = nrow(main_effect_indices)
  num_main = main_effect_indices[V, 2] + 1L # total rows in main-effects matrix
  indicator_index_main = function(i, V) indicator_row_starts(V)[i]

  results = list()
  counter = 0L

  for(v in seq_len(V)) {
    id_idx = indicator_index_main(v, V) # (v,v) position in flattened indicators
    draws_id = indicator_samples[, , id_idx]
    draws_rb = if(!is.null(rb_samples)) rb_samples[, , id_idx] else NULL

    # rows in main-effects matrix belonging to variable v (1-based, inclusive)
    start = main_effect_indices[v, 1] + 1L
    stop = main_effect_indices[v, 2] + 1L

    for(row in start:stop) {
      category = row - start + 1L
      for(h in 1L:(num_groups - 1L)) {
        counter = counter + 1L
        col_index = h * num_main + row # group-major blocks of length num_main
        draws_pw = main_effect_samples[, , col_index]

        pname = if(!is.null(param_names)) {
          # param_names is laid out contrast-major (all rows of contrast 1,
          # then contrast 2, ...).
          param_names[(h - 1L) * num_main + row]
        } else {
          paste0("var", v, " (diff", h, "; ", category, ")")
        }

        results[[counter]] = summarize_mixture_effect(draws_pw, draws_id, pname, draws_rb = draws_rb)
      }
    }
  }

  out = do.call(rbind, results)
  rownames(out) = NULL
  out
}


summarize_pairwise_diff_compare = function(
  fit,
  pairwise_effect_indices,
  num_variables,
  num_groups,
  param_names = NULL
) {
  pairwise_effect_samples = combine_chains(fit, "pairwise_samples")
  indicator_samples = combine_chains(fit, "indicator_samples")
  rb_samples = if(!is.null(fit[[1]][["rb_inclusion_samples"]])) {
    combine_chains(fit, "rb_inclusion_samples")
  } else {
    NULL
  }

  V = num_variables
  num_pair = max(pairwise_effect_indices, na.rm = TRUE) + 1L # total rows in pairwise-effects matrix
  indicator_index_pair = function(i, j, V) indicator_row_starts(V)[i] + (j - i) # (i,j), i<j

  results = list()
  counter = 0L

  for(i in 1L:(V - 1L)) {
    for(j in (i + 1L):V) {
      id_idx = indicator_index_pair(i, j, V) # (i,j) in flattened indicators
      draws_id = indicator_samples[, , id_idx]
      draws_rb = if(!is.null(rb_samples)) rb_samples[, , id_idx] else NULL

      row = pairwise_effect_indices[i, j] + 1L # 1-based row into pairwise-effects matrix
      for(h in 1L:(num_groups - 1L)) {
        counter = counter + 1L
        col_index = h * num_pair + row # group-major blocks of length num_pair
        draws_pw = pairwise_effect_samples[, , col_index]

        pname = if(!is.null(param_names)) {
          # param_names is laid out contrast-major (all edges of contrast 1,
          # then contrast 2, ...).
          param_names[(h - 1L) * num_pair + row]
        } else {
          paste0("V", i, "-", j, " (diff", h, ")")
        }

        results[[counter]] = summarize_mixture_effect(draws_pw, draws_id, pname, draws_rb = draws_rb)
      }
    }
  }

  out = do.call(rbind, results)
  rownames(out) = NULL
  out
}


summarize_fit_compare = function(
  fit,
  main_effect_indices,
  pairwise_effect_indices,
  num_variables,
  num_groups,
  difference_selection = TRUE,
  param_names_main = NULL,
  param_names_pairwise = NULL,
  param_names_main_diff = NULL,
  param_names_pairwise_diff = NULL,
  param_names_indicators = NULL
) {
  count_main = function(main_effect_indices) {
    main_effect_indices[nrow(main_effect_indices), 2] + 1
  }

  count_pairwise = function(pairwise_effect_indices) {
    nr = nrow(pairwise_effect_indices)
    pairwise_effect_indices[nr, nr - 1] + 1
  }


  # --- main baseline
  array3d_main = combine_chains(fit, "main_samples")
  num_main = count_main(main_effect_indices)
  main_baseline = summarize_manual_compare(
    array3d_main[, , 1:num_main, drop = FALSE],
    "main_samples",
    param_names = param_names_main
  )

  # --- pairwise baseline
  array3d_pair = combine_chains(fit, "pairwise_samples")
  num_pair = count_pairwise(pairwise_effect_indices)
  pairwise_baseline = summarize_manual_compare(
    array3d_pair[, , 1:num_pair, drop = FALSE],
    "pairwise_samples",
    param_names = param_names_pairwise
  )

  if(!difference_selection) {
    # --- differences without selection -> treat as plain parameters
    # Drop baseline columns (col 1) and keep group-difference columns
    excl_baseline = 1:num_main
    main_diff_array = array3d_main[, , -excl_baseline, drop = FALSE]
    main_differences = summarize_manual_compare(
      main_diff_array, "main_samples",
      param_names = param_names_main_diff
    )

    excl_baseline = 1:num_pair
    pairwise_diff_array = array3d_pair[, , -excl_baseline, drop = FALSE]
    pairwise_differences = summarize_manual_compare(
      pairwise_diff_array, "pairwise_samples",
      param_names = param_names_pairwise_diff
    )

    indicators = NULL
  } else {
    # --- differences with selection -> use mixture summaries
    main_differences = summarize_main_diff_compare(
      fit, main_effect_indices, num_groups,
      param_names = param_names_main_diff
    )

    pairwise_differences = summarize_pairwise_diff_compare(
      fit, pairwise_effect_indices, num_variables, num_groups,
      param_names = param_names_pairwise_diff
    )

    indicators = summarize_indicator_compare(
      fit, "indicator_samples",
      param_names = param_names_indicators
    )
  }

  list(
    main_baseline        = main_baseline,
    pairwise_baseline    = pairwise_baseline,
    main_differences     = main_differences,
    pairwise_differences = pairwise_differences,
    indicators           = indicators
  )
}
