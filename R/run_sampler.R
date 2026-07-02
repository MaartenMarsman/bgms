# ==============================================================================
# run_sampler: dispatch from bgm_spec to C++ backends
# ==============================================================================
#
# Thin dispatcher that reads a validated bgm_spec and calls the appropriate
# C++ sampling function. Returns the raw per-chain output lists from C++.
# ==============================================================================


# ------------------------------------------------------------------
# bb_between_or_sentinel
# ------------------------------------------------------------------
# Maps NULL to -1.0 (C++ sentinel for "no between-cluster prior").
#
# @param value  Scalar or NULL from the prior spec.
#
# Returns: value unchanged, or -1.0 when NULL.
# ------------------------------------------------------------------
bb_between_or_sentinel = function(value) {
  if(is.null(value)) -1.0 else value
}


# ==============================================================================
# run_sampler()  --- main dispatcher
# ==============================================================================
run_sampler = function(spec) {
  stopifnot(inherits(spec, "bgm_spec"))

  raw = switch(spec$model_type,
    ggm       = run_sampler_ggm(spec),
    omrf      = run_sampler_omrf(spec),
    mixed_mrf = run_sampler_mixed_mrf(spec),
    compare   = run_sampler_compare(spec),
    stop("Unknown model_type: ", spec$model_type)
  )

  # Check for chain-level errors
  chain_errors = vapply(raw, function(ch) isTRUE(ch$error), logical(1L))
  if(all(chain_errors)) {
    msgs = vapply(raw, function(ch) ch$error_msg %||% "unknown error", character(1L))
    stop("All chains failed. First error: ", msgs[1L])
  }
  if(any(chain_errors)) {
    n_fail = sum(chain_errors)
    warning(n_fail, " of ", length(raw), " chain(s) failed and will be dropped.")
    raw = raw[!chain_errors]
  }

  # Check for user interrupt across all chains
  userInterrupt = any(vapply(raw, `[[`, logical(1L), "userInterrupt"))
  attr(raw, "userInterrupt") = userInterrupt
  if(userInterrupt) {
    warning("Stopped sampling after user interrupt, results are likely uninterpretable.")
  }

  raw
}


# ==============================================================================
# run_sampler_ggm()
# ==============================================================================
run_sampler_ggm = function(spec) {
  d = spec$data
  p = spec$prior
  s = spec$sampler
  m = spec$missing

  bb_alpha_between = bb_between_or_sentinel(p$beta_bernoulli_alpha_between)
  bb_beta_between = bb_between_or_sentinel(p$beta_bernoulli_beta_between)

  correction = ggm_edge_prior_correction(p, s, d$num_variables)

  out_raw = sample_ggm(
    inputFromR = list(
      X = d$x,
      pairwise_scale = p$pairwise_scale,
      interaction_prior_type = p$interaction_prior_type,
      interaction_alpha = p$interaction_alpha,
      interaction_beta = p$interaction_beta,
      scale_prior_type = p$scale_prior_type,
      scale_shape = p$scale_shape,
      scale_rate = p$scale_rate
    ),
    prior_inclusion_prob = p$inclusion_probability,
    initial_edge_indicators = matrix(1L,
      nrow = d$num_variables,
      ncol = d$num_variables
    ),
    no_iter = s$iter,
    no_warmup = s$warmup,
    no_chains = s$chains,
    edge_selection = p$edge_selection,
    sampler_type = s$update_method,
    seed = s$seed,
    no_threads = s$cores,
    progress_type = s$progress_type,
    edge_prior = p$edge_prior,
    beta_bernoulli_alpha = p$beta_bernoulli_alpha,
    beta_bernoulli_beta = p$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = bb_alpha_between,
    beta_bernoulli_beta_between = bb_beta_between,
    dirichlet_alpha = p$dirichlet_alpha,
    lambda = p$lambda,
    target_acceptance = s$target_accept,
    max_tree_depth = s$nuts_max_depth,
    na_impute = m$na_impute,
    missing_index_nullable = m$missing_index,
    delta = p$delta,
    edge_prior_correction = correction
  )

  out_raw
}


# ==============================================================================
# run_sampler_omrf()
# ==============================================================================
run_sampler_omrf = function(spec) {
  d = spec$data
  v = spec$variables
  m = spec$missing
  p = spec$prior
  s = spec$sampler

  bb_alpha_between = bb_between_or_sentinel(p$beta_bernoulli_alpha_between)
  bb_beta_between = bb_between_or_sentinel(p$beta_bernoulli_beta_between)

  input_list = list(
    observations           = d$x,
    num_categories         = d$num_categories,
    is_ordinal_variable    = v$is_ordinal,
    baseline_category      = v$baseline_category,
    interaction_prior_type = p$interaction_prior_type,
    pairwise_scale         = p$pairwise_scale,
    interaction_alpha      = p$interaction_alpha,
    interaction_beta       = p$interaction_beta,
    threshold_prior_type   = p$threshold_prior_type,
    main_alpha             = p$main_alpha,
    main_beta              = p$main_beta,
    threshold_scale        = p$threshold_scale
  )

  out_raw = sample_omrf(
    inputFromR = input_list,
    prior_inclusion_prob = p$inclusion_probability,
    initial_edge_indicators = matrix(1L,
      nrow = d$num_variables,
      ncol = d$num_variables
    ),
    no_iter = s$iter,
    no_warmup = s$warmup,
    no_chains = s$chains,
    no_threads = s$cores,
    progress_type = s$progress_type,
    progress_callback = s$progress_callback,
    edge_selection = p$edge_selection,
    sampler_type = s$update_method,
    seed = s$seed,
    edge_prior = p$edge_prior,
    na_impute = m$na_impute,
    missing_index_nullable = m$missing_index,
    beta_bernoulli_alpha = p$beta_bernoulli_alpha,
    beta_bernoulli_beta = p$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = bb_alpha_between,
    beta_bernoulli_beta_between = bb_beta_between,
    dirichlet_alpha = p$dirichlet_alpha,
    lambda = p$lambda,
    target_acceptance = s$target_accept,
    max_tree_depth = s$nuts_max_depth,
    pairwise_scaling_factors_nullable = p$pairwise_scaling_factors
  )

  out_raw
}


# ==============================================================================
# run_sampler_mixed_mrf()
# ==============================================================================
run_sampler_mixed_mrf = function(spec) {
  d = spec$data
  v = spec$variables
  m = spec$missing
  p = spec$prior
  s = spec$sampler

  bb_alpha_between = bb_between_or_sentinel(p$beta_bernoulli_alpha_between)
  bb_beta_between = bb_between_or_sentinel(p$beta_bernoulli_beta_between)

  correction = ggm_edge_prior_correction(
    p, s, d$num_variables, d$num_continuous
  )

  input_list = list(
    discrete_observations   = d$x_discrete,
    continuous_observations = d$x_continuous,
    num_categories          = d$num_categories,
    is_ordinal_variable     = as.integer(v$is_ordinal),
    baseline_category       = v$baseline_category,
    interaction_prior_type  = p$interaction_prior_type,
    pairwise_scale          = p$pairwise_scale,
    interaction_alpha       = p$interaction_alpha,
    interaction_beta        = p$interaction_beta,
    threshold_prior_type    = p$threshold_prior_type,
    main_alpha              = p$main_alpha,
    main_beta               = p$main_beta,
    threshold_scale         = p$threshold_scale,
    means_prior_type        = p$means_prior_type,
    means_scale             = p$means_scale,
    means_alpha             = p$means_alpha,
    means_beta              = p$means_beta,
    scale_prior_type        = p$scale_prior_type,
    scale_shape             = p$scale_shape,
    scale_rate              = p$scale_rate
  )

  out_raw = sample_mixed_mrf(
    inputFromR = input_list,
    prior_inclusion_prob = p$inclusion_probability,
    initial_edge_indicators = matrix(1L,
      nrow = d$num_variables,
      ncol = d$num_variables
    ),
    no_iter = s$iter,
    no_warmup = s$warmup,
    no_chains = s$chains,
    edge_selection = p$edge_selection,
    seed = s$seed,
    no_threads = s$cores,
    progress_type = s$progress_type,
    progress_callback = s$progress_callback,
    edge_prior = p$edge_prior,
    beta_bernoulli_alpha = p$beta_bernoulli_alpha,
    beta_bernoulli_beta = p$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = bb_alpha_between,
    beta_bernoulli_beta_between = bb_beta_between,
    dirichlet_alpha = p$dirichlet_alpha,
    lambda = p$lambda,
    sampler_type = s$update_method,
    target_acceptance = s$target_accept,
    max_tree_depth = s$nuts_max_depth,
    na_impute = m$na_impute,
    missing_index_discrete_nullable = m$missing_index_discrete,
    missing_index_continuous_nullable = m$missing_index_continuous,
    delta = p$delta,
    edge_prior_correction = correction
  )

  out_raw
}


# ==============================================================================
# run_sampler_compare()
# ==============================================================================
run_sampler_compare = function(spec) {
  d = spec$data
  v = spec$variables
  m = spec$missing
  p = spec$prior
  s = spec$sampler
  pc = spec$precomputed

  run_bgmCompare_parallel(
    observations = d$x,
    num_groups = d$num_groups,
    counts_per_category = pc$counts_per_category,
    blume_capel_stats = pc$blume_capel_stats,
    pairwise_stats = pc$pairwise_stats,
    num_categories = d$num_categories,
    main_alpha = p$main_alpha,
    main_beta = p$main_beta,
    pairwise_scale = p$pairwise_scale,
    pairwise_scaling_factors = p$pairwise_scaling_factors,
    difference_scale = p$difference_scale,
    difference_selection_alpha = p$beta_bernoulli_alpha,
    difference_selection_beta = p$beta_bernoulli_beta,
    difference_selection_alpha_between = p$beta_bernoulli_alpha_between,
    difference_selection_beta_between = p$beta_bernoulli_beta_between,
    difference_dirichlet_alpha = p$dirichlet_alpha,
    difference_lambda = p$lambda,
    difference_prior = p$difference_prior,
    iter = s$iter,
    warmup = s$warmup,
    na_impute = m$na_impute,
    missing_data_indices = m$missing_index,
    is_ordinal_variable = v$is_ordinal,
    baseline_category = v$baseline_category,
    difference_selection = p$difference_selection,
    main_difference_selection = p$main_difference_selection,
    main_effect_indices = pc$main_effect_indices,
    pairwise_effect_indices = pc$pairwise_effect_indices,
    target_accept = s$target_accept,
    nuts_max_depth = s$nuts_max_depth,
    learn_mass_matrix = s$learn_mass_matrix,
    projection = d$projection,
    group_membership = sort(d$group) - 1L,
    group_indices = d$group_indices,
    interaction_index_matrix = pc$interaction_index_matrix,
    inclusion_probability = p$inclusion_probability_difference,
    num_chains = s$chains,
    nThreads = s$cores,
    seed = s$seed,
    update_method = s$update_method,
    progress_type = s$progress_type,
    interaction_prior_type_str = p$interaction_prior_type,
    threshold_prior_type_str = p$threshold_prior_type,
    threshold_scale = if(is.na(p$threshold_scale)) 1.0 else p$threshold_scale,
    progress_callback = s$progress_callback
  )
}
