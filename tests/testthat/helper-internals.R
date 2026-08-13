# Make internal functions available when tests run outside R CMD check
# (e.g. via testthat::test_dir or testthat::test_file).
# Under R CMD check / test_local these are already in scope.
if(!exists("bgm_spec", mode = "function")) {
  internals = c(
    "bgm_spec",
    "bgmCompare_test_logp_and_gradient",
    "build_arguments",
    "collapse_categories_across_groups",
    "compute_conditional_ggm",
    "compute_conditional_mixed",
    "compute_conditional_probs",
    "get_explog_switch",
    "ggm_test_forward_map",
    "ggm_test_logp_and_gradient",
    "ggm_test_logp_and_gradient_prior",
    "mixed_test_logp_and_gradient",
    "rcpp_ieee754_exp",
    "rcpp_ieee754_log",
    "reformat_ordinal_data",
    "test_parameter_prior",
    "test_scale_prior",
    "unpack_interaction_prior",
    "unpack_threshold_prior",
    "run_mixed_simulation_parallel",
    "sample_mixed_mrf_gibbs",
    "validate_difference_prior",
    "validate_edge_prior",
    "validate_missing_data",
    "validate_sampler"
  )
  ns = asNamespace("bgms")
  for(fn in internals) {
    if(exists(fn, envir = ns, inherits = FALSE)) {
      assign(fn, get(fn, envir = ns), envir = globalenv())
    }
  }
  rm(internals, fn, ns)
}
