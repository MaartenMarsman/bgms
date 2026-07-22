# ==============================================================================
# Unit tests for build_arguments()
# Phase B.6 of the R scaffolding refactor.
#
# Verifies that build_arguments(spec) produces the same field set and values
# as the arguments lists in prepare_output_bgm / prepare_output_ggm /
# prepare_output_bgmCompare.
# ==============================================================================

# ==============================================================================
# Shared helpers
# ==============================================================================

make_continuous_data = function(n = 20, p = 3) {
  set.seed(42)
  x = matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) = paste0("V", seq_len(p))
  x
}

make_ordinal_data = function(n = 30, p = 3, max_cat = 2) {
  set.seed(99)
  x = matrix(sample(0:max_cat, n * p, replace = TRUE), nrow = n, ncol = p)
  colnames(x) = paste0("V", seq_len(p))
  x
}

spec_ggm = function(...) {
  defaults = list(
    x = make_continuous_data(),
    model_type = "ggm",
    variable_type = "continuous",
    update_method = "adaptive-metropolis",
    na_action = "listwise",
    iter = 500L,
    warmup = 100L,
    chains = 2L,
    cores = 1L,
    seed = 1L,
    display_progress = "none",
    verbose = FALSE
  )
  args = modifyList(defaults, list(...))
  do.call(bgm_spec, args)
}

spec_omrf = function(...) {
  defaults = list(
    x = make_ordinal_data(),
    model_type = "omrf",
    variable_type = "ordinal",
    na_action = "listwise",
    update_method = "nuts",
    iter = 500L,
    warmup = 100L,
    chains = 2L,
    cores = 1L,
    seed = 1L,
    display_progress = "none",
    verbose = FALSE
  )
  args = modifyList(defaults, list(...))
  do.call(bgm_spec, args)
}

spec_compare = function(...) {
  defaults = list(
    x = make_ordinal_data(n = 20, p = 3, max_cat = 2),
    y = make_ordinal_data(n = 25, p = 3, max_cat = 2),
    model_type = "compare",
    variable_type = "ordinal",
    na_action = "listwise",
    update_method = "nuts",
    iter = 500L,
    warmup = 100L,
    chains = 2L,
    cores = 1L,
    seed = 1L,
    display_progress = "none",
    verbose = FALSE
  )
  args = modifyList(defaults, list(...))
  do.call(bgm_spec, args)
}


# ==============================================================================
# 1.  GGM arguments <U+2014> field names match prepare_output_ggm
# ==============================================================================

test_that("GGM build_arguments: all expected field names present", {
  s = spec_ggm()
  a = build_arguments(s)
  expected = c(
    "num_variables", "num_cases", "na_impute", "variable_type",
    "iter", "warmup", "edge_selection", "edge_prior",
    "precision_graph_prior",
    "inclusion_probability", "beta_bernoulli_alpha", "beta_bernoulli_beta",
    "beta_bernoulli_alpha_between", "beta_bernoulli_beta_between",
    "dirichlet_alpha", "lambda", "na_action", "version",
    "update_method", "target_accept", "num_chains",
    "data_columnnames", "no_variables", "column_means", "is_continuous",
    "model_type"
  )
  expect_true(all(expected %in% names(a)),
    info = paste("Missing:", paste(setdiff(expected, names(a)),
      collapse = ", "
    ))
  )
  # no extra fields
  expect_true(all(names(a) %in% expected),
    info = paste("Extra:", paste(setdiff(names(a), expected),
      collapse = ", "
    ))
  )
})

test_that("GGM build_arguments: values are correct", {
  x = make_continuous_data(n = 15, p = 4)
  s = spec_ggm(
    x = x, edge_selection = TRUE, edge_prior = "Bernoulli",
    inclusion_probability = 0.3
  )
  a = build_arguments(s)
  expect_equal(a$num_variables, 4L)
  expect_equal(a$num_cases, 15L)
  expect_equal(a$no_variables, 4L)
  expect_true(a$is_continuous)
  expect_equal(a$update_method, "adaptive-metropolis")
  expect_true(a$edge_selection)
  expect_equal(a$edge_prior, "Bernoulli")
  expect_equal(a$na_action, "listwise")
  expect_equal(a$iter, 500L)
  expect_equal(a$warmup, 100L)
  expect_equal(a$num_chains, 2L)
  expect_equal(a$data_columnnames, colnames(x))
  expect_s3_class(a$version, "package_version")
})

test_that("GGM build_arguments: edge_selection = FALSE", {
  s = spec_ggm(edge_selection = FALSE)
  a = build_arguments(s)
  expect_false(a$edge_selection)
  expect_equal(a$edge_prior, "Not Applicable")
})


# ==============================================================================
# 2.  OMRF arguments <U+2014> field names match prepare_output_bgm
# ==============================================================================

test_that("OMRF build_arguments: all expected field names present", {
  s = spec_omrf()
  a = build_arguments(s)
  expected = c(
    "num_variables", "num_cases", "na_impute", "variable_type",
    "iter", "warmup", "pairwise_scale",
    "main_alpha", "main_beta",
    "edge_selection", "edge_prior", "inclusion_probability",
    "beta_bernoulli_alpha", "beta_bernoulli_beta",
    "beta_bernoulli_alpha_between", "beta_bernoulli_beta_between",
    "dirichlet_alpha", "lambda", "na_action", "version",
    "update_method", "target_accept",
    "nuts_max_depth", "learn_mass_matrix",
    "num_chains", "num_categories", "category_levels", "blume_capel_shift",
    "data_columnnames", "baseline_category",
    "no_variables",
    "model_type"
  )
  expect_true(all(expected %in% names(a)),
    info = paste("Missing:", paste(setdiff(expected, names(a)),
      collapse = ", "
    ))
  )
  expect_true(all(names(a) %in% expected),
    info = paste("Extra:", paste(setdiff(names(a), expected),
      collapse = ", "
    ))
  )
})

test_that("OMRF build_arguments: values are correct", {
  s = spec_omrf(
    pairwise_scale = 3.0, main_alpha = 0.7, main_beta = 0.3
  )
  a = build_arguments(s)
  expect_equal(a$pairwise_scale, 3.0)
  expect_equal(a$main_alpha, 0.7)
  expect_equal(a$main_beta, 0.3)
  expect_equal(length(a$num_categories), a$num_variables)
  expect_equal(length(a$baseline_category), a$num_variables)
  expect_equal(a$no_variables, a$num_variables)
})

test_that("OMRF build_arguments: Beta-Bernoulli prior params preserved", {
  s = spec_omrf(
    edge_selection = TRUE, edge_prior = "Beta-Bernoulli",
    beta_bernoulli_alpha = 2, beta_bernoulli_beta = 3
  )
  a = build_arguments(s)
  expect_equal(a$beta_bernoulli_alpha, 2)
  expect_equal(a$beta_bernoulli_beta, 3)
})


# ==============================================================================
# 3.  Compare arguments <U+2014> field names match prepare_output_bgmCompare
# ==============================================================================

test_that("Compare build_arguments: all expected field names present", {
  s = spec_compare()
  a = build_arguments(s)
  expected = c(
    "num_variables", "num_cases", "iter", "warmup",
    "pairwise_scale", "difference_scale",
    "difference_selection", "main_difference_selection",
    "difference_prior",
    "difference_selection_alpha", "difference_selection_beta",
    "difference_selection_alpha_between", "difference_selection_beta_between",
    "difference_dirichlet_alpha", "difference_lambda",
    "inclusion_probability",
    "version", "update_method", "target_accept",
    "nuts_max_depth", "learn_mass_matrix",
    "num_chains", "num_groups",
    "data_columnnames", "projection",
    "num_categories", "category_levels", "blume_capel_shift",
    "is_ordinal_variable",
    "group",
    "model_type"
  )
  expect_true(all(expected %in% names(a)),
    info = paste("Missing:", paste(setdiff(expected, names(a)),
      collapse = ", "
    ))
  )
  expect_true(all(names(a) %in% expected),
    info = paste("Extra:", paste(setdiff(names(a), expected),
      collapse = ", "
    ))
  )
})

test_that("Compare build_arguments: values are correct", {
  s = spec_compare(difference_scale = 3.0)
  a = build_arguments(s)
  expect_equal(a$num_groups, 2L)
  expect_equal(a$difference_scale, 3.0)
  expect_true(is.matrix(a$projection))
  expect_true(is.matrix(a$inclusion_probability))
  expect_equal(length(a$group), a$num_cases)
  expect_equal(length(a$is_ordinal_variable), a$num_variables)
  expect_equal(a$difference_selection_alpha, 1) # default beta_bernoulli_alpha
  expect_equal(a$difference_selection_beta, 1) # default beta_bernoulli_beta
})

test_that("Compare build_arguments: main_difference_selection preserved", {
  s = spec_compare(main_difference_selection = TRUE)
  a = build_arguments(s)
  expect_true(a$main_difference_selection)
})

test_that("Compare build_arguments: difference_selection = FALSE", {
  s = spec_compare(difference_selection = FALSE)
  a = build_arguments(s)
  expect_false(a$difference_selection)
})

test_that("Compare build_arguments: three groups via group_indicator", {
  x = make_ordinal_data(n = 60, p = 3, max_cat = 2)
  gi = rep(1:3, each = 20)
  s = spec_compare(x = x, y = NULL, group_indicator = gi)
  a = build_arguments(s)
  expect_equal(a$num_groups, 3L)
  expect_equal(ncol(a$projection), 2L)
})


# ==============================================================================
# 4.  Cross-model: version field always present
# ==============================================================================

test_that("build_arguments: version is package_version for all models", {
  for(s in list(spec_ggm(), spec_omrf(), spec_compare())) {
    a = build_arguments(s)
    expect_s3_class(a$version, "package_version")
  }
})
