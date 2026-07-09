# Tests for bgm(precision_graph_prior = "hierarchical"): eligibility validation
# and the end-to-end fit with the Z-ratio engine and alarm suite attached.

hier_test_data = function(q = 10, n = 40, seed = 4) {
  set.seed(seed)
  matrix(rnorm(n * q), n, q)
}

test_that("hierarchical spec eligibility is validated", {
  Y = hier_test_data()

  expect_error(
    bgm(
      x = Y, variable_type = "continuous",
      interaction_prior = beta_prime_prior(),
      precision_graph_prior = "hierarchical",
      update_method = "gibbs", display_progress = "none", verbose = FALSE
    ),
    "normal or Cauchy"
  )
  expect_error(
    bgm(
      x = Y, variable_type = "continuous",
      interaction_prior = normal_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 2, rate = 2),
      precision_graph_prior = "hierarchical",
      update_method = "gibbs", display_progress = "none", verbose = FALSE
    ),
    "shape = 1"
  )
  expect_error(
    bgm(
      x = Y, variable_type = "continuous",
      interaction_prior = normal_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 1, rate = 2),
      edge_selection = FALSE,
      precision_graph_prior = "hierarchical",
      update_method = "gibbs", display_progress = "none", verbose = FALSE
    ),
    "edge_selection"
  )
  expect_error(
    bgm(
      x = matrix(sample(0:3, 200, replace = TRUE), 50, 4),
      interaction_prior = normal_prior(scale = 0.5),
      precision_graph_prior = "hierarchical",
      display_progress = "none", verbose = FALSE
    ),
    "continuous"
  )
})

test_that("the hierarchical spec accepts a Cauchy slab on every update method", {
  skip_on_cran()
  Y = hier_test_data(q = 8)
  for(method in c("nuts", "adaptive-metropolis", "gibbs")) {
    fit = bgm(
      x = Y, variable_type = "continuous",
      iter = 100, warmup = 150,
      interaction_prior = cauchy_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 1, rate = 2),
      precision_graph_prior = "hierarchical", calibration_window = 50,
      update_method = method, chains = 1, cores = 1, seed = 7,
      display_progress = "none", verbose = FALSE
    )
    s = summary(fit)
    expect_true(all(is.finite(s$pairwise$mean)), info = method)
    expect_true(all(s$indicator$mean >= 0 & s$indicator$mean <= 1),
      info = method
    )
    expect_equal(fit@arguments$precision_graph_prior, "hierarchical", info = method)
    expect_false(is.null(fit@zratio_diag), info = method)
  }
})

test_that("bgm fits the hierarchical spec and attaches the trust gauge", {
  skip_on_cran()
  Y = hier_test_data(q = 12)
  fit = bgm(
    x = Y, variable_type = "continuous",
    iter = 150, warmup = 250,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 2),
    edge_prior = beta_bernoulli_prior(2, 4),
    precision_graph_prior = "hierarchical", calibration_window = 100,
    update_method = "gibbs", chains = 2, cores = 2, seed = 11,
    display_progress = "none", verbose = FALSE
  )
  zd = fit@zratio_diag
  expect_false(is.null(zd))
  expect_equal(nrow(zd$per_chain), 2L)
  expect_true(all(is.finite(zd$per_chain$flip_rate)))
  expect_true(all(zd$per_chain$n_ent >= 0))
  expect_false(zd$flagged)
  expect_equal(fit@arguments$precision_graph_prior, "hierarchical")
  # The joint-path hyperparameter correction must not run on this path;
  # inclusion-parameter samples come from the clean conjugate draw.
  expect_equal(length(fit@inclusion_parameter_samples), 2L)
})

test_that("the joint default is unchanged", {
  skip_on_cran()
  Y = hier_test_data(q = 6)
  fit = bgm(
    x = Y, variable_type = "continuous",
    iter = 100, warmup = 150,
    update_method = "gibbs", chains = 1, cores = 1, seed = 3,
    display_progress = "none", verbose = FALSE
  )
  expect_equal(fit@arguments$precision_graph_prior, "joint")
  expect_null(fit@zratio_diag)
})

test_that("mixed data supports the hierarchical spec on the continuous block", {
  skip_on_cran()
  set.seed(9)
  n = 60
  X = cbind(
    matrix(sample(0:2, n * 3, replace = TRUE), n, 3),
    matrix(rnorm(n * 8), n, 8)
  )
  colnames(X) = paste0("V", seq_len(11))
  vt = c(rep("ordinal", 3), rep("continuous", 8))

  expect_error(
    bgm(
      x = X[, 1:4], variable_type = vt[1:4],
      interaction_prior = normal_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 1, rate = 2),
      precision_graph_prior = "hierarchical",
      display_progress = "none", verbose = FALSE
    ),
    "two continuous"
  )

  fit = bgm(
    x = X, variable_type = vt,
    iter = 120, warmup = 200,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 2),
    edge_prior = beta_bernoulli_prior(2, 4),
    precision_graph_prior = "hierarchical", calibration_window = 80,
    update_method = "adaptive-metropolis", chains = 1, cores = 1, seed = 5,
    display_progress = "none", verbose = FALSE
  )
  zd = fit@zratio_diag
  expect_false(is.null(zd))
  expect_true(is.finite(zd$per_chain$flip_rate))
  expect_false(zd$flagged)
  expect_equal(fit@arguments$precision_graph_prior, "hierarchical")
})
