# Tests for bgm(graph_prior_spec = "hierarchical"): eligibility validation
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
      interaction_prior = cauchy_prior(scale = 0.5),
      graph_prior_spec = "hierarchical",
      update_method = "gibbs", display_progress = "none", verbose = FALSE
    ),
    "normal"
  )
  expect_error(
    bgm(
      x = Y, variable_type = "continuous",
      interaction_prior = normal_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 2, rate = 2),
      graph_prior_spec = "hierarchical",
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
      graph_prior_spec = "hierarchical",
      update_method = "gibbs", display_progress = "none", verbose = FALSE
    ),
    "edge_selection"
  )
  expect_error(
    bgm(
      x = matrix(sample(0:3, 200, replace = TRUE), 50, 4),
      interaction_prior = normal_prior(scale = 0.5),
      graph_prior_spec = "hierarchical",
      display_progress = "none", verbose = FALSE
    ),
    "continuous"
  )
})

test_that("bgm fits the hierarchical spec and attaches the alarm suite", {
  skip_on_cran()
  Y = hier_test_data(q = 12)
  fit = bgm(
    x = Y, variable_type = "continuous",
    iter = 150, warmup = 250,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 2),
    edge_prior = beta_bernoulli_prior(2, 4),
    graph_prior_spec = "hierarchical", calibration_window = 100,
    update_method = "gibbs", chains = 2, cores = 2, seed = 11,
    display_progress = "none", verbose = FALSE
  )
  zd = fit@zratio_diag
  expect_false(is.null(zd))
  expect_equal(nrow(zd$per_chain), 2L)
  expect_true(all(zd$per_chain$frozen))
  expect_true(all(zd$per_chain$n_oracle > 0))
  expect_false(zd$verdict_flagged)
  expect_equal(length(zd$audits), 2L)
  expect_equal(fit@arguments$graph_prior_spec, "hierarchical")
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
  expect_equal(fit@arguments$graph_prior_spec, "joint")
  expect_null(fit@zratio_diag)
})

test_that("mixed indicator layout rebuilds the continuous subgraph", {
  q = 5
  G = matrix(0L, q, q)
  G[1, 4] = G[4, 1] = 1L
  G[2, 3] = G[3, 2] = 1L
  diag(G) = 1L
  # Row-major upper triangle without diagonal, as the mixed chain emits it.
  rowmajor = integer(q * (q - 1) / 2)
  e = 1L
  for(i in seq_len(q - 1)) {
    for(j in (i + 1):q) {
      rowmajor[e] = G[i, j]
      e = e + 1L
    }
  }
  expect_identical(
    bgms:::zratio_indicator_graph(rowmajor, q, include_diag = FALSE), G
  )
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
      graph_prior_spec = "hierarchical",
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
    graph_prior_spec = "hierarchical", calibration_window = 80,
    update_method = "adaptive-metropolis", chains = 1, cores = 1, seed = 5,
    display_progress = "none", verbose = FALSE
  )
  zd = fit@zratio_diag
  expect_false(is.null(zd))
  expect_true(zd$per_chain$frozen)
  expect_false(zd$verdict_flagged)
  expect_equal(fit@arguments$graph_prior_spec, "hierarchical")
})
