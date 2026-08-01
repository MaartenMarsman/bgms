# --------------------------------------------------------------------------- #
# RATTLE + edge selection integration tests.
#
# Tests verify that the NUTS sampler correctly uses the RATTLE constrained
# integrator when edge selection is enabled. The RATTLE path operates in
# full Cholesky space with position/momentum projection for zero edges.
#
# Fast tests:
#   5.1  Smoke test: NUTS + edge_selection runs without error
#   5.2  Output structure is correct
#
# Slow tests (gated behind BGMS_RUN_SLOW_TESTS):
#   5.3  Pairwise interaction estimates agree with MH
#   5.4  NUTS diagnostics are reasonable with edge selection
# --------------------------------------------------------------------------- #


# ---- Skip gate ---------------------------------------------------------------

skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run RATTLE edge selection tests"
  )
}


# ---- 5.1  Smoke test --------------------------------------------------------

test_that("NUTS with edge selection runs without error (small problem)", {
  set.seed(51)
  p = 3
  n = 50
  x = matrix(rnorm(n * p), nrow = n)
  colnames(x) = paste0("V", seq_len(p))

  expect_no_error({
    fit = bgm(
      x = x, variable_type = "continuous",
      update_method = "nuts",
      edge_selection = TRUE,
      iter = 100, warmup = 50, chains = 1,
      seed = 510, display_progress = "none"
    )
  })
})


# ---- 5.2  Output structure --------------------------------------------------

test_that("NUTS edge selection output has expected fields", {
  set.seed(52)
  p = 3
  n = 50
  x = matrix(rnorm(n * p), nrow = n)
  colnames(x) = paste0("V", seq_len(p))

  fit = bgm(
    x = x, variable_type = "continuous",
    update_method = "nuts",
    edge_selection = TRUE,
    iter = 100, warmup = 50, chains = 1,
    seed = 520, display_progress = "none"
  )

  # PIPs should exist and be a p x p matrix with values in [0, 1]
  pip = fit$posterior_mean_indicator
  expect_true(is.matrix(pip))
  expect_equal(nrow(pip), p)
  expect_equal(ncol(pip), p)
  expect_true(all(pip >= 0 & pip <= 1))

  # Pairwise interactions should be extractable
  interactions = extract_pairwise_interactions(fit)
  expect_true(is.matrix(interactions))
  expect_equal(ncol(interactions), p * (p - 1) / 2)

  # NUTS diagnostics should exist
  expect_true(!is.null(fit$nuts_diag))
})


# ---- 5.3  Pairwise interaction estimates agree with MH ----------------------

test_that("NUTS RATTLE interaction estimates agree with MH (p=5)", {
  skip_unless_slow()

  p = 5
  n = 300
  set.seed(53)

  # Tridiagonal precision matrix
  K_true = diag(2, p)
  for(i in seq_len(p - 1)) {
    K_true[i, i + 1] = -0.5
    K_true[i + 1, i] = -0.5
  }
  Sigma = solve(K_true)
  x = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(x) = paste0("V", seq_len(p))

  n_iter = 5000
  n_warmup = 2000

  fit_mh = bgm(
    x = x, variable_type = "continuous",
    update_method = "adaptive-metropolis",
    edge_selection = TRUE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 530, display_progress = "none"
  )

  fit_nuts = bgm(
    x = x, variable_type = "continuous",
    update_method = "nuts",
    edge_selection = TRUE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 531, display_progress = "none"
  )

  # Compare posterior mean interactions
  int_mh = colMeans(extract_pairwise_interactions(fit_mh))
  int_nuts = colMeans(extract_pairwise_interactions(fit_nuts))

  # Max absolute difference should be small
  expect_lt(max(abs(int_mh - int_nuts)), 0.15,
    label = "max interaction difference < 0.15"
  )

  # PIPs should also agree
  pip_mh = fit_mh$posterior_mean_indicator[upper.tri(fit_mh$posterior_mean_indicator)]
  pip_nuts = fit_nuts$posterior_mean_indicator[upper.tri(fit_nuts$posterior_mean_indicator)]
  expect_lt(max(abs(pip_mh - pip_nuts)), 0.15,
    label = "max PIP difference < 0.15"
  )
})


# ---- 5.4  NUTS diagnostics with edge selection ------------------------------

test_that("NUTS diagnostics are reasonable with edge selection (p=4)", {
  skip_unless_slow()

  set.seed(54)
  p = 4
  n = 200

  K_true = diag(2, p)
  for(i in seq_len(p - 1)) {
    K_true[i, i + 1] = -0.5
    K_true[i + 1, i] = -0.5
  }
  Sigma = solve(K_true)
  x = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(x) = paste0("V", seq_len(p))

  fit = bgm(
    x = x, variable_type = "continuous",
    update_method = "nuts",
    edge_selection = TRUE,
    iter = 2000, warmup = 1000, chains = 2,
    seed = 540, display_progress = "none"
  )

  expect_true(!is.null(fit$nuts_diag))

  # Divergences should be zero or very few
  expect_lte(fit$nuts_diag$summary$total_divergences, 5)

  # EBFMI should be reasonable
  expect_gt(fit$nuts_diag$summary$min_ebfmi, 0.1)
})
