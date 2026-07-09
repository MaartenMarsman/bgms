# ==============================================================================
# tests/testthat/test-mixed-mrf-simulate-predict.R
# ==============================================================================
# Phase G tests for mixed MRF simulation and prediction.
#
# EXTENDS: test-simulate-predict-regression.R (which handles parameterized
# roundtrip tests via get_bgms_fixtures). This file covers:
#   T25: Gibbs generator sanity (sample_mixed_mrf_gibbs)
#   Mixed-specific structural tests for simulate.bgms / predict.bgms
#   Edge cases: p=1, q=1, binary-only ordinal
#
# PATTERN: Stochastic-robust testing <U+2014> dimensions, ranges, invariants,
# coarse distributional checks. No exact moment matching.
# ==============================================================================


# ==============================================================================
# 1. Gibbs generator sanity (T25)
# ==============================================================================
# Verify that sample_mixed_mrf_gibbs produces structurally correct output
# and coarse statistical properties match known targets.
# ==============================================================================

test_that("sample_mixed_mrf_gibbs returns correct dimensions", {
  p = 2
  q = 2
  n = 100
  nc = c(2L, 2L)
  mux = matrix(c(0, 0.5, -0.3, 0, -0.2, 0.1), nrow = p, ncol = 3, byrow = TRUE)
  pairwise_disc = matrix(c(0, 0.15, 0.15, 0), p, p)
  muy = c(0.5, -0.2)
  pairwise_cont = matrix(c(-0.75, -0.1, -0.1, -0.9), q, q)
  pairwise_cross = matrix(c(0.1, -0.15, 0.2, 0.05), p, q)

  result = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = c("ordinal", "ordinal"),
    baseline_category_r = c(0L, 0L), iter = 200L, seed = 42L
  )

  expect_true(is.list(result))
  expect_equal(dim(result$x), c(n, p))
  expect_equal(dim(result$y), c(n, q))
})

test_that("sample_mixed_mrf_gibbs: discrete values in valid range", {
  p = 3
  q = 1
  n = 500
  nc = c(2L, 3L, 1L) # categories 0-2, 0-3, 0-1 (binary)
  mux = matrix(0, p, 4)
  pairwise_disc = matrix(0, p, p)
  muy = 0
  pairwise_cont = matrix(-1)
  pairwise_cross = matrix(0, p, q)

  result = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = c("ordinal", "ordinal", "ordinal"),
    baseline_category_r = c(0L, 0L, 0L), iter = 100L, seed = 1L
  )

  for(s in 1:p) {
    expect_true(all(result$x[, s] >= 0),
      info = sprintf("var %d has negative values", s)
    )
    expect_true(all(result$x[, s] <= nc[s]),
      info = sprintf("var %d exceeds max category %d", s, nc[s])
    )
  }
})

test_that("sample_mixed_mrf_gibbs: continuous marginal SD matches association scale", {
  # Independent model: pairwise_disc = 0, pairwise_cross = 0, so y ~ N(muy, precision^{-1})
  # precision = -2 * pairwise_cont
  p = 1
  q = 2
  n = 2000
  nc = c(2L)
  mux = matrix(0, p, 3)
  pairwise_disc = matrix(0)
  muy = c(1.0, -0.5)
  pairwise_cont = diag(c(-1.0, -0.25))
  pairwise_cross = matrix(0, p, q)

  result = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = "ordinal", baseline_category_r = 0L,
    iter = 300L, seed = 99L
  )

  for(j in 1:q) {
    expected_sd = 1 / sqrt(-2 * pairwise_cont[j, j])
    empirical_sd = sd(result$y[, j])
    # Loose check: within 30% of expected
    expect_true(
      abs(empirical_sd - expected_sd) / expected_sd < 0.3,
      info = sprintf(
        "y%d SD: expected %.3f, got %.3f",
        j, expected_sd, empirical_sd
      )
    )
  }
})

test_that("sample_mixed_mrf_gibbs: seed reproducibility", {
  p = 2
  q = 1
  n = 50
  nc = c(2L, 2L)
  mux = matrix(c(0, 0.5, -0.3, 0, -0.2, 0.1), nrow = p, ncol = 3, byrow = TRUE)
  pairwise_disc = matrix(c(0, 0.15, 0.15, 0), p, p)
  muy = 0
  pairwise_cont = matrix(-0.5)
  pairwise_cross = matrix(c(0.1, 0.2), p, q)

  args = list(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = c("ordinal", "ordinal"),
    baseline_category_r = c(0L, 0L), iter = 100L, seed = 42L
  )

  r1 = do.call(sample_mixed_mrf_gibbs, args)
  r2 = do.call(sample_mixed_mrf_gibbs, args)

  expect_equal(r1$x, r2$x)
  expect_equal(r1$y, r2$y)
})

test_that("sample_mixed_mrf_gibbs: Blume-Capel variable works", {
  p = 2
  q = 1
  n = 200
  nc = c(2L, 4L)
  # ordinal: mux has num_categories entries; BC: 2 entries (alpha, beta)
  max_cols = max(nc[1], 2)
  mux = matrix(0, p, max_cols)
  mux[2, 1] = 0.5 # alpha
  mux[2, 2] = -0.3 # beta (penalizes distance from reference)

  pairwise_disc = matrix(c(0, 0.1, 0.1, 0), p, p)
  muy = 0
  pairwise_cont = matrix(-0.75)
  pairwise_cross = matrix(c(0.1, -0.1), p, q)

  result = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = c("ordinal", "blume-capel"),
    baseline_category_r = c(0L, 2L), iter = 200L, seed = 7L
  )

  expect_equal(dim(result$x), c(n, p))
  expect_true(all(result$x[, 1] >= 0 & result$x[, 1] <= nc[1]))
  expect_true(all(result$x[, 2] >= 0 & result$x[, 2] <= nc[2]))
})


# ==============================================================================
# 2. Parallel simulation (run_mixed_simulation_parallel)
# ==============================================================================

test_that("run_mixed_simulation_parallel returns correct structure", {
  p = 2L
  q = 1L
  ndraws_total = 3L
  nc = c(2L, 2L)
  mux_s = matrix(rep(c(0, 0.5, -0.3, 0, -0.2, 0.1), each = ndraws_total),
    nrow = ndraws_total
  )
  disc_s = matrix(0.15, nrow = ndraws_total, ncol = 1)
  muy_s = matrix(0, nrow = ndraws_total, ncol = 1)
  cont_s = matrix(-0.75, nrow = ndraws_total, ncol = 1)
  cross_s = matrix(c(0.1, -0.1), nrow = ndraws_total, ncol = 2, byrow = TRUE)

  n_use = 2L
  n_obs = 15L
  res = run_mixed_simulation_parallel(
    mux_samples = mux_s, disc_samples = disc_s,
    muy_samples = muy_s, cont_samples = cont_s, cross_samples = cross_s,
    draw_indices = 1:n_use, num_states = n_obs,
    p = p, q = q, num_categories = nc,
    variable_type_r = c("ordinal", "ordinal"),
    baseline_category = c(0L, 0L),
    iter = 100L, nThreads = 1L, seed = 1L, progress_type = 0L
  )

  expect_true(is.list(res))
  expect_equal(length(res), n_use)

  for(d in seq_len(n_use)) {
    expect_equal(dim(res[[d]]$x), c(n_obs, p))
    expect_equal(dim(res[[d]]$y), c(n_obs, q))
    expect_true(all(res[[d]]$x >= 0 & res[[d]]$x <= 2))
  }
})


# ==============================================================================
# 3. compute_conditional_mixed
# ==============================================================================

test_that("compute_conditional_mixed: discrete probs sum to 1", {
  p = 2
  q = 1
  n = 5
  nc = c(2L, 2L)
  mux = matrix(c(0, 0.5, -0.3, 0, -0.2, 0.1), nrow = p, ncol = 3, byrow = TRUE)
  pairwise_disc = matrix(c(0, 0.15, 0.15, 0), p, p)
  muy = 0
  pairwise_cont = matrix(-0.75)
  pairwise_cross = matrix(c(0.1, 0.2), p, q)

  x_obs = matrix(c(0L, 1L, 2L, 0L, 1L, 1L, 2L, 0L, 1L, 2L),
    nrow = n, ncol = p
  )
  y_obs = matrix(rnorm(n), nrow = n, ncol = q)

  # Predict first discrete variable (0-based index 0)
  preds = compute_conditional_mixed(
    x_observations = x_obs, y_observations = y_obs,
    predict_vars = 0L, pairwise_disc = pairwise_disc, pairwise_cross = pairwise_cross, pairwise_cont = pairwise_cont,
    mux = mux, muy = muy, num_categories = nc,
    variable_type = c("ordinal", "ordinal"),
    baseline_category = c(0L, 0L)
  )

  expect_equal(length(preds), 1)
  expect_equal(nrow(preds[[1]]), n)
  expect_equal(ncol(preds[[1]]), nc[1] + 1) # 3 categories

  row_sums = rowSums(preds[[1]])
  expect_true(all(abs(row_sums - 1) < 1e-8))
  expect_true(all(preds[[1]] >= 0))
})

test_that("compute_conditional_mixed: continuous returns mean and sd", {
  p = 1
  q = 2
  n = 5
  nc = c(2L)
  mux = matrix(c(0, 0.5, -0.3), nrow = 1, ncol = 3)
  pairwise_disc = matrix(0)
  muy = c(1.0, -0.5)
  pairwise_cont = matrix(c(-1, -0.15, -0.15, -0.75), q, q)
  pairwise_cross = matrix(c(0.1, -0.1), 1, q)

  x_obs = matrix(c(0L, 1L, 2L, 1L, 0L), nrow = n, ncol = 1)
  y_obs = matrix(rnorm(n * q), nrow = n, ncol = q)

  # Predict continuous variable at index p=1 (0-based)
  preds = compute_conditional_mixed(
    x_observations = x_obs, y_observations = y_obs,
    predict_vars = as.integer(p), pairwise_disc = pairwise_disc, pairwise_cross = pairwise_cross, pairwise_cont = pairwise_cont,
    mux = mux, muy = muy, num_categories = nc,
    variable_type = "ordinal", baseline_category = 0L
  )

  expect_equal(length(preds), 1)
  expect_equal(nrow(preds[[1]]), n)
  expect_equal(ncol(preds[[1]]), 2) # mean, sd
  expect_true(all(preds[[1]][, 2] > 0)) # sd > 0
})

test_that("compute_conditional_mixed: mixed prediction variables", {
  p = 2
  q = 2
  n = 5
  nc = c(2L, 2L)
  mux = matrix(0, p, 3)
  pairwise_disc = matrix(c(0, 0.15, 0.15, 0), p, p)
  muy = c(0.5, -0.2)
  pairwise_cont = diag(c(-0.75, -0.9))
  pairwise_cross = matrix(0.1, p, q)

  x_obs = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  storage.mode(x_obs) = "integer"
  y_obs = matrix(rnorm(n * q), n, q)

  # Predict one discrete (0) and one continuous (p+0 = 2)
  preds = compute_conditional_mixed(
    x_observations = x_obs, y_observations = y_obs,
    predict_vars = c(0L, 2L), pairwise_disc = pairwise_disc, pairwise_cross = pairwise_cross, pairwise_cont = pairwise_cont,
    mux = mux, muy = muy, num_categories = nc,
    variable_type = c("ordinal", "ordinal"),
    baseline_category = c(0L, 0L)
  )

  expect_equal(length(preds), 2)
  # First is discrete: n x 3
  expect_equal(ncol(preds[[1]]), nc[1] + 1)
  expect_true(all(abs(rowSums(preds[[1]]) - 1) < 1e-8))
  # Second is continuous: n x 2
  expect_equal(ncol(preds[[2]]), 2)
})


# ==============================================================================
# 4. simulate.bgms for mixed MRF (posterior-mean path)
# ==============================================================================

test_that("simulate.bgms works for mixed MRF with posterior-mean", {
  fit = get_bgms_fit_mixed_mrf_no_es()
  args = extract_arguments(fit)

  n_sim = 30
  result = simulate(fit, nsim = n_sim, method = "posterior-mean", seed = 1)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_sim)
  expect_equal(ncol(result), args$num_variables)
  expect_equal(colnames(result), args$data_columnnames)

  # Discrete columns should be non-negative integers
  for(di in args$discrete_indices) {
    vals = result[, di]
    expect_true(all(vals >= 0), info = sprintf("col %d negative", di))
    expect_true(
      all(vals == round(vals)),
      info = sprintf("col %d not integer", di)
    )
  }
})

test_that("simulate.bgms seed reproducibility for mixed MRF", {
  fit = get_bgms_fit_mixed_mrf_no_es()

  r1 = simulate(fit, nsim = 10, method = "posterior-mean", seed = 42)
  r2 = simulate(fit, nsim = 10, method = "posterior-mean", seed = 42)

  expect_equal(r1, r2)
})


# ==============================================================================
# 5. predict.bgms for mixed MRF (posterior-mean path)
# ==============================================================================

test_that("predict.bgms works for mixed MRF with posterior-mean", {
  fit = get_bgms_fit_mixed_mrf_no_es()
  args = extract_arguments(fit)

  newdata = get_prediction_data_mixed(n = 10)
  result = predict(fit, newdata = newdata, type = "probabilities")

  expect_true(is.list(result))
  expect_equal(length(result), args$num_variables)

  for(j in seq_len(args$num_variables)) {
    vname = args$data_columnnames[j]
    expect_equal(nrow(result[[j]]), 10, info = sprintf("var %s nrow", vname))
    expect_false(anyNA(result[[j]]), info = sprintf("var %s has NAs", vname))

    if(args$variable_type[j] %in% c("ordinal", "blume-capel")) {
      row_sums = rowSums(result[[j]])
      expect_true(
        all(abs(row_sums - 1) < 1e-6),
        info = sprintf("var %s probs don't sum to 1", vname)
      )
    } else {
      expect_equal(ncol(result[[j]]), 2, info = sprintf("var %s ncol", vname))
      expect_true(all(result[[j]][, 2] > 0),
        info = sprintf("var %s sd not positive", vname)
      )
    }
  }
})

test_that("predict.bgms response type works for mixed MRF", {
  fit = get_bgms_fit_mixed_mrf_no_es()
  args = extract_arguments(fit)

  newdata = get_prediction_data_mixed(n = 10)
  result = predict(fit, newdata = newdata, type = "response")

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(10, args$num_variables))

  # Discrete columns should be integer-valued
  for(di in args$discrete_indices) {
    expect_true(all(result[, di] == round(result[, di])),
      info = sprintf("response col %d not integer", di)
    )
  }
})


# ==============================================================================
# 6. Edge cases
# ==============================================================================

test_that("sample_mixed_mrf_gibbs: single discrete variable (p=1, q=2)", {
  result = sample_mixed_mrf_gibbs(
    num_states = 50L,
    pairwise_disc_r = matrix(0),
    pairwise_cross_r = matrix(c(0.1, -0.1), 1, 2),
    pairwise_cont_r = diag(c(-0.75, -1.0)),
    mux_r = matrix(c(0, 0.5), 1, 2),
    muy_r = c(0, 0),
    num_categories_r = 1L,
    variable_type_r = "ordinal",
    baseline_category_r = 0L,
    iter = 100L, seed = 3L
  )

  expect_equal(dim(result$x), c(50, 1))
  expect_equal(dim(result$y), c(50, 2))
  expect_true(all(result$x %in% 0:1))
})

test_that("sample_mixed_mrf_gibbs: single continuous variable (p=2, q=1)", {
  result = sample_mixed_mrf_gibbs(
    num_states = 50L,
    pairwise_disc_r = matrix(c(0, 0.1, 0.1, 0), 2, 2),
    pairwise_cross_r = matrix(c(0.1, -0.1), 2, 1),
    pairwise_cont_r = matrix(-1.0),
    mux_r = matrix(c(0, 0.5, -0.3, 0, -0.2, 0.1), 2, 3, byrow = TRUE),
    muy_r = 0.5,
    num_categories_r = c(2L, 2L),
    variable_type_r = c("ordinal", "ordinal"),
    baseline_category_r = c(0L, 0L),
    iter = 100L, seed = 5L
  )

  expect_equal(dim(result$x), c(50, 2))
  expect_equal(dim(result$y), c(50, 1))
})

test_that("sample_mixed_mrf_gibbs: minimal mixed MRF (p=1, q=1)", {
  p = 1L
  q = 1L
  n = 200L
  nc = 2L

  pairwise_disc = matrix(0, p, p)
  pairwise_cross = matrix(0.3, p, q)
  pairwise_cont = matrix(-0.75, q, q)
  mux = matrix(c(0, 0.5, -0.3), p, nc + 1)
  muy = 0.2

  result = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = "ordinal",
    baseline_category_r = 0L, iter = 200L, seed = 77L
  )

  expect_equal(dim(result$x), c(n, p))
  expect_equal(dim(result$y), c(n, q))
  expect_true(all(result$x >= 0 & result$x <= nc))
  expect_true(all(is.finite(result$y)))

  # With non-zero pairwise_cross, discrete and continuous should be associated
  group_means = tapply(result$y[, 1], result$x[, 1], mean)
  expect_true(length(group_means) >= 2,
    info = "at least 2 discrete categories observed"
  )
})


# ==============================================================================
# 7. Edge cases
# ==============================================================================
# Structural smoke tests for boundary configurations: many categories,
# large p / small q, small p / large q, near-singular pairwise_cont, and
# end-to-end bgm() fitting with edge-case data.
# ==============================================================================

# ---- many categories ---------------------------------------------------------
test_that("sample_mixed_mrf_gibbs handles many categories (4+)", {
  p = 2L
  q = 1L
  n = 200L
  nc = c(4L, 3L)

  pairwise_disc = matrix(c(0, 0.075, 0.075, 0), p, p)
  pairwise_cross = matrix(c(0.2, 0.15), p, q)
  pairwise_cont = matrix(-0.75, q, q)
  mux = matrix(0, p, max(nc))
  mux[1, 1:4] = c(-0.5, 0, 0.3, 0.8)
  mux[2, 1:3] = c(-0.3, 0.2, 0.6)
  muy = 0

  result = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = c("ordinal", "ordinal"),
    baseline_category_r = c(0L, 0L), iter = 200L, seed = 33L
  )

  expect_equal(dim(result$x), c(n, p))
  expect_equal(dim(result$y), c(n, q))

  # Discrete values in valid ranges
  expect_true(all(result$x[, 1] >= 0 & result$x[, 1] <= nc[1]))
  expect_true(all(result$x[, 2] >= 0 & result$x[, 2] <= nc[2]))

  # All categories should be observed with enough samples
  expect_true(length(unique(result$x[, 1])) >= 3,
    info = "variable 1: at least 3 of 5 categories observed"
  )
})

# ---- large p, small q -------------------------------------------------------
test_that("sample_mixed_mrf_gibbs handles large p, small q (p=10, q=1)", {
  p = 10L
  q = 1L
  n = 300L
  nc = as.integer(rep(1, p))

  # Sparse pairwise_disc
  pairwise_disc = matrix(0, p, p)
  pairwise_disc[1, 2] = pairwise_disc[2, 1] = 0.1
  pairwise_disc[3, 4] = pairwise_disc[4, 3] = 0.075

  # Sparse pairwise_cross
  pairwise_cross = matrix(0, p, q)
  pairwise_cross[1, 1] = 0.25
  pairwise_cross[5, 1] = 0.2

  pairwise_cont = matrix(-0.75, q, q)
  mux = matrix(0, p, max(nc))
  muy = 0

  result = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = rep("ordinal", p),
    baseline_category_r = rep(0L, p), iter = 200L, seed = 44L
  )

  expect_equal(dim(result$x), c(n, p))
  expect_equal(dim(result$y), c(n, q))
  expect_true(all(is.finite(result$x)))
  expect_true(all(is.finite(result$y)))
  expect_true(all(result$x >= 0 & result$x <= 1))
})

# ---- small p, large q -------------------------------------------------------
test_that("sample_mixed_mrf_gibbs handles small p, large q (p=1, q=5)", {
  p = 1L
  q = 5L
  n = 300L
  nc = 1L

  pairwise_disc = matrix(0, p, p)
  pairwise_cross = matrix(c(0.2, 0.15, 0.1, 0.05, 0.25), p, q)

  # Sparse pairwise_cont (association-scale: negative-definite; precision = -2 * pairwise_cont)
  pairwise_cont = diag(-0.75, q)
  pairwise_cont[1, 2] = pairwise_cont[2, 1] = -0.1
  pairwise_cont[3, 4] = pairwise_cont[4, 3] = -0.075

  mux = matrix(0, p, 1)
  muy = rep(0, q)

  result = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = "ordinal",
    baseline_category_r = 0L, iter = 300L, seed = 55L
  )

  expect_equal(dim(result$x), c(n, p))
  expect_equal(dim(result$y), c(n, q))
  expect_true(all(is.finite(result$y)))

  # pairwise_cont diagonal sets marginal variances: precision = -2 * pairwise_cont, var_j approx 1/precision_jj
  for(j in 1:q) {
    expected_sd = 1 / sqrt(-2 * pairwise_cont[j, j])
    empirical_sd = sd(result$y[, j])
    expect_true(
      abs(empirical_sd - expected_sd) / expected_sd < 0.4,
      info = sprintf("y%d SD: expected %.3f, got %.3f", j, expected_sd, empirical_sd)
    )
  }
})

# ---- near-singular pairwise_cont ------------------------------------------------------
test_that("sample_mixed_mrf_gibbs handles correlated pairwise_cont", {
  p = 1L
  q = 2L
  n = 200L
  nc = 1L

  pairwise_disc = matrix(0, p, p)
  pairwise_cross = matrix(c(0.2, 0.15), p, q)

  # Association-scale pairwise_cont with notable off-diagonal
  pairwise_cont = matrix(c(-0.75, -0.25, -0.25, -0.75), q, q)

  mux = matrix(0, p, 1)
  muy = c(0, 0)

  result = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = "ordinal",
    baseline_category_r = 0L, iter = 300L, seed = 66L
  )

  expect_equal(dim(result$y), c(n, q))
  expect_true(all(is.finite(result$y)))

  # With negative off-diagonal pairwise_cont (positive off-diagonal precision),
  # y1 and y2 should be negatively correlated
  r = cor(result$y[, 1], result$y[, 2])
  expect_true(r < 0.1, info = sprintf("expected negative or near-zero cor, got %.3f", r))
})

# ---- bgm() end-to-end with many categories ----------------------------------
test_that("bgm() fits mixed MRF with many categories", {
  skip_on_cran()
  p = 2L
  q = 1L
  n = 150L
  nc = c(4L, 3L)

  pairwise_disc = matrix(c(0, 0.075, 0.075, 0), p, p)
  pairwise_cross = matrix(c(0.2, 0.15), p, q)
  pairwise_cont = matrix(-0.75, q, q)
  mux = matrix(0, p, max(nc))
  mux[1, 1:4] = c(-0.5, 0, 0.3, 0.8)
  mux[2, 1:3] = c(-0.3, 0.2, 0.6)

  generated = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = 0, num_categories_r = nc,
    variable_type_r = c("ordinal", "ordinal"),
    baseline_category_r = c(0L, 0L), iter = 500L, seed = 77L
  )

  dat = data.frame(generated$x, generated$y)
  vt = c("ordinal", "ordinal", "continuous")

  fit = bgm(dat,
    variable_type = vt, iter = 500, warmup = 250,
    edge_selection = FALSE, verbose = FALSE, display_progress = "none"
  )

  expect_s3_class(fit, "bgms")
  pw = fit$posterior_mean_pairwise
  expect_equal(dim(pw), c(p + q, p + q))
  expect_true(all(is.finite(pw)))
})

# ---- bgm() end-to-end with large p, small q ---------------------------------
test_that("bgm() fits mixed MRF with large p, small q", {
  skip_on_cran()
  p = 6L
  q = 1L
  n = 200L
  nc = as.integer(rep(1, p))

  pairwise_disc = matrix(0, p, p)
  pairwise_disc[1, 2] = pairwise_disc[2, 1] = 0.1
  pairwise_cross = matrix(0, p, q)
  pairwise_cross[1, 1] = 0.15
  pairwise_cont = matrix(-0.75, q, q)
  mux = matrix(0, p, max(nc))

  generated = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = 0, num_categories_r = nc,
    variable_type_r = rep("ordinal", p),
    baseline_category_r = rep(0L, p), iter = 500L, seed = 88L
  )

  dat = data.frame(generated$x, generated$y)
  vt = c(rep("ordinal", p), "continuous")

  fit = bgm(dat,
    variable_type = vt, iter = 500, warmup = 250,
    edge_selection = FALSE, verbose = FALSE, display_progress = "none"
  )

  expect_s3_class(fit, "bgms")
  expect_equal(nrow(fit$posterior_mean_pairwise), p + q)
})

# ---- bgm() end-to-end with small p, large q ---------------------------------
test_that("bgm() fits mixed MRF with small p, large q", {
  skip_on_cran()
  p = 1L
  q = 4L
  n = 200L
  nc = 1L

  pairwise_disc = matrix(0, p, p)
  pairwise_cross = matrix(c(0.2, 0.1, 0.05, 0.15), p, q)
  pairwise_cont = diag(-0.75, q)
  pairwise_cont[1, 2] = pairwise_cont[2, 1] = -0.1
  mux = matrix(0, p, 1)

  generated = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc, pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = rep(0, q), num_categories_r = nc,
    variable_type_r = "ordinal",
    baseline_category_r = 0L, iter = 500L, seed = 99L
  )

  dat = data.frame(generated$x, generated$y)
  vt = c("ordinal", rep("continuous", q))

  fit = bgm(dat,
    variable_type = vt, iter = 500, warmup = 250,
    edge_selection = FALSE, verbose = FALSE, display_progress = "none"
  )

  expect_s3_class(fit, "bgms")
  expect_equal(nrow(fit$posterior_mean_pairwise), p + q)
})
