# ==============================================================================
# Tests for simulate_mrf() - Standalone MRF Simulation
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Range invariants, reproducibility, distributional properties
#
# Tests for the standalone simulation function that generates data from
# a Markov Random Field with user-specified parameters.
#
# Tolerance testing principles applied here:
#   - Range invariants: simulated values in [0, n_cats], integers only
#   - Reproducibility: seed produces identical results
#   - Coarse distributional properties: positive interactions tend to
#     produce positive correlations, thresholds shift marginal distributions
#   - Dimension consistency: output has correct n_states x n_vars shape
#
# Note: We test *tendencies* rather than exact values because simulation
# output is stochastic. See test-tolerance.R for the foundational approach.
# ==============================================================================

# ------------------------------------------------------------------------------
# Basic Functionality Tests
# ------------------------------------------------------------------------------

test_that("simulate_mrf returns matrix of correct dimensions", {
  n_states = 50
  n_vars = 4
  n_cats = 2

  interactions = matrix(0, n_vars, n_vars)
  interactions[1, 2] = interactions[2, 1] = 0.15
  thresholds = matrix(0, n_vars, n_cats)

  result = simulate_mrf(
    num_states = n_states,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 123
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_states)
  expect_equal(ncol(result), n_vars)
})

test_that("simulate_mrf produces values in valid range", {
  n_cats = 3
  n_vars = 5

  interactions = matrix(0, n_vars, n_vars)
  thresholds = matrix(0, n_vars, n_cats)

  result = simulate_mrf(
    num_states = 100,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 42
  )

  # All values should be integers
  expect_true(all(result == round(result)))

  # Values should be in 0 to n_cats range
  expect_true(all(result >= 0))
  expect_true(all(result <= n_cats))
})

test_that("simulate_mrf is reproducible with seed", {
  n_vars = 3
  n_cats = 2

  interactions = matrix(0.1, n_vars, n_vars)
  diag(interactions) = 0
  thresholds = matrix(0.5, n_vars, n_cats)

  result1 = simulate_mrf(
    num_states = 50,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 999
  )

  result2 = simulate_mrf(
    num_states = 50,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 999
  )

  expect_equal(result1, result2)
})


# ------------------------------------------------------------------------------
# Variable Category Tests
# ------------------------------------------------------------------------------

test_that("simulate_mrf handles different categories per variable", {
  n_vars = 4
  n_cats = c(1, 2, 3, 4) # Different number of categories

  interactions = matrix(0, n_vars, n_vars)
  thresholds = matrix(0, n_vars, max(n_cats))

  # Warnings expected: variables with fewer categories have extra threshold columns
  result = suppressWarnings(simulate_mrf(
    num_states = 100,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 123
  ))

  # Check each variable's range
  for(j in 1:n_vars) {
    expect_true(
      all(result[, j] <= n_cats[j]),
      info = sprintf("Variable %d should have max value %d", j, n_cats[j])
    )
  }
})

test_that("simulate_mrf handles binary variables correctly", {
  n_vars = 3
  n_cats = 1 # Binary: 0 or 1

  interactions = matrix(0.15, n_vars, n_vars)
  diag(interactions) = 0
  thresholds = matrix(0, n_vars, n_cats)

  result = simulate_mrf(
    num_states = 100,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 42
  )

  # Only 0 and 1 should appear
  expect_true(all(result %in% c(0, 1)))
})


# ------------------------------------------------------------------------------
# Blume-Capel Model Tests
# ------------------------------------------------------------------------------

test_that("simulate_mrf works with blume-capel variables", {
  n_vars = 3
  n_cats = 4 # Need > 2 categories for Blume-Capel

  interactions = matrix(0, n_vars, n_vars)
  # Blume-Capel thresholds: alpha and beta
  thresholds = matrix(NA, n_vars, n_cats)
  thresholds[, 1] = 0 # alpha
  thresholds[, 2] = -0.5 # beta

  result = simulate_mrf(
    num_states = 100,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    variable_type = "blume-capel",
    baseline_category = 2,
    seed = 123
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 100)
  expect_true(all(result >= 0 & result <= n_cats))
})

test_that("simulate_mrf works with mixed ordinal and blume-capel", {
  n_vars = 4
  n_cats = 4

  interactions = matrix(0, n_vars, n_vars)
  thresholds = matrix(0, n_vars, n_cats)
  # Fill ordinal thresholds for vars 1 and 2
  thresholds[1:2, ] = c(0, 0.3, 0.6, 0.9)
  # Blume-Capel params for vars 3 and 4
  thresholds[3:4, 1:2] = cbind(c(0, 0), c(-0.3, -0.5))

  # Warnings expected: Blume-Capel variables only use 2 threshold columns
  result = suppressWarnings(simulate_mrf(
    num_states = 80,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    variable_type = c("ordinal", "ordinal", "blume-capel", "blume-capel"),
    baseline_category = c(0, 0, 2, 2), # baseline only matters for BC vars
    seed = 42
  ))

  expect_equal(nrow(result), 80)
  expect_equal(ncol(result), n_vars)
})


# ------------------------------------------------------------------------------
# Parameter Effect Tests
# ------------------------------------------------------------------------------

test_that("simulate_mrf: positive interaction tends to align responses", {
  # This is a weaker test that checks the basic behavior of interactions
  # Strong interactions can sometimes collapse variance, so we use moderate values
  n_vars = 2
  n_cats = 2

  # Moderate positive interaction
  pos_int = matrix(c(0, 0.4, 0.4, 0), 2, 2)

  # Use spread-out thresholds to ensure variance in responses
  thresholds = matrix(c(0, 0.5, 0, 0.5), n_vars, n_cats, byrow = TRUE)

  result = simulate_mrf(
    num_states = 1000,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = pos_int,
    main = thresholds,
    iter = 2000,
    seed = 456
  )

  # Check that both variables have variance
  var1 = var(result[, 1])
  var2 = var(result[, 2])

  # With moderate interaction and spread thresholds, should have some variance
  expect_true(var1 > 0 && var2 > 0, info = "Variables should have non-zero variance")

  # Check the correlation is positive (as expected with positive interaction)
  if(var1 > 0 && var2 > 0) {
    cor_val = cor(result[, 1], result[, 2])
    # With positive interaction, correlation should be non-negative
    # (allowing some tolerance for stochastic variation)
    expect_true(
      cor_val > -0.2,
      info = sprintf("Positive interaction should yield non-negative correlation, got: %.3f", cor_val)
    )
  }
})

test_that("simulate_mrf: threshold affects marginal distribution", {
  n_vars = 1
  n_cats = 1 # Binary

  interactions = matrix(0, 1, 1)

  # Positive threshold -> more likely to be 1
  # In MRF parameterization, positive threshold increases log-odds of category 1
  pos_thresh = matrix(3, 1, 1)
  # Negative threshold -> more likely to be 0
  neg_thresh = matrix(-3, 1, 1)

  result_pos = simulate_mrf(
    num_states = 500,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = pos_thresh,
    iter = 1000,
    seed = 42
  )

  result_neg = simulate_mrf(
    num_states = 500,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = neg_thresh,
    iter = 1000,
    seed = 42
  )

  prop_pos = mean(result_pos == 1)
  prop_neg = mean(result_neg == 1)

  # Positive threshold should have higher proportion of 1s
  expect_true(
    prop_pos > prop_neg,
    info = sprintf("Positive threshold prop(1)=%.3f should exceed negative threshold prop(1)=%.3f", prop_pos, prop_neg)
  )
})


# ------------------------------------------------------------------------------
# Deprecated mrfSampler() Test
# ------------------------------------------------------------------------------

test_that("mrfSampler works with deprecation warning", {
  n_vars = 3
  n_cats = 2

  interactions = matrix(0, n_vars, n_vars)
  thresholds = matrix(0, n_vars, n_cats)

  expect_warning(
    result <- mrfSampler(
      num_states = 20,
      num_variables = n_vars,
      num_categories = n_cats,
      pairwise = interactions,
      main = thresholds
    ),
    regexp = "deprecated"
  )

  expect_true(is.matrix(result))
})

test_that("mrfSampler produces identical results to simulate_mrf", {
  args = list(
    num_states = 50,
    num_variables = 4,
    num_categories = 3,
    pairwise = matrix(0.05, 4, 4) - diag(0.05, 4),
    main = matrix(0, 4, 3),
    iter = 100,
    seed = 999
  )

  result_new = do.call(simulate_mrf, args)
  result_old = suppressWarnings(do.call(mrfSampler, args))

  expect_identical(result_new, result_old)
})


# ==============================================================================
# simulate_mrf with continuous (GGM) variables
# ==============================================================================

test_that("simulate_mrf works with variable_type = 'continuous'", {
  p = 4
  n = 100

  # Precision matrix (diagonal dominant for PD)
  omega = diag(p)
  omega[1, 2] = omega[2, 1] = 0.3
  omega[3, 4] = omega[4, 3] = -0.2

  result = simulate_mrf(
    num_states = n,
    num_variables = p,
    pairwise = omega,
    variable_type = "continuous",
    seed = 42
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), p)
  expect_true(all(is.finite(result)))
  expect_true(is.numeric(result))
})

test_that("simulate_mrf continuous is reproducible with seed", {
  p = 3
  omega = diag(p)
  omega[1, 2] = omega[2, 1] = 0.5

  sim1 = simulate_mrf(
    num_states = 50,
    num_variables = p,
    pairwise = omega,
    variable_type = "continuous",
    seed = 123
  )
  sim2 = simulate_mrf(
    num_states = 50,
    num_variables = p,
    pairwise = omega,
    variable_type = "continuous",
    seed = 123
  )

  expect_equal(sim1, sim2)
})

test_that("simulate_mrf continuous: sample covariance approaches true covariance", {
  p = 3
  n = 5000

  # Known precision matrix
  omega = matrix(c(
    2.0,  0.5, 0.0,
    0.5,  1.5, 0.3,
    0.0,  0.3, 1.0
  ), nrow = p, byrow = TRUE)

  result = simulate_mrf(
    num_states = n,
    num_variables = p,
    pairwise = omega,
    variable_type = "continuous",
    seed = 42
  )

  # True covariance
  sigma_true = solve(omega)

  # Sample covariance should be close
  sigma_hat = cov(result)

  # Off-diagonal elements should correlate highly
  expect_true(
    cor(sigma_true[lower.tri(sigma_true)], sigma_hat[lower.tri(sigma_hat)]) > 0.95,
    info = "Sample covariance should track true covariance"
  )

  # Diagonal elements should be close
  for(j in 1:p) {
    expect_equal(sigma_hat[j, j], sigma_true[j, j],
      tolerance = 0.1,
      info = sprintf("Variance of variable %d", j)
    )
  }
})

test_that("simulate_mrf continuous rejects non-positive diagonal", {
  p = 3
  omega = diag(c(1, -1, 1)) # negative diagonal element

  expect_error(
    simulate_mrf(
      num_states = 10,
      num_variables = p,
      pairwise = omega,
      variable_type = "continuous"
    ),
    "positive"
  )
})

test_that("simulate_mrf rejects mixed continuous and ordinal", {
  expect_error(
    simulate_mrf(
      num_states = 10,
      num_variables = 3,
      num_categories = c(2, 2, 2),
      pairwise = matrix(0, 3, 3),
      main = matrix(0, 3, 2),
      variable_type = c("continuous", "ordinal", "ordinal")
    ),
    "all variables must be of type"
  )
})


# ------------------------------------------------------------------------------
# Numerical stability of the full conditional
# ------------------------------------------------------------------------------
# The full conditional subtracts the largest category exponent before
# exponentiating. Without that, a large rest score sends exp() to Inf, the
# cumulative weights become Inf or NaN, and the category search degenerates.

test_that("large parameters do not overflow the full conditional", {
  num_categories = c(2L, 2L, 3L)
  pairwise = matrix(0, 3, 3)
  pairwise[lower.tri(pairwise)] = c(80, -60, 50)
  pairwise = pairwise + t(pairwise)
  main = matrix(NA_real_, 3, max(num_categories))
  main[1, 1:2] = c(60, -40)
  main[2, 1:2] = c(-100, 80)
  main[3, 1:3] = c(20, 40, -80)

  result = simulate_mrf(
    num_states = 40, num_variables = 3, num_categories = num_categories,
    pairwise = pairwise, main = main,
    variable_type = rep("ordinal", 3), baseline_category = rep(0L, 3),
    iter = 50, seed = 3
  )

  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))
  expect_true(all(result <= rep(num_categories, each = nrow(result))))
})

test_that("the stabilized ordinal conditional samples the intended law", {
  skip_on_cran()
  # Three variables small enough to enumerate the joint exactly, so the
  # simulated marginals have something other than themselves to agree with.
  num_categories = c(2L, 2L, 3L)
  pairwise = matrix(0, 3, 3)
  pairwise[lower.tri(pairwise)] = c(0.4, -0.3, 0.25)
  pairwise = pairwise + t(pairwise)
  main = matrix(NA_real_, 3, max(num_categories))
  main[1, 1:2] = c(0.3, -0.2)
  main[2, 1:2] = c(-0.5, 0.4)
  main[3, 1:3] = c(0.1, 0.2, -0.4)

  simulated = simulate_mrf(
    num_states = 50000, num_variables = 3, num_categories = num_categories,
    pairwise = pairwise, main = main,
    variable_type = rep("ordinal", 3), baseline_category = rep(0L, 3),
    iter = 100, seed = 99
  )

  grid = as.matrix(expand.grid(lapply(num_categories, function(k) 0:k)))
  log_weight = apply(grid, 1, function(z) {
    total = 0
    for(v in seq_len(3)) if(z[v] > 0) total = total + main[v, z[v]]
    for(i in 1:2) {
      for(j in (i + 1):3) total = total + 2 * z[i] * z[j] * pairwise[i, j]
    }
    total
  })
  weight = exp(log_weight - max(log_weight))
  weight = weight / sum(weight)

  for(v in seq_len(3)) {
    exact = as.numeric(tapply(weight, grid[, v], sum))
    empirical = as.numeric(
      table(factor(simulated[, v], levels = 0:num_categories[v]))
    ) / nrow(simulated)
    expect_lt(max(abs(empirical - exact)), 0.01)
  }
})
