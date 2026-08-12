# ==============================================================================
# Unit tests for validate_difference_prior()
# Phase A.4 of the R scaffolding refactor.
# ==============================================================================

# ==============================================================================
# 1. difference_selection = FALSE  <U+2192>  bypass
# ==============================================================================

test_that("difference_selection = FALSE returns Not applicable with 1x1 matrix", {
  result = validate_difference_prior(
    difference_selection = FALSE,
    num_variables        = 5
  )
  expect_false(result$difference_selection)
  expect_equal(result$difference_prior, "Not applicable")
  expect_equal(result$inclusion_probability_difference, matrix(0.5, 1, 1))
})

test_that("difference_selection NA errors", {
  expect_error(
    validate_difference_prior(difference_selection = NA, num_variables = 3),
    "difference_selection needs to be TRUE or FALSE"
  )
})

# ==============================================================================
# 2. Bernoulli <U+2014> scalar difference_probability
# ==============================================================================

test_that("Bernoulli with scalar 0.5 builds p x p matrix", {
  result = validate_difference_prior(
    difference_selection = TRUE,
    difference_prior = "Bernoulli",
    difference_probability = 0.5,
    num_variables = 4
  )
  expect_true(result$difference_selection)
  expect_equal(result$difference_prior, "Bernoulli")
  expect_equal(result$inclusion_probability_difference, matrix(0.5, 4, 4))
})

test_that("Bernoulli scalar: NA errors", {
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Bernoulli",
      difference_probability = NA, num_variables = 3
    ),
    "no value specified"
  )
})

test_that("Bernoulli scalar: zero errors", {
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Bernoulli",
      difference_probability = 0, num_variables = 3
    ),
    "needs to be positive"
  )
})

test_that("Bernoulli scalar: >= 1 errors", {
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Bernoulli",
      difference_probability = 1, num_variables = 3
    ),
    "cannot equal or exceed"
  )
})

# ==============================================================================
# 3. Bernoulli <U+2014> matrix difference_probability
# ==============================================================================

test_that("Bernoulli with symmetric matrix accepted", {
  mat = matrix(0.3, 3, 3)
  diag(mat) = 0.5
  result = validate_difference_prior(
    difference_selection = TRUE,
    difference_prior = "Bernoulli",
    difference_probability = mat,
    num_variables = 3
  )
  expect_equal(result$inclusion_probability_difference, mat)
})

test_that("Bernoulli matrix: non-symmetric errors", {
  mat = matrix(c(0.5, 0.3, 0.4, 0.5), 2, 2)
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Bernoulli",
      difference_probability = mat, num_variables = 2
    ),
    "symmetric"
  )
})

test_that("Bernoulli matrix: wrong dimension errors", {
  mat = matrix(0.3, 2, 2)
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Bernoulli",
      difference_probability = mat, num_variables = 4
    ),
    "as many rows"
  )
})

test_that("Bernoulli matrix: NA in lower tri (diag=TRUE) errors", {
  mat = matrix(0.3, 3, 3)
  mat[2, 2] = NA
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Bernoulli",
      difference_probability = mat, num_variables = 3
    ),
    "not specified"
  )
})

test_that("Bernoulli matrix: zero in lower tri errors", {
  mat = matrix(0.3, 3, 3)
  mat[2, 1] = 0
  mat[1, 2] = 0
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Bernoulli",
      difference_probability = mat, num_variables = 3
    ),
    "negative or zero"
  )
})

test_that("Bernoulli matrix: value >= 1 in lower tri errors", {
  mat = matrix(0.3, 3, 3)
  mat[3, 3] = 1.0
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Bernoulli",
      difference_probability = mat, num_variables = 3
    ),
    "one or larger"
  )
})

test_that("Bernoulli with vector (not matrix) errors", {
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Bernoulli",
      difference_probability = c(0.3, 0.4, 0.5), num_variables = 3
    ),
    "single number, matrix, or dataframe"
  )
})

# ==============================================================================
# 4. Beta-Bernoulli
# ==============================================================================

test_that("Beta-Bernoulli with valid params returns 0.5 matrix", {
  result = validate_difference_prior(
    difference_selection = TRUE,
    difference_prior     = "Beta-Bernoulli",
    num_variables        = 4,
    beta_bernoulli_alpha = 2,
    beta_bernoulli_beta  = 2
  )
  expect_equal(result$difference_prior, "Beta-Bernoulli")
  expect_equal(result$inclusion_probability_difference, matrix(0.5, 4, 4))
})

test_that("Beta-Bernoulli: NA alpha errors (check ordering)", {
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Beta-Bernoulli",
      num_variables = 3, beta_bernoulli_alpha = NA, beta_bernoulli_beta = 1
    ),
    "need to be specified"
  )
})

test_that("Beta-Bernoulli: non-positive beta errors", {
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Beta-Bernoulli",
      num_variables = 3, beta_bernoulli_alpha = 1, beta_bernoulli_beta = 0
    ),
    "need to be positive"
  )
})

test_that("Beta-Bernoulli: infinite alpha errors", {
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Beta-Bernoulli",
      num_variables = 3, beta_bernoulli_alpha = Inf, beta_bernoulli_beta = 1
    ),
    "need to be finite"
  )
})

# ==============================================================================
# 5. Invalid difference_prior
# ==============================================================================

test_that("Invalid difference_prior string errors", {
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Not-A-Real-Prior",
      num_variables = 3
    ),
    "arg"
  )
})

test_that("Stochastic-Block: returns SBM hyperparameters", {
  out = validate_difference_prior(
    difference_selection = TRUE, difference_prior = "Stochastic-Block",
    num_variables = 4,
    beta_bernoulli_alpha = 1, beta_bernoulli_beta = 1,
    beta_bernoulli_alpha_between = 2, beta_bernoulli_beta_between = 3,
    dirichlet_alpha = 0.5, lambda = 1
  )
  expect_equal(out$difference_prior, "Stochastic-Block")
  expect_equal(out$beta_bernoulli_alpha_between, 2)
  expect_equal(out$beta_bernoulli_beta_between, 3)
  expect_equal(out$dirichlet_alpha, 0.5)
  expect_equal(out$lambda, 1)
  expect_equal(dim(out$inclusion_probability_difference), c(4, 4))
})

test_that("Stochastic-Block: rejects non-positive hyperparameters", {
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Stochastic-Block",
      num_variables = 3, dirichlet_alpha = 0
    ),
    "positive"
  )
  expect_error(
    validate_difference_prior(
      difference_selection = TRUE, difference_prior = "Stochastic-Block",
      num_variables = 3, lambda = -1
    ),
    "positive"
  )
})
