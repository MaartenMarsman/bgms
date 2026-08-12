# ==============================================================================
# Unit tests for validate_edge_prior()
# Phase A.3 of the R scaffolding refactor.
# ==============================================================================

# ==============================================================================
# 1. edge_selection = FALSE  <U+2192>  bypass
# ==============================================================================

test_that("edge_selection = FALSE returns Not Applicable with 1x1 theta", {
  result = validate_edge_prior(
    edge_selection = FALSE,
    num_variables  = 5
  )
  expect_false(result$edge_selection)
  expect_equal(result$edge_prior, "Not Applicable")
  expect_equal(result$inclusion_probability, matrix(0.5, 1, 1))
})

test_that("edge_selection coerced to logical; NA errors", {
  expect_error(
    validate_edge_prior(edge_selection = NA, num_variables = 3),
    "edge_selection needs to be TRUE or FALSE"
  )
})

# ==============================================================================
# 2. Bernoulli <U+2014> scalar inclusion_probability
# ==============================================================================

test_that("Bernoulli with scalar 0.5 builds p x p matrix", {
  result = validate_edge_prior(
    edge_selection        = TRUE,
    edge_prior            = "Bernoulli",
    inclusion_probability = 0.5,
    num_variables         = 4
  )
  expect_true(result$edge_selection)
  expect_equal(result$edge_prior, "Bernoulli")
  expect_equal(result$inclusion_probability, matrix(0.5, 4, 4))
})

test_that("Bernoulli scalar: NA errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = NA, num_variables = 3
    ),
    "no value specified"
  )
})

test_that("Bernoulli scalar: zero errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = 0, num_variables = 3
    ),
    "needs to be positive"
  )
})

test_that("Bernoulli scalar: negative errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = -0.1, num_variables = 3
    ),
    "needs to be positive"
  )
})

test_that("Bernoulli scalar: > 1 errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = 1.5, num_variables = 3
    ),
    "cannot equal or exceed the value one"
  )
})

test_that("Bernoulli scalar: exactly 1 errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = 1, num_variables = 3
    ),
    "cannot equal or exceed the value one"
  )
})

# ==============================================================================
# 3. Bernoulli <U+2014> matrix inclusion_probability
# ==============================================================================

test_that("Bernoulli with symmetric matrix accepted", {
  mat = matrix(0.3, 3, 3)
  diag(mat) = 0.5
  result = validate_edge_prior(
    edge_selection        = TRUE,
    edge_prior            = "Bernoulli",
    inclusion_probability = mat,
    num_variables         = 3
  )
  expect_equal(result$inclusion_probability, mat)
})

test_that("Bernoulli with data.frame coerced", {
  mat = matrix(0.3, 3, 3)
  diag(mat) = 0.5
  rownames(mat) = colnames(mat) = paste0("V", 1:3)
  df = as.data.frame(mat)
  result = validate_edge_prior(
    edge_selection        = TRUE,
    edge_prior            = "Bernoulli",
    inclusion_probability = df,
    num_variables         = 3
  )
  expect_equal(result$inclusion_probability, data.matrix(df))
})

test_that("Bernoulli matrix: non-symmetric errors", {
  mat = matrix(c(0.5, 0.3, 0.4, 0.5), 2, 2)
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = mat, num_variables = 2
    ),
    "symmetric"
  )
})

test_that("Bernoulli matrix: wrong dimension errors", {
  mat = matrix(0.3, 2, 2)
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = mat, num_variables = 4
    ),
    "as many rows"
  )
})

test_that("Bernoulli matrix: NA in lower triangle errors", {
  mat = matrix(0.3, 3, 3)
  mat[3, 1] = NA
  mat[1, 3] = NA
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = mat, num_variables = 3
    ),
    "not specified"
  )
})

test_that("Bernoulli matrix: zero in lower triangle errors", {
  mat = matrix(0.3, 3, 3)
  mat[2, 1] = 0
  mat[1, 2] = 0
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = mat, num_variables = 3
    ),
    "negative or zero"
  )
})

test_that("Bernoulli matrix: value >= 1 in lower triangle errors", {
  mat = matrix(0.3, 3, 3)
  mat[2, 1] = 1.0
  mat[1, 2] = 1.0
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = mat, num_variables = 3
    ),
    "greater than or equal to one"
  )
})

test_that("Bernoulli with vector (not matrix) errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Bernoulli",
      inclusion_probability = c(0.3, 0.4, 0.5), num_variables = 3
    ),
    "single number, matrix, or dataframe"
  )
})

# ==============================================================================
# 4. Beta-Bernoulli
# ==============================================================================

test_that("Beta-Bernoulli with valid params returns 0.5 matrix", {
  result = validate_edge_prior(
    edge_selection       = TRUE,
    edge_prior           = "Beta-Bernoulli",
    num_variables        = 4,
    beta_bernoulli_alpha = 1,
    beta_bernoulli_beta  = 1
  )
  expect_equal(result$edge_prior, "Beta-Bernoulli")
  expect_equal(result$inclusion_probability, matrix(0.5, 4, 4))
})

test_that("Beta-Bernoulli: non-positive alpha errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Beta-Bernoulli",
      num_variables = 3, beta_bernoulli_alpha = 0, beta_bernoulli_beta = 1
    ),
    "need to be positive"
  )
})

test_that("Beta-Bernoulli: infinite beta errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Beta-Bernoulli",
      num_variables = 3, beta_bernoulli_alpha = 1, beta_bernoulli_beta = Inf
    ),
    "need to be finite"
  )
})

test_that("Beta-Bernoulli: NA alpha errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Beta-Bernoulli",
      num_variables = 3, beta_bernoulli_alpha = NA, beta_bernoulli_beta = 1
    ),
    "need to be specified"
  )
})

# ==============================================================================
# 5. Stochastic-Block
# ==============================================================================

test_that("Stochastic-Block with valid params succeeds", {
  result = validate_edge_prior(
    edge_selection               = TRUE,
    edge_prior                   = "Stochastic-Block",
    num_variables                = 3,
    beta_bernoulli_alpha         = 1,
    beta_bernoulli_beta          = 1,
    beta_bernoulli_alpha_between = 1,
    beta_bernoulli_beta_between  = 1,
    dirichlet_alpha              = 1,
    lambda                       = 1
  )
  expect_equal(result$edge_prior, "Stochastic-Block")
  expect_equal(result$inclusion_probability, matrix(0.5, 3, 3))
})

test_that("Stochastic-Block: NULL between-params errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Stochastic-Block",
      num_variables = 3,
      beta_bernoulli_alpha = 1, beta_bernoulli_beta = 1,
      beta_bernoulli_alpha_between = NULL, beta_bernoulli_beta_between = 1,
      dirichlet_alpha = 1, lambda = 1
    ),
    "requires all four beta parameters"
  )
})

test_that("Stochastic-Block: non-positive dirichlet_alpha errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Stochastic-Block",
      num_variables = 3,
      beta_bernoulli_alpha = 1, beta_bernoulli_beta = 1,
      beta_bernoulli_alpha_between = 1, beta_bernoulli_beta_between = 1,
      dirichlet_alpha = 0, lambda = 1
    ),
    "need to be positive"
  )
})

test_that("Stochastic-Block: non-positive lambda errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Stochastic-Block",
      num_variables = 3,
      beta_bernoulli_alpha = 1, beta_bernoulli_beta = 1,
      beta_bernoulli_alpha_between = 1, beta_bernoulli_beta_between = 1,
      dirichlet_alpha = 1, lambda = -1
    ),
    "need to be positive"
  )
})

test_that("Stochastic-Block: infinite between-beta errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Stochastic-Block",
      num_variables = 3,
      beta_bernoulli_alpha = 1, beta_bernoulli_beta = 1,
      beta_bernoulli_alpha_between = Inf, beta_bernoulli_beta_between = 1,
      dirichlet_alpha = 1, lambda = 1
    ),
    "need to be finite"
  )
})

test_that("Stochastic-Block: NA parameters error", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "Stochastic-Block",
      num_variables = 3,
      beta_bernoulli_alpha = 1, beta_bernoulli_beta = 1,
      beta_bernoulli_alpha_between = 1, beta_bernoulli_beta_between = NA,
      dirichlet_alpha = 1, lambda = 1
    ),
    "cannot be NA"
  )
})

# ==============================================================================
# 6. edge_prior match.arg
# ==============================================================================

test_that("Invalid edge_prior string errors", {
  expect_error(
    validate_edge_prior(
      edge_selection = TRUE, edge_prior = "InvalidPrior",
      num_variables = 3
    ),
    "arg"
  )
})
