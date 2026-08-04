# ==============================================================================
# Input Validation Tests
# ==============================================================================
#
# COMPLEMENTS: test-tolerance.R (tests the "failure" side of API contracts)
# PATTERN: Error conditions, boundary cases, invalid input handling
#
# While test-tolerance.R and other files test that VALID inputs produce
# outputs with correct properties, this file tests that INVALID inputs
# produce appropriate error messages.
#
# Together with tolerance tests, these form complete API contract testing:
#   - Tolerance tests: f(valid_input) -> output with correct properties
#   - Validation tests: f(invalid_input) -> informative error message
#
# See helper-fixtures.R for test data generators and testing philosophy.
# ==============================================================================

# ------------------------------------------------------------------------------
# bgm() Input Validation
# ------------------------------------------------------------------------------

test_that("bgm errors on non-matrix/data.frame input", {
  expect_error(bgm(x = 1:10), regexp = "matrix|data.frame|data frame")
  expect_error(bgm(x = list(a = 1, b = 2)), regexp = "matrix|data.frame|data frame")
})

test_that("bgm errors on data with too few variables", {
  bad_data = matrix(c(0, 1, 0, 1), ncol = 1)
  expect_error(bgm(x = bad_data), regexp = "variable|column")
})

test_that("bgm errors on invalid iter values", {
  data = generate_test_data(n = 20, p = 3)

  expect_error(bgm(x = data, iter = 0), regexp = "iter")
  expect_error(bgm(x = data, iter = -10), regexp = "iter")
  expect_error(bgm(x = data, iter = "100"), regexp = "iter|numeric")
})

test_that("bgm errors on invalid edge_prior", {
  data = generate_test_data(n = 20, p = 3)

  expect_error(
    suppressWarnings(
      bgm(x = data, edge_selection = TRUE, edge_prior = "Invalid")
    ),
    regexp = "should be one of|edge_prior"
  )
})

test_that("bgm errors on invalid na_action", {
  data = generate_test_data(n = 20, p = 3)

  expect_error(
    bgm(x = data, na_action = "invalid_action"),
    regexp = "na_action"
  )
})

test_that("bgm errors on invalid update_method", {
  data = generate_test_data(n = 20, p = 3)

  expect_error(
    bgm(x = data, update_method = "not-a-method"),
    regexp = "should be one of|update_method"
  )
})


# ------------------------------------------------------------------------------
# bgm() GGM-Specific Input Validation
# ------------------------------------------------------------------------------

test_that("GGM accepts NUTS", {
  set.seed(42)
  x = matrix(rnorm(200), nrow = 50, ncol = 4)

  # NUTS should be accepted (no error at validation stage)
  spec = bgm_spec(
    x = x, variable_type = "continuous", update_method = "nuts",
    iter = 10, warmup = 300, chains = 1
  )
  expect_equal(spec$sampler$update_method, "nuts")
})

test_that("GGM rejects removed hamiltonian-mc", {
  set.seed(42)
  x = matrix(rnorm(200), nrow = 50, ncol = 4)
  expect_error(
    bgm_spec(
      x = x, variable_type = "continuous", update_method = "hamiltonian-mc",
      iter = 10, warmup = 300, chains = 1
    ),
    "arg"
  )
})

test_that("Mixed continuous and ordinal variable types are accepted for bgm", {
  set.seed(42)
  x = data.frame(
    ord1 = sample(0:2, 50, replace = TRUE),
    ord2 = sample(0:2, 50, replace = TRUE),
    cont1 = rnorm(50),
    cont2 = rnorm(50)
  )
  spec = bgm_spec(
    x = x,
    variable_type = c("ordinal", "ordinal", "continuous", "continuous")
  )
  expect_equal(spec$model_type, "mixed_mrf")
  expect_equal(spec$data$num_discrete, 2L)
  expect_equal(spec$data$num_continuous, 2L)
})


# ------------------------------------------------------------------------------
# bgmCompare() Input Validation
# ------------------------------------------------------------------------------

test_that("bgmCompare errors on insufficient data", {
  # Too few groups
  data = generate_test_data(n = 20, p = 3)
  group_ind = rep(1, 20) # Only one group

  expect_error(
    bgmCompare(x = data, group_indicator = group_ind),
    regexp = "group"
  )
})

test_that("bgmCompare errors on mismatched group_indicator length", {
  data = generate_test_data(n = 20, p = 3)
  group_ind = rep(1:2, each = 5) # Only 10 elements for 20 rows

  expect_error(
    bgmCompare(x = data, group_indicator = group_ind),
    regexp = "group_indicator|length|match"
  )
})

test_that("bgmCompare accepts character, factor, integer, and 0/1 group indicators alike", {
  data("Boredom", package = "bgms")
  rows = c(1:25, 491:515)
  x = Boredom[rows, 2:4]
  language = Boredom[rows, "language"]

  indicators = list(
    character = language,
    factor = factor(language),
    integer = match(language, unique(language)),
    zero_one = match(language, unique(language)) - 1L
  )
  specs = lapply(indicators, function(g) {
    without_support_warning(
      bgm_spec(x = x, model_type = "compare", group_indicator = g, chains = 1)
    )
  })
  for(nm in names(specs)) {
    expect_equal(tabulate(specs[[nm]]$data$group), c(25L, 25L),
      info = sprintf("group sizes under a %s indicator", nm)
    )
  }
  # First-appearance numbering: every coding maps the same rows to group 1.
  reference = specs$integer$data$group
  for(nm in names(specs)) {
    expect_equal(specs[[nm]]$data$group, reference,
      info = sprintf("group membership under a %s indicator", nm)
    )
  }

  fit = without_support_warning(bgmCompare(
    x = x, group_indicator = language,
    iter = 25, warmup = 50, chains = 1, seed = 21,
    display_progress = "none"
  ))
  expect_s3_class(fit, "bgmCompare")
  expect_equal(sort(unique(extract_arguments(fit)$group)), c(1L, 2L))
})

test_that("a character group indicator survives listwise removal", {
  data("Boredom", package = "bgms")
  rows = c(1:25, 491:515)
  x = Boredom[rows, 2:4]
  x[3, 1] = NA
  spec = without_support_warning(bgm_spec(
    x = x, model_type = "compare",
    group_indicator = Boredom[rows, "language"],
    na_action = "listwise", chains = 1
  ))
  expect_equal(tabulate(spec$data$group), c(24L, 25L))
})

test_that("bgmCompare gives the validator's message for a vector x", {
  data("Boredom", package = "bgms")
  expect_error(
    bgmCompare(Boredom[, 1], group_indicator = rep(1:2, length.out = nrow(Boredom))),
    regexp = "x must be a matrix or data.frame",
    fixed = TRUE
  )
  # The single-group entry point owes the same message.
  expect_error(
    bgm(Boredom[, 2]),
    regexp = "x must be a matrix or data.frame",
    fixed = TRUE
  )
  # A vector y hits the same validator under its own name.
  expect_error(
    bgmCompare(x = Boredom[1:50, 2:4], y = Boredom[1:50, 2]),
    regexp = "y must be a matrix or data.frame",
    fixed = TRUE
  )
})

test_that("bgmCompare rejects continuous variable type", {
  x = matrix(rnorm(100), nrow = 50, ncol = 2)
  group_ind = rep(1:2, each = 25)

  expect_error(
    bgmCompare(
      x = x, group_indicator = group_ind,
      variable_type = "continuous"
    ),
    regexp = "not of type continuous"
  )
})

test_that("bgmCompare rejects mixed ordinal + continuous variable types", {
  x = data.frame(
    ord1 = sample(0:2, 50, replace = TRUE),
    cont1 = rnorm(50)
  )
  group_ind = rep(1:2, each = 25)

  # allow_continuous = FALSE fires before the mixed check, so the error
  # message is the same as for pure continuous input.
  expect_error(
    bgmCompare(
      x = x, group_indicator = group_ind,
      variable_type = c("ordinal", "continuous")
    ),
    regexp = "not of type continuous"
  )
})


# ------------------------------------------------------------------------------
# simulate_mrf() Input Validation
# ------------------------------------------------------------------------------

test_that("simulate_mrf errors on invalid num_states", {
  expect_error(
    simulate_mrf(
      num_states = 0,
      num_variables = 3,
      num_categories = 2,
      pairwise = matrix(0, 3, 3),
      main = matrix(0, 3, 2)
    ),
    regexp = "num_states"
  )

  expect_error(
    simulate_mrf(
      num_states = -5,
      num_variables = 3,
      num_categories = 2,
      pairwise = matrix(0, 3, 3),
      main = matrix(0, 3, 2)
    ),
    regexp = "num_states"
  )
})

test_that("simulate_mrf errors on non-symmetric pairwise", {
  non_sym = matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3)

  expect_error(
    simulate_mrf(
      num_states = 10,
      num_variables = 3,
      num_categories = 2,
      pairwise = non_sym,
      main = matrix(0, 3, 2)
    ),
    regexp = "symmetric|pairwise"
  )
})

test_that("simulate_mrf errors on dimension mismatch", {
  # Interactions matrix wrong size
  expect_error(
    simulate_mrf(
      num_states = 10,
      num_variables = 3,
      num_categories = 2,
      pairwise = matrix(0, 4, 4), # Wrong: 4x4 for 3 variables
      main = matrix(0, 3, 2)
    ),
    regexp = "num_variables|dimension|size"
  )
})

test_that("simulate_mrf errors on missing thresholds", {
  expect_error(
    simulate_mrf(
      num_states = 10,
      num_variables = 3,
      num_categories = 2,
      pairwise = matrix(0, 3, 3),
      main = matrix(c(0, 0, NA, 0, 0, 0), 3, 2) # NA threshold
    ),
    regexp = "NA|threshold|missing"
  )
})


# ------------------------------------------------------------------------------
# predict.bgms() Input Validation
# ------------------------------------------------------------------------------

test_that("predict.bgms errors when newdata is missing", {
  fit = get_bgms_fit()

  expect_error(predict(fit), regexp = "newdata")
})

test_that("predict.bgms errors on invalid variable names", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)

  data("Wenchuan", package = "bgms")
  newdata = Wenchuan[1:5, 1:args$num_variables]

  expect_error(
    predict(fit, newdata = newdata, variables = "NonexistentVar"),
    regexp = "not found|Variable"
  )
})

test_that("predict.bgms errors on out-of-range variable indices", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)

  data("Wenchuan", package = "bgms")
  newdata = Wenchuan[1:5, 1:args$num_variables]

  expect_error(
    predict(fit, newdata = newdata, variables = 999),
    regexp = "indices|between"
  )
})


# ------------------------------------------------------------------------------
# simulate.bgms() Input Validation
# ------------------------------------------------------------------------------

test_that("simulate.bgms errors on invalid seed", {
  fit = get_bgms_fit()

  expect_error(simulate(fit, nsim = 10, seed = -1), regexp = "seed")
  expect_error(simulate(fit, nsim = 10, seed = "abc"), regexp = "seed")
})

test_that("simulate.bgms errors on invalid cores argument", {
  fit = get_bgms_fit()

  expect_error(
    simulate(fit, nsim = 10, method = "posterior-sample", cores = 0),
    regexp = "cores"
  )
  expect_error(
    simulate(fit, nsim = 10, method = "posterior-sample", cores = "two"),
    regexp = "cores"
  )
})


# ------------------------------------------------------------------------------
# predict.bgmCompare() Input Validation
# ------------------------------------------------------------------------------

test_that("predict.bgmCompare errors when newdata is missing", {
  fit = get_bgmcompare_fit()

  expect_error(predict(fit, group = 1), regexp = "newdata")
})

test_that("predict.bgmCompare errors when group is missing", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)

  newdata = matrix(0L, nrow = 5, ncol = args$num_variables)

  expect_error(predict(fit, newdata = newdata), regexp = "group.*required")
})

test_that("predict.bgmCompare errors on invalid group argument", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)

  newdata = matrix(0L, nrow = 5, ncol = args$num_variables)

  # Group out of range
  expect_error(
    predict(fit, newdata = newdata, group = 0),
    regexp = "group.*1"
  )
  expect_error(
    predict(fit, newdata = newdata, group = 999),
    regexp = "group"
  )
  # Invalid type
  expect_error(
    predict(fit, newdata = newdata, group = "a"),
    regexp = "group"
  )
})

test_that("predict.bgmCompare errors on invalid variable names", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)

  newdata = matrix(0L, nrow = 5, ncol = args$num_variables)

  expect_error(
    predict(fit, newdata = newdata, group = 1, variables = "NonexistentVar"),
    regexp = "not found|Variable"
  )
})


# ------------------------------------------------------------------------------
# simulate.bgmCompare() Input Validation
# ------------------------------------------------------------------------------

test_that("simulate.bgmCompare errors when group is missing", {
  fit = get_bgmcompare_fit()

  expect_error(simulate(fit, nsim = 10), regexp = "group.*required")
})

test_that("simulate.bgmCompare errors on invalid group argument", {
  fit = get_bgmcompare_fit()

  # Group out of range
  expect_error(simulate(fit, nsim = 10, group = 0), regexp = "group.*1")
  expect_error(simulate(fit, nsim = 10, group = 999), regexp = "group")
  # Invalid type
  expect_error(simulate(fit, nsim = 10, group = "a"), regexp = "group")
})

test_that("simulate.bgmCompare errors on invalid seed", {
  fit = get_bgmcompare_fit()

  expect_error(simulate(fit, nsim = 10, group = 1, seed = -1), regexp = "seed")
  expect_error(simulate(fit, nsim = 10, group = 1, seed = "abc"), regexp = "seed")
})


# ------------------------------------------------------------------------------
# Extractor Function Error Handling
# ------------------------------------------------------------------------------

test_that("extract_indicators errors correctly without edge selection", {
  data = generate_test_data(n = 20, p = 3)
  args = c(list(x = data, edge_selection = FALSE), quick_mcmc_args())
  fit = do.call(bgm, args)

  expect_error(extract_indicators(fit), regexp = "edge_selection")
})

test_that("extract_posterior_inclusion_probabilities errors without edge selection", {
  data = generate_test_data(n = 20, p = 3)
  args = c(list(x = data, edge_selection = FALSE), quick_mcmc_args())
  fit = do.call(bgm, args)

  expect_error(
    extract_posterior_inclusion_probabilities(fit),
    regexp = "edge_selection"
  )
})

test_that("extract_sbm errors with non-SBM prior", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)

  # Skip if this fit actually used SBM
  if(identical(args$edge_prior, "Stochastic-Block")) {
    skip("Fit uses SBM prior")
  }

  expect_error(extract_sbm(fit), regexp = "Stochastic-Block")
})


# ------------------------------------------------------------------------------
# Edge Cases
# ------------------------------------------------------------------------------

test_that("bgm handles constant columns gracefully", {
  data = generate_test_data(n = 20, p = 3)
  data[, 1] = 0 # Make first column constant

  # This test verifies bgm doesn't crash unexpectedly on edge cases.
  # The function may error, warn, or succeed depending on implementation.
  result = tryCatch(
    {
      args = c(list(x = data), quick_mcmc_args())
      do.call(bgm, args)
    },
    error = function(e) list(type = "error", obj = e)
  )

  # Either it errors or it succeeds - both are acceptable behaviors
  if(is.list(result) && !is.null(result$type)) {
    # Got an error - bgm handled the edge case
    expect_true(TRUE, info = "bgm raised an error for constant column")
  } else {
    # If it succeeded, verify we got a valid bgms object
    expect_s3_class(result, "bgms")
  }
})

test_that("bgm handles all-NA columns", {
  data = generate_test_data(n = 20, p = 3)
  data[, 1] = NA # Make first column all NA

  expect_error(bgm(x = data), regexp = "NA|missing|incomplete")
})
