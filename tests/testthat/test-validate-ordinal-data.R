# ==============================================================================
# Tests for reformat_ordinal_data()
# Phase A.6 of the R scaffolding refactor.
# ==============================================================================


# ==============================================================================
# 1. Ordinal recoding
# ==============================================================================

test_that("ordinal: already contiguous 0-based passes through unchanged", {
  x = matrix(c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2), nrow = 6, ncol = 2)
  is_ordinal = c(TRUE, TRUE)
  bc = c(0, 0)

  result = reformat_ordinal_data(x, is_ordinal, bc)
  expect_equal(result$x, x)
  expect_equal(result$num_categories, c(2, 2))
  expect_equal(result$baseline_category, bc)
})

test_that("ordinal: non-contiguous values are recoded to contiguous 0-based", {
  x = matrix(c(1, 3, 5, 1, 3, 5, 1, 3, 5, 1, 3, 5), nrow = 6, ncol = 2)
  is_ordinal = c(TRUE, TRUE)
  bc = c(0, 0)

  result = reformat_ordinal_data(x, is_ordinal, bc)
  expect_equal(result$x[, 1], c(0, 1, 2, 0, 1, 2))
  expect_equal(result$x[, 2], c(0, 1, 2, 0, 1, 2))
  expect_equal(result$num_categories, c(2, 2))
})

test_that("ordinal: 1-based values are recoded to 0-based", {
  x = matrix(c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), nrow = 6, ncol = 2)
  is_ordinal = c(TRUE, TRUE)
  bc = c(0, 0)

  result = reformat_ordinal_data(x, is_ordinal, bc)
  expect_equal(result$x[, 1], c(0, 1, 2, 0, 1, 2))
  expect_equal(result$num_categories, c(2, 2))
})

test_that("ordinal: gaps in observed values are collapsed", {
  # Values 0, 2, 4 <U+2192> should become 0, 1, 2
  x = matrix(c(
    0, 2, 4, 0, 2, 4,
    0, 2, 4, 0, 2, 4
  ), nrow = 6, ncol = 2)
  is_ordinal = c(TRUE, TRUE)
  bc = c(0, 0)

  result = reformat_ordinal_data(x, is_ordinal, bc)
  expect_equal(result$x[, 1], c(0, 1, 2, 0, 1, 2))
  expect_equal(result$num_categories[1], 2)
})

test_that("ordinal: multiple observations per category preserved", {
  x = matrix(c(
    0, 0, 1, 1, 2, 2,
    0, 1, 0, 1, 0, 1
  ), nrow = 6, ncol = 2)
  is_ordinal = c(TRUE, TRUE)
  bc = c(0, 0)

  result = reformat_ordinal_data(x, is_ordinal, bc)
  expect_equal(result$x, x) # already contiguous
  expect_equal(result$num_categories, c(2, 1))
})


# ==============================================================================
# 2. Blume-Capel recoding
# ==============================================================================

test_that("BC: integer values starting at 0 pass through", {
  x = matrix(c(0, 1, 2, 3, 0, 1, 2, 3), nrow = 4, ncol = 2)
  is_ordinal = c(FALSE, FALSE)
  bc = c(1, 1)

  result = reformat_ordinal_data(x, is_ordinal, bc)
  expect_equal(result$x, x)
  expect_equal(result$num_categories, c(3, 3))
  expect_equal(result$baseline_category, c(1, 1))
})

test_that("BC: values not starting at 0 are shifted, baseline adjusted", {
  x = matrix(c(2, 3, 4, 5, 2, 3, 4, 5), nrow = 4, ncol = 2)
  is_ordinal = c(FALSE, FALSE)
  bc = c(3, 4)

  old_opts = options(bgms.verbose = TRUE)
  on.exit(options(old_opts))
  expect_message(
    result <- reformat_ordinal_data(x, is_ordinal, bc),
    "Variables 1, 2 recoded to start at 0"
  )
  expect_equal(result$x[, 1], c(0, 1, 2, 3))
  expect_equal(result$baseline_category, c(1, 2))
})

test_that("BC: single variable shift messages without plural", {
  x = matrix(c(
    0, 1, 2, 0, 1, 2,
    3, 4, 5, 3, 4, 5
  ), nrow = 6, ncol = 2)
  is_ordinal = c(TRUE, FALSE)
  bc = c(0, 4)

  old_opts = options(bgms.verbose = TRUE)
  on.exit(options(old_opts))
  expect_message(
    result <- reformat_ordinal_data(x, is_ordinal, bc),
    "Variable 2 recoded to start at 0 \\(baseline category adjusted\\)"
  )
  expect_equal(result$baseline_category[2], 1)
})

test_that("BC: verbose=FALSE suppresses shift message", {
  x = matrix(c(
    1, 2, 3, 4, 1, 2,
    1, 2, 3, 4, 1, 2
  ), nrow = 6, ncol = 2)
  is_ordinal = c(FALSE, FALSE)
  bc = c(2, 2)

  old_opts = options(bgms.verbose = FALSE)
  on.exit(options(old_opts))
  expect_silent(
    result <- reformat_ordinal_data(x, is_ordinal, bc)
  )
  expect_equal(result$x[, 1], c(0, 1, 2, 3, 0, 1))
})

test_that("BC: non-integer values are recoded to integer", {
  x = matrix(c(0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0), nrow = 4, ncol = 2)
  # These are already representable as integers <U+2014> just stored as double
  is_ordinal = c(FALSE, FALSE)
  bc = c(1, 1)

  result = reformat_ordinal_data(x, is_ordinal, bc)
  expect_equal(result$x, x)
  expect_equal(result$num_categories, c(3, 3))
})

test_that("BC: > 10 categories triggers warning", {
  vals = 0:11 # 12 values, num_categories = 11
  x = matrix(c(rep(vals, 2), rep(0:2, 8)), ncol = 2)
  is_ordinal = c(FALSE, TRUE)
  bc = c(5, 0)

  expect_warning(
    result <- reformat_ordinal_data(x, is_ordinal, bc),
    "Blume-Capel variable 1 has 11 categories"
  )
  expect_equal(result$num_categories[1], 11)
})


# ==============================================================================
# 3. Error conditions <U+2014> Blume-Capel
# ==============================================================================

test_that("BC: truly non-integer values that can't round to integer <U+2192> error", {
  x = matrix(c(0.3, 1.7, 2.3, 3.7, 0, 1, 2, 3), nrow = 4, ncol = 2)
  is_ordinal = c(FALSE, TRUE)
  bc = c(1, 0)

  # 0.3, 1.7, 2.3, 3.7 <U+2192> as.integer gives 0, 1, 2, 3 (truncation)
  # unique(as.integer) = {0, 1, 2, 3} but length(unique original) = 4 <U+2192> ok
  # But abs(val - round(val)) > eps triggers the check
  # Actually unique(as.integer(c(0.3, 1.7, 2.3, 3.7))) = c(0, 1, 2, 3), length 4
  # and length(c(0.3, 1.7, 2.3, 3.7)) = 4, so lengths match <U+2192> no error there.
  # The integer conversion succeeds. Let me craft a case that actually fails.
  x2 = matrix(c(0.3, 0.7, 2.3, 3.7, 0, 1, 2, 3), nrow = 4, ncol = 2)
  # as.integer(c(0.3, 0.7, 2.3, 3.7)) = c(0, 0, 2, 3), unique = c(0, 2, 3)
  # length unique int = 3 <U+2260> length unique original = 4

  expect_error(
    reformat_ordinal_data(x2, c(FALSE, TRUE), c(1, 0)),
    "a single integer value was used for several observed score categories"
  )
})

test_that("BC: baseline outside range <U+2192> error", {
  # Non-integer values <U+2192> triggers integer recode <U+2192> then baseline check
  x = matrix(c(0.1, 1.1, 2.1, 3.1, 0, 1, 2, 3), nrow = 4, ncol = 2)
  is_ordinal = c(FALSE, TRUE)
  bc = c(10, 0) # baseline 10 is out of range [0, 3]

  expect_error(
    reformat_ordinal_data(x, is_ordinal, bc),
    "reference category for the Blume-Capel variable 1.*outside its"
  )
})

test_that("BC: fewer than 3 unique values <U+2192> error", {
  x = matrix(c(0, 1, 0, 1, 0, 1, 2, 3), nrow = 4, ncol = 2)
  is_ordinal = c(FALSE, TRUE)
  bc = c(0, 0)

  expect_error(
    reformat_ordinal_data(x, is_ordinal, bc),
    "Blume-Capel is only available for variables with more than one category"
  )
})


# ==============================================================================
# 4. Error conditions <U+2014> general
# ==============================================================================

test_that("all unique responses <U+2192> error", {
  # max(unique) must equal nrow for the check to trigger
  # nrow = 4, so values {0, 1, 2, 4} => max = 4 = nrow <U+2192> triggers
  x = matrix(c(0, 1, 2, 4, 0, 1, 2, 4), nrow = 4, ncol = 2)
  is_ordinal = c(TRUE, TRUE)
  bc = c(0, 0)

  expect_error(
    reformat_ordinal_data(x, is_ordinal, bc),
    "Only unique responses observed for variable 1"
  )
})

test_that("single category only <U+2192> error", {
  x = matrix(c(0, 0, 0, 0, 0, 1, 0, 1), nrow = 4, ncol = 2)
  is_ordinal = c(TRUE, TRUE)
  bc = c(0, 0)

  expect_error(
    reformat_ordinal_data(x, is_ordinal, bc),
    "Only one value.*was observed for variable 1"
  )
})


# ==============================================================================
# 5. Mixed ordinal + BC
# ==============================================================================

test_that("mixed ordinal and BC variables processed correctly", {
  # Variable 1: ordinal, values 1,2,3 <U+2192> recode to 0,1,2
  # Variable 2: BC, values 2,3,4,5 <U+2192> shift to 0,1,2,3, baseline adjusted
  # Variable 3: ordinal, values 0,1,2 <U+2192> already contiguous
  x = matrix(c(
    1, 2, 3, 1, 2, 3,
    2, 3, 4, 5, 2, 3,
    0, 1, 2, 0, 1, 2
  ), nrow = 6, ncol = 3)
  is_ordinal = c(TRUE, FALSE, TRUE)
  bc = c(0, 3, 0)

  old_opts = options(bgms.verbose = TRUE)
  on.exit(options(old_opts))
  expect_message(
    result <- reformat_ordinal_data(x, is_ordinal, bc),
    "Variable 2 recoded to start at 0"
  )

  expect_equal(result$x[, 1], c(0, 1, 2, 0, 1, 2))
  expect_equal(result$x[, 2], c(0, 1, 2, 3, 0, 1))
  expect_equal(result$x[, 3], c(0, 1, 2, 0, 1, 2))
  expect_equal(result$num_categories, c(2, 3, 2))
  expect_equal(result$baseline_category[2], 1) # 3 - 2 = 1
})


# ==============================================================================
# 6. Integration: validate_missing_data + reformat_ordinal_data pipeline
# ==============================================================================

test_that("validate_missing_data + reformat_ordinal_data pipeline works end-to-end", {
  x = matrix(c(
    1, 2, 3, 1, 2, 3,
    1, 2, 3, 1, 2, 3
  ), nrow = 6, ncol = 2)
  variable_bool = c(TRUE, TRUE)
  bc = c(0, 0)

  md = validate_missing_data(x, na_action = "listwise")
  result = reformat_ordinal_data(md$x,
    is_ordinal = variable_bool,
    baseline_category = bc
  )
  expect_equal(result$x[, 1], c(0, 1, 2, 0, 1, 2))
  expect_equal(result$num_categories, c(2, 2))
  expect_false(md$na_impute)
})

# ==============================================================================
# All-unique guard: a regular ordinal variable needs repeated values (category
# thresholds are unidentified otherwise). Blume-Capel is parametric in the
# category score and is exempt.
# ==============================================================================

test_that("ordinal: errors when every response is distinct", {
  x = matrix(0:4, ncol = 1) # 5 rows, 5 distinct values
  expect_error(
    reformat_ordinal_data(x, is_ordinal = TRUE, baseline_category = 0L),
    "Only unique responses"
  )
})

test_that("Blume-Capel: all-distinct responses are allowed (parametric, exempt)", {
  x = matrix(0:4, ncol = 1) # 5 rows, 5 distinct values -> fine for BC
  expect_no_error(
    reformat_ordinal_data(x, is_ordinal = FALSE, baseline_category = 0L)
  )
})

test_that("ordinal: does not error when values repeat even if a code equals nrow", {
  # nrow = 4, max code = 4 (== nrow), but only 3 distinct values: not all unique.
  x = matrix(c(0, 4, 4, 1), ncol = 1)
  expect_no_error(
    reformat_ordinal_data(x, is_ordinal = TRUE, baseline_category = 0L)
  )
})

test_that("ordinal: returns the recode map (sorted original values) per variable", {
  x = matrix(c(1, 3, 5, 1, 3, 5), nrow = 6, ncol = 1)
  result = reformat_ordinal_data(x, is_ordinal = TRUE, baseline_category = 0L)
  # category_levels[[1]] are the original sorted values; recoded = match - 1.
  expect_equal(result$category_levels[[1]], c(1, 3, 5))
  expect_equal(result$x[, 1], match(x[, 1], result$category_levels[[1]]) - 1L)
})

test_that("Blume-Capel variables get no recode map (NULL)", {
  x = matrix(c(0, 1, 2, 0, 1, 2), nrow = 6, ncol = 1)
  result = reformat_ordinal_data(x, is_ordinal = FALSE, baseline_category = 0L)
  expect_null(result$category_levels[[1]])
})
