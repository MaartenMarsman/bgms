# ==============================================================================
# Unit tests for validate_missing_data()
# Phase A.5 of the R scaffolding refactor.
# ==============================================================================

# ==============================================================================
# 1. Listwise <U+2014> no NAs
# ==============================================================================

test_that("listwise with no NAs returns data unchanged", {
  x = matrix(1:12, nrow = 4, ncol = 3)
  result = validate_missing_data(x, na_action = "listwise")
  expect_equal(result$x, x)
  expect_false(result$na_impute)
  expect_equal(result$missing_index, matrix(NA, 1, 1))
  expect_equal(result$n_removed, 0L)
})

# ==============================================================================
# 2. Listwise <U+2014> some NAs removed
# ==============================================================================

test_that("listwise removes rows with NAs and messages", {
  x = matrix(1:12, nrow = 4, ncol = 3)
  x[2, 1] = NA
  x[4, 3] = NA

  old_opts = options(bgms.verbose = TRUE)
  on.exit(options(old_opts))
  expect_message(
    result <- validate_missing_data(x, na_action = "listwise"),
    "2 rows with missing values excluded"
  )
  expect_equal(nrow(result$x), 2)
  expect_false(result$na_impute)
  expect_equal(result$n_removed, 2L)
})

test_that("listwise suppresses message when bgms.verbose is FALSE", {
  x = matrix(1:12, nrow = 4, ncol = 3)
  x[2, 1] = NA
  old_opts = options(bgms.verbose = FALSE)
  on.exit(options(old_opts))
  expect_silent(
    result <- validate_missing_data(x, na_action = "listwise")
  )
  expect_equal(nrow(result$x), 3)
})

# ==============================================================================
# 3. Listwise <U+2014> all rows have NAs
# ==============================================================================

test_that("listwise errors when all rows have NAs", {
  x = matrix(NA_real_, nrow = 3, ncol = 2)
  expect_error(
    validate_missing_data(x, na_action = "listwise"),
    "All rows in x contain at least one missing response"
  )
})

# ==============================================================================
# 4. Listwise <U+2014> too few rows after removal
# ==============================================================================

test_that("listwise errors when < 2 rows remain", {
  x = matrix(1:6, nrow = 3, ncol = 2)
  x[1, 1] = NA
  x[2, 2] = NA
  old_opts = options(bgms.verbose = FALSE)
  on.exit(options(old_opts))
  expect_error(
    validate_missing_data(x, na_action = "listwise"),
    "less than two rows"
  )
})

# ==============================================================================
# 5. Impute <U+2014> no NAs
# ==============================================================================

test_that("impute with no NAs returns data unchanged", {
  x = matrix(1:12, nrow = 4, ncol = 3)
  result = validate_missing_data(x, na_action = "impute")
  expect_equal(result$x, x)
  expect_false(result$na_impute)
  expect_equal(result$missing_index, matrix(NA, 1, 1))
})

# ==============================================================================
# 6. Impute <U+2014> with NAs
# ==============================================================================

test_that("impute fills NAs and builds missing_index", {
  set.seed(42)
  x = matrix(c(1, 2, 3, 4, 5, NA, 7, 8, 9), nrow = 3, ncol = 3)
  result = validate_missing_data(x, na_action = "impute")

  expect_true(result$na_impute)
  expect_false(anyNA(result$x))
  expect_equal(nrow(result$missing_index), 1)
  # 0-based indices: row 2 (0-based: 1), col 1 (0-based: 1)
  expect_equal(result$missing_index[1, 1], 2) # row index (0-based)
  expect_equal(result$missing_index[1, 2], 1) # col index (0-based)
})

test_that("impute builds correct missing_index for multiple NAs", {
  set.seed(123)
  x = matrix(c(1, NA, 3, NA, 5, 6, 7, 8, NA), nrow = 3, ncol = 3)
  colnames(x) = paste0("V", 1:3)
  result = validate_missing_data(x, na_action = "impute")

  expect_true(result$na_impute)
  expect_equal(nrow(result$missing_index), 3)
  expect_false(anyNA(result$x))
})

# ==============================================================================
# 7. GGM + impute now supported
# ==============================================================================

test_that("GGM + impute works", {
  x = matrix(rnorm(12), nrow = 4, ncol = 3)
  colnames(x) = paste0("V", 1:3)
  x[1, 2] = NA
  result = validate_missing_data(x, na_action = "impute", is_continuous = TRUE)
  expect_true(result$na_impute)
  expect_equal(nrow(result$missing_index), 1L)
})

test_that("GGM + impute: entire-column-missing gives clear error", {
  x = matrix(rnorm(12), nrow = 4, ncol = 3)
  colnames(x) = paste0("V", 1:3)
  x[, 2] = NA
  expect_error(
    validate_missing_data(x, na_action = "impute", is_continuous = TRUE),
    "no observed values"
  )
})

test_that("GGM + listwise works", {
  x = matrix(rnorm(12), nrow = 4, ncol = 3)
  x[1, 2] = NA
  old_opts = options(bgms.verbose = FALSE)
  on.exit(options(old_opts))
  result = validate_missing_data(x, na_action = "listwise", is_continuous = TRUE)
  expect_equal(nrow(result$x), 3)
  expect_false(result$na_impute)
})

# ==============================================================================
# 8. Group vector filtering (bgmCompare path)
# ==============================================================================

test_that("listwise with group filters group vector", {
  x = matrix(1:12, nrow = 4, ncol = 3)
  x[2, 1] = NA
  group = c(1, 1, 2, 2)

  old_opts = options(bgms.verbose = FALSE)
  on.exit(options(old_opts))
  result = validate_missing_data(x, na_action = "listwise", group = group)
  expect_equal(length(result$group), nrow(result$x))
  expect_equal(result$group, c(1, 2, 2))
})

test_that("impute with group preserves group vector", {
  x = matrix(1:12, nrow = 4, ncol = 3)
  group = c(1, 1, 2, 2)
  result = validate_missing_data(x, na_action = "impute", group = group)
  expect_equal(result$group, group)
})

test_that("group is absent from result when not provided", {
  x = matrix(1:12, nrow = 4, ncol = 3)
  result = validate_missing_data(x, na_action = "listwise")
  expect_null(result$group)
})

# ==============================================================================
# Imputation draws from the observed values, indexing explicitly so a column
# with a single observed value is filled with that value (not sample(v, 1),
# which would treat a length-1 numeric v as the range 1:v).
# ==============================================================================

test_that("impute fills a single-observed-value column with that value", {
  x = matrix(c(
    7, NA, NA, NA, NA, NA,
    1, 0, 1, 0, 1, 0
  ), nrow = 6, ncol = 2)
  set.seed(1)
  result = validate_missing_data(x, na_action = "impute")
  expect_true(all(result$x[, 1] == 7))
})
