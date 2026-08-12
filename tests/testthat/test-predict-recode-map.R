# ==============================================================================
# recode_data_for_prediction(): newdata must be recoded to the 0-based
# categories the model was fitted on, using the stored training recode map
# (category_levels). Subtracting the per-column minimum (the legacy fallback)
# is wrong when training categories were non-contiguous or newdata spans a
# different range.
# ==============================================================================

recode = bgms:::recode_data_for_prediction

test_that("map recodes non-contiguous training categories correctly", {
  # Training values {1,3,5} -> categories {0,1,2}. Legacy subtract-min gives
  # {0,2,4} (wrong); the map gives {0,1,2}.
  levels = list(c(1, 3, 5))
  out = recode(matrix(c(1, 3, 5), ncol = 1), TRUE, category_levels = levels)
  expect_equal(out[, 1], c(0, 1, 2))
})

test_that("map is absolute, not relative to newdata's range", {
  # Training {1,2,3} -> {0,1,2}; value 2 is category 1 regardless of whether
  # the lowest category appears in newdata.
  levels = list(c(1, 2, 3))
  out = recode(matrix(c(2, 3), ncol = 1), TRUE, category_levels = levels)
  expect_equal(out[, 1], c(1, 2))
})

test_that("values not observed in training give NA with a warning", {
  levels = list(c(1, 2, 3))
  expect_warning(
    out <- recode(matrix(c(1, 9), ncol = 1), TRUE, category_levels = levels),
    "not\\s+observed in the training data"
  )
  expect_true(is.na(out[2, 1]))
  expect_equal(out[1, 1], 0)
})

test_that("an ordinal column with a NULL level entry is an error, not a guess", {
  # The old fallback shifted by the newdata column minimum, which is the
  # training offset only by coincidence. Every current fit carries a map.
  levels = list(NULL)
  expect_error(
    recode(matrix(c(2, 3, 4), ncol = 1), TRUE, category_levels = levels),
    "predates the recode map"
  )
})

test_that("no map at all (old fit) is an error naming the refit", {
  expect_error(
    recode(matrix(c(2, 3, 4), ncol = 1), TRUE, category_levels = NULL),
    "predates the recode map"
  )
})

test_that("continuous/non-ordinal columns are left unchanged", {
  levels = list(NULL, c(0, 1, 2))
  x = matrix(c(1.5, 2.5, 3.5, 0, 1, 2), ncol = 2)
  out = recode(x, c(FALSE, TRUE), category_levels = levels)
  expect_equal(out[, 1], c(1.5, 2.5, 3.5))
  expect_equal(out[, 2], c(0, 1, 2))
})

# ------------------------------------------------------------------------------
# bgmCompare uses a NAMED lookup (names = original values, values = final
# collapsed categories) which may be many-to-one when cross-group collapsing
# merged categories.
# ------------------------------------------------------------------------------

test_that("named lookup maps original values to final categories (bijective)", {
  lk = c(0L, 1L, 2L)
  names(lk) = c("1", "3", "5")
  out = recode(matrix(c(1, 3, 5), ncol = 1), TRUE, category_levels = list(lk))
  expect_equal(out[, 1], c(0, 1, 2))
})

test_that("named lookup handles many-to-one (merged categories)", {
  # original 1 and 3 both map to category 0 (merged); 5 -> 1.
  lk = c(0L, 0L, 1L)
  names(lk) = c("1", "3", "5")
  out = recode(matrix(c(1, 3, 5, 1, 5), ncol = 1), TRUE, category_levels = list(lk))
  expect_equal(out[, 1], c(0, 0, 1, 0, 1))
})

test_that("named lookup warns + NA on values absent from the lookup", {
  lk = c(0L, 1L)
  names(lk) = c("1", "2")
  expect_warning(
    out <- recode(matrix(c(1, 9), ncol = 1), TRUE, category_levels = list(lk)),
    "not\\s+observed in the training data"
  )
  expect_true(is.na(out[2, 1]))
})
