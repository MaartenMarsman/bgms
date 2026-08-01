# --------------------------------------------------------------------------- #
# Tests for the rank-1 Cholesky downdate failure signal.
#
# A downdate that would leave the matrix non-positive-definite must report
# failure so the model classes rebuild the factor from the source matrix
# instead of carrying a partially updated one.
# --------------------------------------------------------------------------- #

test_that("successful downdate updates the factor and reports ok", {
  set.seed(3)
  A = crossprod(matrix(rnorm(25), 5, 5)) + diag(5)
  R = chol(A)
  u = rnorm(5, sd = 0.1)

  out = test_cholesky_downdate(R, u)

  expect_true(out$ok)
  expect_equal(out$R[upper.tri(out$R, diag = TRUE)],
    chol(A - tcrossprod(u))[upper.tri(diag(5), diag = TRUE)],
    tolerance = 1e-10
  )
})

test_that("non-positive-definite downdate reports failure", {
  R = diag(3)
  u = c(1.2, 0, 0)

  out = test_cholesky_downdate(R, u)

  expect_false(out$ok)
})

test_that("1 x 1 downdate handles both outcomes", {
  ok_case = test_cholesky_downdate(matrix(2), 1)
  expect_true(ok_case$ok)
  expect_equal(ok_case$R[1, 1], sqrt(3))

  fail_case = test_cholesky_downdate(matrix(2), 3)
  expect_false(fail_case$ok)
})
