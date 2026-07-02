# --------------------------------------------------------------------------- #
# Tests for the mixed MRF edge-indicator bookkeeping.
#
# The stochastic block edge prior reads full columns of the internal
# (p+q) x (p+q) indicator matrix, so both triangles must stay in sync under
# the discrete, continuous, and cross edge moves.
# --------------------------------------------------------------------------- #

test_that("mixed edge-indicator matrix stays symmetric under indicator sweeps", {
  set.seed(1)
  n = 100
  p = 2
  q = 3
  X = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  Y = matrix(rnorm(n * q), n, q)

  G = test_mixed_edge_indicator_matrix(X, Y, c(3L, 3L), 200L, 42L)

  expect_identical(dim(G), as.integer(c(p + q, p + q)))
  expect_true(all(G == t(G)))

  # The sweeps must actually have moved cross edges for the check to bite:
  # with all-ones initialization, at least one cross indicator should have
  # left its initial state after 200 sweeps.
  cross_block = G[1:p, (p + 1):(p + q), drop = FALSE]
  expect_true(any(cross_block == 0L))
})
