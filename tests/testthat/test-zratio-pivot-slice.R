# The block oracle's diagonal pivot update. At a non-unit Gamma shape the row
# is two blocks -- an accept/reject on the off-diagonals at a held pivot, then
# an exact slice update of the pivot itself -- replacing a single joint
# independence-Metropolis step whose weight carried K_ii^(alpha - 1) and froze
# the whole row as the shape grew.

dense_block = function(q) {
  g = matrix(0L, q, q)
  for(a in 1:(q - 1)) for(b in (a + 1):q) {
    g[a, b] = 1L
    g[b, a] = 1L
  }
  g
}

oracle = function(g, alpha, sweeps, burn, seed, delta = 0.5 * log(12), eta = 2) {
  zc = bgms:::zratio_constants(delta, eta, alpha = alpha)
  zratio_test_gold_moments(
    g, 1, 2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    zc$delta, zc$eta, as.integer(sweeps), as.integer(burn), as.integer(seed),
    FALSE, alpha
  )
}

test_that("the exponential shape is bit-identical to the pre-split kernel", {
  skip_on_cran()
  # The guard for every validated cell: at alpha = 1 there is no coupling, so
  # the row stays a direct Gibbs draw and the RNG call order is unchanged
  # (nq normals, then the conjugate Gamma, no uniform). These values were
  # produced by the kernel as it stood before the pivot split, same seed, and
  # are pinned to every digit a decimal literal round-trips: a changed draw
  # stream moves them in the first few, not the last.
  r = oracle(dense_block(8), 1, 2000L, 500L, 11L)
  expect_equal(r$S1, 0.1798251775467119, tolerance = 1e-15)
  expect_equal(r$logR, 0.60844358005728083, tolerance = 1e-15)
  # No accept step and no slice at the exponential shape.
  expect_equal(r$im_proposed, 0)
  expect_equal(r$n_slice_cap, 0)
})

test_that("the slice path reproduces the conjugate law it replaces", {
  skip_on_cran()
  # Runs every time rather than behind the slow gate: it is the only check that
  # scores the split branch against ground truth instead of against another
  # approximation, and it costs a few seconds.
  #
  # A shape a hair off 1 takes the split branch, so the off-diagonal step
  # always accepts and the pivot conditional collapses to the same
  # Gamma(delta + 1, beta) the conjugate path draws directly. The two paths
  # must therefore agree on the block moments: this is the slice sampler and
  # the two-block restructure scored against the bit-identical kernel above.
  g = dense_block(12)
  seeds = c(11L, 23L, 37L)
  draw = function(alpha) {
    vapply(seeds, function(s) {
      r = oracle(g, alpha, 2500L, 500L, s)
      c(r$S1, r$S2, r$logR)
    }, numeric(3))
  }
  conj = draw(1)
  slice = draw(1 + 1e-9)
  se = sqrt(apply(conj, 1, var) / length(seeds) +
              apply(slice, 1, var) / length(seeds))
  expect_lt(max(abs(rowMeans(slice) - rowMeans(conj)) / se), 3)
})

test_that("the pivot slice never exhausts its shrinkage budget", {
  skip_on_cran()
  # The cap is a backstop, not a working part: a hit leaves the pivot unmoved
  # and biases the chain. Zero across the deployed shape range, including well
  # past it, is the reading that keeps it a backstop.
  for(alpha in c(0.5, 2, 5, 10)) {
    r = suppressWarnings(oracle(dense_block(10), alpha, 600L, 150L, 5L))
    expect_equal(r$n_slice_cap, 0)
    expect_gt(r$im_accepted / r$im_proposed, 0.1)
  }
})

test_that("the pivot slice handles a non-log-concave conditional", {
  skip_on_cran()
  # The conditional is xi^delta e^(-beta xi) (xi + q)^(alpha - 1), with
  # curvature -delta/xi^2 - (alpha-1)/(xi+q)^2. Below shape 1 the second term
  # is positive, so concavity needs delta >= 1 - alpha. At delta = 0.2,
  # alpha = 0.5 the curvature flips at xi = 1.72 with 3.2% of the target
  # beyond it -- a region where adaptive rejection's guarantee genuinely
  # fails. Slice does not need concavity; this cell is pinned so that any
  # future proposal that quietly assumes it (an ARS step, a Laplace proposal)
  # fails against the case that proves it cannot.
  curvature = function(xi, delta, alpha, q) {
    -delta / xi^2 - (alpha - 1) / (xi + q)^2
  }
  expect_lt(curvature(0.3, 0.2, 0.5, 1), 0)   # concave near the origin
  expect_gt(curvature(3.0, 0.2, 0.5, 1), 0)   # and not, past the flip
  expect_lt(curvature(3.0, 0.5 * log(12), 0.5, 1), 0)  # default stays concave

  r = oracle(dense_block(10), 0.5, 800L, 200L, 9L, delta = 0.2)
  expect_true(isTRUE(r$valid))
  expect_equal(r$n_slice_cap, 0)
  expect_true(is.finite(r$S1) && r$S1 > 0)
})
