# Trust-gauge reference route (block_reference_logR): the block-local exact
# Monte-Carlo reference log R_e that the in-chain gauge compares against the
# deployed log J. These drive the engine reference in isolation.

test_that("the block reference is finite and lands near the deployed ratio", {
  skip_on_cran()
  p = 10
  delta = 0.5 * log(p)
  eta = 2
  zc = bgms:::zratio_constants(delta, eta = eta)
  set.seed(3)
  G = matrix(0L, p, p)
  ut = upper.tri(G)
  G[ut] = rbinom(sum(ut), 1, 0.6)
  G = G + t(G)
  diag(G) = 1L

  # First edge whose mediating block is non-trivial (m >= 2).
  edge = NULL
  for(i in seq_len(p - 1)) {
    for(j in (i + 1):p) {
      r = bgms:::zratio_test_reference(
        G, i, j, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
        delta, eta, 120L, 30L, 7L, FALSE
      )
      if(isTRUE(r$valid) && r$m >= 2) {
        edge = c(i, j)
        break
      }
    }
    if(!is.null(edge)) break
  }
  expect_false(is.null(edge))

  r = bgms:::zratio_test_reference(
    G, edge[1], edge[2], zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    delta, eta, 240L, 30L, 7L, FALSE
  )
  expect_true(r$ok)
  expect_true(is.finite(r$logR))
  expect_true(is.finite(r$mcse) && r$mcse > 0)
  # The deployed corrected ratio is a validated approximation of the exact
  # reference on an entangled edge: they agree to within a few MCSE.
  expect_lt(abs(r$logR - r$log_zratio), 0.05)
})

test_that("the reference MCSE shrinks with the draw count", {
  skip_on_cran()
  p = 10
  delta = 0.5 * log(p)
  eta = 2
  zc = bgms:::zratio_constants(delta, eta = eta)
  set.seed(3)
  G = matrix(0L, p, p)
  ut = upper.tri(G)
  G[ut] = rbinom(sum(ut), 1, 0.6)
  G = G + t(G)
  diag(G) = 1L

  # A high-degree endpoint pair has an entangled block.
  small = bgms:::zratio_test_reference(
    G, 9, 10, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    delta, eta, 60L, 30L, 11L, FALSE
  )
  large = bgms:::zratio_test_reference(
    G, 9, 10, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    delta, eta, 960L, 30L, 11L, FALSE
  )
  expect_true(small$valid && large$valid)
  # Two estimates from independent streams agree within combined MCSE.
  expect_lt(abs(small$logR - large$logR),
            5 * sqrt(small$mcse^2 + large$mcse^2))
  expect_lt(large$mcse, small$mcse)
})
