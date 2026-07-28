# Unit tests for the anchored-curve reweighting engine (anchor_curve.R):
# per-draw log weights against direct density ratios, the identity at the
# anchor's own scale, the ESS honesty mask, and best-anchor assembly.

test_that("normal and cauchy log weights match direct slab-density ratios", {
  set.seed(1)
  n_iter = 50L
  n_edge = 4L
  theta = matrix(rnorm(n_iter * n_edge), n_iter, n_edge)
  gamma = matrix(rbinom(n_iter * n_edge, 1, 0.6), n_iter, n_edge)
  s_a = 1.2
  s_grid = c(0.8, 1.2, 2.0)

  direct = function(density) {
    vapply(s_grid, function(s) {
      rowSums(gamma * (density(theta, s) - density(theta, s_a)))
    }, numeric(n_iter))
  }

  lw = anchor_log_weights(theta, gamma, "normal", s_a, s_grid)
  ref = direct(function(x, s) stats::dnorm(x, 0, s, log = TRUE))
  expect_equal(lw, ref, tolerance = 1e-12)

  lw_c = anchor_log_weights(theta, gamma, "cauchy", s_a, s_grid)
  ref_c = direct(function(x, s) stats::dcauchy(x, 0, s, log = TRUE))
  expect_equal(lw_c, ref_c, tolerance = 1e-12)
})


test_that("reweighting at the anchor's own scale is the identity", {
  set.seed(2)
  n_iter = 200L
  n_edge = 3L
  draws = list(
    theta = list(
      matrix(rnorm(n_iter * n_edge), n_iter, n_edge),
      matrix(rnorm(n_iter * n_edge), n_iter, n_edge)
    ),
    gamma = list(
      matrix(rbinom(n_iter * n_edge, 1, 0.5), n_iter, n_edge),
      matrix(rbinom(n_iter * n_edge, 1, 0.5), n_iter, n_edge)
    ),
    family = "normal"
  )
  rw = anchor_reweight(draws, s_a = 1, s_grid = c(1, 1.3))

  pooled_gamma = rbind(draws$gamma[[1]], draws$gamma[[2]])
  expect_equal(rw$pip[1, ], colMeans(pooled_gamma), tolerance = 1e-12)
  expect_equal(rw$ess[1], 2 * n_iter, tolerance = 1e-9)
  expect_equal(rw$chain_pip[[1]][1, ], colMeans(draws$gamma[[1]]), tolerance = 1e-12)
  expect_equal(rw$chain_ess[1, 1], n_iter, tolerance = 1e-9)
  # Away from the anchor the ESS strictly decays.
  expect_lt(rw$ess[2], 2 * n_iter)
})


test_that("the ESS floor masks far scales instead of reporting a value", {
  # Many included edges make the weights degenerate at a far scale: the point
  # must mask NA, never report a number (regression against the phase-B
  # failure mode).
  set.seed(3)
  n_iter = 400L
  n_edge = 50L
  draws = list(
    theta = list(matrix(rnorm(n_iter * n_edge), n_iter, n_edge)),
    gamma = list(matrix(1L, n_iter, n_edge)),
    family = "normal"
  )
  s_grid = c(1, 2.5)
  rw = anchor_reweight(draws, s_a = 1, s_grid = s_grid)
  expect_lt(rw$ess[2], 10) # importance weights degenerate at 2.5x
  curve = assemble_curve(list(rw), usable = TRUE, ess_floor = 400)
  expect_true(is.na(curve$anchor_used[2]))
  expect_true(all(is.na(curve$pip[2, ])))
  # The anchor's own scale stays exact.
  expect_false(is.na(curve$anchor_used[1]))
  expect_equal(curve$pip[1, ], rep(1, n_edge))
})


test_that("assembly picks the highest-ESS usable anchor per point", {
  set.seed(4)
  n_iter = 300L
  n_edge = 5L
  mk = function() {
    list(
      theta = list(matrix(rnorm(n_iter * n_edge, sd = 0.5), n_iter, n_edge)),
      gamma = list(matrix(rbinom(n_iter * n_edge, 1, 0.7), n_iter, n_edge)),
      family = "normal"
    )
  }
  s_grid = c(0.8, 1.0, 1.6, 2.0)
  rw_lo = anchor_reweight(mk(), s_a = 1.0, s_grid = s_grid)
  rw_hi = anchor_reweight(mk(), s_a = 2.0, s_grid = s_grid)
  curve = assemble_curve(list(rw_lo, rw_hi), usable = c(TRUE, TRUE), ess_floor = 10)
  expect_equal(curve$anchor_used[2], 1L) # own scale of anchor 1
  expect_equal(curve$anchor_used[4], 2L) # own scale of anchor 2
  # An unusable anchor never contributes, even at its own scale.
  curve2 = assemble_curve(list(rw_lo, rw_hi), usable = c(TRUE, FALSE), ess_floor = 10)
  expect_true(all(curve2$anchor_used == 1L, na.rm = TRUE))
})


test_that("precision-weighted pooling blends anchors by inverse variance", {
  # Two anchors reach a common target scale. The pooled PIP is the
  # inverse-variance weighted mean (weight = ESS / (p(1-p))), not either
  # anchor alone; the dominant (highest-ESS) anchor is tagged in anchor_used.
  mk = function(p_true, n) {
    list(
      theta = list(matrix(0, n, 1L)), # slab value irrelevant at s = s_a
      gamma = list(matrix(rbinom(n, 1, p_true), n, 1L)),
      family = "normal"
    )
  }
  set.seed(5)
  a1 = mk(0.4, 4000L)
  a2 = mk(0.6, 1000L)
  s_grid = 1.0 # single point, both anchors at their own scale (identity)
  rw1 = anchor_reweight(a1, s_a = 1.0, s_grid = s_grid)
  rw2 = anchor_reweight(a2, s_a = 1.0, s_grid = s_grid)
  p1 = rw1$pip[1, 1]
  p2 = rw2$pip[1, 1]
  w1 = rw1$ess[1] / (p1 * (1 - p1))
  w2 = rw2$ess[1] / (p2 * (1 - p2))
  expected = (w1 * p1 + w2 * p2) / (w1 + w2)
  curve = assemble_curve(list(rw1, rw2), usable = c(TRUE, TRUE), ess_floor = 10)
  expect_equal(curve$pip[1, 1], expected, tolerance = 1e-12)
  # pooled estimate lies strictly between the two anchor estimates
  expect_gt(curve$pip[1, 1], min(p1, p2))
  expect_lt(curve$pip[1, 1], max(p1, p2))
  # the higher-ESS anchor (a1) dominates the tag
  expect_equal(curve$anchor_used[1], 1L)
})
