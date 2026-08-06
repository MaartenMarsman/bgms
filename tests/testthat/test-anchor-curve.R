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
  # inverse-variance weighted mean, not either anchor alone; the dominant
  # (highest-ESS) anchor is tagged in anchor_used. The variance plug-in is the
  # smoothed proportion (ESS * p + 0.5) / (ESS + 1), so a saturated anchor
  # cannot claim near-zero variance; the pooled estimate itself is unsmoothed.
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
  smoothed = function(p, ess) (ess * p + 0.5) / (ess + 1)
  w1 = rw1$ess[1] / (smoothed(p1, rw1$ess[1]) * (1 - smoothed(p1, rw1$ess[1])))
  w2 = rw2$ess[1] / (smoothed(p2, rw2$ess[1]) * (1 - smoothed(p2, rw2$ess[1])))
  expected = (w1 * p1 + w2 * p2) / (w1 + w2)
  curve = assemble_curve(list(rw1, rw2), usable = c(TRUE, TRUE), ess_floor = 10)
  expect_equal(curve$pip[1, 1], expected, tolerance = 1e-12)
  # pooled estimate lies strictly between the two anchor estimates
  expect_gt(curve$pip[1, 1], min(p1, p2))
  expect_lt(curve$pip[1, 1], max(p1, p2))
  # the higher-ESS anchor (a1) dominates the tag
  expect_equal(curve$anchor_used[1], 1L)
})


test_that("the diagonal prior ratio enters the weights when the sweep moves it", {
  # Under vary = "slab-and-diagonal" the refits rewrite the raw diagonal rate
  # as eta / s, so the anchor and the target differ in the Gamma prior on the
  # precision diagonal as well as in the slab. The weight must carry both
  # ratios; the diagonal one is the direct Gamma log-density ratio, summed over
  # the p diagonal values the sampler evaluates the prior at (K_jj / 2).
  set.seed(11)
  n_iter = 40L
  n_edge = 3L
  p = 4L
  theta = matrix(rnorm(n_iter * n_edge, sd = 0.2), n_iter, n_edge)
  gamma = matrix(rbinom(n_iter * n_edge, 1, 0.6), n_iter, n_edge)
  k_diag = matrix(stats::rgamma(n_iter * p, shape = 2, rate = 1), n_iter, p)
  shape = 1.3
  diagonal = list(sum = 0.5 * rowSums(k_diag), n = p, shape = shape)
  eta = 0.7
  s_a = 0.4
  s_grid = c(0.3, 0.4, 0.55)

  base = anchor_log_weights(theta, gamma, "cauchy", s_a, s_grid)
  moved = anchor_log_weights(theta, gamma, "cauchy", s_a, s_grid,
    diagonal = diagonal, eta = eta
  )
  direct = vapply(s_grid, function(s) {
    rowSums(vapply(seq_len(p), function(j) {
      x = 0.5 * k_diag[, j]
      stats::dgamma(x, shape = shape, rate = eta / s, log = TRUE) -
        stats::dgamma(x, shape = shape, rate = eta / s_a, log = TRUE)
    }, numeric(n_iter)))
  }, numeric(n_iter))
  expect_equal(moved - base, direct, tolerance = 1e-10)

  # At the anchor's own scale the diagonal ratio is zero, as the slab one is.
  expect_equal(
    anchor_log_weights(theta, gamma, "cauchy", s_a, s_a,
      diagonal = diagonal, eta = eta
    ),
    matrix(0, n_iter, 1L)
  )

  # A sweep that leaves the diagonal prior alone carries no such term.
  expect_equal(
    anchor_log_weights(theta, gamma, "cauchy", s_a, s_grid,
      diagonal = diagonal, eta = NULL
    ),
    base
  )
})


test_that("anchor_diagonal reads the precision diagonal a model actually has", {
  # GGM: every main column is a precision diagonal; the statistic is the sum
  # of the values the prior is evaluated at, K_jj / 2.
  set.seed(12)
  main = list(matrix(stats::rgamma(30, 2, 1), 10, 3), matrix(stats::rgamma(30, 2, 1), 10, 3))
  spec_ggm = list(model_type = "ggm", prior = list(scale_shape = 2.5))
  raw = list(main = main, parameter_names = list(main = paste0("v", 1:3, " (precision)")))
  d = anchor_diagonal(spec_ggm, raw)
  expect_equal(d$n, 3L)
  expect_equal(d$shape, 2.5)
  expect_equal(d$sum[[1]], 0.5 * rowSums(main[[1]]))
  expect_equal(d$sum[[2]], 0.5 * rowSums(main[[2]]))

  # Mixed MRF: only the columns that are precision diagonals.
  raw_mixed = list(
    main = list(matrix(1:20, 5, 4)),
    parameter_names = list(main = c(
      "a (1)", "b (mean)", "b (precision diag)", "c (precision diag)"
    ))
  )
  spec_mixed = list(model_type = "mixed_mrf", prior = list(scale_shape = 1))
  dm = anchor_diagonal(spec_mixed, raw_mixed)
  expect_equal(dm$n, 2L)
  expect_equal(dm$sum[[1]], 0.5 * rowSums(matrix(1:20, 5, 4)[, 3:4]))

  # A discrete model has no precision diagonal at all.
  expect_null(anchor_diagonal(
    list(model_type = "omrf", prior = list()),
    list(main = list(matrix(0, 5, 3)), parameter_names = list(main = letters[1:3]))
  ))
})


test_that("a slab family the reweighting identity is not written for stops", {
  # The two-way if/else handed every non-normal family the Cauchy formula, so
  # a beta-prime fit got a curve rather than an error.
  theta = matrix(rnorm(20), 10, 2)
  gamma = matrix(1L, 10, 2)
  expect_error(
    anchor_log_weights(theta, gamma, "beta-prime", 1, c(0.5, 1, 2)),
    "normal or Cauchy slab"
  )
  expect_silent(anchor_log_weights(theta, gamma, "normal", 1, c(0.5, 1, 2)))
  expect_silent(anchor_log_weights(theta, gamma, "cauchy", 1, c(0.5, 1, 2)))
})


test_that("a saturated low-ESS anchor no longer overrides a high-ESS one", {
  # BEHAVIOR-CHECK GATE for the smoothed variance plug-in. The raw plug-in
  # p(1-p) collapses to the 1e-6 floor when an anchor's reweighted PIP
  # saturates, giving that anchor a weight of ESS / 1e-6 -- so a distant,
  # low-ESS, saturated anchor buried a near, high-ESS, well-measured one.
  mk = function(pip, ess) {
    list(
      pip = matrix(pip, 1L, 1L), ess = ess,
      chain_pip = list(matrix(pip, 1L, 1L))
    )
  }
  old_pool = function(pips, esss) {
    w = esss / pmax(pips * (1 - pips), 1e-6)
    sum(w * pips) / sum(w)
  }

  # (i) high-ESS non-saturated anchor vs distant saturated low-ESS anchor
  pips = c(0.60, 1.00)
  esss = c(2000, 25)
  before = old_pool(pips, esss)
  after = assemble_curve(
    list(mk(pips[1], esss[1]), mk(pips[2], esss[2])),
    usable = c(TRUE, TRUE), ess_floor = 10
  )$pip[1, 1]
  # Before, the saturated anchor set the pooled value outright.
  expect_gt(before, 0.99)
  expect_lt(abs(before - pips[2]), 0.01)
  # After, the well-measured anchor leads and the pooled value sits nearer it.
  expect_lt(after, before)
  expect_lt(abs(after - pips[1]), abs(after - pips[2]))
  # It is still a genuine pooled value, not a switch to winner-take-all.
  expect_gt(after, pips[1])
  expect_lt(after, pips[2])

  # (ii) away from saturation the smoothing is a formality: < 1% relative move
  configs = list(
    list(pips = c(0.30, 0.50), esss = c(2000, 1500)),
    list(pips = c(0.12, 0.22), esss = c(900, 400)),
    list(pips = c(0.55, 0.80), esss = c(1200, 600)),
    list(pips = c(0.05, 0.10), esss = c(800, 800))
  )
  rel = vapply(configs, function(cf) {
    b = old_pool(cf$pips, cf$esss)
    a = assemble_curve(
      list(mk(cf$pips[1], cf$esss[1]), mk(cf$pips[2], cf$esss[2])),
      usable = c(TRUE, TRUE), ess_floor = 10
    )$pip[1, 1]
    abs(a - b) / b
  }, numeric(1))
  expect_true(all(rel < 0.01))
})
