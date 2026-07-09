# Cauchy-slab Z-ratio constants. The hierarchical spec with a Cauchy
# interaction prior normalizes by the marginal-Cauchy per-graph constant
# Z_C(Gamma) = E_omega[Z_N(Gamma; sigma sqrt(omega))]: the omega updates
# stay conjugate, the between-move slab factor stays omega-conditional, and
# the Z-ratio becomes a fixed function of Gamma with slab-family-specific
# channel constants. Each quadrature channel is checked here against a
# direct Monte Carlo of the same F-measure moment ratio; the graph-law
# identity itself is gated in test-hier-zratio-identity.R.

test_that("Cauchy node channel matches direct Monte Carlo", {
  skip_on_cran()
  delta = 0.5 * log(6)
  sigma = 1
  beta = 1
  t2 = 2 * beta * sigma^2
  w_quad = bgms:::zratio_node_channel(delta, sigma, beta, slab = "cauchy")
  set.seed(101)
  n = 1e6
  x = rgamma(n, delta + 2, rate = beta)
  oa = 1 / rnorm(n)^2
  ob = 1 / rnorm(n)^2
  da = x + t2 * oa
  db = x + t2 * ob
  den = mean((da * db)^-0.5)
  w1 = sigma^4 * mean(oa * ob * (da * db)^-1.5) / den
  w2 = sigma^8 * mean((oa * ob)^2 * (da * db)^-2.5) / den
  expect_lt(abs(w1 / w_quad[1] - 1), 0.01)
  expect_lt(abs(w2 / w_quad[2] - 1), 0.02)
})

test_that("Cauchy bridge channel matches direct Monte Carlo", {
  skip_on_cran()
  delta = 0.5 * log(6)
  sigma = 1
  beta = 1
  t2 = 2 * beta * sigma^2
  cb_quad = bgms:::zratio_bridge_channel(delta, sigma, beta, slab = "cauchy")
  set.seed(202)
  n = 2e6
  kaa = rexp(n, beta)
  kbb = rexp(n, beta)
  u = rcauchy(n, 0, sigma)
  oa = 1 / rnorm(n)^2
  ob = 1 / rnorm(n)^2
  d = kaa * kbb - u^2
  pd = d > 0
  d = d[pd]
  da = d + t2 * oa[pd] * kbb[pd]
  db = d + t2 * ob[pd] * kaa[pd]
  base = d^delta * d
  den = mean(base * (da * db)^-0.5)
  cb1 = mean(base * sigma^4 * oa[pd] * ob[pd] * u[pd]^2 * (da * db)^-1.5) / den
  cb2 = mean(
    base * sigma^8 * (oa[pd] * ob[pd])^2 * u[pd]^4 * (da * db)^-2.5
  ) / den
  expect_lt(abs(cb1 / cb_quad[1] - 1), 0.02)
  expect_lt(abs(cb2 / cb_quad[2] - 1), 0.04)
})

test_that("Cauchy isolated-edge ratio psi0 matches direct Monte Carlo", {
  skip_on_cran()
  delta = 0.5 * log(6)
  pair = bgms:::zratio_pair_integrals(delta,
    sigma = 1, beta = 1,
    slab = "cauchy"
  )
  psi0 = pair$ispike(0) / pair$g(0)
  set.seed(303)
  n = 2e6
  x1 = rexp(n, 1)
  x2 = rexp(n, 1)
  k = rcauchy(n, 0, 1)
  i0 = mean((x1 * x2)^delta)
  g0 = mean(ifelse(k^2 < x1 * x2, (x1 * x2 - k^2)^delta, 0))
  expect_lt(abs((i0 / g0) / psi0 - 1), 0.01)
})

test_that("Cauchy constants live in their own cache cell", {
  skip_on_cran()
  d = 0.5 * log(6)
  zn = bgms:::zratio_cell_constants(d, 0.5, 2)
  zc = bgms:::zratio_cell_constants(d, 0.5, 2, slab = "cauchy")
  expect_identical(zn$slab, "normal")
  expect_identical(zc$slab, "cauchy")
  expect_false(isTRUE(all.equal(zn$addc, zc$addc)))
  expect_false(isTRUE(all.equal(zn$psi0, zc$psi0)))
  # Scale invariance holds per family: same eta, different frames.
  zc2 = bgms:::zratio_cell_constants(d, 0.25, 4, slab = "cauchy")
  expect_identical(zc, zc2)
})

test_that("Cauchy exact Monte-Carlo evaluation tracks the additive prediction", {
  skip_on_cran()
  # Two common neighbours, no CN-CN edge: the additive form is the exact
  # no-edge baseline (kappa_2 = m * w1), so the oracle correction
  # log(saddle on oracle moments) - log(additive) must sit at MC noise.
  # Exercises the omega-augmented sweep and the leg-dressed moments.
  q = 4L
  G = matrix(0L, q, q)
  diag(G) = 1L
  for(e in list(c(1, 3), c(2, 3), c(1, 4), c(2, 4))) {
    G[e[1], e[2]] = 1L
    G[e[2], e[1]] = 1L
  }
  for(slab in c("normal", "cauchy")) {
    zc = bgms:::zratio_cell_constants(0.5 * log(6), 0.5, 2, slab = slab)
    # Block-local exact reference on the (common-neighbour, no-bridge) block:
    # the deployed additive ratio is accurate here, so the reference must sit
    # near it. Exercises the omega-augmented sweep and the leg-dressed
    # endpoint transform under the Cauchy slab.
    r = bgms:::zratio_test_reference(
      G, 1L, 2L, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
      zc$delta, zc$eta, 600L, 50L, 42L, identical(slab, "cauchy")
    )
    expect_true(isTRUE(r$valid) && isTRUE(r$ok), label = slab)
    expect_true(is.finite(r$logR), label = slab)
    expect_lt(abs(r$logR - r$log_zratio), 0.06, label = slab)
  }
})
