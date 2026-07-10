# Gamma-shape (alpha != 1) Z-ratio constants. The channel builders carry
# the diagonal prior K_ii ~ Gamma(alpha, beta) with the prior weight on the
# raw diagonals: the spike/slab pair integrals gain a (s1 s2)^(alpha - 1)
# factor (generalized Gauss-Laguerre axes), the node channel folds the
# weight into the tilt exponent (delta + alpha), and the clique-2 sweep
# corrects the alpha = 1 conjugate proposal by an independence-Metropolis
# ratio (K_ii_new / K_ii_old)^(alpha - 1). Each generalized channel is
# checked against a direct Monte Carlo of the same moment ratio, and the
# alpha = 1 path against the closed forms it must reduce to.

test_that("generalized Gauss-Laguerre rule reproduces gamma moments", {
  for(a in c(-0.5, 1, 2.5)) {
    q = bgms:::zratio_gauss_quad(48, "laguerre", glag_a = a)
    for(k in 0:4) {
      expect_lt(
        abs(sum(q$weights * q$nodes^k) / gamma(a + 1 + k) - 1),
        1e-12,
        label = paste0("a = ", a, ", moment ", k)
      )
    }
  }
  # a = 0 keeps the plain rule bit-exactly.
  q0 = bgms:::zratio_gauss_quad(48, "laguerre", glag_a = 0)
  qp = bgms:::zratio_gauss_quad(48, "laguerre")
  expect_identical(q0, qp)
})

test_that("generalized spike integral matches quadrature-free reference", {
  # The alpha != 1 Laguerre evaluation against stats::integrate on the
  # shifted w-form, and the alpha -> 1 limit against the closed Bessel form.
  delta = 0.5 * log(6)
  beta = 1.5
  cs = c(0.05, 0.5, 2, 8)
  for(alpha in c(0.5, 2)) {
    quad = bgms:::zratio_ispike(cs, delta, beta, alpha)
    ref = vapply(cs, function(cc) {
      4 * exp(-2 * beta * cc) * integrate(
        function(s) {
          w = cc + s
          w^(2 * alpha - 1) * (s * (s + 2 * cc))^delta *
            besselK(2 * beta * w, 0, expon.scaled = TRUE) *
            exp(-2 * beta * s)
        },
        0, Inf,
        rel.tol = 1e-12
      )$value
    }, 0.0)
    expect_lt(max(abs(quad / ref - 1)), 1e-8, label = paste("alpha", alpha))
  }
  near_one = bgms:::zratio_ispike(cs, delta, beta, 1 + 1e-9)
  closed = bgms:::zratio_ispike(cs, delta, beta, 1)
  expect_lt(max(abs(near_one / closed - 1)), 1e-6)
})

test_that("spike ratio at alpha != 1 matches direct Monte Carlo", {
  skip_on_cran()
  delta = 0.5 * log(6)
  beta = 1
  for(alpha in c(0.5, 2)) {
    set.seed(404)
    n = 2e6
    s1 = rgamma(n, alpha, beta)
    s2 = rgamma(n, alpha, beta)
    for(cc in c(0.4, 1.2)) {
      ratio_mc = mean(pmax(s1 * s2 - cc^2, 0)^delta) / mean((s1 * s2)^delta)
      ratio_quad = bgms:::zratio_ispike(cc, delta, beta, alpha) /
        bgms:::zratio_ispike(0, delta, beta, alpha)
      expect_lt(
        abs(ratio_mc / ratio_quad - 1), 0.02,
        label = paste0("alpha = ", alpha, ", c = ", cc)
      )
    }
  }
})

test_that("isolated-edge ratio psi0 at alpha != 1 matches direct Monte Carlo", {
  skip_on_cran()
  delta = 0.5 * log(6)
  beta = 1
  for(alpha in c(0.5, 2)) {
    for(slab in c("normal", "cauchy")) {
      pair = bgms:::zratio_pair_integrals(delta,
        sigma = 1, beta = beta,
        slab = slab, alpha = alpha
      )
      psi0 = pair$ispike(0) / pair$g(0)
      set.seed(505)
      n = 2e6
      x1 = rgamma(n, alpha, beta)
      x2 = rgamma(n, alpha, beta)
      k = if(identical(slab, "cauchy")) rcauchy(n, 0, 1) else rnorm(n, 0, 1)
      i0 = mean((x1 * x2)^delta)
      g0 = mean(ifelse(k^2 < x1 * x2, (x1 * x2 - k^2)^delta, 0))
      expect_lt(
        abs((i0 / g0) / psi0 - 1), 0.02,
        label = paste0("alpha = ", alpha, ", ", slab)
      )
    }
  }
})

test_that("node channel at alpha != 1 matches direct Monte Carlo", {
  skip_on_cran()
  delta = 0.5 * log(6)
  sigma = 1
  beta = 1
  t2 = 2 * beta * sigma^2
  for(alpha in c(0.5, 2)) {
    w_quad = bgms:::zratio_node_channel(delta, sigma, beta, alpha = alpha)
    set.seed(606)
    n = 1e6
    x = rgamma(n, delta + alpha + 1, rate = beta)
    den = mean((x + t2)^-1)
    w1 = sigma^4 * mean((x + t2)^-3) / den
    w2 = sigma^8 * mean((x + t2)^-5) / den
    expect_lt(abs(w1 / w_quad[1] - 1), 0.01, label = paste("alpha", alpha))
    expect_lt(abs(w2 / w_quad[2] - 1), 0.02, label = paste("alpha", alpha))
  }
})

test_that("bridge channel at alpha != 1 matches direct Monte Carlo", {
  skip_on_cran()
  delta = 0.5 * log(6)
  sigma = 1
  beta = 1
  t2 = 2 * beta * sigma^2
  alpha = 2
  cb_quad = bgms:::zratio_bridge_channel(delta, sigma, beta, alpha = alpha)
  set.seed(707)
  n = 2e6
  kaa = rgamma(n, alpha, beta)
  kbb = rgamma(n, alpha, beta)
  kab = rnorm(n, 0, sigma)
  d = kaa * kbb - kab^2
  pd = d > 0
  d = d[pd]
  da = d + t2 * kbb[pd]
  db = d + t2 * kaa[pd]
  w = d^delta * d / sqrt(da * db)
  den = mean(w)
  cb1 = mean(w * sigma^4 * kab[pd]^2 / (da * db)) / den
  cb2 = mean(w * sigma^8 * kab[pd]^4 / (da * db)^2) / den
  expect_lt(abs(cb1 / cb_quad[1] - 1), 0.02)
  expect_lt(abs(cb2 / cb_quad[2] - 1), 0.04)
})

test_that("clique-2 moments at alpha != 1 match an importance-sampled reference", {
  skip_on_cran()
  delta = 0.5 * log(6)
  sigma = 1
  beta = 1
  t2 = 2 * beta * sigma^2
  alpha = 2
  e2 = bgms:::zratio_clique2_moments(delta, sigma, beta,
    alpha = alpha,
    n_mc = 60000
  )
  # Importance sampling from the untilted prior: k11, k22 ~ Gamma(alpha,
  # beta), k12 ~ N(0, sigma^2), tilt weight d^delta on the PD region, then
  # the same F-measure moments as the sweep.
  set.seed(808)
  n = 2e6
  k11 = rgamma(n, alpha, beta)
  k22 = rgamma(n, alpha, beta)
  k12 = rnorm(n, 0, sigma)
  d = k11 * k22 - k12^2
  pd = d > 0
  k11 = k11[pd]
  k22 = k22[pd]
  k12 = k12[pd]
  d = d[pd]
  tilt = d^delta
  d2 = (k11 + t2) * (k22 + t2) - k12^2
  w = d / d2
  # tr(Ri Ri) and tr((Ri Ri)^2) for Ri = (K + t2 I)^{-1} in closed 2x2 form.
  tr_ri2 = ((k22 + t2)^2 + (k11 + t2)^2 + 2 * k12^2) / d2^2
  m11 = ((k22 + t2)^2 + k12^2) / d2^2
  m22 = ((k11 + t2)^2 + k12^2) / d2^2
  m12 = -k12 * (k11 + k22 + 2 * t2) / d2^2
  tr_ri4 = m11^2 + m22^2 + 2 * m12^2
  sw = sum(tilt * w)
  p1_ref = sigma^4 * sum(tilt * w * tr_ri2) / sw
  p2_ref = sigma^8 * sum(tilt * w * tr_ri4) / sw
  expect_lt(abs(e2[1] / p1_ref - 1), 0.03)
  expect_lt(abs(e2[2] / p2_ref - 1), 0.05)
})

test_that("gamma-shape constants live in their own cache cell", {
  skip_on_cran()
  d = 0.5 * log(6)
  z1 = bgms:::zratio_cell_constants(d, 0.5, 2)
  z2 = bgms:::zratio_cell_constants(d, 0.5, 2, scale_shape = 2)
  expect_identical(z1$alpha, 1)
  expect_identical(z2$alpha, 2)
  expect_false(isTRUE(all.equal(z1$addc, z2$addc)))
  expect_false(isTRUE(all.equal(z1$psi0, z2$psi0)))
  # Standardized-cell invariance holds per shape: same eta, different
  # frames.
  z3 = bgms:::zratio_cell_constants(d, 0.25, 4, scale_shape = 2)
  expect_identical(z2, z3)
})

test_that("calibrated engine tracks the block oracle at alpha = 2", {
  skip_on_cran()
  # Additive tables + warm-up calibration vs the long-run block-Gibbs
  # oracle, both at the gamma-shape cell: internal consistency of the
  # generalized constants and the independence-MH oracle sweep.
  q = 10L
  dlt = 0.5 * log(q)
  et = 1
  alpha = 2
  zc = bgms:::zratio_constants(dlt, et, alpha = alpha)

  set.seed(11)
  graphs = list()
  edges = NULL
  while(is.null(edges) || nrow(edges) < 10) {
    a = matrix(0L, q, q)
    ut = which(upper.tri(a))
    a[ut] = rbinom(length(ut), 1L, 0.3)
    a = a + t(a)
    diag(a) = 1L
    pr = which(upper.tri(matrix(0, q, q)), arr.ind = TRUE)
    for(r in sample(nrow(pr))) {
      i = pr[r, 1]
      j = pr[r, 2]
      oth = setdiff(1:q, c(i, j))
      sio = oth[a[i, oth] == 1 & a[j, oth] == 0]
      sjo = oth[a[j, oth] == 1 & a[i, oth] == 0]
      if(!length(sio) || !length(sjo)) next
      bd = max(c(
        0L, vapply(sio, function(x) sum(a[x, sjo]), 0L),
        vapply(sjo, function(x) sum(a[sio, x]), 0L)
      ))
      if(bd >= 2) {
        graphs[[length(graphs) + 1]] = a
        edges = rbind(edges, c(i, j))
      }
    }
  }
  n = nrow(edges)

  truth = vapply(seq_len(n), function(e) {
    bgms:::zratio_test_calibrated_eval(
      graphs[e], edges[e, , drop = FALSE], zc$addc, zc$tg, zc$ihat,
      zc$ghat, zc$wt, zc$psi0, dlt, et,
      seed = 900L + e, n_sweep = 2000L, burn = 50L, freeze_after = 0L,
      alpha = alpha
    )$log_zratio[1]
  }, 0.0)

  fz = ceiling(0.6 * n)
  res = bgms:::zratio_test_calibrated_eval(
    graphs, edges, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    dlt, et,
    seed = 13L, n_sweep = 300L, burn = 30L, freeze_after = fz,
    alpha = alpha
  )
  add = vapply(seq_len(n), function(e) {
    as.numeric(bgms:::zratio_test_eval(
      graphs[[e]], edges[e, , drop = FALSE], zc$addc, zc$tg, zc$ihat,
      zc$ghat, zc$wt, zc$psi0
    )$log_zratio)
  }, 0.0)

  est = as.numeric(res$log_zratio)
  expect_true(res$frozen)
  rmse = function(x) sqrt(mean((x - truth)^2))
  # Both routes must track the oracle; at these block sizes the additive
  # tables are already near-exact, so no beats-additive ordering is
  # asserted.
  expect_lt(rmse(est), 0.01)
  expect_lt(rmse(add), 0.01)
  expect_lt(max(abs(est - truth)), 0.05)
})
