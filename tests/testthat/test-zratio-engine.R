# Tests for the hierarchical-spec per-edge Z-ratio engine: parity against the
# reference implementation (fixture generated from the z_graph_prior deployed
# kernel), the fit-time constant builders, and engine invariants.

fixture_path = testthat::test_path("fixtures", "zratio_reference.rds")

test_that("engine reproduces the reference log Z-ratios (all cells/variants)", {
  skip_if_not(file.exists(fixture_path), "zratio_reference.rds not generated")
  fx = readRDS(fixture_path)
  for(cell in fx) {
    addc13 = c(cell$addc6, cell$fc_coef, 1)
    addc23 = c(addc13, cell$hull)
    for(gn in names(cell$graphs)) {
      g = cell$graphs[[gn]]
      eu = which(upper.tri(g), arr.ind = TRUE)
      edges = cbind(eu[, 1], eu[, 2])
      for(variant in c("base", "direct", "clamp")) {
        addc = switch(variant,
          base = cell$addc6,
          direct = addc13,
          clamp = addc23
        )
        res = zratio_test_eval(
          g, edges, addc, cell$tg, cell$ihat, cell$ghat, cell$wt, cell$psi0
        )
        # The hull box (addc[13..22]) is inert since the deployed correction
        # now applies to every maxbd >= 2 block: the "clamp" variant collapses
        # onto the "direct" (correction-everywhere) reference.
        ref_variant = if(variant == "clamp") "direct" else variant
        ref = cell$evals[[paste(gn, ref_variant, sep = "_")]]
        expect_equal(
          as.numeric(res$log_zratio), unname(ref),
          tolerance = 1e-12,
          label = sprintf(
            "delta=%.3f sigma=%g beta=%g %s %s",
            cell$delta, cell$sigma, cell$beta, gn, variant
          )
        )
      }
    }
  }
})

test_that("fit-time constant builders match the reference builders", {
  skip_if_not(file.exists(fixture_path), "zratio_reference.rds not generated")
  fx = readRDS(fixture_path)
  # The builders run in the standardized cell (sigma = 1), so match against
  # the two sigma = 1 reference cells; eta = sigma * beta = beta there.
  for(cell in fx[c(1, 4)]) {
    zc = bgms:::zratio_constants(cell$delta, cell$sigma * cell$beta)
    expect_equal(zc$addc, cell$addc6, tolerance = 1e-8)
    expect_equal(zc$psi0, cell$psi0, tolerance = 1e-8)
    expect_equal(zc$tg, cell$tg, tolerance = 1e-12)
    expect_equal(zc$wt, cell$wt, tolerance = 1e-12)
    expect_equal(zc$ihat, cell$ihat, tolerance = 1e-8)
    expect_equal(zc$ghat, cell$ghat, tolerance = 1e-8)
  }
})

test_that("Gauss quadrature nodes integrate known moments exactly", {
  gl = bgms:::zratio_gauss_quad(24, "laguerre")
  expect_equal(sum(gl$weights), 1, tolerance = 1e-12) # int e^-x
  expect_equal(sum(gl$weights * gl$nodes^3), 6, tolerance = 1e-9) # Gamma(4)
  gh = bgms:::zratio_gauss_quad(24, "hermite")
  expect_equal(sum(gh$weights), sqrt(pi), tolerance = 1e-12)
  expect_equal(sum(gh$weights * gh$nodes^2), sqrt(pi) / 2, tolerance = 1e-9)
  lg = bgms:::zratio_gauss_quad(24, "legendre")
  expect_equal(sum(lg$weights), 2, tolerance = 1e-12)
  expect_equal(sum(lg$weights * lg$nodes^2), 2 / 3, tolerance = 1e-12)
})

test_that("I_spike matches its delta = 0 Bessel identity", {
  cg = seq(0.2, 8, length.out = 12)
  # At delta = 0, beta = 0.5: I_spike(c) = 4 c K_1(c).
  expect_equal(
    bgms:::zratio_ispike(cg, delta = 0, beta = 0.5),
    4 * cg * besselK(cg, 1),
    tolerance = 1e-10
  )
})

test_that("engine invariants: state-invariance, isolated edge, precompute", {
  skip_if_not(file.exists(fixture_path), "zratio_reference.rds not generated")
  cell = readRDS(fixture_path)[[2]]
  g = cell$graphs$er3
  edges = cbind(1L, 2L)

  # The ratio is invariant to the toggled edge's own state.
  g_on = g
  g_on[1, 2] = g_on[2, 1] = 1L
  g_off = g
  g_off[1, 2] = g_off[2, 1] = 0L
  r_on = zratio_test_eval(
    g_on, edges, cell$addc6, cell$tg, cell$ihat, cell$ghat, cell$wt, cell$psi0
  )
  r_off = zratio_test_eval(
    g_off, edges, cell$addc6, cell$tg, cell$ihat, cell$ghat, cell$wt, cell$psi0
  )
  expect_identical(r_on$log_zratio, r_off$log_zratio)

  # No mediating structure: the isolated-edge ratio psi0.
  g_empty = diag(1L, 10)
  r_iso = zratio_test_eval(
    g_empty, edges, cell$addc6, cell$tg, cell$ihat, cell$ghat, cell$wt,
    cell$psi0
  )
  expect_equal(r_iso$log_zratio[1], log(cell$psi0), tolerance = 1e-14)

  # Preloading the count table removes every cache miss and preserves values.
  eu = which(upper.tri(g), arr.ind = TRUE)
  all_edges = cbind(eu[, 1], eu[, 2])
  lazy = zratio_test_eval(
    g, all_edges, cell$addc6, cell$tg, cell$ihat, cell$ghat, cell$wt,
    cell$psi0
  )
  pre = zratio_test_precompute(
    g, all_edges, cell$addc6, cell$tg, cell$ihat, cell$ghat, cell$wt,
    cell$psi0, 10L, 30L
  )
  expect_equal(as.numeric(pre$log_zratio), as.numeric(lazy$log_zratio), tolerance = 1e-14)
  expect_true(pre$n_miss == 0)
})

test_that("online calibrator: oracle taper, freeze pack, accuracy gain", {
  skip_on_cran()
  skip_if_not(file.exists(fixture_path), "zratio_reference.rds not generated")
  q = 14L
  dlt = 0.5 * log(30)
  et = 1
  zc = bgms:::zratio_constants(dlt, et)

  set.seed(4)
  graphs = list()
  edges = NULL
  while(is.null(edges) || nrow(edges) < 24) {
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
      seed = 500L + e, n_sweep = 2000L, burn = 50L, freeze_after = 0L
    )$log_zratio[1]
  }, 0.0)

  fz = ceiling(0.6 * n)
  res = bgms:::zratio_test_calibrated_eval(
    graphs, edges, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    dlt, et,
    seed = 7L, n_sweep = 300L, burn = 30L, freeze_after = fz
  )
  add = vapply(seq_len(n), function(e) {
    as.numeric(bgms:::zratio_test_eval(
      graphs[[e]], edges[e, , drop = FALSE], zc$addc, zc$tg, zc$ihat,
      zc$ghat, zc$wt, zc$psi0
    )$log_zratio)
  }, 0.0)

  est = as.numeric(res$log_zratio)
  expect_true(res$frozen)
  expect_length(res$addc, 23) # coef + hull packed on freeze
  expect_lte(res$n_oracle, fz) # oracle only during warm-up
  rmse = function(x) sqrt(mean((x - truth)^2))
  expect_lt(rmse(est), rmse(add)) # correction beats additive-only
  expect_lt(max(abs(est - truth)), 0.02)
})
