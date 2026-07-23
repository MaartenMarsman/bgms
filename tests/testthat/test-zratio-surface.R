# Option-B absolute-moment surface deploy in the Z-ratio engine: mediating-block
# decomposition into disjoint CN clusters + bipartite bridge structures,
# per-component raw-poly surface evaluation (or additive below size_min), the
# summed two-moment saddle, and the alpha != 1 slab-family fence. The surface is
# SYNTHETIC (fixed coefficients), isolating the C++ deploy from any Monte-Carlo
# build. Decomposition/eval are checked against an independent R reimplementation
# on hand-built blocks with a known component structure.

zc = bgms:::zratio_constants(0.5 * log(50), 3)

# Synthetic family: log-S1 ~ 0 (exp ~ 1, unclamped), log-S2 ~ -4 (unclamped);
# wide hulls so clamps are inert except where a test narrows them.
mkfam = function(seed) {
  set.seed(seed)
  list(
    c1 = round(rnorm(9) * 0.1, 4),
    c2 = c(-4, round(rnorm(8) * 0.02, 4)),
    size_lo = 3, size_hi = 40, dens_lo = 0.1, dens_hi = 1.0,
    l1_lo = -10, l1_hi = 10, l2_lo = -10, l2_hi = 10, size_min = 3
  )
}
surface = list(cn = mkfam(11), bip = mkfam(22))

# R mirror of surface_eval_ (clamp size/dens to hull, raw quadratic in
# (log size, dens), clamp the log-moment to its range +/- 0.1).
eval_surf = function(f, s2, size, dens) {
  n = min(max(size, f$size_lo), f$size_hi)
  d = min(max(dens, f$dens_lo), f$dens_hi)
  L = log(n)
  x = c(1, L, L^2, d, d^2, L * d, L^2 * d, L * d^2, L^2 * d^2)
  cc = if(s2) f$c2 else f$c1
  lo = (if(s2) f$l2_lo else f$l1_lo) - 0.1
  hi = (if(s2) f$l2_hi else f$l1_hi) + 0.1
  exp(min(max(sum(cc * x), lo), hi))
}

# CN K4 (size 4, dens 1) + bipartite bridge (na 2, nb 3, e 4) on edge (1, 2).
make_graph1 = function() {
  q = 12
  G = matrix(0L, q, q)
  ei = function(a, b) {
    G[a, b] <<- 1L
    G[b, a] <<- 1L
  }
  for(v in 3:6) {                              # CN nodes 3..6 adjacent to both
    ei(1, v)
    ei(2, v)
  }
  for(a in 3:5) for(b in (a + 1):6) ei(a, b)   # CN-CN complete K4
  for(b in c(7, 8)) ei(1, b)                   # A-side adjacent to i only
  for(b in c(9, 10, 11)) ei(2, b)              # B-side adjacent to j only
  for(e in list(c(7, 9), c(7, 10), c(8, 10), c(8, 11))) ei(e[1], e[2])  # bridges
  G
}

surf_eval = function(G, i, j, surf, alpha = 1) {
  zratio_test_surface_eval(
    G, i, j, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    surf, zc$delta, zc$eta, slab_cauchy = FALSE, alpha = alpha
  )
}

test_that("block decomposes into CN clusters + bipartite bridges", {
  r = surf_eval(make_graph1(), 1, 2, surface)
  cn = r$comp[r$comp$family == 0, ]
  bp = r$comp[r$comp$family == 1, ]
  expect_equal(nrow(cn), 1)
  expect_equal(cn$size, 4)
  expect_equal(cn$e, 6)
  expect_equal(cn$dens, 1)
  expect_equal(cn$used_surface, 1)
  expect_equal(nrow(bp), 1)
  expect_equal(bp$size, 5)
  expect_equal(c(bp$na, bp$nb, bp$e), c(2, 3, 4))
  expect_equal(bp$dens, 4 / 6)
  expect_equal(bp$used_surface, 1)
})

test_that("per-component moments match the raw-poly surface, summed logR matches the saddle", {
  r = surf_eval(make_graph1(), 1, 2, surface)
  cn = r$comp[r$comp$family == 0, ]
  bp = r$comp[r$comp$family == 1, ]
  s1_cn = eval_surf(surface$cn, FALSE, 4, 1)
  s2_cn = eval_surf(surface$cn, TRUE, 4, 1)
  s1_bp = eval_surf(surface$bip, FALSE, 5, 4 / 6)
  s2_bp = eval_surf(surface$bip, TRUE, 5, 4 / 6)
  expect_equal(cn$s1, s1_cn, tolerance = 1e-12)
  expect_equal(cn$s2, s2_cn, tolerance = 1e-12)
  expect_equal(bp$s1, s1_bp, tolerance = 1e-12)
  expect_equal(bp$s2, s2_bp, tolerance = 1e-12)
  # summed moments + saddle
  S1 = s1_cn + s1_bp
  S2 = s2_cn + s2_bp
  ratio = zratio_test_saddle(S1, S2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0)
  expect_equal(r$S1, S1, tolerance = 1e-12)
  expect_equal(r$S2, S2, tolerance = 1e-12)
  expect_equal(r$logR, log(ratio), tolerance = 1e-12)
  # log_zratio routes through the surface branch under the alpha = 1 cell
  expect_equal(r$log_zratio, r$logR, tolerance = 1e-12)
})

test_that("alpha != 1 fences log_zratio to the additive path, bypassing the surface", {
  G = make_graph1()
  r = surf_eval(G, 1, 2, surface, alpha = 2)
  add = zratio_test_eval(
    G, matrix(c(1, 2), 1, 2), zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0
  )
  expect_equal(r$log_zratio, as.numeric(add$log_zratio[1]), tolerance = 1e-12)
  expect_gt(abs(r$log_zratio - r$logR), 1e-6)
})

test_that("components below size_min fall back to additive (exact through pairwise overlap)", {
  q = 8
  G = matrix(0L, q, q)
  ei = function(a, b) {
    G[a, b] <<- 1L
    G[b, a] <<- 1L
  }
  for(e in list(c(1, 3), c(2, 3), c(1, 4), c(2, 4), c(3, 4))) ei(e[1], e[2])  # CN pair {3,4}, e 1
  for(e in list(c(1, 5), c(2, 6), c(5, 6))) ei(e[1], e[2])                    # single bridge 5-6
  r = surf_eval(G, 1, 2, surface)
  cn = r$comp[r$comp$family == 0, ]
  bp = r$comp[r$comp$family == 1, ]
  expect_equal(cn$size, 2)
  expect_equal(cn$used_surface, 0)
  expect_equal(cn$s1, 2 * zc$addc[1] + 1 * zc$addc[3], tolerance = 1e-12)  # 2*kcn + kcc
  expect_equal(cn$s2, 2 * zc$addc[2] + 1 * zc$addc[4], tolerance = 1e-12)
  expect_equal(bp$size, 2)
  expect_equal(bp$e, 1)
  expect_equal(bp$used_surface, 0)
  expect_equal(bp$s1, zc$addc[5], tolerance = 1e-12)                       # kbr
  expect_equal(bp$s2, zc$addc[6], tolerance = 1e-12)
})

test_that("gold reference collapses to additive when every component is trivial", {
  # CN pair + single bridge: all components are trivial (size < 3), so gold
  # (per-component oracle), additive, and the surface deploy all reduce to the
  # exact additive moment -> identical logR, deterministically (no Monte Carlo).
  q = 8
  G = matrix(0L, q, q)
  ei = function(a, b) {
    G[a, b] <<- 1L
    G[b, a] <<- 1L
  }
  for(e in list(c(1, 3), c(2, 3), c(1, 4), c(2, 4), c(3, 4))) ei(e[1], e[2])
  for(e in list(c(1, 5), c(2, 6), c(5, 6))) ei(e[1], e[2])
  gold = zratio_test_gold_moments(
    G, 1, 2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    zc$delta, zc$eta, 1000L, 200L, 1L, FALSE, 1
  )
  addl = zratio_test_eval(
    G, matrix(c(1, 2), 1, 2), zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0
  )$log_zratio[1]
  surf = surf_eval(G, 1, 2, surface)
  expect_true(gold$valid)
  expect_equal(gold$logR, addl, tolerance = 1e-12)
  expect_equal(surf$logR, addl, tolerance = 1e-12)
})

test_that("gold reference is finite on a non-trivial block", {
  r = zratio_test_gold_moments(
    make_graph1(), 1, 2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    zc$delta, zc$eta, 800L, 150L, 7L, FALSE, 1
  )
  expect_true(r$valid)
  expect_true(is.finite(r$logR))
  expect_gt(r$S1, 0)
  expect_gt(r$S2, 0)
})

test_that("size / density / log-moment clamps match the R predictor", {
  # Narrow hulls so every clamp fires: size_hi 4 (bip size 5 -> 4), dens_hi 0.5
  # (dens 0.667 -> 0.5), and a degenerate log-S2 range -> exp(-3.1).
  clamped = mkfam(22)
  clamped$size_hi = 4
  clamped$dens_hi = 0.5
  clamped$l2_lo = -3
  clamped$l2_hi = -3
  surf = list(cn = surface$cn, bip = clamped)
  r = surf_eval(make_graph1(), 1, 2, surf)
  bp = r$comp[r$comp$family == 1, ]
  expect_equal(bp$s1, eval_surf(clamped, FALSE, 5, 4 / 6), tolerance = 1e-12)
  expect_equal(bp$s2, exp(-3.1), tolerance = 1e-12)   # log-moment clamp
})

test_that("the extrapolation counter fires only beyond the trained hull", {
  # The synthetic surface has size_hi = 40 (mkfam). A CN clique of size 30 stays
  # within the hull; size 50 exceeds it, so every edge whose CN block is the
  # clique deploys an extrapolated (clamped) moment.
  cn_clique = function(k) {
    q = k + 2
    G = matrix(0L, q, q)
    ei = function(a, b) {
      G[a, b] <<- 1L
      G[b, a] <<- 1L
    }
    for(v in 3:(k + 2)) {
      ei(1, v)
      ei(2, v)
    }
    for(a in 3:(k + 1)) for(b in (a + 1):(k + 2)) ei(a, b)
    G
  }
  run = function(G) {
    ed = which(upper.tri(G) & G == 1, arr.ind = TRUE)
    storage.mode(ed) = "integer"
    zratio_test_surface_batch(G, ed, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt,
      zc$psi0, surface, zc$delta, zc$eta, FALSE, 1)
  }
  r30 = run(cn_clique(30))
  expect_equal(r30$n_extrap, 0)
  expect_equal(r30$max_extrap_size, 0)

  r50 = run(cn_clique(50))
  expect_gt(r50$n_extrap, 0)
  expect_equal(r50$max_extrap_size, 50)
})
