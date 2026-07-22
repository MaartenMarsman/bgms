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
  ei = function(a, b) { G[a, b] <<- 1L; G[b, a] <<- 1L }
  for(v in 3:6) { ei(1, v); ei(2, v) }         # CN nodes 3..6 adjacent to both
  for(a in 3:5) for(b in (a + 1):6) ei(a, b)   # CN-CN complete K4
  ei(1, 7); ei(1, 8)                           # A-side adjacent to i only
  ei(2, 9); ei(2, 10); ei(2, 11)               # B-side adjacent to j only
  ei(7, 9); ei(7, 10); ei(8, 10); ei(8, 11)    # bridges
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
  ei = function(a, b) { G[a, b] <<- 1L; G[b, a] <<- 1L }
  ei(1, 3); ei(2, 3); ei(1, 4); ei(2, 4); ei(3, 4)   # CN pair {3,4}, e 1
  ei(1, 5); ei(2, 6); ei(5, 6)                        # single bridge 5-6
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
