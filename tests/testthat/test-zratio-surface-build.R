# End-to-end accuracy gate for the Option-B surface: build the absolute-moment
# surface from real block-Gibbs anchors (build_surfaces_allmc) and check that its
# deployed per-edge log-ratio tracks the per-component gold oracle far more
# tightly than the additive kernel it replaces. This is the in-suite counterpart
# of the release cert (which runs the same comparison at q = 50); here a small
# size cap keeps the build well under a second. Caches are disabled so the build
# is fresh and never touches the user cache directory.

test_that("the built surface tracks the gold oracle far tighter than additive", {
  skip_on_cran()
  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE
  )
  zc = bgms:::zratio_constants(0.5 * log(12), 3)

  # A dense K8 common-neighbour cluster (nodes 3..10 adjacent to both 1 and 2 and
  # complete among themselves) is exactly where summing pairwise overlaps biases
  # the additive kernel; the surface is trained on such CN blocks.
  q = 10
  G = matrix(0L, q, q)
  ei = function(a, b) {
    G[a, b] <<- 1L
    G[b, a] <<- 1L
  }
  for(v in 3:10) {
    ei(1, v)
    ei(2, v)
  }
  for(a in 3:9) for(b in (a + 1):10) ei(a, b)

  surf = bgms:::build_surfaces_allmc(zc, max_size = 10L, cores = 1L)
  expect_false(is.null(surf))

  sv = zratio_test_surface_eval(
    G, 1, 2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    surf, zc$delta, zc$eta, FALSE, 1
  )
  gold = zratio_test_gold_moments(
    G, 1, 2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    zc$delta, zc$eta, 5000L, 1000L, 11L, FALSE, 1
  )
  add = zratio_test_eval(
    G, matrix(c(1, 2), 1, 2), zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0
  )$log_zratio[1]

  expect_true(gold$valid)
  err_surf = abs(sv$logR - gold$logR)
  err_add = abs(add - gold$logR)
  # Surface is within a few thousandths of gold; additive is off by ~0.01, so the
  # surface is at least twice as tight. Thresholds carry headroom over the
  # measured margin (err_surf ~ 4e-4, err_add ~ 1e-2) for MC slack.
  expect_lt(err_surf, 0.004)
  expect_gt(err_add, 0.006)
  expect_lt(err_surf, 0.4 * err_add)
})

test_that("the build fences a non-exponential precision diagonal to NULL", {
  skip_on_cran()
  zc_gamma = bgms:::zratio_constants(0.5 * log(12), 3, alpha = 2)
  expect_null(bgms:::build_surfaces_allmc(zc_gamma, max_size = 10L, cores = 1L))
})

test_that("the surface build is invariant to the core count", {
  skip_on_cran()
  skip_on_os("windows")
  # Anchors self-seed per job and rows reassemble in grid order, so the merged
  # cost-sorted dynamic schedule must return bit-identical surfaces at any
  # core count.
  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE
  )
  zc = bgms:::zratio_constants(0.5 * log(12), 3)
  s1 = bgms:::build_surfaces_allmc(zc, max_size = 8L, cores = 1L)
  s2 = bgms:::build_surfaces_allmc(zc, max_size = 8L, cores = 2L)
  expect_identical(s1, s2)
})
