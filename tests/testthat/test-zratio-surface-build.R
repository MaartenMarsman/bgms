# End-to-end accuracy gate for the Option-B surface: build the absolute-moment
# surface from real block-Gibbs anchors (zratio_build_surfaces) and check that its
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

  surf = bgms:::zratio_build_surfaces(zc, max_size = 10L, cores = 1L)
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
  expect_null(bgms:::zratio_build_surfaces(zc_gamma, max_size = 10L, cores = 1L))
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
  s1 = bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 1L)
  s2 = bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 2L)
  expect_identical(s1, s2)
})

test_that("the socket-cluster build path matches the serial build", {
  skip_on_cran()
  # The Windows branch (PSOCK workers instead of forks), forced cross-platform
  # via the override option. Workers load the INSTALLED bgms namespace, which
  # on a dev machine can be a different binary than the loaded dev build (a
  # last-ULP anchor difference amplifies to ~1e-11 in the fitted
  # coefficients), so this compares at a tight tolerance rather than
  # expect_identical; under R CMD check both sides are the same installed
  # binary. A scheduling or seeding bug would differ at O(1), not O(1e-8).
  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE
  )
  zc = bgms:::zratio_constants(0.5 * log(12), 3)
  s1 = bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 1L)
  withr::local_options(bgms.zratio_surface_psock = TRUE)
  s2 = bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 2L)
  expect_equal(s1, s2, tolerance = 1e-8)
})

test_that("the Cauchy slab builds and deploys its own surface cell", {
  skip_on_cran()
  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE
  )
  zc_n = bgms:::zratio_constants(0.5 * log(12), 3)
  zc_c = bgms:::zratio_constants(0.5 * log(12), 3, slab = "cauchy")
  s_n = bgms:::zratio_build_surfaces(zc_n, max_size = 10L, cores = 2L)
  s_c = bgms:::zratio_build_surfaces(zc_c, max_size = 10L, cores = 2L)
  expect_false(is.null(s_c))
  # The cells differ: a Cauchy fit must not be served the Normal surface.
  expect_gt(max(abs(s_c$cn$c1 - s_n$cn$c1)), 0.01)

  # Deployed Cauchy value on the dense K8 CN block against a pooled gold
  # reference. The Cauchy block oracle mixes slowly (gold spread ~0.014
  # across seeds at this budget), so this is a wiring gate at loose
  # tolerance, not an accuracy cert; the tight Cauchy gates (eta = 1
  # surface 0.0033 nats vs additive 0.0094) are in the migration record.
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
  sv = zratio_test_surface_eval(
    G, 1, 2, zc_c$addc, zc_c$tg, zc_c$ihat, zc_c$ghat, zc_c$wt, zc_c$psi0,
    s_c, zc_c$delta, zc_c$eta, TRUE, 1
  )
  gold = vapply(c(7L, 11L, 21L), function(sd) {
    zratio_test_gold_moments(
      G, 1, 2, zc_c$addc, zc_c$tg, zc_c$ihat, zc_c$ghat, zc_c$wt, zc_c$psi0,
      zc_c$delta, zc_c$eta, 5000L, 1000L, sd, TRUE, 1
    )$logR
  }, numeric(1))
  expect_lt(abs(sv$logR - mean(gold)), 0.04)
})
