# End-to-end accuracy gate for the Option-B surface: build the absolute-moment
# surface from real block-Gibbs anchors (zratio_build_surfaces) and check that its
# deployed per-edge log-ratio tracks the per-component gold oracle far more
# tightly than the additive kernel it replaces. This is the in-suite counterpart
# of the release cert (which runs the same comparison at q = 50); here a small
# size cap keeps the build well under a second. Caches are disabled so the build
# is fresh and never touches the user cache directory.
#
# The serial Normal builds are session-cached across tests: the core-count
# invariance test proves serial and parallel builds bit-identical, so a single
# serial build per max_size serves every comparison.

surf_cache = new.env(parent = emptyenv())

normal_surface = function(max_size) {
  key = paste0("n", max_size)
  if(is.null(surf_cache[[key]])) {
    surf_cache[[key]] = withr::with_options(
      list(
        bgms.zratio_surface_cache = FALSE,
        bgms.correction_table_cache = FALSE
      ),
      bgms:::zratio_build_surfaces(
        bgms:::zratio_constants(0.5 * log(12), 3),
        max_size = max_size, cores = 1L
      )
    )
  }
  surf_cache[[key]]
}

test_that("the built surface tracks the gold oracle far tighter than additive", {
  skip_on_cran()
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to run the surface-vs-gold accuracy cert"
  )
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

  surf = normal_surface(10L)
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

test_that("the deployed route serves the surface at a non-unit shape", {
  skip_on_cran()
  # Deployment probe. Every accuracy figure in the validation program was
  # scored on `logR`, the surface path; a fit takes `log_zratio`, the hot path.
  # Those were not the same function: the hot path carried its own alpha == 1
  # test, so a surface built and attached at a non-unit shape was silently
  # ignored and the additive kernel served instead. Component-level scoring
  # cannot certify deployment, so this asserts the deployed value IS the
  # surface value -- not merely that it differs from the additive one, which
  # would still pass under any future half-wired state that perturbs additive.
  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE
  )
  q = 12
  G = matrix(0L, q, q)
  for(v in 3:q) G[1, v] = G[v, 1] = G[2, v] = G[v, 2] = 1L
  for(a in 3:(q - 1)) for(b in (a + 1):q) G[a, b] = G[b, a] = 1L

  for(alpha in c(0.5, 2)) {
    zc = bgms:::zratio_constants(0.5 * log(12), 2, alpha = alpha)
    surf = bgms:::zratio_build_surfaces(zc, max_size = 12L, cores = 1L)
    expect_false(is.null(surf), label = paste("surface built at shape", alpha))
    res = zratio_test_surface_eval(
      G, 1, 2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
      surf, zc$delta, zc$eta, FALSE, alpha
    )
    expect_equal(
      res$log_zratio, res$logR,
      tolerance = 1e-12,
      label = sprintf("deployed route at shape %g", alpha)
    )
  }
})

test_that("anchors are drawn at the cell's own Gamma shape", {
  skip_on_cran()
  # The oracle's shape argument defaults to the exponential, so an anchor helper
  # that forgets to pass it draws at shape 1 while the constants carry the
  # cell's own. The alpha != 1 fence hides that today; it would surface as an
  # unattributable surface error the moment the fence comes off. Same block,
  # same seed, different shape: the moments have to move.
  n = 6L
  a1 = bgms:::zratio_anchor_cn(
    n, 0.9, bgms:::zratio_constants(0.5 * log(12), 2, alpha = 1),
    200L, 50L, 5L
  )
  a2 = bgms:::zratio_anchor_cn(
    n, 0.9, bgms:::zratio_constants(0.5 * log(12), 2, alpha = 2),
    200L, 50L, 5L
  )
  expect_false(is.null(a1))
  expect_false(is.null(a2))
  expect_true(abs(a1$S1 - a2$S1) / a1$S1 > 1e-6)
  expect_true(abs(a1$S2 - a2$S2) / a1$S2 > 1e-6)
})

test_that("the build fences shapes outside the validated range", {
  skip_on_cran()
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to build the gamma-shape constants cells"
  )
  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE
  )
  # Above the range the anchor oracle's independence-Metropolis step stops
  # mixing (about 1% acceptance at shape 5), so no surface is built.
  for(shape in c(2.5, 5)) {
    zc = bgms:::zratio_constants(0.5 * log(12), 2, alpha = shape)
    expect_null(bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 1L))
  }
  # Inside it the build proceeds at both validated endpoints.
  for(shape in c(0.5, 2)) {
    zc = bgms:::zratio_constants(0.5 * log(12), 2, alpha = shape)
    expect_false(is.null(bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 1L)))
  }
})

test_that("a non-unit shape gets a raised anchor budget", {
  # The independence-Metropolis row update means the same nominal budget buys
  # fewer effective sweeps, so the anchors are run longer off shape 1. The
  # multiplier was resolved by matching measured anchor spread to the shape-1
  # reference, not by 1 / acceptance.
  expect_equal(bgms:::zratio_anchor_shape_multiplier(1), 1L)
  expect_equal(bgms:::zratio_anchor_shape_multiplier(2), 2L)
  expect_equal(bgms:::zratio_anchor_shape_multiplier(0.5), 4L)
})

test_that("the fence message names the validated shapes and the reason", {
  withr::local_options(bgms.verbose = TRUE)
  zc = suppressWarnings(bgms:::zratio_constants(0.5 * log(12), 2, alpha = 15))
  # The claim is the scored points, not the interval they span: the message
  # must not read as if every shape in between had been measured.
  expect_message(
    bgms:::zratio_surface_fence_message(zc),
    "shapes 0.5, 1, 2, 3 and 5"
  )
  expect_message(bgms:::zratio_surface_fence_message(zc), "interpolated, not measured")
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
  s1 = normal_surface(8L)
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
  s1 = normal_surface(8L)
  withr::local_options(bgms.zratio_surface_psock = TRUE)
  s2 = bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 2L)
  expect_equal(s1, s2, tolerance = 1e-8)
})

test_that("the Cauchy slab builds and deploys its own surface cell", {
  skip_on_cran()
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to run the Cauchy surface deploy cert"
  )
  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE
  )
  zc_c = bgms:::zratio_constants(0.5 * log(12), 3, slab = "cauchy")
  s_n = normal_surface(10L)
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
