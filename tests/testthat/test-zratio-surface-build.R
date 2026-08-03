# End-to-end accuracy gate for the Option-B surface: build the absolute-moment
# surface from real block-Gibbs anchors (zratio_build_surfaces) and check that its
# deployed per-edge log-ratio tracks the per-component gold oracle far more
# tightly than the additive kernel it replaces. This is the in-suite counterpart
# of the release cert (which runs the same comparison at q = 50); here a small
# size cap keeps the build well under a second. Caches are disabled so the build
# is fresh and never touches the user cache directory.
#
# The three gated blocks -- surface-vs-gold, the shape fence, and the Cauchy
# build-and-deploy -- are the heavy end-to-end build machinery and run in the
# weekly certification tier (T2, BGMS_RUN_CERTIFICATION). The tier contract
# keeps surface-vs-gold SINGLE CELLS nightly; these are builds, not cells.
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
  skip_unless_certification()
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

  # The claim is that two routes over one surface return the same value, which
  # does not depend on how well that surface is anchored: an eight-variable hull
  # leaves the ten-node mediating block extrapolating past its boundary, and the
  # two routes have to agree there too. Anchoring the full block instead costs
  # three times as long and asserts nothing further.
  for(alpha in c(0.5, 2)) {
    zc = bgms:::zratio_constants(0.5 * log(12), 2, alpha = alpha)
    surf = bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 1L)
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
  skip_unless_certification()
  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE
  )
  # The deployment range is [.zratio_surface_shape_lo, .zratio_surface_shape_hi]
  # and zratio_build_surfaces is its single owner, so the contract is read off
  # those two constants rather than restated as literals.
  lo = bgms:::.zratio_surface_shape_lo
  hi = bgms:::.zratio_surface_shape_hi
  expect_equal(c(lo, hi), c(0.5, 10))

  # Inside the range the build proceeds -- at both endpoints and at interior
  # shapes, which are interpolated rather than separately scored.
  for(shape in c(lo, 2.5, 5, hi)) {
    zc = bgms:::zratio_constants(0.5 * log(12), 2, alpha = shape)
    expect_false(
      is.null(bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 1L)),
      label = sprintf("surface built at shape %g", shape)
    )
  }
  # Outside it no surface is built and the engine keeps the route the fence
  # assigns: the isolated-edge ratio above the range, the additive path below.
  for(shape in c(0.25, 12)) {
    zc = suppressWarnings(
      bgms:::zratio_constants(0.5 * log(12), 2, alpha = shape)
    )
    expect_null(
      bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 1L),
      label = sprintf("surface at shape %g", shape)
    )
  }
})

test_that("an analysis too small to anchor either family builds no surface", {
  # A bipartite bridge needs 2 + 2 nodes, so the bipartite anchor grid starts at
  # size 4 and filters to nothing at a cap of 3 or less. Assigning the family
  # tag into that empty job table used to abort the fit ("replacement has 1 row,
  # data has 0"), which made every hierarchical fit at 2 or 3 variables an error.
  expect_true(bgms:::zratio_anchor_grids_empty(2L))
  expect_true(bgms:::zratio_anchor_grids_empty(3L))
  expect_false(bgms:::zratio_anchor_grids_empty(4L))
  expect_equal(nrow(bgms:::zratio_anchor_grids(3L)$bip), 0L)

  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE
  )
  # The guard returns before any anchor runs, and only the cell's shape is read
  # on the way there, so it is supplied directly; building the cell's constants
  # would cost seconds and nothing else is used.
  zc = list(alpha = 1)
  for(cap in c(2L, 3L)) {
    expect_null(bgms:::zratio_build_surfaces(zc, max_size = cap, cores = 1L))
  }
  # The route message says the surface is unnecessary here, not that its build
  # failed.
  expect_message(
    bgms:::zratio_surface_fence_message(list(alpha = 1, eta = 1), size = 3L),
    "none is needed"
  )
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
  # Below the range, where the additive kernel is what actually serves. A shape
  # ABOVE the range no longer reaches the additive path at all -- it routes to
  # the isolated-edge ratio, and its wording is pinned in
  # test-zratio-isolated-edge-routing.R.
  # The wording is a function of the cell's shape and rate alone, so those two
  # fields are supplied directly; building the cell's constants would cost
  # seconds and none of them are read.
  zc = list(alpha = 0.25, eta = 2)
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
  # The workers load the INSTALLED bgms namespace, so this block needs bgms on
  # a library path -- which is the R CMD check situation the comment above
  # describes, and R-CMD-check.yaml runs that on five platforms on every push
  # and PR. Under a bare devtools::test() on a machine that has never installed
  # bgms the workers cannot load it and the block has nothing to compare, so it
  # skips rather than erroring. Asking a worker is the exact question; anything
  # read in this session is confounded by pkgload's shims.
  probe = parallel::makePSOCKcluster(1L)
  on.exit(parallel::stopCluster(probe), add = TRUE)
  installed = isTRUE(unlist(parallel::clusterEvalQ(
    probe, requireNamespace("bgms", quietly = TRUE)
  )))
  skip_if(!installed, "bgms is not installed; PSOCK workers cannot load it")

  zc = bgms:::zratio_constants(0.5 * log(12), 3)
  s1 = normal_surface(8L)
  withr::local_options(bgms.zratio_surface_psock = TRUE)
  s2 = bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 2L)
  expect_equal(s1, s2, tolerance = 1e-8)
})

test_that("the Cauchy slab builds and deploys its own surface cell", {
  skip_on_cran()
  skip_unless_certification()
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


# ---- The build runs at the width it was asked for (F-104) --------------------

test_that("a worker count that is not one usable number collapses to one", {
  # parallel::detectCores() is documented to return NA when it cannot tell, and
  # min(k, NA) is NA -- which would then be handed to a cluster constructor,
  # where it is not a width but an error or a default. Every shape that is not
  # a single usable number resolves to 1 rather than travelling on.
  expect_identical(bgms:::normalize_parallel_cores(NA_integer_), 1L)
  expect_identical(bgms:::normalize_parallel_cores(integer(0)), 1L)
  expect_identical(bgms:::normalize_parallel_cores(0L), 1L)
  expect_identical(bgms:::normalize_parallel_cores(-3L), 1L)
  expect_identical(bgms:::normalize_parallel_cores("nonsense"), 1L)
  expect_identical(bgms:::normalize_parallel_cores(1L), 1L)

  local_mocked_bindings(
    detectCores = function(...) NA_integer_, .package = "parallel"
  )
  expect_identical(bgms:::normalize_parallel_cores(2L), 2L)
})

test_that("the surface build refuses a cluster that is not the width asked for", {
  skip_on_cran()
  # The constructor's own default is getOption("mc.cores", 2L), so a worker
  # count that failed to reach it does not announce itself: the build simply
  # runs at a width nobody asked for. The builder reads the cluster back, and
  # this is the read-back firing -- a cluster of the wrong length, however it
  # got that way, stops the build with both numbers named. A stand-in cluster
  # is enough: the check is on its length, and no job is ever dispatched.
  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE,
    bgms.zratio_surface_psock = TRUE
  )
  stand_in = structure(list("node"), class = c("SOCKcluster", "cluster"))
  local_mocked_bindings(
    makePSOCKcluster = function(...) stand_in,
    stopCluster = function(...) invisible(NULL),
    .package = "parallel"
  )
  zc = bgms:::zratio_constants(0.5 * log(12), 3)
  expect_error(
    bgms:::zratio_build_surfaces(zc, max_size = 8L, cores = 2L),
    "asked for 2 PSOCK worker\\(s\\) and got 1"
  )
})
