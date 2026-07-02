# --------------------------------------------------------------------------- #
# Tests for the hierarchical-graph-prior correction tables.
#
# The table math has closed forms when the per-edge tilt slope is constant c:
#   C(theta) = (1 - theta + theta * exp(c))^E,
# so the per-pair curve is f(theta) = log(1 - theta + theta * exp(c)) up to
# an additive constant, the implied local slope is c everywhere, and the
# fed(theta, d) read-off equals f(theta) in every density column. c = 0 is
# the untilted case where every curve is identically zero.
# --------------------------------------------------------------------------- #

test_that("theta grid is dense at both ends and stays inside (0, 1)", {
  g = ggm_correction_theta_grid(120L)

  expect_false(is.unsorted(g))
  expect_true(all(g > 0 & g < 1))
  expect_gte(sum(g < 0.1), 25)
  expect_gte(sum(g > 0.9), 25)
})

test_that("untilted sweep gives identically zero curves", {
  theta = ggm_correction_theta_grid(80L)
  tab = correction_table_from_edens(theta,
    edens_raw = theta,
    num_pairs = 1e6
  )

  expect_equal(tab$f, rep(0, length(theta)))
  expect_equal(tab$logC, rep(0, length(theta)))
  expect_equal(tab$fprime, rep(0, length(tab$fprime_density)),
    tolerance = 1e-12
  )
  expect_equal(max(abs(tab$fed)), 0, tolerance = 1e-12)
  expect_equal(tab$num_repaired, 0)
})

test_that("constant-slope tilt recovers the closed form", {
  c0 = -0.3
  theta = ggm_correction_theta_grid(120L)
  edens = theta * exp(c0) / (1 - theta + theta * exp(c0))
  tab = correction_table_from_edens(theta, edens, num_pairs = 1e6)

  expect_equal(tab$fprime, rep(c0, length(tab$fprime_density)),
    tolerance = 1e-10
  )
  expect_lte(length(tab$fprime_density), sum(theta <= 0.90))

  f_exact = log(1 - theta + theta * exp(c0))
  expect_equal(tab$f, f_exact - f_exact[1], tolerance = 2e-4)

  fed_exact = log(1 - tab$fed_theta + tab$fed_theta * exp(c0))
  for(j in seq_along(tab$fed_density)) {
    expect_equal(tab$fed[, j], fed_exact, tolerance = 1e-10)
  }
})

test_that("isotonic repair lifts a crash dip and reports it", {
  theta = ggm_correction_theta_grid(60L)
  edens = theta * exp(-0.3) / (1 - theta + theta * exp(-0.3))
  edens_dipped = edens
  edens_dipped[50] = 0.01

  tab = correction_table_from_edens(theta, edens_dipped, num_pairs = 1e6)

  expect_gte(tab$num_repaired, 1)
  expect_false(is.unsorted(tab$edens))
})

test_that("table build with cache round-trips and reuses the file", {
  cache_dir = file.path(tempdir(), "bgms-ctable-test")
  unlink(cache_dir, recursive = TRUE)
  old = options(bgms.correction_cache_dir = cache_dir)
  on.exit(options(old), add = TRUE)

  tab1 = ggm_correction_table(
    p = 4, n_grid = 12L, n_samples = 100L, n_warmup = 100L, n_seeds = 1L,
    update_method = "gibbs"
  )
  files = list.files(cache_dir)
  expect_length(files, 1)

  tab2 = ggm_correction_table(
    p = 4, n_grid = 12L, n_samples = 100L, n_warmup = 100L, n_seeds = 1L,
    update_method = "gibbs"
  )
  expect_identical(tab1, tab2)

  expect_identical(tab1$cell$q, 4L)
  expect_identical(tab1$cell$slab_family, "cauchy")
  expect_equal(tab1$cell$eta, 1)
  expect_equal(tab1$cell$delta, 0.5 * log(4))
})

test_that("gibbs and adaptive-metropolis sweeps agree on edge density", {
  skip_on_cran()

  theta = c(0.2, 0.5, 0.8)
  args = list(
    p = 6, theta = theta, delta = 0.5 * log(6),
    interaction_prior = cauchy_prior(scale = 2.5),
    precision_scale_prior = gamma_prior(shape = 1, eta = 1),
    n_samples = 3000L, n_warmup = 500L, n_seeds = 2L
  )
  sweep_gibbs = do.call(
    sweep_prior_edge_density,
    c(args, update_method = "gibbs")
  )
  sweep_am = do.call(
    sweep_prior_edge_density,
    c(args, update_method = "adaptive-metropolis")
  )

  expect_lt(max(abs(sweep_gibbs$edens_raw - sweep_am$edens_raw)), 0.02)
})
