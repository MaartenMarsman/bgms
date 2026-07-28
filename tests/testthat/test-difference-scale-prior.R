# Tests for the random difference-slab-scale hyperprior in bgmCompare().
#
# The difference scale s = s0 * u is sampled through a one-dimensional
# MH-within-Gibbs step whose target is the mean-1 hyperprior pi(u) times the
# slab densities of the currently included difference parameters (main-effect
# and pairwise). These tests cover the storage/summary/extractor plumbing, the
# scope and mean-1 validation, and prior-only recovery of the hyperprior.

make_two_group_ordinal = function(n = 180, p = 5, seed = 1) {
  set.seed(seed)
  list(
    x = matrix(sample(0:2, n * p, TRUE), n, p),
    y = matrix(sample(0:2, n * p, TRUE), n, p)
  )
}


test_that("the sampled difference scale is stored, summarized, and extractable", {
  d = make_two_group_ordinal(seed = 1)
  fit = bgmCompare(
    d$x, d$y,
    difference_family = "Normal",
    difference_scale_prior = gamma_prior(2, 2),
    iter = 800, warmup = 500, chains = 2, seed = 3,
    update_method = "adaptive-metropolis", display_progress = "none"
  )

  raw = fit$raw_samples
  expect_false(is.null(raw$difference_scale))
  expect_equal(length(raw$difference_scale), 2L)
  expect_equal(length(raw$difference_scale[[1]]), 800L)
  expect_true(all(unlist(raw$difference_scale) > 0))

  draws = extract_scale_draws(fit)
  expect_equal(ncol(draws), 1L)
  expect_equal(colnames(draws), "difference_scale")
  expect_true(all(draws > 0))

  s = summary(fit)
  expect_false(is.null(s$difference_scale))
  expect_true(all(c("2.5%", "97.5%", "Rhat", "n_eff") %in%
    colnames(s$difference_scale)))
})


test_that("fixed-scale compare fits expose no scale draws and error helpfully", {
  d = make_two_group_ordinal(seed = 2)
  fit = bgmCompare(
    d$x, d$y,
    iter = 400, warmup = 400, chains = 1, seed = 2,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  expect_null(fit$raw_samples$difference_scale)
  expect_null(summary(fit)$difference_scale)
  expect_error(extract_scale_draws(fit), "fixed difference slab scale")
})


test_that("difference_scale_prior validates the mean-1 and normal-slab scope", {
  d = make_two_group_ordinal(seed = 3)

  # A cauchy difference slab cannot be randomized yet.
  expect_error(
    bgmCompare(d$x, d$y,
      difference_family = "Cauchy",
      difference_scale_prior = gamma_prior(2, 2),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"
    ),
    "normal difference slab"
  )
  # Non-mean-1 gamma is rejected.
  expect_error(
    bgmCompare(d$x, d$y,
      difference_family = "Normal",
      difference_scale_prior = gamma_prior(2, 3),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"
    ),
    "mean 1"
  )
  # Standardized (eta) frame has no meaning for the multiplier; rejected.
  expect_error(
    bgmCompare(d$x, d$y,
      difference_family = "Normal",
      difference_scale_prior = gamma_prior(shape = 2),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"
    ),
    "raw frame"
  )
  # Exponential(rate = 1) has mean 1 and is accepted.
  expect_no_error(
    bgmCompare(d$x, d$y,
      difference_family = "Normal",
      difference_scale_prior = exponential_prior(rate = 1),
      iter = 50, warmup = 50, chains = 1, seed = 4,
      update_method = "adaptive-metropolis", display_progress = "none"
    )
  )
})


test_that("with the differences held out the scale recovers the mean-1 hyperprior", {
  skip_on_cran()
  # A near-zero inclusion probability for both main-effect and pairwise
  # differences keeps every difference excluded, so the scale full conditional
  # is exactly the mean-1 hyperprior throughout.
  d = make_two_group_ordinal(n = 200, p = 4, seed = 42)
  fit = bgmCompare(
    d$x, d$y,
    difference_family = "Normal",
    difference_scale_prior = gamma_prior(2, 2),
    difference_prior = bernoulli_prior(1e-3),
    main_difference_selection = TRUE,
    iter = 6000, warmup = 1500, chains = 2, seed = 7,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  ind = do.call(rbind, fit$raw_samples$indicator)
  s = do.call(c, fit$raw_samples$difference_scale)

  # The differences are held out essentially always.
  expect_lt(mean(ind), 0.02)
  # gamma(shape = 2, rate = 2): mean 1, variance 0.5, central fourth moment
  # 1.5. The scale chain is a slow-mixing MH walk (integrated autocorrelation
  # time ~ 200, so ~65 effective draws), so both moment checks bound the
  # error by 4 MCSE at the chain's own effective sample size; the target is
  # gross wiring bugs, not fine calibration.
  ac = stats::acf(s, lag.max = 500, plot = FALSE)$acf[-1]
  ess = length(s) / (1 + 2 * sum(ac[cumsum(ac < 0.01) == 0]))
  expect_lt(abs(mean(s) - 1), 4 * sqrt(0.5 / ess))
  expect_lt(abs(stats::var(s) - 0.5), 4 * sqrt((1.5 - 0.5^2) / ess))
})


test_that("the random difference scale works under NUTS", {
  d = make_two_group_ordinal(n = 150, p = 4, seed = 5)
  fit = bgmCompare(
    d$x, d$y,
    difference_family = "Normal",
    difference_scale_prior = gamma_prior(2, 2),
    iter = 300, warmup = 300, chains = 2, seed = 8,
    update_method = "nuts", display_progress = "none"
  )
  sc = fit@difference_scale_samples
  expect_equal(length(sc), 2L)
  expect_true(all(unlist(sc) > 0))
  expect_gt(length(unique(round(sc[[1]], 6))), 10L)
})
