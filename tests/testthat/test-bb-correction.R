# --------------------------------------------------------------------------- #
# Tests for the corrected beta-bernoulli theta draw.
#
# With logC constant the corrected conditional is the conjugate Beta, so the
# inverse-CDF grid draw must reproduce rbeta up to grid discretization. With
# the constant-slope tilt logC(theta) = E * log(1 - theta + theta * exp(c)),
# the corrected conditional has no standard form and the draw is checked
# against numerical quadrature.
# --------------------------------------------------------------------------- #

test_that("logC interpolation is linear inside and extends linearly outside", {
  theta_grid = c(0.1, 0.2, 0.4, 0.8)
  logC = c(0, -1, -1.5, -3.5)

  eval = c(0.05, 0.1, 0.15, 0.3, 0.4, 0.9)
  out = test_correction_logC_interp(theta_grid, logC, eval)

  expect_equal(out[2], 0)
  expect_equal(out[3], -0.5)
  expect_equal(out[4], -1.25)
  expect_equal(out[5], -1.5)
  # extrapolation with boundary slopes: below with slope -10, above with -5
  expect_equal(out[1], 0 + (-1 - 0) / 0.1 * (0.05 - 0.1))
  expect_equal(out[6], -3.5 + (-3.5 - (-1.5)) / 0.4 * (0.9 - 0.8))
})

test_that("corrected draw with flat logC matches the conjugate Beta", {
  theta_grid = seq(0.001, 0.999, length.out = 50)
  logC = rep(0, 50)

  draws = test_corrected_bb_theta_draw(
    a_post = 5, b_post = 10,
    theta_grid = theta_grid, logC = logC,
    n_draws = 6000, seed = 11, theta_init = 0.3
  )

  probs = seq(0.05, 0.95, by = 0.05)
  expect_lt(max(abs(quantile(draws, probs) - qbeta(probs, 5, 10))), 0.012)
  expect_lt(abs(mean(draws) - 5 / 15), 0.005)
})

test_that("corrected draw matches quadrature under a constant-slope tilt", {
  E = 190
  c0 = -0.3
  theta_grid = seq(0.001, 0.999, length.out = 200)
  logC = E * log(1 - theta_grid + theta_grid * exp(c0))

  a_post = 20
  b_post = 170
  draws = test_corrected_bb_theta_draw(
    a_post = a_post, b_post = b_post,
    theta_grid = theta_grid, logC = logC,
    n_draws = 6000, seed = 12, theta_init = 0.1
  )

  # numerical quadrature of the corrected conditional
  grid = seq(1e-6, 1 - 1e-6, length.out = 20000)
  logd = (a_post - 1) * log(grid) + (b_post - 1) * log1p(-grid) -
    E * log(1 - grid + grid * exp(c0))
  d = exp(logd - max(logd))
  d = d / sum(d)
  mean_exact = sum(grid * d)
  cdf = cumsum(d)
  q_exact = sapply(c(0.25, 0.5, 0.75), function(p) grid[which(cdf >= p)[1]])

  expect_lt(abs(mean(draws) - mean_exact), 0.003)
  expect_lt(
    max(abs(quantile(draws, c(0.25, 0.5, 0.75)) - q_exact)), 0.005
  )

  # the tilt must actually move the draw: the conjugate mean is biased here
  expect_gt(abs(mean(draws) - a_post / (a_post + b_post)), 0.005)
})

# --------------------------------------------------------------------------- #
# Prior-only identity: summing the corrected joint over all graphs returns
# the Beta hyperprior for theta exactly, because 1/C(theta) cancels. The
# uncorrected chain leaves the theta marginal proportional to
# p(theta) C(theta), which is biased sparse. The identity tests share one
# correction-table cache so the deployed-default cell at p = 5 builds once.
# --------------------------------------------------------------------------- #

ct_identity_cache = file.path(tempdir(), "bgms-ctable-identity")

test_that("corrected prior-only chain returns the Beta hyperprior", {
  skip_on_cran()
  old = options(bgms.correction_cache_dir = ct_identity_cache)
  on.exit(options(old), add = TRUE)

  draws = sample_ggm_prior(
    p = 5, n_samples = 4000, n_warmup = 500,
    seed = 21, verbose = FALSE, spec = "joint",
    update_method = "gibbs", edge_prior = "beta-bernoulli"
  )

  expect_length(draws$theta, 4000)
  expect_lt(abs(mean(draws$theta) - 0.5), 0.02)
  ks = suppressWarnings(stats::ks.test(draws$theta, "punif"))
  expect_lt(unname(ks$statistic), 0.05)
})

test_that("uncorrected prior-only chain misses the hyperprior", {
  skip_on_cran()

  draws = sample_ggm_prior(
    p = 5, n_samples = 4000, n_warmup = 500,
    seed = 22, verbose = FALSE, spec = "joint",
    update_method = "gibbs", edge_prior = "beta-bernoulli",
    apply_correction = FALSE
  )

  expect_lt(mean(draws$theta), 0.46)
  ks = suppressWarnings(stats::ks.test(draws$theta, "punif"))
  expect_gt(unname(ks$statistic), 0.08)
})

test_that("bgm() with a beta-bernoulli prior applies the correction and stores theta", {
  skip_on_cran()
  old = options(bgms.correction_cache_dir = ct_identity_cache)
  on.exit(options(old), add = TRUE)

  set.seed(5)
  x = matrix(rnorm(60 * 5), 60, 5)
  fit = suppressMessages(bgm(x,
    variable_type = "continuous",
    edge_prior = beta_bernoulli_prior(1, 1),
    update_method = "gibbs",
    iter = 400, warmup = 300, chains = 1, cores = 1,
    display_progress = "none", verbose = FALSE
  ))

  theta = fit$inclusion_parameter_samples
  expect_length(theta, 1)
  expect_length(theta[[1]], 400)
  expect_true(all(theta[[1]] > 0 & theta[[1]] < 1))
})
