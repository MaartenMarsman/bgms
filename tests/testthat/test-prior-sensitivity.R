# Tests for the random interaction-slab-scale hyperprior and
# prior_sensitivity_check().
#
# The scale s = s0 * u is sampled through a one-dimensional MH-within-Gibbs
# step whose target is the prior-only full conditional pi(u) times the included
# slab densities. These tests cover the storage/summary/extractor plumbing, the
# mean-1 and scope validation, prior-only recovery of the hyperprior, the
# marginalized-inclusion definition, the prior-odds convention, and agreement
# of the reweighted fixed-scale curve with brute-force refits.

test_that("the sampled scale is stored, summarized, and extractable", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:6],
    interaction_scale_prior = gamma_prior(2, 2),
    iter = 1000, warmup = 500, chains = 2, seed = 1,
    update_method = "adaptive-metropolis", display_progress = "none"
  )

  raw = fit$raw_samples
  expect_false(is.null(raw$interaction_scale))
  expect_equal(length(raw$interaction_scale), 2L)
  expect_equal(length(raw$interaction_scale[[1]]), 1000L)
  expect_true(all(unlist(raw$interaction_scale) > 0))

  draws = extract_scale_draws(fit)
  expect_equal(ncol(draws), 1L)
  expect_equal(colnames(draws), "interaction_scale")
  expect_true(all(draws > 0))

  s = summary(fit)
  expect_false(is.null(s$interaction_scale))
  expect_true(all(c("2.5%", "97.5%", "Rhat", "n_eff") %in%
    colnames(s$interaction_scale)))
})


test_that("fixed-scale fits expose no scale draws and error helpfully", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:5],
    iter = 500, warmup = 500, chains = 1, seed = 2,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  expect_null(fit$raw_samples$interaction_scale)
  expect_null(summary(fit)$interaction_scale)
  expect_error(extract_scale_draws(fit), "fixed interaction slab scale")
  expect_error(prior_sensitivity_check(fit), "fixed interaction slab scale")
})


test_that("interaction_scale_prior validates the mean-1 and normal/omrf scope", {
  data("Wenchuan", package = "bgms")
  x = Wenchuan[, 1:4]

  # Non-mean-1 gamma is rejected.
  expect_error(
    bgm(x, interaction_scale_prior = gamma_prior(2, 3),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"),
    "mean 1"
  )
  # Standardized (eta) frame has no meaning for the multiplier; rejected.
  expect_error(
    bgm(x, interaction_scale_prior = gamma_prior(shape = 2),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"),
    "raw frame"
  )
  # A cauchy slab cannot be randomized yet.
  expect_error(
    bgm(x, interaction_prior = cauchy_prior(2.5),
      interaction_scale_prior = gamma_prior(2, 2),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"),
    "normal slab"
  )
  # Continuous data routes to GGM, which is out of scope for now.
  set.seed(10)
  xc = matrix(rnorm(200 * 3), 200, 3)
  expect_error(
    bgm(xc, variable_type = "continuous",
      interaction_scale_prior = gamma_prior(2, 2),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"),
    "discrete"
  )
  # Exponential(rate = 1) has mean 1 and is accepted.
  expect_no_error(
    bgm(x, interaction_scale_prior = exponential_prior(rate = 1),
      iter = 50, warmup = 50, chains = 1, seed = 4,
      update_method = "adaptive-metropolis", display_progress = "none")
  )
})


test_that("with the edge held out the scale draws recover the mean-1 hyperprior", {
  # Two independent ordinal variables with a near-zero prior inclusion
  # probability: the edge stays excluded for the whole chain, so the scale
  # full conditional is exactly the hyperprior throughout (no transient
  # contamination from included periods).
  set.seed(42)
  n = 400
  x = cbind(sample(0:2, n, TRUE), sample(0:2, n, TRUE))
  fit = bgm(
    x, interaction_scale_prior = gamma_prior(2, 2),
    edge_prior = bernoulli_prior(1e-3),
    iter = 6000, warmup = 1500, chains = 2, seed = 3,
    edge_selection = TRUE,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  raw = fit$raw_samples
  s = do.call(c, raw$interaction_scale)
  g = do.call(rbind, raw$indicator)[, 1]

  # The edge is held out essentially always.
  expect_lt(mean(g), 0.05)
  u = s[g == 0] # s0 = 1, so u = s
  # gamma(shape = 2, rate = 2): mean 1, variance shape / rate^2 = 0.5.
  expect_equal(mean(u), 1, tolerance = 0.1)
  expect_equal(stats::var(u), 0.5, tolerance = 0.2)
})


test_that("the marginalized inclusion probability equals the indicator-draw mean", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:6],
    interaction_scale_prior = gamma_prior(2, 2),
    iter = 1500, warmup = 500, chains = 2, seed = 6,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  ps = prior_sensitivity_check(fit)
  gpool = do.call(rbind, fit$raw_samples$indicator)
  expect_equal(
    ps$edges$marginalized_inclusion_probability,
    unname(colMeans(gpool)),
    tolerance = 1e-12
  )
})


test_that("the Bayes factor divides by the per-edge prior odds, not 1/2", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:5],
    interaction_scale_prior = gamma_prior(2, 2),
    edge_prior = bernoulli_prior(0.2),
    iter = 1200, warmup = 500, chains = 2, seed = 5,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  ps = prior_sensitivity_check(fit)
  expect_true(all(abs(ps$edges$prior_inclusion_probability - 0.2) < 1e-8))

  pip = ps$edges$marginalized_inclusion_probability
  expected = (pip / (1 - pip)) / (0.2 / 0.8)
  expect_equal(ps$edges$marginalized_bf, expected, tolerance = 1e-8)
})


test_that("the reweighted fixed-scale curve matches brute-force refits", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  x = Wenchuan[, 1:8]

  fit = bgm(
    x, interaction_scale_prior = gamma_prior(2, 2),
    iter = 6000, warmup = 1500, chains = 4, seed = 11,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  # Compare inside the reliable window, at the scale posterior median.
  sx = round(stats::median(do.call(c, fit$raw_samples$interaction_scale)), 2)
  ps = prior_sensitivity_check(fit, grid = c(sx * 0.8, sx, sx * 1.25))
  gi = which.min(abs(ps$grid$scale - sx))
  expect_gte(ps$grid$ess[gi], 400)

  # Brute-force refit at the fixed scale; inclusion Bayes factor from the
  # indicator-draw mean (prior odds 1 at the default 0.5 edge prior).
  f = bgm(
    x, interaction_prior = normal_prior(scale = sx),
    iter = 6000, warmup = 1500, chains = 4, seed = 22,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  pip_fixed = colMeans(do.call(rbind, f$raw_samples$indicator))
  b = log10(pip_fixed / (1 - pip_fixed))
  a = ps$log10_bf[gi, ]

  keep = !ps$edges$saturated & is.finite(a) & is.finite(b)
  err = a[keep] - b[keep]
  # Tutorial validation reports a bulk error near 0.04. The log-odds transform
  # inflates the error for near-saturated edges, so gate median error and bias
  # on all non-saturated edges, and the RMSE on well-identified bulk edges
  # (inclusion probability away from 0 and 1) where log10 BF is stable.
  expect_lt(abs(mean(err)), 0.06)
  expect_lt(stats::median(abs(err)), 0.1)
  bulk = keep & pip_fixed > 0.05 & pip_fixed < 0.95
  bulk_err = a[bulk] - b[bulk]
  expect_lt(sqrt(mean(bulk_err^2)), 0.15)
})
