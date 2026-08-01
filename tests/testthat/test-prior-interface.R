# ==============================================================================
# Tests for Prior Class Interface
# ==============================================================================
#
# Verifies that all prior constructors work correctly with bgm() and
# bgmCompare() across all model types (ordinal, continuous, mixed,
# blume-capel).
#
# Tests cover:
#   1. Prior constructor output structure
#   2. bgm() runs without error for each prior x model_type combination
#   3. bgmCompare() runs without error for each prior combination
#   4. Backward compatibility of deprecated scalar parameters
# ==============================================================================


# ==============================================================================
# 1. Prior Constructor Tests
# ==============================================================================

test_that("cauchy_prior creates valid prior object (single class)", {
  p = cauchy_prior(scale = 2.5)
  expect_identical(class(p), "bgms_parameter_prior")
  expect_equal(p$family, "cauchy")
  expect_equal(p$hyper.parameters$scale, 2.5)
})

test_that("normal_prior creates valid prior object (single class)", {
  p = normal_prior(scale = 0.5)
  expect_identical(class(p), "bgms_parameter_prior")
  expect_equal(p$family, "normal")
  expect_equal(p$hyper.parameters$scale, 0.5)
})

test_that("beta_prime_prior creates valid prior object (single class)", {
  p = beta_prime_prior(alpha = 1, beta = 1)
  expect_identical(class(p), "bgms_parameter_prior")
  expect_equal(p$family, "beta-prime")
  expect_equal(p$hyper.parameters$alpha, 1)
  expect_equal(p$hyper.parameters$beta, 1)
})

test_that("normal_prior works for threshold_prior argument", {
  p = normal_prior(scale = 2)
  expect_s3_class(p, "bgms_parameter_prior")
  expect_equal(p$family, "normal")
  expect_equal(p$hyper.parameters$scale, 2)
  # Can be unpacked as threshold prior
  tp = unpack_threshold_prior(p)
  expect_equal(tp$threshold_prior_type, "normal")
  expect_equal(tp$threshold_scale, 2)
})

test_that("all parameter priors are interchangeable for interaction_prior and threshold_prior", {
  # cauchy_prior works for threshold_prior
  pp = unpack_threshold_prior(cauchy_prior(scale = 2))
  expect_equal(pp$threshold_prior_type, "cauchy")
  expect_equal(pp$threshold_scale, 2)

  # beta_prime_prior works for interaction_prior
  pp = unpack_interaction_prior(beta_prime_prior(1, 1))
  expect_equal(pp$interaction_prior_type, "beta-prime")
  expect_equal(pp$interaction_alpha, 1)
  expect_equal(pp$interaction_beta, 1)
})

test_that("bernoulli_prior creates valid prior object", {
  p = bernoulli_prior(0.3)
  expect_s3_class(p, "bgms_indicator_prior")
  expect_equal(p$family, "Bernoulli")
  expect_equal(p$hyper.parameters$inclusion_probability, 0.3)
})

test_that("beta_bernoulli_prior creates valid prior object", {
  p = beta_bernoulli_prior(alpha = 2, beta = 5)
  expect_s3_class(p, "bgms_indicator_prior")
  expect_equal(p$family, "Beta-Bernoulli")
  expect_equal(p$hyper.parameters$alpha, 2)
  expect_equal(p$hyper.parameters$beta, 5)
})

test_that("sbm_prior creates valid prior object", {
  p = sbm_prior()
  expect_s3_class(p, "bgms_indicator_prior")
  expect_equal(p$family, "Stochastic-Block")
})

test_that("gamma_prior creates valid prior object", {
  p = gamma_prior(shape = 2, rate = 0.5)
  expect_s3_class(p, "bgms_scale_prior")
  expect_equal(p$family, "gamma")
  expect_equal(p$hyper.parameters$shape, 2)
  expect_equal(p$hyper.parameters$rate, 0.5)
})

test_that("exponential_prior creates valid prior object", {
  p = exponential_prior(rate = 2)
  expect_s3_class(p, "bgms_scale_prior")
  expect_equal(p$family, "exponential")
  expect_equal(p$hyper.parameters$rate, 2)
})

test_that("gamma_prior supports the standardized frame via eta", {
  p = gamma_prior(shape = 2, eta = 0.5)
  expect_equal(p$hyper.parameters$eta, 0.5)
  expect_true(is.na(p$hyper.parameters$rate))

  # Default is the standardized frame with eta = 1
  d = gamma_prior()
  expect_equal(d$hyper.parameters$eta, 1)
  expect_true(is.na(d$hyper.parameters$rate))

  # Raw frame stores the rate and leaves eta unset
  r = gamma_prior(shape = 2, rate = 0.5)
  expect_equal(r$hyper.parameters$rate, 0.5)
  expect_true(is.na(r$hyper.parameters$eta))

  expect_error(gamma_prior(shape = 1, rate = 1, eta = 1), "not both")
  expect_error(gamma_prior(eta = -1), "positive")
})

test_that("exponential_prior supports the standardized frame via eta", {
  p = exponential_prior(eta = 2)
  expect_equal(p$hyper.parameters$eta, 2)
  expect_true(is.na(p$hyper.parameters$rate))

  d = exponential_prior()
  expect_equal(d$hyper.parameters$eta, 1)

  expect_error(exponential_prior(rate = 1, eta = 1), "not both")
  expect_error(exponential_prior(eta = Inf), "positive and finite")
})

test_that("resolve_scale_rate derives the raw rate and guards scale-free slabs", {
  # Raw frame passes through untouched
  expect_equal(bgms:::resolve_scale_rate(0.7, NA_real_, 4), 0.7)
  # Standardized frame: rate = eta / slab scale
  expect_equal(bgms:::resolve_scale_rate(NA_real_, 1, 4), 0.25)
  # No slab scale (e.g. beta-prime interaction prior): eta cannot resolve
  expect_error(
    bgms:::resolve_scale_rate(NA_real_, 1, NA_real_),
    "scale parameter"
  )
})


# ==============================================================================
# 2. bgm() with Prior Objects <U+2014> Ordinal
# ==============================================================================

test_that("bgm ordinal works with cauchy_prior", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    interaction_prior = cauchy_prior(scale = 2.5),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  expect_true(all(is.finite(do.call(rbind, fit$raw_samples$pairwise))))
})

test_that("bgm ordinal works with normal_prior", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    interaction_prior = normal_prior(scale = 1),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  expect_true(all(is.finite(do.call(rbind, fit$raw_samples$pairwise))))
})

test_that("bgm ordinal works with normal_prior", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    threshold_prior = normal_prior(scale = 1),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  expect_true(all(is.finite(do.call(rbind, fit$raw_samples$main))))
})

test_that("bgm ordinal works with both normal priors", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    interaction_prior = normal_prior(scale = 0.5),
    threshold_prior = normal_prior(scale = 0.5),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
})

test_that("bgm ordinal works with normal priors + edge selection", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    interaction_prior = normal_prior(scale = 0.5),
    threshold_prior = normal_prior(scale = 0.5),
    edge_prior = bernoulli_prior(0.5),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
})

test_that("bgm ordinal works with normal priors + beta_bernoulli edge prior", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    interaction_prior = normal_prior(scale = 0.5),
    threshold_prior = normal_prior(scale = 0.5),
    edge_prior = beta_bernoulli_prior(1, 1),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
})


# ==============================================================================
# 3. bgm() with Prior Objects <U+2014> GGM (Continuous)
# ==============================================================================

test_that("bgm ggm works with normal_prior", {
  set.seed(42)
  Y = as.data.frame(matrix(rnorm(200), nrow = 50, ncol = 4))
  fit = bgm(Y,
    variable_type = "continuous",
    interaction_prior = normal_prior(scale = 1),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  expect_true(all(is.finite(do.call(rbind, fit$raw_samples$pairwise))))
})

test_that("bgm ggm works with cauchy_prior", {
  set.seed(42)
  Y = as.data.frame(matrix(rnorm(200), nrow = 50, ncol = 4))
  fit = bgm(Y,
    variable_type = "continuous",
    interaction_prior = cauchy_prior(scale = 2.5),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
})

test_that("bgm ggm works with precision_scale_prior", {
  set.seed(42)
  Y = as.data.frame(matrix(rnorm(200), nrow = 50, ncol = 4))
  fit = bgm(Y,
    variable_type = "continuous",
    precision_scale_prior = exponential_prior(rate = 2),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  spec = fit$.bgm_spec
  expect_equal(spec$prior$scale_prior_type, "exponential")
  expect_equal(spec$prior$scale_shape, 1)
  expect_equal(spec$prior$scale_rate, 2)
})

test_that("bgm ggm resolves eta to the raw rate eta / slab scale", {
  set.seed(42)
  Y = as.data.frame(matrix(rnorm(200), nrow = 50, ncol = 4))
  common = list(
    variable_type = "continuous",
    interaction_prior = cauchy_prior(scale = 4),
    update_method = "adaptive-metropolis",
    iter = 25, warmup = 50, chains = 1, seed = 3,
    display_progress = "none", verbose = FALSE
  )
  fit_eta = do.call(bgm, c(
    list(Y, precision_scale_prior = gamma_prior(shape = 1, eta = 1)), common
  ))
  spec = fit_eta$.bgm_spec
  expect_equal(spec$prior$scale_rate, 0.25)
  expect_equal(spec$prior$scale_eta, 1)

  # The eta spec and its pre-derived raw-rate spec give identical chains
  fit_rate = do.call(bgm, c(
    list(Y, precision_scale_prior = gamma_prior(shape = 1, rate = 0.25)),
    common
  ))
  raw_eta = S7::prop(fit_eta, "raw_samples")
  raw_rate = S7::prop(fit_rate, "raw_samples")
  expect_equal(raw_eta$main, raw_rate$main)
  expect_equal(raw_eta$pairwise, raw_rate$pairwise)
})


# ==============================================================================
# 4. bgm() with Prior Objects <U+2014> Mixed MRF
# ==============================================================================

test_that("bgm mixed works with default priors", {
  set.seed(42)
  data("Wenchuan", package = "bgms")
  dat = data.frame(Wenchuan[1:100, 1:4], V5 = rnorm(100), V6 = rnorm(100))
  fit = bgm(dat,
    variable_type = c(rep("ordinal", 4), rep("continuous", 2)),
    iter = 25, warmup = 100, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  expect_true(all(is.finite(do.call(rbind, fit$raw_samples$pairwise))))
})

test_that("bgm mixed works with normal_prior", {
  set.seed(42)
  data("Wenchuan", package = "bgms")
  dat = data.frame(Wenchuan[1:100, 1:4], V5 = rnorm(100), V6 = rnorm(100))
  fit = bgm(dat,
    variable_type = c(rep("ordinal", 4), rep("continuous", 2)),
    interaction_prior = normal_prior(scale = 0.5),
    iter = 25, warmup = 100, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  expect_true(all(is.finite(do.call(rbind, fit$raw_samples$pairwise))))
})

test_that("bgm mixed works with normal_prior", {
  set.seed(42)
  data("Wenchuan", package = "bgms")
  dat = data.frame(Wenchuan[1:100, 1:4], V5 = rnorm(100), V6 = rnorm(100))
  fit = bgm(dat,
    variable_type = c(rep("ordinal", 4), rep("continuous", 2)),
    threshold_prior = normal_prior(scale = 1),
    iter = 25, warmup = 100, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  expect_true(all(is.finite(do.call(rbind, fit$raw_samples$main))))
})

test_that("bgm mixed works with both normal priors", {
  set.seed(42)
  data("Wenchuan", package = "bgms")
  dat = data.frame(Wenchuan[1:100, 1:4], V5 = rnorm(100), V6 = rnorm(100))
  fit = bgm(dat,
    variable_type = c(rep("ordinal", 4), rep("continuous", 2)),
    interaction_prior = normal_prior(scale = 0.5),
    threshold_prior = normal_prior(scale = 0.5),
    iter = 25, warmup = 100, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
})

test_that("bgm mixed works with means_prior", {
  set.seed(42)
  data("Wenchuan", package = "bgms")
  dat = data.frame(Wenchuan[1:100, 1:4], V5 = rnorm(100), V6 = rnorm(100))
  fit = bgm(dat,
    variable_type = c(rep("ordinal", 4), rep("continuous", 2)),
    means_prior = cauchy_prior(scale = 2),
    iter = 25, warmup = 100, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  # Verify means_prior is stored in spec
  expect_equal(fit$.bgm_spec$prior$means_prior_type, "cauchy")
  expect_equal(fit$.bgm_spec$prior$means_scale, 2)
})


# ==============================================================================
# 5. bgm() with Prior Objects <U+2014> Blume-Capel
# ==============================================================================

test_that("bgm blume-capel works with normal priors", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    variable_type = "blume-capel",
    baseline_category = 2,
    interaction_prior = normal_prior(scale = 0.5),
    threshold_prior = normal_prior(scale = 0.5),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  main_names = fit$raw_samples$parameter_names$main
  expect_true(any(grepl("linear", main_names)))
  expect_true(any(grepl("quadratic", main_names)))
})


# ==============================================================================
# 6. bgmCompare() with Prior Objects
# ==============================================================================

test_that("bgmCompare works with normal priors", {
  data("Wenchuan", package = "bgms")
  x = Wenchuan[1:25, 1:4]
  y = Wenchuan[26:50, 1:4]
  fit = bgmCompare(
    x = x, y = y,
    interaction_prior = normal_prior(scale = 0.5),
    threshold_prior = normal_prior(scale = 0.5),
    difference_selection = FALSE,
    iter = 25, warmup = 100, chains = 1,
    update_method = "adaptive-metropolis",
    display_progress = "none"
  )
  expect_s3_class(fit, "bgmCompare")
})

test_that("bgmCompare works with normal priors + difference selection", {
  data("Wenchuan", package = "bgms")
  x = Wenchuan[1:25, 1:4]
  y = Wenchuan[26:50, 1:4]
  fit = bgmCompare(
    x = x, y = y,
    interaction_prior = normal_prior(scale = 0.5),
    threshold_prior = normal_prior(scale = 0.5),
    difference_selection = TRUE,
    iter = 25, warmup = 100, chains = 1,
    update_method = "adaptive-metropolis",
    display_progress = "none"
  )
  expect_s3_class(fit, "bgmCompare")
})


# ==============================================================================
# 7. Backward Compatibility <U+2014> Deprecated Scalar Parameters
# ==============================================================================

test_that("deprecated pairwise_scale still works", {
  data("Wenchuan", package = "bgms")
  expect_warning(
    fit <- bgm(Wenchuan[1:50, 1:4],
      pairwise_scale = 2.5,
      iter = 25, warmup = 50, chains = 1,
      display_progress = "none"
    ),
    "pairwise_scale"
  )
  expect_s3_class(fit, "bgms")
})

test_that("deprecated main_alpha/main_beta still works", {
  data("Wenchuan", package = "bgms")
  expect_warning(
    fit <- bgm(Wenchuan[1:50, 1:4],
      main_alpha = 1, main_beta = 1,
      iter = 25, warmup = 50, chains = 1,
      display_progress = "none"
    ),
    "main_alpha"
  )
  expect_s3_class(fit, "bgms")
})

test_that("deprecated standardize = FALSE warns and proceeds", {
  data("Wenchuan", package = "bgms")
  expect_warning(
    fit <- bgm(Wenchuan[1:50, 1:4],
      standardize = FALSE,
      iter = 25, warmup = 50, chains = 1,
      display_progress = "none"
    ),
    "standardize"
  )
  expect_s3_class(fit, "bgms")
})

test_that("deprecated standardize = TRUE errors with the manual alternative", {
  data("Wenchuan", package = "bgms")
  # The removed per-pair adjustment has no replacement, so this must fail
  # through lifecycle rather than as an unused-argument error.
  expect_error(
    bgm(Wenchuan[1:50, 1:4],
      standardize = TRUE,
      iter = 25, warmup = 50, chains = 1,
      display_progress = "none"
    ),
    "interaction_prior"
  )
  expect_error(
    bgm(Wenchuan[1:50, 1:4],
      standardize = TRUE,
      iter = 25, warmup = 50, chains = 1,
      display_progress = "none"
    ),
    class = "defunctError"
  )
})

test_that("bgmCompare handles the deprecated standardize on both paths", {
  data("Wenchuan", package = "bgms")
  x = Wenchuan[1:40, 1:4]
  y = Wenchuan[41:80, 1:4]

  expect_warning(
    fit <- bgmCompare(x, y,
      standardize = FALSE,
      iter = 25, warmup = 50, chains = 1,
      display_progress = "none"
    ),
    "standardize"
  )
  expect_s3_class(fit, "bgmCompare")

  expect_error(
    bgmCompare(x, y,
      standardize = TRUE,
      iter = 25, warmup = 50, chains = 1,
      display_progress = "none"
    ),
    "difference_scale"
  )
})


# ==============================================================================
# 8. Backward Compatibility <U+2014> Deprecated String Edge Priors
# ==============================================================================

test_that("deprecated string edge_prior 'Bernoulli' warns", {
  data("Wenchuan", package = "bgms")
  expect_warning(
    fit <- bgm(Wenchuan[1:50, 1:4],
      edge_prior = "Bernoulli",
      iter = 25, warmup = 50, chains = 1,
      display_progress = "none"
    ),
    "edge_prior"
  )
  expect_s3_class(fit, "bgms")
})

test_that("deprecated string edge_prior 'Beta-Bernoulli' warns", {
  data("Wenchuan", package = "bgms")
  expect_warning(
    fit <- bgm(Wenchuan[1:50, 1:4],
      edge_prior = "Beta-Bernoulli",
      beta_bernoulli_alpha = 1, beta_bernoulli_beta = 1,
      iter = 25, warmup = 50, chains = 1,
      display_progress = "none"
    ),
    "edge_prior"
  )
  expect_s3_class(fit, "bgms")
})

test_that("deprecated string edge_prior 'Stochastic-Block' warns", {
  data("Wenchuan", package = "bgms")
  expect_warning(
    fit <- bgm(Wenchuan[1:50, 1:4],
      edge_prior = "Stochastic-Block",
      iter = 25, warmup = 50, chains = 1,
      display_progress = "none"
    ),
    "edge_prior"
  )
  expect_s3_class(fit, "bgms")
})


# ==============================================================================
# 9. Edge Prior Objects <U+2014> Correct Plumbing
# ==============================================================================

test_that("bernoulli_prior object is correctly plumbed", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    edge_prior = bernoulli_prior(0.3),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  spec = fit$.bgm_spec
  expect_equal(spec$prior$edge_prior, "Bernoulli")
  expect_true(all(spec$prior$inclusion_probability[
    upper.tri(spec$prior$inclusion_probability)
  ] == 0.3))
})

test_that("beta_bernoulli_prior object is correctly plumbed", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    edge_prior = beta_bernoulli_prior(alpha = 2, beta = 5),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  spec = fit$.bgm_spec
  expect_equal(spec$prior$edge_prior, "Beta-Bernoulli")
  expect_equal(spec$prior$beta_bernoulli_alpha, 2)
  expect_equal(spec$prior$beta_bernoulli_beta, 5)
})

test_that("sbm_prior object is correctly plumbed", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    edge_prior = sbm_prior(
      alpha = 2, beta = 3,
      alpha_between = 1, beta_between = 4,
      dirichlet_alpha = 0.5, lambda = 2
    ),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
  spec = fit$.bgm_spec
  expect_equal(spec$prior$edge_prior, "Stochastic-Block")
  expect_equal(spec$prior$beta_bernoulli_alpha, 2)
  expect_equal(spec$prior$beta_bernoulli_beta, 3)
  expect_equal(spec$prior$beta_bernoulli_alpha_between, 1)
  expect_equal(spec$prior$beta_bernoulli_beta_between, 4)
  expect_equal(spec$prior$dirichlet_alpha, 0.5)
  expect_equal(spec$prior$lambda, 2)
})


# ==============================================================================
# 10. Prior Object Stored in Fit
# ==============================================================================

test_that("prior info is stored in .bgm_spec", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[1:50, 1:4],
    interaction_prior = normal_prior(scale = 0.5),
    threshold_prior = normal_prior(scale = 0.5),
    iter = 25, warmup = 50, chains = 1,
    display_progress = "none"
  )
  spec = fit$.bgm_spec
  expect_equal(spec$prior$interaction_prior_type, "normal")
  expect_equal(spec$prior$pairwise_scale, 0.5)
  expect_equal(spec$prior$threshold_prior_type, "normal")
  expect_equal(spec$prior$threshold_scale, 0.5)
})

test_that("the exponential default matches gamma_prior(shape = 1, eta = 1) exactly", {
  set.seed(42)
  Y = as.data.frame(matrix(rnorm(200), nrow = 50, ncol = 4))

  fit_default = bgm(Y,
    variable_type = "continuous",
    iter = 25, warmup = 50, chains = 1, seed = 7,
    display_progress = "none"
  )
  fit_gamma = bgm(Y,
    variable_type = "continuous",
    precision_scale_prior = gamma_prior(shape = 1, eta = 1),
    iter = 25, warmup = 50, chains = 1, seed = 7,
    display_progress = "none"
  )

  p_default = fit_default$.bgm_spec$prior
  p_gamma = fit_gamma$.bgm_spec$prior
  expect_equal(p_default$scale_prior_type, "exponential")
  expect_equal(p_gamma$scale_prior_type, "gamma")
  expect_identical(p_default$scale_shape, p_gamma$scale_shape)
  expect_identical(p_default$scale_rate, p_gamma$scale_rate)
  expect_identical(p_default$scale_eta, p_gamma$scale_eta)

  expect_identical(fit_default$raw_samples, fit_gamma$raw_samples)
})
