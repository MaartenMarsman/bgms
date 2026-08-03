# --------------------------------------------------------------------------- #
# Mixed-MRF normalizing-constant correction.
#
# The determinant tilt acts on the continuous precision block, so the
# correction table is built for the continuous variables and the block-
# structure corrections read continuous-continuous edges only. The prior-only
# identities certify the fit-path wiring end to end: driving sample_mixed_mrf
# with n = 0 targets the joint prior, under which the corrected
# hyperparameter updates must return their hyperpriors exactly (theta ~
# Beta(a, b); number of blocks ~ MFM partition prior), because 1/C cancels
# when the graph is summed out. The uncorrected controls show the bias the
# correction removes.
# --------------------------------------------------------------------------- #

mixed_correction_cache = function() {
  options(
    bgms.correction_cache_dir = file.path(tempdir(), "bgms-ctable-mixed")
  )
}

run_mixed_prior_chain = function(p_disc, q_cont, n_samples, n_warmup, seed,
                                 edge_prior_obj, apply_correction, delta) {
  q_tot = p_disc + q_cont
  pairwise_scale = 2.5
  ep = unpack_indicator_prior(edge_prior_obj, num_variables = q_tot)
  correction = NULL
  if(apply_correction && !identical(ep$edge_prior, "Bernoulli")) {
    table = ggm_correction_table(
      p = q_cont, delta = delta,
      interaction_prior = cauchy_prior(scale = pairwise_scale),
      precision_scale_prior = gamma_prior(shape = 1, eta = 1),
      update_method = "gibbs"
    )
    correction = correction_list_from_table(
      table, ep$edge_prior,
      is_continuous = c(rep(0L, p_disc), rep(1L, q_cont))
    )
  }
  input = list(
    discrete_observations = matrix(0L, 0, p_disc),
    continuous_observations = matrix(0, 0, q_cont),
    num_categories = rep(1L, p_disc),
    is_ordinal_variable = rep(1L, p_disc),
    baseline_category = rep(0L, p_disc),
    pairwise_scale = pairwise_scale,
    interaction_prior_type = "cauchy",
    scale_prior_type = "gamma",
    scale_shape = 1.0,
    scale_rate = 1 / pairwise_scale
  )
  res = sample_mixed_mrf(
    inputFromR = input,
    prior_inclusion_prob = ep$inclusion_probability,
    initial_edge_indicators = matrix(1L, q_tot, q_tot),
    no_iter = as.integer(n_samples),
    no_warmup = as.integer(n_warmup),
    no_chains = 1L,
    edge_selection = TRUE,
    seed = as.integer(seed),
    no_threads = 1L,
    progress_type = 0L,
    edge_prior = ep$edge_prior,
    beta_bernoulli_alpha = ep$beta_bernoulli_alpha,
    beta_bernoulli_beta = ep$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = bb_between_or_sentinel(ep$beta_bernoulli_alpha_between),
    beta_bernoulli_beta_between = bb_between_or_sentinel(ep$beta_bernoulli_beta_between),
    dirichlet_alpha = ep$dirichlet_alpha,
    lambda = ep$lambda,
    sampler_type = "adaptive-metropolis",
    delta = delta,
    edge_prior_correction = correction
  )
  expect_length(res, 1L)
  expect_false(isTRUE(res[[1L]]$error))
  res[[1L]]
}

r_mfm_prior_num_blocks = function(n, q, lambda, dirichlet_alpha) {
  K = rpois(n, lambda) + 1L
  vapply(K, function(k) {
    pi = rgamma(k, dirichlet_alpha)
    z = sample.int(k, q, replace = TRUE, prob = pi)
    length(unique(z))
  }, integer(1))
}

nc_tv = function(nc_a, nc_b, q) {
  pa = tabulate(nc_a, nbins = q) / length(nc_a)
  pb = tabulate(nc_b, nbins = q) / length(nc_b)
  0.5 * sum(abs(pa - pb))
}


# --------------------------------------------------------------------------- #
# Resolver: table keyed on the continuous block, mask in model order
# --------------------------------------------------------------------------- #

test_that("the resolver keys the mixed correction on the continuous block", {
  skip_on_cran()
  old = mixed_correction_cache()
  on.exit(options(old), add = TRUE)

  prior = list(
    edge_selection = TRUE,
    edge_prior = "Stochastic-Block",
    interaction_prior_type = "cauchy",
    pairwise_scale = 2.5,
    scale_prior_type = "gamma",
    scale_shape = 1,
    scale_rate = 1 / 2.5,
    delta = 0.5 * log(3)
  )
  sampler = list(verbose = FALSE, cores = 1L)

  correction = ggm_edge_prior_correction(
    prior, sampler,
    num_variables = 5L, num_continuous = 3L
  )
  expect_identical(correction$is_continuous, c(0L, 0L, 1L, 1L, 1L))
  expect_true(all(
    c("theta", "logC", "fprime_density", "fprime", "quad_theta", "quad_f")
    %in% names(correction)
  ))

  # Same cell, all-continuous call: identical curves, no mask.
  ggm_correction = ggm_edge_prior_correction(
    prior, sampler,
    num_variables = 3L
  )
  expect_null(ggm_correction$is_continuous)
  expect_identical(ggm_correction$logC, correction$logC)

  # The beta-bernoulli list needs no mask: its draw reads logC alone.
  prior$edge_prior = "Beta-Bernoulli"
  bb_correction = ggm_edge_prior_correction(
    prior, sampler,
    num_variables = 5L, num_continuous = 3L
  )
  expect_null(bb_correction$is_continuous)
  expect_identical(bb_correction$logC, correction$logC)

  # Fewer than two continuous variables: nothing is tilted.
  expect_null(ggm_edge_prior_correction(
    prior, sampler,
    num_variables = 3L, num_continuous = 1L
  ))
})


# --------------------------------------------------------------------------- #
# Prior-only identities
# --------------------------------------------------------------------------- #

test_that("corrected mixed prior chain returns the Beta hyperprior on theta", {
  skip_on_cran()
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to run the prior-chain identity certifications"
  )
  old = mixed_correction_cache()
  on.exit(options(old), add = TRUE)

  ch = run_mixed_prior_chain(
    p_disc = 2L, q_cont = 3L, n_samples = 8000, n_warmup = 1000, seed = 42,
    edge_prior_obj = beta_bernoulli_prior(1, 1),
    apply_correction = TRUE, delta = 0.5 * log(3)
  )
  th = as.numeric(ch$inclusion_parameter_samples)
  ks = suppressWarnings(ks.test(th, "punif"))

  expect_lt(abs(mean(th) - 0.5), 0.08)
  expect_lt(unname(ks$statistic), 0.12)
})

test_that("uncorrected mixed prior chain biases theta toward sparsity", {
  skip_on_cran()
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to run the prior-chain identity certifications"
  )

  ch = run_mixed_prior_chain(
    p_disc = 2L, q_cont = 3L, n_samples = 8000, n_warmup = 1000, seed = 42,
    edge_prior_obj = beta_bernoulli_prior(1, 1),
    apply_correction = FALSE, delta = 0.5 * log(3)
  )
  th = as.numeric(ch$inclusion_parameter_samples)
  ks = suppressWarnings(ks.test(th, "punif"))

  expect_lt(mean(th), 0.42)
  expect_gt(unname(ks$statistic), 0.15)
})

test_that("corrected mixed prior chain returns the MFM partition prior", {
  skip_on_cran()
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to run the prior-chain identity certifications"
  )
  old = mixed_correction_cache()
  on.exit(options(old), add = TRUE)

  q_tot = 7L
  ch = run_mixed_prior_chain(
    p_disc = 2L, q_cont = 5L, n_samples = 12000, n_warmup = 1000, seed = 77,
    edge_prior_obj = sbm_prior(),
    apply_correction = TRUE, delta = 2
  )
  nc = apply(ch$allocation_samples, 2, function(z) length(unique(z)))

  set.seed(99)
  nc_prior = r_mfm_prior_num_blocks(2e5, q_tot, lambda = 1, dirichlet_alpha = 1)

  expect_lt(nc_tv(nc, nc_prior, q_tot), 0.04)
})

test_that("uncorrected mixed prior chain under-segments the partition", {
  skip_on_cran()
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to run the prior-chain identity certifications"
  )

  q_tot = 7L
  ch = run_mixed_prior_chain(
    p_disc = 2L, q_cont = 5L, n_samples = 12000, n_warmup = 1000, seed = 77,
    edge_prior_obj = sbm_prior(),
    apply_correction = FALSE, delta = 2
  )
  nc = apply(ch$allocation_samples, 2, function(z) length(unique(z)))

  set.seed(99)
  nc_prior = r_mfm_prior_num_blocks(2e5, q_tot, lambda = 1, dirichlet_alpha = 1)

  expect_gt(nc_tv(nc, nc_prior, q_tot), 0.07)
  expect_lt(mean(nc), mean(nc_prior))
})


# --------------------------------------------------------------------------- #
# bgm() end to end
#
# The normalizing-constant correction is a joint-specification artifact: the
# hierarchical path tracks Z(Gamma) in the edge moves and builds no table at
# all. Since F-010 the default is hierarchical, so every fit below names the
# joint specification -- these tests are about the correction, not about what a
# default fit does.
# --------------------------------------------------------------------------- #

mixed_smoke_data = function() {
  set.seed(7)
  n = 60
  x = cbind(
    matrix(sample(0:2, n * 2, replace = TRUE), n, 2),
    matrix(rnorm(n * 3), n, 3)
  )
  list(
    x = x,
    variable_type = c("ordinal", "ordinal", rep("continuous", 3))
  )
}

test_that("bgm() applies the correction to a mixed beta-bernoulli fit", {
  skip_on_cran()
  old = mixed_correction_cache()
  on.exit(options(old), add = TRUE)
  # A fresh cache forces an actual build, the only case that announces
  # itself; a cache hit is silent.
  unlink(getOption("bgms.correction_cache_dir"), recursive = TRUE)

  d = mixed_smoke_data()
  # Capture stdout too, so an interactive test run does not leak the
  # build progress bar; the announcement is a message (stderr).
  msgs = NULL
  capture.output(
    msgs <- capture.output(
      fit <- bgm(d$x,
        variable_type = d$variable_type,
        edge_prior = beta_bernoulli_prior(),
        precision_graph_prior = "joint",
        iter = 300, warmup = 300, chains = 1, cores = 1,
        display_progress = "none", verbose = TRUE
      ),
      type = "message"
    ),
    type = "output"
  )

  expect_true(any(grepl("correction table", msgs)))
  expect_s3_class(fit, "bgms")
  th = fit$inclusion_parameter_samples
  expect_length(th, 1L)
  expect_length(th[[1L]], 300L)
  expect_true(all(th[[1L]] > 0 & th[[1L]] < 1))
})

test_that("bgm() applies the correction to a mixed sbm fit", {
  skip_on_cran()
  old = mixed_correction_cache()
  on.exit(options(old), add = TRUE)

  d = mixed_smoke_data()
  fit = suppressMessages(bgm(d$x,
    variable_type = d$variable_type,
    edge_prior = sbm_prior(),
    precision_graph_prior = "joint",
    iter = 300, warmup = 300, chains = 1, cores = 1,
    display_progress = "none"
  ))

  expect_s3_class(fit, "bgms")
  expect_false(is.null(fit$posterior_num_blocks))
  expect_false(is.null(fit$posterior_mean_indicator))
})

test_that("bgm() skips the correction with one continuous variable", {
  skip_on_cran()

  d = mixed_smoke_data()
  msgs = capture.output(
    fit <- bgm(d$x[, 1:3],
      variable_type = d$variable_type[1:3],
      edge_prior = beta_bernoulli_prior(),
      precision_graph_prior = "joint",
      iter = 200, warmup = 300, chains = 1, cores = 1,
      display_progress = "none", verbose = TRUE
    ),
    type = "message"
  )

  expect_false(any(grepl("correction table", msgs)))
  expect_s3_class(fit, "bgms")
})

test_that("two continuous variables: beta-bernoulli corrects, sbm falls back", {
  skip_on_cran()
  old = mixed_correction_cache()
  on.exit(options(old), add = TRUE)

  # One tilted pair: the logC curve exists (the beta-bernoulli draw needs
  # nothing else), but the slope curve is not resolvable, so the block
  # model warns and keeps the plain conjugate updates.
  d = mixed_smoke_data()
  x2 = d$x[, 1:4]
  vt2 = d$variable_type[1:4]

  fit_bb = suppressMessages(bgm(x2,
    variable_type = vt2,
    edge_prior = beta_bernoulli_prior(),
    precision_graph_prior = "joint",
    iter = 200, warmup = 300, chains = 1, cores = 1,
    display_progress = "none"
  ))
  expect_s3_class(fit_bb, "bgms")
  th = fit_bb$inclusion_parameter_samples
  expect_length(th, 1L)
  expect_true(all(th[[1L]] > 0 & th[[1L]] < 1))

  expect_warning(
    fit_sbm <- suppressMessages(bgm(x2,
      variable_type = vt2,
      edge_prior = sbm_prior(),
      precision_graph_prior = "joint",
      iter = 200, warmup = 300, chains = 1, cores = 1,
      display_progress = "none"
    )),
    "slope curve is not resolvable"
  )
  expect_s3_class(fit_sbm, "bgms")
  expect_false(is.null(fit_sbm$posterior_num_blocks))
})
