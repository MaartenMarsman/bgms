# --------------------------------------------------------------------------- #
# Prior inclusion probabilities.
#
# Under the joint spike-and-slab prior on a continuous block the graph
# marginal is reweighted by the per-graph normalizer (positive-definite-cone
# mass shaped by the determinant tilt), so continuous-continuous edges carry
# the joint-block prior edge density instead of the edge-prior marginal — at
# any delta, including zero. The extractor reads that density off the
# correction table (Bernoulli, Beta-Bernoulli) or estimates it with a
# prior-only chain (Stochastic-Block); the consistency gates check the two
# routes against each other.
# --------------------------------------------------------------------------- #

prior_pip_ggm_data = function() {
  set.seed(7)
  x = matrix(rnorm(80 * 3), 80, 3)
  colnames(x) = c("c1", "c2", "c3")
  x
}

# Interleaved variable order so the class placement in the output matrix is
# exercised through the reordering, not just the block layout.
prior_pip_mixed_data = function() {
  set.seed(9)
  n = 80
  x = cbind(
    sample(0:2, n, replace = TRUE),
    rnorm(n),
    sample(0:2, n, replace = TRUE),
    rnorm(n),
    rnorm(n)
  )
  colnames(x) = c("d1", "c1", "d2", "c2", "c3")
  list(
    x = x,
    variable_type = c(
      "ordinal", "continuous", "ordinal", "continuous", "continuous"
    )
  )
}

small_fit = function(x, variable_type, edge_prior, ...) {
  suppressMessages(bgm(x,
    variable_type = variable_type,
    interaction_prior = cauchy_prior(scale = 2.5),
    edge_prior = edge_prior,
    iter = 100, warmup = 200, chains = 1, cores = 1,
    display_progress = "none", ...
  ))
}


# --------------------------------------------------------------------------- #
# Analytic pieces
# --------------------------------------------------------------------------- #

test_that("the factorized SBM marginal matches partition simulation", {
  prior = list(
    dirichlet_alpha = 1, lambda = 1,
    beta_bernoulli_alpha = 1, beta_bernoulli_beta = 1,
    beta_bernoulli_alpha_between = 2, beta_bernoulli_beta_between = 3
  )
  analytic = sbm_factorized_edge_marginal(prior)

  set.seed(1)
  K = rpois(1e5, prior$lambda) + 1L
  same = vapply(K, function(k) {
    w = rgamma(k, prior$dirichlet_alpha)
    z = sample.int(k, 2, replace = TRUE, prob = w)
    z[1L] == z[2L]
  }, logical(1))
  simulated = mean(same) * 0.5 + (1 - mean(same)) * 0.4

  expect_equal(analytic, simulated, tolerance = 0.015)

  # Shared between-block hyperparameters collapse the mixture.
  prior$beta_bernoulli_alpha_between = NULL
  expect_equal(sbm_factorized_edge_marginal(prior), 0.5, tolerance = 1e-12)
})

test_that("the Beta quadrature reduces to the Beta mean without reweighting", {
  theta = ggm_correction_theta_grid(120L)
  table = list(theta = theta, edens = theta)
  expect_equal(
    tilted_bb_edge_marginal(table, 1, 1), 0.5,
    tolerance = 5e-3
  )
  expect_equal(
    tilted_bb_edge_marginal(table, 2, 6), 0.25,
    tolerance = 5e-3
  )
})


# --------------------------------------------------------------------------- #
# GGM: table route vs prior-only chain
# --------------------------------------------------------------------------- #

test_that("GGM beta-bernoulli prior PIPs match the prior-only chain", {
  skip_on_cran()

  fit = small_fit(prior_pip_ggm_data(), "continuous",
    beta_bernoulli_prior(),
    update_method = "gibbs"
  )
  pip = extract_prior_inclusion_probabilities(fit)

  expect_identical(dim(pip), c(3L, 3L))
  expect_identical(rownames(pip), c("c1", "c2", "c3"))
  expect_true(isTRUE(all.equal(pip, t(pip))))
  expect_identical(diag(pip), c(c1 = 0, c2 = 0, c3 = 0))
  offdiag = pip[upper.tri(pip)]
  expect_length(unique(offdiag), 1L)
  expect_true(offdiag[1L] > 0 && offdiag[1L] < 0.5)

  chain = prior_only_chain_pips(get_fit_spec(fit), iter = 6000, warmup = 1000)
  expect_lt(abs(offdiag[1L] - chain), 0.025)
})

test_that("the joint block reweights a fixed Bernoulli prior at delta = 0", {
  skip_on_cran()

  fit = small_fit(prior_pip_ggm_data(), "continuous",
    bernoulli_prior(0.5),
    update_method = "gibbs", delta = 0
  )
  pip = extract_prior_inclusion_probabilities(fit)[1, 2]

  # The positive-definite-cone mass favors sparser patterns, so the prior
  # edge probability sits well below the nominal 0.5 even untilted.
  expect_lt(pip, 0.35)
  chain = prior_only_chain_pips(get_fit_spec(fit), iter = 6000, warmup = 1000)
  expect_lt(abs(pip - chain), 0.02)
})

test_that("GGM SBM prior PIPs come from a cached deterministic chain", {
  skip_on_cran()

  fit = small_fit(prior_pip_ggm_data(), "continuous",
    sbm_prior(),
    update_method = "gibbs"
  )
  pip1 = extract_prior_inclusion_probabilities(fit, iter = 3000, warmup = 500)
  pip2 = extract_prior_inclusion_probabilities(fit)
  pip3 = extract_prior_inclusion_probabilities(fit,
    iter = 3000, warmup = 500,
    recompute = TRUE
  )

  offdiag = pip1[upper.tri(pip1)]
  expect_length(unique(offdiag), 1L)
  expect_true(offdiag[1L] > 0 && offdiag[1L] < 1)
  expect_identical(pip1, pip2)
  expect_identical(pip1, pip3)
})


# --------------------------------------------------------------------------- #
# Mixed MRF: per-class values in user variable order
# --------------------------------------------------------------------------- #

test_that("mixed beta-bernoulli prior PIPs split by edge class", {
  skip_on_cran()

  d = prior_pip_mixed_data()
  fit = small_fit(d$x, d$variable_type, beta_bernoulli_prior())
  pip = extract_prior_inclusion_probabilities(fit)

  expect_true(isTRUE(all.equal(pip, t(pip))))
  # Discrete-discrete and cross edges keep the Beta mean exactly.
  expect_identical(pip["d1", "d2"], 0.5)
  expect_identical(pip["d1", "c1"], 0.5)
  expect_identical(pip["d2", "c3"], 0.5)
  # Continuous-continuous edges carry the joint-block density.
  cc = c(pip["c1", "c2"], pip["c1", "c3"], pip["c2", "c3"])
  expect_length(unique(cc), 1L)
  expect_true(cc[1L] < 0.5)

  chain = prior_only_chain_pips(get_fit_spec(fit), iter = 8000, warmup = 1000)
  expect_lt(abs(cc[1L] - chain[["cc"]]), 0.03)
  expect_lt(abs(chain[["dd"]] - 0.5), 0.05)
  expect_lt(abs(chain[["cross"]] - 0.5), 0.05)
})


# --------------------------------------------------------------------------- #
# Factorized models: analytic values
# --------------------------------------------------------------------------- #

test_that("ordinal fits use the analytic edge-prior marginals", {
  skip_on_cran()

  set.seed(3)
  x = matrix(sample(0:2, 60 * 4, replace = TRUE), 60, 4)
  colnames(x) = paste0("o", 1:4)

  fit_bb = suppressMessages(bgm(x,
    edge_prior = beta_bernoulli_prior(alpha = 2, beta = 6),
    iter = 100, warmup = 200, chains = 1, cores = 1,
    display_progress = "none"
  ))
  pip_bb = extract_prior_inclusion_probabilities(fit_bb)
  expect_identical(unique(pip_bb[upper.tri(pip_bb)]), 0.25)

  fit_sbm = suppressMessages(bgm(x,
    edge_prior = sbm_prior(),
    iter = 100, warmup = 200, chains = 1, cores = 1,
    display_progress = "none"
  ))
  pip_sbm = extract_prior_inclusion_probabilities(fit_sbm)
  expected = sbm_factorized_edge_marginal(get_fit_spec(fit_sbm)$prior)
  expect_identical(unique(pip_sbm[upper.tri(pip_sbm)]), expected)
})


# --------------------------------------------------------------------------- #
# Errors and posterior-extractor regression
# --------------------------------------------------------------------------- #

test_that("extraction requires edge selection", {
  skip_on_cran()

  fit = suppressMessages(bgm(prior_pip_ggm_data(),
    variable_type = "continuous",
    edge_selection = FALSE, update_method = "gibbs",
    iter = 100, warmup = 200, chains = 1, cores = 1,
    display_progress = "none"
  ))
  expect_error(
    extract_prior_inclusion_probabilities(fit),
    "edge_selection = TRUE"
  )
})

test_that("posterior extractor maps mixed indicators through block order", {
  skip_on_cran()

  d = prior_pip_mixed_data()
  fit = small_fit(d$x, d$variable_type, beta_bernoulli_prior())

  expect_identical(
    extract_posterior_inclusion_probabilities(fit),
    fit$posterior_mean_indicator
  )
})
