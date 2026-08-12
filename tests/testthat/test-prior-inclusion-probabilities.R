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

# Tiers. The acceptance-target wiring, the analytic pieces, the extractor
# contracts and the hierarchical pass-through stay local. The cells that
# check one route against another (table vs prior-only chain, the mixed
# per-class split, the hierarchical edge-prior identity) are settled-numerics
# agreement, so they are T1.

skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run the prior-PIP route agreement cells"
  )
}

prior_pip_ggm_data = function() {
  set.seed(7)
  x = matrix(rnorm(80 * 3), 80, 3)
  colnames(x) = c("c1", "c2", "c3")
  x
}

# The hierarchical specification builds a Z-ratio surface at fit time, whose
# bipartite anchor grid starts at four nodes, so its fits need a wider block.
prior_pip_hier_data = function() {
  set.seed(7)
  x = matrix(rnorm(80 * 5), 80, 5)
  colnames(x) = paste0("c", 1:5)
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
# The prior-only chain tunes the way the method it names does
# --------------------------------------------------------------------------- #
#
# The chain describes the prior the fit targeted, so it has to carry the same
# acceptance target the fit resolved for its update method. The C++ default is
# 0.80 (the NUTS target), so omitting the argument silently retuned an
# adaptive-metropolis chain. resolve_target_acceptance() is the single source
# shared with validate_sampler() and sample_ggm_prior().

test_that("the GGM prior-only chain carries its method's acceptance target", {
  skip_on_cran()

  fit = small_fit(prior_pip_ggm_data(), "continuous", bernoulli_prior(0.5))

  captured = NULL
  local_mocked_bindings(
    sample_ggm = function(...) {
      captured <<- list(...)
      stop("captured-call")
    },
    .package = "bgms"
  )
  expect_error(
    prior_only_chain_pips(get_fit_spec(fit), iter = 10L, warmup = 10L),
    "captured-call"
  )
  expect_identical(
    captured$target_acceptance,
    resolve_target_acceptance(captured$sampler_type)
  )
})

test_that("the mixed prior-only chain carries its method's acceptance target", {
  skip_on_cran()

  d = prior_pip_mixed_data()
  fit = small_fit(d$x, d$variable_type, bernoulli_prior(0.5))

  captured = NULL
  local_mocked_bindings(
    sample_mixed_mrf = function(...) {
      captured <<- list(...)
      stop("captured-call")
    },
    .package = "bgms"
  )
  expect_error(
    prior_only_chain_pips(get_fit_spec(fit), iter = 10L, warmup = 10L),
    "captured-call"
  )
  expect_identical(captured$sampler_type, "adaptive-metropolis")
  expect_identical(
    captured$target_acceptance,
    resolve_target_acceptance("adaptive-metropolis")
  )
})


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
  skip_unless_slow()
  skip_on_cran()

  fit = small_fit(prior_pip_ggm_data(), "continuous",
    beta_bernoulli_prior(),
    update_method = "gibbs",
    # The tilted table and the prior-only chain are the joint route; the
    # default is hierarchical since F-010, and it returns the edge prior
    # itself (asserted in its own section below).
    precision_graph_prior = "joint"
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
    update_method = "gibbs", delta = 0,
    precision_graph_prior = "joint"
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
    update_method = "gibbs",
    precision_graph_prior = "joint"
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
  skip_unless_slow()
  skip_on_cran()

  d = prior_pip_mixed_data()
  fit = small_fit(d$x, d$variable_type, beta_bernoulli_prior(),
    precision_graph_prior = "joint"
  )
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

  # A single prior-only chain is a noisy estimator here: the Beta-Bernoulli
  # inclusion parameter mixes slowly, and across seeds the class proportions
  # have SD ~= 0.03 even at 40000 iterations. This gate is a fixed-seed spot
  # check, not a tolerance derived from an MCSE. The chain was previously run
  # for 8000 iterations, which left it ~0.036 short of the table value once it
  # was retuned to its method's acceptance target (0.44, not the C++ default
  # 0.80); from 20000 on the deviations settle below 0.02 at this seed.
  chain = prior_only_chain_pips(get_fit_spec(fit), iter = 40000, warmup = 4000)
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


# --------------------------------------------------------------------------- #
# The hierarchical specification
#
# p(K | Gamma) is normalized per graph there, so integrating K out returns
# pi(Gamma) exactly: the prior inclusion probability of a
# continuous-continuous edge is the edge prior's own marginal, with no tilt to
# read off a table and no prior-only chain to run.
# --------------------------------------------------------------------------- #

test_that("a hierarchical fit's prior inclusion probability is the edge prior", {
  skip_unless_slow()
  skip_on_cran()

  x = prior_pip_hier_data()
  hier = small_fit(x, "continuous", bernoulli_prior(0.4),
    precision_graph_prior = "hierarchical"
  )
  joint = small_fit(x, "continuous", bernoulli_prior(0.4),
    precision_graph_prior = "joint"
  )

  hier_pip = extract_prior_inclusion_probabilities(hier)
  expect_equal(hier_pip[upper.tri(hier_pip)], rep(0.4, 10L))

  # The joint specification's tilt is what the two differ by, and it is real.
  joint_pip = extract_prior_inclusion_probabilities(joint)
  expect_false(isTRUE(all.equal(joint_pip[1, 2], 0.4)))

  bb = small_fit(x, "continuous", beta_bernoulli_prior(alpha = 2, beta = 3),
    precision_graph_prior = "hierarchical"
  )
  bb_pip = extract_prior_inclusion_probabilities(bb)
  expect_equal(bb_pip[upper.tri(bb_pip)], rep(2 / 5, 10L))
})

test_that("a hierarchical fit needs neither a correction table nor a chain", {
  skip_on_cran()

  # An empty cache directory with table building disabled: a route that reaches
  # ggm_correction_table() would have to run the sweep, which takes minutes.
  # The closed form returns at once.
  withr::local_options(
    bgms.correction_cache_dir = withr::local_tempdir()
  )
  x = prior_pip_hier_data()
  fit = small_fit(x, "continuous", bernoulli_prior(0.5),
    precision_graph_prior = "hierarchical"
  )

  elapsed = system.time(
    pip <- extract_prior_inclusion_probabilities(fit, recompute = TRUE)
  )[["elapsed"]]
  expect_lt(elapsed, 5)
  expect_equal(pip[1, 2], 0.5)
})

test_that("hierarchical per-edge Bernoulli probabilities pass through", {
  skip_on_cran()

  x = prior_pip_hier_data()
  probability = matrix(0.3, 5L, 5L)
  probability[1, 2] = probability[2, 1] = 0.7
  fit = small_fit(x, "continuous", bernoulli_prior(probability),
    precision_graph_prior = "hierarchical"
  )

  pip = extract_prior_inclusion_probabilities(fit)
  expect_equal(pip[1, 2], 0.7)
  expect_equal(pip[1, 3], 0.3)
  expect_equal(unname(diag(pip)), rep(0, 5L))
})

test_that("the class values are cached on the fit and survive a round trip", {
  skip_on_cran()

  x = prior_pip_hier_data()
  fit = small_fit(x, "continuous", bernoulli_prior(0.4),
    precision_graph_prior = "hierarchical"
  )
  cache = get_fit_cache(fit)
  expect_null(cache$prior_inclusion_class_values)

  first = extract_prior_inclusion_probabilities(fit)
  expect_false(is.null(cache$prior_inclusion_class_values))

  path = withr::local_tempfile(fileext = ".rds")
  saveRDS(fit, path)
  restored = readRDS(path)
  expect_equal(
    get_fit_cache(restored)$prior_inclusion_class_values,
    cache$prior_inclusion_class_values
  )
  expect_equal(extract_prior_inclusion_probabilities(restored), first)
})
