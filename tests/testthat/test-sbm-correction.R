# --------------------------------------------------------------------------- #
# Tests for the stochastic-block normalizing-constant corrections.
#
# With a CONSTANT slope curve f'(d) = c0 the ported pieces have closed forms:
# the self-consistent per-edge slopes are all c0; the mini-TI label-move
# correction integrates exactly to sum_j [F(th_new_j) - F(th_old_j)] with
# F(th) = log(1 - th + th * exp(c0)); and with a constant per-pair curve
# f(theta) = f0 the corrected new-cluster marginal equals the clean
# Beta-Bernoulli marginal minus f0 * (number of candidate edges).
# --------------------------------------------------------------------------- #

const_curves = function(c0, f0 = 0) {
  list(
    fprime_density = c(0.01, 0.99),
    fprime = c(c0, c0),
    quad_theta = seq(0.0025, 0.9975, length.out = 200),
    quad_f = rep(f0, 200)
  )
}

test_that("self-consistent slopes reduce to the constant slope", {
  cv = const_curves(-0.4)
  z = c(1L, 1L, 2L, 2L, 1L)
  bp = matrix(c(0.7, 0.2, 0.2, 0.5), 2, 2)

  ce = test_sbm_compute_ce(
    z, bp, cv$fprime_density, cv$fprime,
    cv$quad_theta, cv$quad_f
  )

  off = ce[upper.tri(ce)]
  expect_equal(off, rep(-0.4, length(off)), tolerance = 1e-12)
})

test_that("mini-TI matches the constant-slope closed form", {
  c0 = -0.35
  cv = const_curves(c0)
  z = c(1L, 1L, 2L, 2L, 3L, 3L)
  bp = matrix(c(
    0.75, 0.20, 0.10,
    0.20, 0.60, 0.30,
    0.10, 0.30, 0.55
  ), 3, 3)

  dlogC = test_sbm_miniti_node(
    1L, z, bp, 1L, 3L,
    cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f
  )

  Ff = function(th) log(1 - th + th * exp(c0))
  exact = 0
  for(j in 2:6) {
    exact = exact + Ff(bp[3, z[j]]) - Ff(bp[1, z[j]])
  }
  expect_equal(dlogC, exact, tolerance = 5e-4)
})

test_that("removal morph matches the constant-slope closed form", {
  c0 = -0.35
  cv = const_curves(c0)
  z = c(1L, 1L, 2L, 2L, 3L, 3L)
  bp = matrix(c(
    0.75, 0.20, 0.10,
    0.20, 0.60, 0.30,
    0.10, 0.30, 0.55
  ), 3, 3)

  removal = test_sbm_miniti_removal(
    1L, z, bp, 1L,
    cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f
  )

  # F(0) = 0, so the removal is -sum_j F(theta_old_j)
  Ff = function(th) log(1 - th + th * exp(c0))
  exact = -sum(vapply(2:6, function(j) Ff(bp[1, z[j]]), numeric(1)))
  expect_equal(removal, exact, tolerance = 5e-4)
})

test_that("corrected new-cluster marginal reduces to the clean marginal", {
  set.seed(8)
  q = 6
  z = c(1L, 1L, 2L, 2L, 2L, 1L)
  G = matrix(0L, q, q)
  for(i in 1:(q - 1)) {
    for(j in (i + 1):q) {
      G[i, j] = G[j, i] = rbinom(1, 1, 0.5)
    }
  }
  diag(G) = 1L
  a = 1.5
  b = 2.5
  node = 3L

  clean_exact = 0
  for(r in 1:2) {
    members = setdiff(which(z == r), node)
    nr = length(members)
    mr = sum(G[node, members])
    clean_exact = clean_exact + lbeta(a + mr, b + nr - mr) - lbeta(a, b)
  }

  cv0 = const_curves(0, f0 = 0)
  out0 = test_sbm_corrected_log_marginal(
    node, z, G, a, b,
    cv0$fprime_density, cv0$fprime, cv0$quad_theta, cv0$quad_f
  )
  expect_equal(out0, clean_exact, tolerance = 2e-3)

  f0 = -0.25
  cvf = const_curves(0, f0 = f0)
  outf = test_sbm_corrected_log_marginal(
    node, z, G, a, b,
    cvf$fprime_density, cvf$fprime, cvf$quad_theta, cvf$quad_f
  )
  n_candidates = q - 1
  expect_equal(outf, clean_exact - f0 * n_candidates, tolerance = 2e-3)
})

# --------------------------------------------------------------------------- #
# Mixed-MRF masking: with an is_continuous mask only continuous-continuous
# pairs carry the tilt. A full mask must reproduce the unmasked results
# exactly; under a mixed mask the constant-slope closed forms restrict their
# sums to continuous neighbours, discrete nodes contribute nothing, and the
# corrected marginal reduces to the clean conjugate marginal for a discrete
# node.
# --------------------------------------------------------------------------- #

test_that("a full mask reproduces the unmasked results exactly", {
  c0 = -0.3
  cv = const_curves(c0, f0 = -0.2)
  z = c(1L, 1L, 2L, 2L, 3L, 3L)
  bp = matrix(c(
    0.75, 0.20, 0.10,
    0.20, 0.60, 0.30,
    0.10, 0.30, 0.55
  ), 3, 3)
  full = rep(1L, 6)

  expect_identical(
    test_sbm_compute_ce(
      z, bp, cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
      is_continuous = full
    ),
    test_sbm_compute_ce(
      z, bp, cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f
    )
  )
  expect_identical(
    test_sbm_miniti_node(
      1L, z, bp, 1L, 3L,
      cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
      is_continuous = full
    ),
    test_sbm_miniti_node(
      1L, z, bp, 1L, 3L,
      cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f
    )
  )
  expect_identical(
    test_sbm_miniti_removal(
      1L, z, bp, 1L,
      cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
      is_continuous = full
    ),
    test_sbm_miniti_removal(
      1L, z, bp, 1L,
      cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f
    )
  )

  set.seed(11)
  G = matrix(rbinom(36, 1, 0.5), 6, 6)
  G[lower.tri(G)] = t(G)[lower.tri(G)]
  diag(G) = 1L
  expect_identical(
    test_sbm_corrected_log_marginal(
      3L, z, G, 1.5, 2.5,
      cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
      is_continuous = full
    ),
    test_sbm_corrected_log_marginal(
      3L, z, G, 1.5, 2.5,
      cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f
    )
  )
})

test_that("masked slopes vanish on pairs with a discrete endpoint", {
  cv = const_curves(-0.4)
  z = c(1L, 1L, 2L, 2L, 1L, 2L)
  bp = matrix(c(0.7, 0.2, 0.2, 0.5), 2, 2)
  mask = c(0L, 0L, 1L, 1L, 1L, 1L)

  ce = test_sbm_compute_ce(
    z, bp, cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
    is_continuous = mask
  )

  for(i in 1:5) {
    for(j in (i + 1):6) {
      expected = if(mask[i] == 1L && mask[j] == 1L) -0.4 else 0
      expect_equal(ce[i, j], expected, tolerance = 1e-12)
    }
  }
})

test_that("masked mini-TI sums over continuous neighbours only", {
  c0 = -0.35
  cv = const_curves(c0)
  z = c(1L, 1L, 2L, 2L, 3L, 3L)
  bp = matrix(c(
    0.75, 0.20, 0.10,
    0.20, 0.60, 0.30,
    0.10, 0.30, 0.55
  ), 3, 3)
  mask = c(0L, 0L, 1L, 1L, 1L, 1L)
  Ff = function(th) log(1 - th + th * exp(c0))

  # Continuous node 3 (cluster 2), morph cluster 2 -> 1: continuous
  # neighbours are 4, 5, 6.
  dlogC = test_sbm_miniti_node(
    3L, z, bp, 2L, 1L,
    cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
    is_continuous = mask
  )
  exact = 0
  for(j in 4:6) {
    exact = exact + Ff(bp[1, z[j]]) - Ff(bp[2, z[j]])
  }
  expect_equal(dlogC, exact, tolerance = 5e-4)

  removal = test_sbm_miniti_removal(
    3L, z, bp, 2L,
    cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
    is_continuous = mask
  )
  exact_removal = -sum(vapply(4:6, function(j) Ff(bp[2, z[j]]), numeric(1)))
  expect_equal(removal, exact_removal, tolerance = 5e-4)

  # Discrete node 1: no tilted pairs, so both corrections are exactly zero.
  expect_identical(
    test_sbm_miniti_node(
      1L, z, bp, 1L, 3L,
      cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
      is_continuous = mask
    ),
    0
  )
  expect_identical(
    test_sbm_miniti_removal(
      1L, z, bp, 1L,
      cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
      is_continuous = mask
    ),
    0
  )
})

test_that("masked corrected marginal tilts continuous candidates only", {
  set.seed(12)
  q = 6
  z = c(1L, 1L, 2L, 2L, 2L, 1L)
  G = matrix(0L, q, q)
  for(i in 1:(q - 1)) {
    for(j in (i + 1):q) {
      G[i, j] = G[j, i] = rbinom(1, 1, 0.5)
    }
  }
  diag(G) = 1L
  a = 1.5
  b = 2.5
  mask = c(0L, 0L, 1L, 1L, 1L, 1L)

  clean_for_node = function(node) {
    out = 0
    for(r in 1:2) {
      members = setdiff(which(z == r), node)
      nr = length(members)
      mr = sum(G[node, members])
      out = out + lbeta(a + mr, b + nr - mr) - lbeta(a, b)
    }
    out
  }

  # Continuous node 3: constant f0 exposes one tilt factor per continuous
  # candidate (nodes 4, 5, 6).
  f0 = -0.25
  cvf = const_curves(0, f0 = f0)
  outf = test_sbm_corrected_log_marginal(
    3L, z, G, a, b,
    cvf$fprime_density, cvf$fprime, cvf$quad_theta, cvf$quad_f,
    is_continuous = mask
  )
  expect_equal(outf, clean_for_node(3L) - f0 * 3, tolerance = 2e-3)

  # Discrete node 1: no tilted pairs, so the marginal is the exact clean
  # conjugate value (analytic branch, no quadrature error).
  outd = test_sbm_corrected_log_marginal(
    1L, z, G, a, b,
    cvf$fprime_density, cvf$fprime, cvf$quad_theta, cvf$quad_f,
    is_continuous = mask
  )
  expect_equal(outd, clean_for_node(1L), tolerance = 1e-10)
})

test_that("corrections vanish with fewer than two continuous nodes", {
  cv = const_curves(-0.5, f0 = -0.3)
  z = c(1L, 1L, 2L, 2L, 3L, 3L)
  bp = matrix(c(
    0.75, 0.20, 0.10,
    0.20, 0.60, 0.30,
    0.10, 0.30, 0.55
  ), 3, 3)
  mask = c(0L, 0L, 0L, 0L, 0L, 1L)

  ce = test_sbm_compute_ce(
    z, bp, cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
    is_continuous = mask
  )
  expect_true(all(ce == 0))

  expect_identical(
    test_sbm_miniti_node(
      6L, z, bp, 3L, 1L,
      cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
      is_continuous = mask
    ),
    0
  )
  expect_identical(
    test_sbm_miniti_removal(
      6L, z, bp, 3L,
      cv$fprime_density, cv$fprime, cv$quad_theta, cv$quad_f,
      is_continuous = mask
    ),
    0
  )
})


# --------------------------------------------------------------------------- #
# Prior-only identity: the corrected prior chain's number-of-blocks
# distribution must return the MFM partition prior (K - 1 ~ Poisson(lambda)
# components, symmetric Dirichlet allocation), because 1/C(z, theta) cancels
# when the graph is summed out. The uncorrected chain under-segments. Shares
# the correction-table cache cell with the beta-bernoulli identity tests.
# --------------------------------------------------------------------------- #

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

test_that("corrected prior-only chain returns the MFM partition prior", {
  skip_on_cran()
  old = options(
    bgms.correction_cache_dir = file.path(tempdir(), "bgms-ctable-identity")
  )
  on.exit(options(old), add = TRUE)

  q = 5
  draws = sample_ggm_prior(
    p = q, n_samples = 6000, n_warmup = 500,
    seed = 31, verbose = FALSE, spec = "joint",
    update_method = "gibbs", edge_prior = sbm_prior()
  )
  nc_chain = apply(draws$allocations, 1, function(z) length(unique(z)))

  set.seed(32)
  nc_prior = r_mfm_prior_num_blocks(2e5, q, lambda = 1, dirichlet_alpha = 1)

  expect_lt(nc_tv(nc_chain, nc_prior, q), 0.06)
})

test_that("uncorrected prior-only chain misses the partition prior", {
  skip_on_cran()

  q = 5
  draws = sample_ggm_prior(
    p = q, n_samples = 6000, n_warmup = 500,
    seed = 33, verbose = FALSE, spec = "joint",
    update_method = "gibbs", edge_prior = sbm_prior(),
    apply_correction = FALSE
  )
  nc_chain = apply(draws$allocations, 1, function(z) length(unique(z)))

  set.seed(34)
  nc_prior = r_mfm_prior_num_blocks(2e5, q, lambda = 1, dirichlet_alpha = 1)

  expect_gt(nc_tv(nc_chain, nc_prior, q), 0.10)
  expect_lt(mean(nc_chain), mean(nc_prior))
})

test_that("bgm() with an sbm prior applies the correction end to end", {
  skip_on_cran()
  old = options(
    bgms.correction_cache_dir = file.path(tempdir(), "bgms-ctable-identity")
  )
  on.exit(options(old), add = TRUE)

  set.seed(9)
  x = matrix(rnorm(60 * 5), 60, 5)
  fit = suppressMessages(bgm(x,
    variable_type = "continuous",
    edge_prior = sbm_prior(),
    update_method = "gibbs",
    iter = 300, warmup = 300, chains = 1, cores = 1,
    display_progress = "none", verbose = FALSE
  ))

  expect_s3_class(fit, "bgms")
  expect_false(is.null(fit$posterior_num_blocks))
  expect_false(is.null(fit$posterior_mean_indicator))
})
