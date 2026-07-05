# Tests for sample_ggm_prior(): structure, prior plumbing, edge constraints,
# and input validation. Sampling correctness (uniform SBC ranks) is covered
# in test-sbc-ggm.R.

short_run = function(p, n_samples = 50L, n_warmup = 50L, seed = 1L, ...) {
  sample_ggm_prior(
    p = p, n_samples = n_samples, n_warmup = n_warmup,
    seed = seed, verbose = FALSE, ...
  )
}


# ---- Return value structure --------------------------------------------------

test_that("sample_ggm_prior returns the documented list shape", {
  p = 4L
  n_samples = 60L
  draws = short_run(p = p, n_samples = n_samples)

  expect_named(
    draws,
    c(
      "K_offdiag", "K_diag", "offdiag_names", "diag_names",
      "step_size", "edge_indicators"
    )
  )

  n_edges = p * (p - 1L) / 2L
  expect_equal(dim(draws$K_offdiag), c(n_samples, n_edges))
  expect_equal(dim(draws$K_diag), c(n_samples, p))
  expect_length(draws$offdiag_names, n_edges)
  expect_length(draws$diag_names, p)
  expect_equal(draws$diag_names, paste0("K_", seq_len(p), "_", seq_len(p)))
  expect_true(all(grepl("^K_[0-9]+_[0-9]+$", draws$offdiag_names)))

  expect_true(all(is.finite(draws$K_offdiag)))
  expect_true(all(is.finite(draws$K_diag)))
  expect_true(all(draws$K_diag > 0)) # diagonals are positive (Cholesky)

  expect_equal(dim(draws$edge_indicators), c(p, p))
  expect_true(all(draws$edge_indicators == 1L))
})


# ---- edge_indicators constraint ----------------------------------------------

test_that("excluded edges are exactly zero in K_offdiag", {
  p = 4L
  E = matrix(1L, p, p)
  # Drop edges (1, 4) and (2, 3) on both triangles.
  E[1, 4] = E[4, 1] = 0L
  E[2, 3] = E[3, 2] = 0L

  draws = short_run(p = p, n_samples = 80L, edge_indicators = E)

  colnames(draws$K_offdiag) = draws$offdiag_names
  # Excluded edges are zeroed out by the constraint structure / RATTLE
  # projection; with the Cholesky parameterization "structural" zeros
  # (Phi[i, j] = 0 directly) are exact, while quadratic constraints are
  # held to numerical tolerance.
  expect_true(all(draws$K_offdiag[, "K_1_4"] == 0))
  expect_true(max(abs(draws$K_offdiag[, "K_2_3"])) < 1e-8)

  # Included edges should be (almost surely) non-zero somewhere.
  expect_true(any(draws$K_offdiag[, "K_1_2"] != 0))
  expect_true(any(draws$K_offdiag[, "K_3_4"] != 0))

  # The returned edge_indicators should round-trip.
  expect_equal(draws$edge_indicators, E)
})


# ---- Scale prior plumbing ----------------------------------------------------

test_that("gamma_prior(shape, rate) shifts diagonal mean toward 2*shape/rate", {
  # The diagonal prior is on the partial-association diagonal
  # -K_yy_{ii} = K_{ii}/2, so gamma_prior(shape, rate) implies
  #   K_ii ~ 2 * Gamma(shape, rate),  mean = 2 * shape/rate.
  # Gamma(4, 2) -> mean(K_ii) = 4.0; Gamma(1, 1) -> mean(K_ii) = 2.0.
  draws_unit = short_run(
    p = 3L, n_samples = 400L, n_warmup = 200L,
    precision_scale_prior = gamma_prior(shape = 1, rate = 1)
  )
  draws_heavy = short_run(
    p = 3L, n_samples = 400L, n_warmup = 200L,
    precision_scale_prior = gamma_prior(shape = 4, rate = 2)
  )

  expect_lt(mean(draws_unit$K_diag), mean(draws_heavy$K_diag))
  expect_equal(mean(draws_heavy$K_diag), 4.0, tolerance = 0.6)
})


test_that("gamma_prior(eta) resolves against the interaction-prior scale", {
  # eta = 1 at cauchy scale 2.5 resolves to rate = 0.4; the same raw-rate
  # spec with the same seed gives identical draws.
  d_eta = short_run(
    p = 3L, n_samples = 100L, n_warmup = 100L,
    precision_scale_prior = gamma_prior(shape = 1, eta = 1)
  )
  d_rate = short_run(
    p = 3L, n_samples = 100L, n_warmup = 100L,
    precision_scale_prior = gamma_prior(shape = 1, rate = 1 / 2.5)
  )
  expect_equal(d_eta$K_diag, d_rate$K_diag)
  expect_equal(d_eta$K_offdiag, d_rate$K_offdiag)
})


test_that("exponential_prior(rate) is equivalent to gamma_prior(1, rate)", {
  d_exp = short_run(
    p = 3L, n_samples = 200L, n_warmup = 100L,
    precision_scale_prior = exponential_prior(rate = 2)
  )
  d_gamma = short_run(
    p = 3L, n_samples = 200L, n_warmup = 100L,
    precision_scale_prior = gamma_prior(shape = 1, rate = 2)
  )
  # Same seed + same prior parameters => identical draws.
  expect_equal(d_exp$K_diag, d_gamma$K_diag)
  expect_equal(d_exp$K_offdiag, d_gamma$K_offdiag)
})


# ---- Interaction prior plumbing ---------------------------------------------

test_that("normal_prior shrinks off-diagonals more than the default Cauchy", {
  draws_cauchy = short_run(
    p = 3L, n_samples = 300L, n_warmup = 200L,
    interaction_prior = cauchy_prior(scale = 2.5)
  )
  draws_normal = short_run(
    p = 3L, n_samples = 300L, n_warmup = 200L,
    interaction_prior = normal_prior(scale = 0.1)
  )

  expect_lt(
    sd(as.numeric(draws_normal$K_offdiag)),
    sd(as.numeric(draws_cauchy$K_offdiag))
  )
})


# ---- Error paths -------------------------------------------------------------

test_that("invalid scalar arguments error early with informative messages", {
  expect_error(short_run(p = 1L, n_samples = 10L), "'p' must be >= 2")
  expect_error(short_run(p = 3.5, n_samples = 10L), "'p' must be an integer")
  expect_error(short_run(p = 3L, n_samples = 0L), "'n_samples' must be >= 1")
  expect_error(
    short_run(p = 3L, n_samples = 10L, step_size = -0.1),
    "'step_size' must be positive"
  )
  expect_error(
    short_run(p = 3L, n_samples = 10L, max_depth = 0L),
    "'max_depth' must be >= 1"
  )
  expect_error(
    sample_ggm_prior(p = 3L, n_samples = 10L, n_warmup = 10L, verbose = NA),
    "'verbose' must be TRUE or FALSE"
  )
})


test_that("beta_prime_prior is rejected for interaction_prior", {
  expect_error(
    short_run(
      p = 3L, n_samples = 10L,
      interaction_prior = beta_prime_prior(1, 1)
    ),
    "beta_prime_prior\\(\\) is not supported"
  )
})


test_that("scale prior must be a bgms_scale_prior", {
  expect_error(
    short_run(
      p = 3L, n_samples = 10L,
      precision_scale_prior = cauchy_prior(1)
    ),
    "bgms_scale_prior"
  )
})


test_that("interaction prior must be a bgms_parameter_prior", {
  expect_error(
    short_run(
      p = 3L, n_samples = 10L,
      interaction_prior = gamma_prior(1, 1)
    ),
    "bgms_parameter_prior"
  )
})


test_that("malformed edge_indicators are rejected", {
  expect_error(
    short_run(
      p = 3L, n_samples = 10L,
      edge_indicators = matrix(1L, 2, 2)
    ),
    "p x p matrix"
  )
  E_asym = matrix(c(
    1, 1, 0,
    0, 1, 1,
    0, 1, 1
  ), nrow = 3, byrow = TRUE)
  expect_error(
    short_run(p = 3L, n_samples = 10L, edge_indicators = E_asym),
    "symmetric"
  )
  E_baddiag = matrix(1L, 3, 3)
  E_baddiag[1, 1] = 0L
  expect_error(
    short_run(p = 3L, n_samples = 10L, edge_indicators = E_baddiag),
    "diagonal"
  )
  E_bad = matrix(1L, 3, 3)
  E_bad[1, 2] = E_bad[2, 1] = 2L
  expect_error(
    short_run(p = 3L, n_samples = 10L, edge_indicators = E_bad),
    "0 or 1"
  )
})


# ---- Ancestral initialization (spec = "joint") -------------------------------

test_that("ancestral indicator draws are seed-keyed and restore the RNG", {
  ep = bgms:::unpack_indicator_prior(bernoulli_prior(0.2), num_variables = 8L)
  set.seed(123)
  before = .Random.seed
  g1 = bgms:::ggm_prior_ancestral_indicators(8L, ep, seed = 7L)
  expect_identical(.Random.seed, before)
  g2 = bgms:::ggm_prior_ancestral_indicators(8L, ep, seed = 7L)
  expect_identical(g1, g2)

  expect_true(isSymmetric(g1))
  expect_true(all(diag(g1) == 1L))
  expect_true(all(g1 %in% c(0L, 1L)))
})

test_that("ancestral indicator draws track the edge prior", {
  p = 20L
  offdiag_mean = function(g) {
    mean(g[upper.tri(g)])
  }

  ep_dense = bgms:::unpack_indicator_prior(
    bernoulli_prior(0.9),
    num_variables = p
  )
  ep_sparse = bgms:::unpack_indicator_prior(
    bernoulli_prior(0.1),
    num_variables = p
  )
  g_dense = bgms:::ggm_prior_ancestral_indicators(p, ep_dense, seed = 11L)
  g_sparse = bgms:::ggm_prior_ancestral_indicators(p, ep_sparse, seed = 11L)
  expect_gt(offdiag_mean(g_dense), offdiag_mean(g_sparse))

  # Beta-Bernoulli: a concentrated hyperprior pins the drawn theta.
  ep_bb = bgms:::unpack_indicator_prior(
    beta_bernoulli_prior(alpha = 200, beta = 1),
    num_variables = p
  )
  g_bb = bgms:::ggm_prior_ancestral_indicators(p, ep_bb, seed = 11L)
  expect_gt(offdiag_mean(g_bb), 0.8)

  # SBM: valid draw with block-pair probabilities in (0, 1).
  ep_sbm = bgms:::unpack_indicator_prior(
    sbm_prior(),
    num_variables = p
  )
  g_sbm = bgms:::ggm_prior_ancestral_indicators(p, ep_sbm, seed = 11L)
  expect_true(isSymmetric(g_sbm))
  expect_true(all(diag(g_sbm) == 1L))
  expect_true(all(g_sbm %in% c(0L, 1L)))
})

test_that("joint-spec chains run from the ancestral start for each prior", {
  skip_on_cran()
  for(prior in list(
    bernoulli_prior(0.3),
    beta_bernoulli_prior(alpha = 2, beta = 4),
    sbm_prior()
  )) {
    draws = sample_ggm_prior(
      p = 6L, n_samples = 25L, n_warmup = 50L,
      spec = "joint", edge_prior = prior, apply_correction = FALSE,
      seed = 3L, verbose = FALSE
    )
    expect_equal(dim(draws$edge_indicators), c(25L, 15L))
    expect_true(all(is.finite(draws$K_diag)))
  }
})
