# --------------------------------------------------------------------------- #
# Row-block Gibbs within-step for the Gaussian graphical model
# (update_method = "gibbs", fixed graph).
#
# Two kinds of checks:
#   1. Closed-form anchors. The C++ entry ggm_test_gibbs_sweep drives the row
#      sweep directly. On an empty graph at n = 0 every row collapses to its
#      diagonal, whose marginal is a known Gamma, so the draw is checked against
#      that Gamma. This also exercises the Gamma(shape != 1) independent-MH,
#      whose target is Gamma(shape, rate/2).
#   2. Agreement with adaptive-Metropolis. On a fixed complete graph both
#      samplers target the same K | graph posterior, so their posterior means
#      must agree within MCMC error. This catches a wrong target that a
#      "runs without error" check would miss.
# --------------------------------------------------------------------------- #


# ---- Closed-form anchors ---------------------------------------------------- #

test_that("row sweep draws K_ii from Gamma(1, rate/2) on an empty graph at n = 0", {
  p = 4L
  edges = matrix(0L, p, p)
  S = matrix(0, p, p) # n = 0: sufficient statistics are irrelevant
  beta0 = 2.0
  n_sweeps = 5000L

  out = ggm_test_gibbs_sweep(
    suf_stat = S, n = 0L, edge_indicators = edges,
    pairwise_scale = 1.0, gamma_shape = 1.0, gamma_rate = beta0,
    n_sweeps = n_sweeps, seed = 42L
  )

  burn = 201:n_sweeps # forget the identity init
  K_ii = sapply(seq_len(p), function(i) out$K_samples[i, i, burn])

  # Gamma(shape = 1, rate = beta0 / 2): mean 2/beta0, var 4/beta0^2.
  theory_mean = 2.0 / beta0
  tol_mean = 4 * sqrt((4.0 / beta0^2) / length(burn))
  expect_true(all(abs(colMeans(K_ii) - theory_mean) < tol_mean))

  ks_p = vapply(seq_len(p), function(i) {
    suppressWarnings(
      stats::ks.test(K_ii[, i], "pgamma", shape = 1, rate = beta0 / 2)$p.value
    )
  }, numeric(1L))
  expect_gt(min(ks_p), 0.01 / p)
})


test_that("Gamma shape != 1 targets Gamma(shape, rate/2) via the independent-MH", {
  # The shape != 1 path is an independent-MH chain, so its draws are
  # autocorrelated (rejections repeat values). A KS test would be invalid;
  # check the first two moments instead, pooled across rows.
  p = 4L
  edges = matrix(0L, p, p)
  S = matrix(0, p, p)
  alpha = 2.0
  beta0 = 2.0
  n_sweeps = 8000L

  out = ggm_test_gibbs_sweep(
    suf_stat = S, n = 0L, edge_indicators = edges,
    pairwise_scale = 1.0, gamma_shape = alpha, gamma_rate = beta0,
    n_sweeps = n_sweeps, seed = 7L
  )

  burn = 1001:n_sweeps
  K_ii = as.vector(sapply(seq_len(p), function(i) out$K_samples[i, i, burn]))

  # Gamma(shape = alpha, rate = beta0 / 2): mean 2 alpha / beta0, var 4 alpha / beta0^2.
  theory_mean = 2.0 * alpha / beta0
  theory_var = 4.0 * alpha / beta0^2
  expect_lt(abs(mean(K_ii) - theory_mean) / theory_mean, 0.05)
  expect_lt(abs(var(K_ii) - theory_var) / theory_var, 0.15)
})


# ---- Gibbs-vs-AM agreement helpers ------------------------------------------ #

# Posterior means of K (diagonal + off-diagonal) from a fixed-graph bgm() fit.
ggm_K_means = function(Y, update_method, interaction_prior, alpha, delta,
                       iter, warmup, seed = 1L) {
  fit = bgm(
    Y,
    variable_type = "continuous",
    interaction_prior = interaction_prior,
    precision_scale_prior = gamma_prior(shape = alpha, rate = 1),
    delta = delta,
    edge_selection = FALSE,
    iter = iter, warmup = warmup,
    update_method = update_method,
    chains = 1L, cores = 1L, seed = seed,
    display_progress = "none", verbose = FALSE
  )
  raw = S7::prop(fit, "raw_samples")
  list(
    diag = colMeans(raw$main[[1L]]),
    off  = colMeans(raw$pairwise[[1L]])
  )
}

# A small dense GGM, centered as bgm() expects.
ggm_agreement_data = function(seed = 7L, p = 6L, n = 200L) {
  set.seed(seed)
  K_true = diag(p) + 0.3 * (abs(row(diag(p)) - col(diag(p))) == 1)
  Y = MASS::mvrnorm(n, mu = rep(0, p), Sigma = solve(K_true))
  Y = scale(Y, center = TRUE, scale = FALSE)
  colnames(Y) = paste0("V", seq_len(p))
  Y
}

expect_K_agreement = function(gibbs, am, tol_abs, tol_rel) {
  d_diag = abs(gibbs$diag - am$diag)
  d_off = abs(gibbs$off - am$off)
  expect_lt(max(d_diag, d_off), tol_abs)
  expect_lt(mean(c(d_diag, d_off)) / mean(abs(c(am$diag, am$off))), tol_rel)
}


# ---- Agreement with adaptive-Metropolis ------------------------------------- #

test_that("gibbs routes through bgm() and agrees with AM (Normal, alpha=1, delta=0)", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  Y = ggm_agreement_data()
  args = list(
    Y = Y, interaction_prior = normal_prior(scale = 1),
    alpha = 1, delta = 0, iter = 2500L, warmup = 600L
  )
  am = do.call(ggm_K_means, c(args, update_method = "adaptive-metropolis"))
  gibb = do.call(ggm_K_means, c(args, update_method = "gibbs"))
  expect_K_agreement(gibb, am, tol_abs = 0.06, tol_rel = 0.025)
})


test_that("gibbs agrees with AM at delta = 0.5 (xi shape shift)", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  Y = ggm_agreement_data()
  args = list(
    Y = Y, interaction_prior = normal_prior(scale = 1),
    alpha = 1, delta = 0.5, iter = 2500L, warmup = 600L
  )
  am = do.call(ggm_K_means, c(args, update_method = "adaptive-metropolis"))
  gibb = do.call(ggm_K_means, c(args, update_method = "gibbs"))
  expect_K_agreement(gibb, am, tol_abs = 0.06, tol_rel = 0.025)
})


test_that("gibbs agrees with AM at alpha = 2 (independent-MH on the diagonal)", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  Y = ggm_agreement_data()
  args = list(
    Y = Y, interaction_prior = normal_prior(scale = 1),
    alpha = 2, delta = 0, iter = 2500L, warmup = 600L
  )
  am = do.call(ggm_K_means, c(args, update_method = "adaptive-metropolis"))
  gibb = do.call(ggm_K_means, c(args, update_method = "gibbs"))
  expect_K_agreement(gibb, am, tol_abs = 0.07, tol_rel = 0.035)
})


test_that("gibbs agrees with AM under a Cauchy slab (scale mixture of normals)", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  Y = ggm_agreement_data()
  # Cauchy mixing is slower than Normal: longer chain, looser tolerance.
  args = list(
    Y = Y, interaction_prior = cauchy_prior(scale = 1),
    alpha = 1, delta = 0, iter = 4000L, warmup = 1000L
  )
  am = do.call(ggm_K_means, c(args, update_method = "adaptive-metropolis"))
  gibb = do.call(ggm_K_means, c(args, update_method = "gibbs"))
  expect_K_agreement(gibb, am, tol_abs = 0.10, tol_rel = 0.05)
})


# ---- Edge selection: between-model agreement with NUTS ---------------------- #

# Posterior inclusion probabilities from an edge-selection fit.
ggm_pips = function(Y, update_method, iter, warmup, shape = 1, seed = 7L) {
  fit = bgm(
    Y,
    variable_type = "continuous",
    interaction_prior = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = shape, rate = 1),
    edge_selection = TRUE,
    iter = iter, warmup = warmup,
    update_method = update_method,
    chains = 1L, cores = 1L, seed = seed,
    display_progress = "none", verbose = FALSE
  )
  colMeans(S7::prop(fit, "raw_samples")$indicator[[1L]])
}

test_that("gibbs edge selection recovers the same inclusion probabilities as NUTS", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  set.seed(11)
  p = 8L
  n = 250L
  # Sparse truth: a chain graph.
  K = diag(p)
  for(i in seq_len(p - 1L)) K[i, i + 1L] = K[i + 1L, i] = -0.35
  diag(K) = 2
  Y = MASS::mvrnorm(n, mu = rep(0, p), Sigma = solve(K))
  Y = scale(Y, center = TRUE, scale = FALSE)
  colnames(Y) = paste0("V", seq_len(p))

  pip_gibbs = ggm_pips(Y, "gibbs", iter = 4000L, warmup = 1500L)
  pip_nuts = ggm_pips(Y, "nuts", iter = 4000L, warmup = 1500L)

  expect_lt(max(abs(pip_gibbs - pip_nuts)), 0.10)
  expect_lt(mean(abs(pip_gibbs - pip_nuts)), 0.03)
})

test_that("gibbs edge selection agrees with NUTS at Gamma shape = 2", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  set.seed(11)
  p = 8L
  n = 250L
  K = diag(p)
  for(i in seq_len(p - 1L)) K[i, i + 1L] = K[i + 1L, i] = -0.35
  diag(K) = 2
  Y = MASS::mvrnorm(n, mu = rep(0, p), Sigma = solve(K))
  Y = scale(Y, center = TRUE, scale = FALSE)
  colnames(Y) = paste0("V", seq_len(p))

  pip_gibbs = ggm_pips(Y, "gibbs", iter = 4000L, warmup = 1500L, shape = 2)
  pip_nuts = ggm_pips(Y, "nuts", iter = 4000L, warmup = 1500L, shape = 2)

  expect_lt(max(abs(pip_gibbs - pip_nuts)), 0.10)
  expect_lt(mean(abs(pip_gibbs - pip_nuts)), 0.03)
})
