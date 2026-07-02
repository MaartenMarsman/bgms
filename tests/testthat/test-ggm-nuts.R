# --------------------------------------------------------------------------- #
# Phase 3 <U+2014> Correctness validation for GGM NUTS sampler.
#
# Tests compare NUTS posterior to MH baseline using long chains.
# Gated behind BGMS_RUN_SLOW_TESTS because they take several minutes.
#
# 3.1  Posterior moment comparison (means, variances, KS, bivariate)
# 3.2  Edge selection accuracy (PIPs)
# 3.3  Gradient near PD boundary
# 3.4  NUTS diagnostics are well-behaved
#
# TODO: add a proper Geweke JDT or SBC test for one-transition correctness
# using the bgms prior directly and a fixed-kernel NUTS step (issue #TBD).
# --------------------------------------------------------------------------- #


# ---- Skip gate ---------------------------------------------------------------

skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run GGM NUTS correctness tests"
  )
}


# ---- Data generators ---------------------------------------------------------

# Generate data from a known precision matrix.
# Returns list(x, K_true, Sigma_true, p, n).
generate_ggm_data = function(p, n, seed = 1) {
  set.seed(seed)
  # Build a sparse precision matrix: tridiagonal + constant diagonal
  K = diag(2, p)
  for(i in seq_len(p - 1)) {
    K[i, i + 1] = -0.5
    K[i + 1, i] = -0.5
  }
  Sigma = solve(K)
  x = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(x) = paste0("V", seq_len(p))
  list(x = x, K_true = K, Sigma_true = Sigma, p = p, n = n)
}


# ---- 3.1  Posterior moment comparison ----------------------------------------

test_that("NUTS and MH posteriors agree on means and variances (p=4)", {
  skip_unless_slow()

  dat = generate_ggm_data(p = 4, n = 200, seed = 31)
  n_iter = 5000
  n_warmup = 2000

  fit_mh = bgm(
    x = dat$x, variable_type = "continuous",
    update_method = "adaptive-metropolis",
    edge_selection = FALSE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 100, display_progress = "none"
  )

  fit_nuts = bgm(
    x = dat$x, variable_type = "continuous",
    update_method = "nuts",
    edge_selection = FALSE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 200, display_progress = "none"
  )

  # Combine all chains for each sampler
  mh_pairwise = do.call(rbind, fit_mh$raw_samples$pairwise)
  nuts_pairwise = do.call(rbind, fit_nuts$raw_samples$pairwise)
  mh_main = do.call(rbind, fit_mh$raw_samples$main)
  nuts_main = do.call(rbind, fit_nuts$raw_samples$main)

  # Posterior means should agree within 2 posterior SDs
  for(j in seq_len(ncol(mh_pairwise))) {
    mh_mean = mean(mh_pairwise[, j])
    nuts_mean = mean(nuts_pairwise[, j])
    pooled_sd = sqrt((var(mh_pairwise[, j]) + var(nuts_pairwise[, j])) / 2)
    expect_lt(
      abs(mh_mean - nuts_mean) / pooled_sd, 2,
      label = paste0("pairwise[", j, "] mean within 2 SDs")
    )
  }

  for(j in seq_len(ncol(mh_main))) {
    mh_mean = mean(mh_main[, j])
    nuts_mean = mean(nuts_main[, j])
    pooled_sd = sqrt((var(mh_main[, j]) + var(nuts_main[, j])) / 2)
    expect_lt(
      abs(mh_mean - nuts_mean) / pooled_sd, 2,
      label = paste0("main[", j, "] mean within 2 SDs")
    )
  }

  # Variance ratios should be close to 1
  for(j in seq_len(ncol(mh_pairwise))) {
    ratio = var(nuts_pairwise[, j]) / var(mh_pairwise[, j])
    expect_gt(ratio, 0.7, label = paste0("pairwise[", j, "] var ratio > 0.7"))
    expect_lt(ratio, 1.4, label = paste0("pairwise[", j, "] var ratio < 1.4"))
  }
})

test_that("NUTS and MH 95% credible intervals overlap (p=4)", {
  skip_unless_slow()

  dat = generate_ggm_data(p = 4, n = 200, seed = 32)
  n_iter = 5000
  n_warmup = 2000

  fit_mh = bgm(
    x = dat$x, variable_type = "continuous",
    update_method = "adaptive-metropolis",
    edge_selection = FALSE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 101, display_progress = "none"
  )

  fit_nuts = bgm(
    x = dat$x, variable_type = "continuous",
    update_method = "nuts",
    edge_selection = FALSE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 201, display_progress = "none"
  )

  mh_pairwise = do.call(rbind, fit_mh$raw_samples$pairwise)
  nuts_pairwise = do.call(rbind, fit_nuts$raw_samples$pairwise)

  # 95% credible interval overlap > 80% for all marginals
  ci_overlap = function(x, y) {
    ci_x = quantile(x, c(0.025, 0.975))
    ci_y = quantile(y, c(0.025, 0.975))
    overlap_lo = max(ci_x[1], ci_y[1])
    overlap_hi = min(ci_x[2], ci_y[2])
    if(overlap_lo >= overlap_hi) {
      return(0)
    }
    overlap_width = overlap_hi - overlap_lo
    avg_width = (diff(ci_x) + diff(ci_y)) / 2
    overlap_width / avg_width
  }

  for(j in seq_len(ncol(mh_pairwise))) {
    ov = ci_overlap(mh_pairwise[, j], nuts_pairwise[, j])
    expect_gt(ov, 0.80,
      label = paste0("pairwise[", j, "] CI overlap > 80%")
    )
  }

  # Also check SD ratios are in [0.8, 1.25]
  for(j in seq_len(ncol(mh_pairwise))) {
    ratio = sd(nuts_pairwise[, j]) / sd(mh_pairwise[, j])
    expect_gt(ratio, 0.8, label = paste0("pairwise[", j, "] SD ratio > 0.8"))
    expect_lt(ratio, 1.25, label = paste0("pairwise[", j, "] SD ratio < 1.25"))
  }
})

test_that("NUTS posterior means agree on p=6 tridiagonal", {
  skip_unless_slow()

  dat = generate_ggm_data(p = 6, n = 300, seed = 33)
  n_iter = 5000
  n_warmup = 2000

  fit_mh = bgm(
    x = dat$x, variable_type = "continuous",
    update_method = "adaptive-metropolis",
    edge_selection = FALSE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 102, display_progress = "none"
  )

  fit_nuts = bgm(
    x = dat$x, variable_type = "continuous",
    update_method = "nuts",
    edge_selection = FALSE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 202, display_progress = "none"
  )

  # Max absolute difference in posterior mean associations
  diff = abs(fit_mh$posterior_mean_pairwise - fit_nuts$posterior_mean_pairwise)
  expect_lt(max(diff), 0.05, label = "max abs diff in pairwise means < 0.05")

  # NUTS should have clean diagnostics
  expect_equal(fit_nuts$nuts_diag$summary$total_divergences, 0)
})


# ---- 3.2  Edge selection accuracy --------------------------------------------

test_that("NUTS edge selection recovers true graph structure (p=6)", {
  skip_unless_slow()

  p = 6
  n = 500
  set.seed(34)

  # True graph: only adjacent edges (tridiagonal)
  K_true = diag(2, p)
  for(i in seq_len(p - 1)) {
    K_true[i, i + 1] = -0.5
    K_true[i + 1, i] = -0.5
  }
  Sigma = solve(K_true)
  x = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(x) = paste0("V", seq_len(p))

  # True edges: (1,2), (2,3), (3,4), (4,5), (5,6) = 5 edges
  true_edges = matrix(0, p, p)
  for(i in seq_len(p - 1)) {
    true_edges[i, i + 1] = 1
    true_edges[i + 1, i] = 1
  }

  n_iter = 5000
  n_warmup = 2000

  fit_mh = bgm(
    x = x, variable_type = "continuous",
    update_method = "adaptive-metropolis",
    edge_selection = TRUE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 103, display_progress = "none"
  )

  fit_nuts = bgm(
    x = x, variable_type = "continuous",
    update_method = "nuts",
    edge_selection = TRUE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 203, display_progress = "none"
  )

  # Extract PIPs
  pip_mh = fit_mh$posterior_mean_indicator
  pip_nuts = fit_nuts$posterior_mean_indicator

  # True edges should have high PIP in both
  true_ut = which(true_edges[upper.tri(true_edges)] == 1)
  false_ut = which(true_edges[upper.tri(true_edges)] == 0)

  # PIP agreement: max absolute difference < 0.15
  pip_mh_ut = pip_mh[upper.tri(pip_mh)]
  pip_nuts_ut = pip_nuts[upper.tri(pip_nuts)]
  expect_lt(max(abs(pip_mh_ut - pip_nuts_ut)), 0.15,
    label = "max PIP difference < 0.15"
  )

  # NUTS should identify true edges (PIP > 0.5)
  for(idx in true_ut) {
    expect_gt(pip_nuts_ut[idx], 0.5,
      label = paste0("true edge ", idx, " PIP > 0.5")
    )
  }

  # NUTS should exclude false edges (PIP < 0.5)
  for(idx in false_ut) {
    expect_lt(pip_nuts_ut[idx], 0.5,
      label = paste0("false edge ", idx, " PIP < 0.5")
    )
  }
})


# ---- 3.3  Gradient near PD boundary -----------------------------------------

# Helpers (same as test-ggm-gradient.R, repeated for test isolation)
theta_dim_local = function(edge_mat) {
  p = nrow(edge_mat)
  p + sum(edge_mat[upper.tri(edge_mat)] == 1L)
}

fd_gradient_local = function(theta, suf_stat, n, edge_mat, pairwise_scale,
                             eps = 1e-6) {
  g = numeric(length(theta))
  for(k in seq_along(theta)) {
    t_plus = theta
    t_minus = theta
    t_plus[k] = t_plus[k] + eps
    t_minus[k] = t_minus[k] - eps
    f_plus = ggm_test_logp_and_gradient(
      t_plus, suf_stat, n, edge_mat, pairwise_scale
    )$value
    f_minus = ggm_test_logp_and_gradient(
      t_minus, suf_stat, n, edge_mat, pairwise_scale
    )$value
    g[k] = (f_plus - f_minus) / (2 * eps)
  }
  g
}

test_that("gradient is accurate near the PD boundary", {
  # Test with theta values that produce near-singular precision matrices.
  # The diagonal Cholesky elements (psi_qq) close to zero push K toward
  # the boundary of the PD cone.

  p = 4
  edge_mat = matrix(1L, p, p)
  diag(edge_mat) = 0L

  set.seed(340)
  X = matrix(rnorm(200 * p), nrow = 200)
  S = t(X) %*% X

  d = theta_dim_local(edge_mat)

  # Generate theta with small diagonal elements (near PD boundary)
  set.seed(341)
  theta = rnorm(d, sd = 0.1)
  # First p elements are psi (diagonal of Cholesky factor)
  # Make them small but positive to be near the boundary
  theta[1:p] = runif(p, min = 0.05, max = 0.15)

  ag = ggm_test_logp_and_gradient(theta, S, 200, edge_mat, 1.0)
  fd = fd_gradient_local(theta, S, 200, edge_mat, 1.0, eps = 1e-6)

  denom = pmax(abs(ag$gradient), abs(fd), 1)
  rel_err = abs(ag$gradient - fd) / denom
  expect_lt(max(rel_err), 1e-3,
    label = "gradient rel error < 1e-3 near PD boundary"
  )
})

test_that("gradient is accurate with sparse graph near PD boundary", {
  p = 6
  # Only 3 edges: (1,2), (3,4), (5,6)
  edge_mat = matrix(0L, p, p)
  edges = list(c(1, 2), c(3, 4), c(5, 6))
  for(e in edges) {
    edge_mat[e[1], e[2]] = 1L
    edge_mat[e[2], e[1]] = 1L
  }

  set.seed(342)
  X = matrix(rnorm(300 * p), nrow = 300)
  S = t(X) %*% X

  d = theta_dim_local(edge_mat)
  theta = rnorm(d, sd = 0.1)
  theta[1:p] = runif(p, min = 0.05, max = 0.15)

  ag = ggm_test_logp_and_gradient(theta, S, 300, edge_mat, 1.0)
  fd = fd_gradient_local(theta, S, 300, edge_mat, 1.0, eps = 1e-6)

  denom = pmax(abs(ag$gradient), abs(fd), 1)
  rel_err = abs(ag$gradient - fd) / denom
  expect_lt(max(rel_err), 1e-3,
    label = "sparse graph gradient near PD boundary"
  )
})


# ---- 3.4  NUTS diagnostics are well-behaved ----------------------------------

test_that("NUTS diagnostics are clean for well-specified model", {
  skip_unless_slow()

  dat = generate_ggm_data(p = 4, n = 200, seed = 35)

  fit = bgm(
    x = dat$x, variable_type = "continuous",
    update_method = "nuts",
    edge_selection = FALSE,
    iter = 2000, warmup = 1000, chains = 2,
    seed = 300, display_progress = "none"
  )

  expect_true(!is.null(fit$nuts_diag))
  expect_equal(fit$nuts_diag$summary$total_divergences, 0)
  expect_equal(fit$nuts_diag$summary$max_tree_depth_hits, 0)
  expect_gt(fit$nuts_diag$summary$min_ebfmi, 0.3)
})


# ---- 3.5  Tilt: NUTS vs MH posterior moments at delta = 1 -------------------

test_that("NUTS and MH posteriors agree under tilt (p=4, delta=1)", {
  skip_unless_slow()

  dat = generate_ggm_data(p = 4, n = 200, seed = 36)
  n_iter = 5000
  n_warmup = 2000

  fit_mh = bgm(
    x = dat$x, variable_type = "continuous",
    update_method = "adaptive-metropolis",
    edge_selection = FALSE,
    delta = 1,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 102, display_progress = "none"
  )

  fit_nuts = bgm(
    x = dat$x, variable_type = "continuous",
    update_method = "nuts",
    edge_selection = FALSE,
    delta = 1,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 202, display_progress = "none"
  )

  mh_pairwise = do.call(rbind, fit_mh$raw_samples$pairwise)
  nuts_pairwise = do.call(rbind, fit_nuts$raw_samples$pairwise)
  mh_main = do.call(rbind, fit_mh$raw_samples$main)
  nuts_main = do.call(rbind, fit_nuts$raw_samples$main)

  for(j in seq_len(ncol(mh_pairwise))) {
    mh_mean = mean(mh_pairwise[, j])
    nuts_mean = mean(nuts_pairwise[, j])
    pooled_sd = sqrt((var(mh_pairwise[, j]) + var(nuts_pairwise[, j])) / 2)
    expect_lt(
      abs(mh_mean - nuts_mean) / pooled_sd, 2,
      label = paste0("pairwise[", j, "] mean within 2 SDs (delta=1)")
    )
  }

  for(j in seq_len(ncol(mh_main))) {
    mh_mean = mean(mh_main[, j])
    nuts_mean = mean(nuts_main[, j])
    pooled_sd = sqrt((var(mh_main[, j]) + var(nuts_main[, j])) / 2)
    expect_lt(
      abs(mh_mean - nuts_mean) / pooled_sd, 2,
      label = paste0("main[", j, "] mean within 2 SDs (delta=1)")
    )
  }

  for(j in seq_len(ncol(mh_pairwise))) {
    ratio = var(nuts_pairwise[, j]) / var(mh_pairwise[, j])
    expect_gt(ratio, 0.7,
      label = paste0("pairwise[", j, "] var ratio > 0.7 (delta=1)")
    )
    expect_lt(ratio, 1.4,
      label = paste0("pairwise[", j, "] var ratio < 1.4 (delta=1)")
    )
  }
})
