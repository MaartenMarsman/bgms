# --------------------------------------------------------------------------- #
# SBC for GGM (no edge selection)
#
# Simulation-based calibration checks that posterior ranks are uniform
# when data is generated from the prior predictive distribution. This
# validates the NUTS and Adaptive-Metropolis samplers for the GGM with
#
# p = 3, R = 200 replications, L = 999 posterior draws per replication.
# Prior: Cauchy(0, 2.5) off-diagonal, Gamma(1, 1) diagonal. PD
# constraint enforced by rejection sampling.
#
# Edge selection is tested separately via parameter recovery (PR.2).
# Omitting edge selection here avoids the tied-rank issue inherent
# in spike-and-slab SBC (discrete indicators + point-mass at zero
# make standard uniformity tests unreliable).
#
# Uniformity tested by KS test at alpha = 0.01 per parameter.
# Global chi-squared test as a fallback.
#
# Weekly certification tier (T2, BGMS_RUN_CERTIFICATION): SBC replicate suites.
# --------------------------------------------------------------------------- #


# ---- Skip gate ---------------------------------------------------------------

# ---- Prior sampler -----------------------------------------------------------

# Draw a symmetric PD precision matrix K from the GGM prior (no edge
# selection). The prior is specified on the partial-association scale
# Omega = -K/2 to match the bgms sampler convention:
#   Omega_ij    ~ Cauchy(0, scale)         (off-diagonal, i < j, symmetrised)
#   -Omega_ii   ~ Gamma(shape, rate)       (diagonal, positive)
# Implies on K:
#   K_ij = -2 * Omega_ij ~ Cauchy(0, 2*scale)
#   K_ii =  2 * (-Omega_ii) ~ 2 * Gamma(shape, rate)
# Rejection-samples until K is positive definite.
draw_prior_K = function(p, scale = 2.5, shape = 1, rate = 1,
                        max_tries = 10000) {
  for(attempt in seq_len(max_tries)) {
    K = matrix(0, p, p)

    # Off-diagonal (upper triangle): Omega ~ Cauchy(0, scale), K = -2 * Omega
    for(i in seq_len(p - 1)) {
      for(j in (i + 1):p) {
        omega_ij = rcauchy(1, 0, scale)
        K[i, j] = -2 * omega_ij
        K[j, i] = K[i, j]
      }
    }

    # Diagonal: -Omega_ii ~ Gamma(shape, rate), K_ii = 2 * (-Omega_ii)
    for(i in seq_len(p)) {
      K[i, i] = 2 * rgamma(1, shape = shape, rate = rate)
    }

    # Check positive definiteness
    ev = eigen(K, symmetric = TRUE, only.values = TRUE)$values
    if(all(ev > 1e-8)) {
      return(K)
    }
  }
  stop("draw_prior_K: failed to draw PD matrix in ", max_tries, " attempts")
}


# ---- Rank computation --------------------------------------------------------

# Compute SBC rank: number of posterior draws less than the true value.
# Returns a named vector of ranks (one per parameter).
compute_sbc_ranks = function(K_true, p, fit, thin_idx = NULL) {
  # Off-diagonal precision entries (raw precision scale)
  pw_samples = do.call(rbind, fit$raw_samples$pairwise)
  if(!is.null(thin_idx)) pw_samples = pw_samples[thin_idx, , drop = FALSE]

  # Diagonal precision entries
  main_samples = do.call(rbind, fit$raw_samples$main)
  if(!is.null(thin_idx)) main_samples = main_samples[thin_idx, , drop = FALSE]

  ranks = numeric(0)
  names_out = character(0)

  # Off-diagonal K entries
  col_idx = 0
  for(i in seq_len(p - 1)) {
    for(j in (i + 1):p) {
      col_idx = col_idx + 1
      true_k = K_true[i, j]
      ranks = c(ranks, sum(pw_samples[, col_idx] < true_k))
      names_out = c(names_out, paste0("K_", i, j))
    }
  }

  # Diagonal K entries
  for(i in seq_len(p)) {
    true_k = K_true[i, i]
    ranks = c(ranks, sum(main_samples[, i] < true_k))
    names_out = c(names_out, paste0("K_", i, i))
  }

  names(ranks) = names_out
  ranks
}


# ---- SBC test ----------------------------------------------------------------

test_that("SBC: GGM NUTS produces uniform ranks (p=3, no edge selection)", {
  skip_unless_certification()

  p = 3
  n = 100
  R = 200
  L = 999
  scale = 2.5

  set.seed(2026)

  # Pre-draw all prior K matrices to separate randomness
  prior_draws = vector("list", R)
  for(r in seq_len(R)) {
    prior_draws[[r]] = draw_prior_K(p, scale)
  }

  # Storage: 3 off-diagonal + 3 diagonal = 6 parameters
  n_off = p * (p - 1) / 2
  n_params = n_off + p
  ranks = matrix(NA_real_, nrow = R, ncol = n_params)

  for(r in seq_len(R)) {
    K_true = prior_draws[[r]]
    Sigma = solve(K_true)
    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = bgm(dat,
      variable_type = "continuous",
      iter = L, warmup = 1000, chains = 1,
      edge_selection = FALSE, update_method = "nuts",
      interaction_prior = cauchy_prior(scale = scale), delta = 0,
      precision_scale_prior = gamma_prior(shape = 1, rate = 1),
      display_progress = "none", seed = 2026L + r
    )

    ranks[r, ] = compute_sbc_ranks(K_true, p, fit)
    if(r == 1) colnames(ranks) = names(ranks[1, ])
  }

  # Per-parameter KS test: ranks / (L + 1) ~ Uniform(0, 1)
  n_fail_ks = 0
  for(j in seq_len(ncol(ranks))) {
    u = ranks[, j] / (L + 1)
    p_val = suppressWarnings(ks.test(u, "punif")$p.value)
    if(p_val <= 0.01) n_fail_ks = n_fail_ks + 1
  }

  # At alpha=0.01 with 6 parameters, allow at most 1 false positive
  max_fail = max(1, ceiling(n_params * 0.01 * 2))
  expect_true(n_fail_ks <= max_fail,
    info = sprintf(
      "SBC KS: %d/%d parameters failed (limit %d)",
      n_fail_ks, n_params, max_fail
    )
  )

  # Global chi-squared on aggregated ranks (20 bins)
  all_ranks = as.vector(ranks)
  bins = cut(all_ranks / (L + 1),
    breaks = seq(0, 1, length.out = 21),
    include.lowest = TRUE
  )
  counts = tabulate(bins, nbins = 20)
  chisq_p = chisq.test(counts)$p.value
  expect_true(chisq_p > 0.001,
    info = sprintf("SBC global chi-squared p=%.4f", chisq_p)
  )
})


# ---- SBC test: Adaptive Metropolis ------------------------------------------

test_that("SBC: GGM MH produces uniform ranks (p=3, no edge selection)", {
  skip_unless_certification()

  p = 3
  n = 100
  R = 200
  L = 999
  thin = 5
  L_raw = L * thin
  scale = 2.5

  set.seed(2027)

  prior_draws = vector("list", R)
  for(r in seq_len(R)) {
    prior_draws[[r]] = draw_prior_K(p, scale)
  }

  n_off = p * (p - 1) / 2
  n_params = n_off + p
  ranks = matrix(NA_real_, nrow = R, ncol = n_params)

  for(r in seq_len(R)) {
    K_true = prior_draws[[r]]
    Sigma = solve(K_true)
    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = bgm(dat,
      variable_type = "continuous",
      iter = L_raw, warmup = 5000, chains = 1,
      edge_selection = FALSE, update_method = "adaptive-metropolis",
      interaction_prior = cauchy_prior(scale = scale), delta = 0,
      precision_scale_prior = gamma_prior(shape = 1, rate = 1),
      display_progress = "none", seed = 2027L + r
    )

    # Thin to reduce autocorrelation in MH chains
    thin_idx = seq(1, L_raw, by = thin)
    ranks[r, ] = compute_sbc_ranks(K_true, p, fit, thin_idx = thin_idx)
    if(r == 1) colnames(ranks) = names(ranks[1, ])
  }

  n_fail_ks = 0
  for(j in seq_len(ncol(ranks))) {
    u = ranks[, j] / (L + 1)
    p_val = suppressWarnings(ks.test(u, "punif")$p.value)
    if(p_val <= 0.01) n_fail_ks = n_fail_ks + 1
  }

  max_fail = max(1, ceiling(n_params * 0.01 * 2))
  expect_true(n_fail_ks <= max_fail,
    info = sprintf(
      "SBC KS: %d/%d parameters failed (limit %d)",
      n_fail_ks, n_params, max_fail
    )
  )

  all_ranks = as.vector(ranks)
  bins = cut(all_ranks / (L + 1),
    breaks = seq(0, 1, length.out = 21),
    include.lowest = TRUE
  )
  counts = tabulate(bins, nbins = 20)
  chisq_p = chisq.test(counts)$p.value
  expect_true(chisq_p > 0.001,
    info = sprintf("SBC global chi-squared (MH) p=%.4f", chisq_p)
  )
})


# ---- Prior sampler with edge selection ---------------------------------------

# Draw a precision matrix K from the spike-and-slab GGM prior, where the
# prior is specified on the partial-association scale Omega = -K/2:
#   gamma_ij ~ Bernoulli(0.5)
#   Omega_ij | gamma_ij=1 ~ Cauchy(0, scale); Omega_ij | gamma_ij=0 = 0
#   -Omega_ii ~ Gamma(1, 1)
# Implies on K:
#   K_ij | gamma_ij=1 ~ Cauchy(0, 2*scale); K_ij | gamma_ij=0 = 0
#   K_ii ~ 2 * Gamma(1, 1)
# Rejection-samples until K is positive definite.
draw_prior_K_es = function(p, scale = 2.5, inclusion_prob = 0.5,
                           max_tries = 50000) {
  for(attempt in seq_len(max_tries)) {
    K = matrix(0, p, p)
    gamma = matrix(0L, p, p)

    for(i in seq_len(p - 1)) {
      for(j in (i + 1):p) {
        if(runif(1) < inclusion_prob) {
          gamma[i, j] = 1L
          gamma[j, i] = 1L
          omega_ij = rcauchy(1, 0, scale)
          K[i, j] = -2 * omega_ij
          K[j, i] = K[i, j]
        }
      }
    }

    for(i in seq_len(p)) {
      K[i, i] = 2 * rgamma(1, shape = 1, rate = 1)
    }

    ev = eigen(K, symmetric = TRUE, only.values = TRUE)$values
    if(all(ev > 1e-8)) {
      return(list(K = K, gamma = gamma))
    }
  }
  stop("draw_prior_K_es: failed to draw PD matrix in ", max_tries, " attempts")
}


# ---- Diagonal rank computation for edge-selection SBC ------------------------

# Compute SBC ranks for diagonal K_ii only (avoids tied-rank issues
# from the spike-and-slab on off-diagonals).
compute_sbc_ranks_diag = function(K_true, p, fit, thin_idx = NULL) {
  main_samples = do.call(rbind, fit$raw_samples$main)
  if(!is.null(thin_idx)) main_samples = main_samples[thin_idx, , drop = FALSE]

  ranks = numeric(p)
  for(i in seq_len(p)) {
    ranks[i] = sum(main_samples[, i] < K_true[i, i])
  }
  names(ranks) = paste0("K_", seq_len(p), seq_len(p))
  ranks
}


# ---- SBC test: Edge selection (MH, diagonal elements) ------------------------

test_that("SBC: GGM MH produces uniform diagonal ranks (p=3, edge selection)", {
  skip_unless_certification()

  p = 3
  n = 100
  R = 200
  L = 999
  thin = 5
  L_raw = L * thin
  scale = 2.5

  set.seed(2028)

  prior_draws = vector("list", R)
  for(r in seq_len(R)) {
    prior_draws[[r]] = draw_prior_K_es(p, scale)
  }

  n_params = p
  ranks = matrix(NA_real_, nrow = R, ncol = n_params)

  for(r in seq_len(R)) {
    draw = prior_draws[[r]]
    K_true = draw$K
    Sigma = solve(K_true)
    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = bgm(dat,
      variable_type = "continuous",
      iter = L_raw, warmup = 5000, chains = 1,
      edge_selection = TRUE, update_method = "adaptive-metropolis",
      interaction_prior = cauchy_prior(scale = scale), delta = 0,
      precision_scale_prior = gamma_prior(shape = 1, rate = 1),
      display_progress = "none", seed = 2028L + r
    )

    # Thin to reduce autocorrelation
    thin_idx = seq(1, L_raw, by = thin)
    ranks[r, ] = compute_sbc_ranks_diag(K_true, p, fit, thin_idx = thin_idx)
    if(r == 1) colnames(ranks) = names(ranks[1, ])
  }

  n_fail_ks = 0
  for(j in seq_len(ncol(ranks))) {
    u = ranks[, j] / (L + 1)
    p_val = suppressWarnings(ks.test(u, "punif")$p.value)
    if(p_val <= 0.01) n_fail_ks = n_fail_ks + 1
  }

  max_fail = max(1, ceiling(n_params * 0.01 * 2))
  expect_true(n_fail_ks <= max_fail,
    info = sprintf(
      "SBC KS (edge selection, diag): %d/%d parameters failed (limit %d)",
      n_fail_ks, n_params, max_fail
    )
  )

  all_ranks = as.vector(ranks)
  bins = cut(all_ranks / (L + 1),
    breaks = seq(0, 1, length.out = 21),
    include.lowest = TRUE
  )
  counts = tabulate(bins, nbins = 20)
  chisq_p = chisq.test(counts)$p.value
  expect_true(chisq_p > 0.001,
    info = sprintf("SBC global chi-squared p=%.4f (edge selection, diag)", chisq_p)
  )
})


# ---- Tilted prior sampler (delta > 0) ----------------------------------------

# Draw K from the determinant-tilted GGM prior:
#   p(K) propto |K|^delta * slab(K_ij) * diag(K_ii) * 1{K in M+}
# using the rejection scheme from the spikeslab manuscript sec:tilt-empirics.
# Proposal shifts the diagonal Gamma shape up by delta on the K_ii scale;
# acceptance probability r = (|K| / prod(K_ii))^delta lies in [0, 1] on PD
# support by Hadamard's inequality. Off-diagonals are unchanged from the
# untilted Cauchy-slab proposal.
draw_prior_K_tilted = function(p, scale = 2.5, delta = 1,
                               max_tries = 100000) {
  for(attempt in seq_len(max_tries)) {
    K = matrix(0, p, p)

    # Off-diagonal: Cauchy slab, K_ij = -2 * Omega_ij as in untilted sampler.
    for(i in seq_len(p - 1)) {
      for(j in (i + 1):p) {
        omega_ij = rcauchy(1, 0, scale)
        K[i, j] = -2 * omega_ij
        K[j, i] = K[i, j]
      }
    }

    # Diagonal: shifted-shape Gamma proposal. Untilted draws K_ii ~
    # 2 * Gamma(1, 1) = Gamma(1, 1/2) on K_ii; tilted proposal draws K_ii ~
    # 2 * Gamma(1 + delta, 1) = Gamma(1 + delta, 1/2). The factor K_ii^delta
    # from the bound is absorbed into the proposal density.
    for(i in seq_len(p)) {
      K[i, i] = 2 * rgamma(1, shape = 1 + delta, rate = 1)
    }

    ev = eigen(K, symmetric = TRUE, only.values = TRUE)$values
    if(!all(ev > 1e-8)) next

    # Tilt acceptance: log r = delta * (sum log eigenvalues - sum log diag).
    log_r = delta * (sum(log(ev)) - sum(log(diag(K))))
    if(log(runif(1)) < log_r) {
      return(K)
    }
  }
  stop("draw_prior_K_tilted: failed in ", max_tries, " attempts")
}


# ---- SBC test: NUTS under determinant tilt ----------------------------------

test_that("SBC: GGM NUTS produces uniform ranks under tilt (p=3, delta=1)", {
  skip_unless_certification()

  p = 3
  n = 100
  R = 200
  L = 999
  scale = 2.5
  delta = 1

  set.seed(2029)

  prior_draws = vector("list", R)
  for(r in seq_len(R)) {
    prior_draws[[r]] = draw_prior_K_tilted(p, scale, delta)
  }

  n_off = p * (p - 1) / 2
  n_params = n_off + p
  ranks = matrix(NA_real_, nrow = R, ncol = n_params)

  for(r in seq_len(R)) {
    K_true = prior_draws[[r]]
    Sigma = solve(K_true)
    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = bgm(dat,
      variable_type = "continuous",
      iter = L, warmup = 1000, chains = 1,
      edge_selection = FALSE, update_method = "nuts",
      interaction_prior = cauchy_prior(scale = scale), delta = delta,
      precision_scale_prior = gamma_prior(shape = 1, rate = 1),
      display_progress = "none", seed = 2029L + r
    )

    ranks[r, ] = compute_sbc_ranks(K_true, p, fit)
    if(r == 1) colnames(ranks) = names(ranks[1, ])
  }

  n_fail_ks = 0
  for(j in seq_len(ncol(ranks))) {
    u = ranks[, j] / (L + 1)
    p_val = suppressWarnings(ks.test(u, "punif")$p.value)
    if(p_val <= 0.01) n_fail_ks = n_fail_ks + 1
  }

  max_fail = max(1, ceiling(n_params * 0.01 * 2))
  expect_true(n_fail_ks <= max_fail,
    info = sprintf(
      "SBC KS (delta=1): %d/%d parameters failed (limit %d)",
      n_fail_ks, n_params, max_fail
    )
  )

  all_ranks = as.vector(ranks)
  bins = cut(all_ranks / (L + 1),
    breaks = seq(0, 1, length.out = 21),
    include.lowest = TRUE
  )
  counts = tabulate(bins, nbins = 20)
  chisq_p = chisq.test(counts)$p.value
  expect_true(chisq_p > 0.001,
    info = sprintf("SBC global chi-squared (delta=1) p=%.4f", chisq_p)
  )
})


# ---- SBC test: Joint specification (edge selection on) ----------------------

# Reconstruct a p x p symmetric K matrix from the (off-diag, diag) row pair
# returned by sample_ggm_prior(spec = "joint"). Column order in off matches
# the (i < j) upper-triangle visit order.
reconstruct_K = function(off_row, diag_row, p) {
  K = matrix(0, p, p)
  idx = 1L
  for(i in seq_len(p - 1)) {
    for(j in (i + 1):p) {
      K[i, j] = off_row[idx]
      K[j, i] = off_row[idx]
      idx = idx + 1L
    }
  }
  diag(K) = diag_row
  K
}

test_that("SBC: GGM joint-spec produces uniform ranks (p=5, edge selection)", {
  skip_unless_certification()

  # Joint-specification SBC: draw (K_true, Gamma_true) from the un-normalised
  # joint prior via sample_ggm_prior(spec = "joint"), simulate Y | K_true, and
  # check that bgm()'s default-MH posterior reproduces the prior in rank
  # histograms. Generator and fit both use the auto-default tilt (delta =
  # 0.5 * log(p)), and both target the joint specification, so the SBC ranks
  # should be uniform.
  #
  # Test functions are the continuous K diagonals and log|K|. Off-diagonal
  # K_ij entries are skipped here because the spike-and-slab point mass at
  # zero induces tied ranks that break standard uniformity tests (the same
  # rationale as the existing edge-selection SBC test, which only checks the
  # diagonal).
  p = 5
  n = 100
  R = 300
  L = 999
  thin = 10

  set.seed(2030)

  # Step 1: a single thinned joint-prior chain produces R quasi-independent
  # (K_true, Gamma_true) draws. Auto-default delta on the prior generator.
  joint = sample_ggm_prior(
    p = p, n_samples = as.integer(R * thin), n_warmup = 2000L,
    spec = "joint", delta = NULL, seed = 2030L, verbose = FALSE
  )
  thin_idx = seq(thin, R * thin, by = thin)
  K_off_true = joint$K_offdiag[thin_idx, , drop = FALSE]
  K_diag_true = joint$K_diag[thin_idx, , drop = FALSE]

  n_params = p + 1L # K_11, ..., K_pp, log|K|
  ranks = matrix(NA_real_, nrow = R, ncol = n_params)

  for(r in seq_len(R)) {
    K_true = reconstruct_K(K_off_true[r, ], K_diag_true[r, ], p)
    Sigma = solve(K_true)
    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = bgm(dat,
      variable_type = "continuous",
      iter = L, warmup = 1000, chains = 1,
      edge_selection = TRUE, update_method = "adaptive-metropolis",
      delta = NULL, # auto-default, matches generator
      display_progress = "none", seed = 2030L + r
    )

    main_samples = do.call(rbind, fit$raw_samples$main)
    pw_samples = do.call(rbind, fit$raw_samples$pairwise)

    # K diagonal ranks
    for(i in seq_len(p)) {
      ranks[r, i] = sum(main_samples[, i] < K_true[i, i])
    }

    # log|K| rank: rebuild K_post per iteration, compute log|K|.
    log_det_K_true = determinant(K_true, logarithm = TRUE)$modulus
    log_det_K_post = vapply(seq_len(nrow(main_samples)), function(it) {
      K_it = reconstruct_K(pw_samples[it, ], main_samples[it, ], p)
      determinant(K_it, logarithm = TRUE)$modulus
    }, numeric(1))
    ranks[r, p + 1L] = sum(log_det_K_post < log_det_K_true)
  }

  # Per-parameter KS test, allowing 1 false positive across (p + 1)
  # parameters at alpha = 0.01.
  n_fail_ks = 0
  for(j in seq_len(ncol(ranks))) {
    u = ranks[, j] / (L + 1)
    p_val = suppressWarnings(ks.test(u, "punif")$p.value)
    if(p_val <= 0.01) n_fail_ks = n_fail_ks + 1
  }
  max_fail = max(1, ceiling(n_params * 0.01 * 2))
  expect_true(n_fail_ks <= max_fail,
    info = sprintf(
      "SBC KS (joint): %d/%d parameters failed (limit %d)",
      n_fail_ks, n_params, max_fail
    )
  )

  # Global chi-squared on aggregated ranks.
  all_ranks = as.vector(ranks)
  bins = cut(all_ranks / (L + 1),
    breaks = seq(0, 1, length.out = 21),
    include.lowest = TRUE
  )
  counts = tabulate(bins, nbins = 20)
  chisq_p = chisq.test(counts)$p.value
  expect_true(chisq_p > 0.001,
    info = sprintf("SBC global chi-squared (joint) p=%.4f", chisq_p)
  )
})


# ---- SBC test: gamma-shape diagonal (Gibbs, no edge selection) ----------------

test_that("SBC: GGM Gibbs produces uniform ranks at a gamma-shape diagonal", {
  skip_unless_certification()

  # The row-block Gibbs corrects its shape-1 conjugate row proposal by an
  # independence-Metropolis accept at shape != 1; uniform ranks under a
  # Gamma(2, 1) diagonal certify the corrected sweep against the prior
  # predictive.
  p = 3
  n = 100
  R = 200
  L = 999
  thin = 5
  L_raw = L * thin
  scale = 2.5
  shape = 2

  set.seed(2029)

  prior_draws = vector("list", R)
  for(r in seq_len(R)) {
    prior_draws[[r]] = draw_prior_K(p, scale, shape = shape)
  }

  n_off = p * (p - 1) / 2
  n_params = n_off + p
  ranks = matrix(NA_real_, nrow = R, ncol = n_params)

  for(r in seq_len(R)) {
    K_true = prior_draws[[r]]
    Sigma = solve(K_true)
    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = bgm(dat,
      variable_type = "continuous",
      iter = L_raw, warmup = 5000, chains = 1,
      edge_selection = FALSE, update_method = "gibbs",
      interaction_prior = cauchy_prior(scale = scale), delta = 0,
      precision_scale_prior = gamma_prior(shape = shape, rate = 1),
      display_progress = "none", seed = 2029L + r
    )

    thin_idx = seq(1, L_raw, by = thin)
    ranks[r, ] = compute_sbc_ranks(K_true, p, fit, thin_idx = thin_idx)
    if(r == 1) colnames(ranks) = names(ranks[1, ])
  }

  n_fail_ks = 0
  for(j in seq_len(ncol(ranks))) {
    u = ranks[, j] / (L + 1)
    p_val = suppressWarnings(ks.test(u, "punif")$p.value)
    if(p_val <= 0.01) n_fail_ks = n_fail_ks + 1
  }

  max_fail = max(1, ceiling(n_params * 0.01 * 2))
  expect_true(n_fail_ks <= max_fail,
    info = sprintf(
      "SBC KS (gibbs, shape 2): %d/%d parameters failed (limit %d)",
      n_fail_ks, n_params, max_fail
    )
  )

  all_ranks = as.vector(ranks)
  bins = cut(all_ranks / (L + 1),
    breaks = seq(0, 1, length.out = 21),
    include.lowest = TRUE
  )
  counts = tabulate(bins, nbins = 20)
  chisq_p = chisq.test(counts)$p.value
  expect_true(chisq_p > 0.001,
    info = sprintf("SBC global chi-squared (gibbs, shape 2) p=%.4f", chisq_p)
  )
})
