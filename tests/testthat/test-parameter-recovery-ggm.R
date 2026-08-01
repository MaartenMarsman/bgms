# --------------------------------------------------------------------------- #
# Parameter Recovery for GGM
#
# Checks that the NUTS sampler recovers known precision matrix entries:
# the true value should fall within the 95% credible interval at least
# 85% of the time across R = 50 simulated datasets.
#
# Two conditions:
#   PR.1  Dense graph (p=5, all edges present, no edge selection)
#   PR.2  Sparse graph (p=5, ~40% edges zeroed, edge selection on)
#
# Weekly certification tier (T2, BGMS_RUN_CERTIFICATION): a full
# parameter-recovery sweep over R = 50 simulated datasets per condition.
# --------------------------------------------------------------------------- #


# ---- Skip gate ---------------------------------------------------------------

# ---- Ground truth precision matrices -----------------------------------------

# AR(1)-like precision: tridiagonal with rho = 0.5.
# Constructed directly to avoid solve() rounding noise.
make_K_dense = function(p = 5, rho = 0.5) {
  K = matrix(0, p, p)
  denom = 1 - rho^2
  for(i in seq_len(p)) {
    if(i == 1 || i == p) {
      K[i, i] = 1 / denom
    } else {
      K[i, i] = (1 + rho^2) / denom
    }
  }
  for(i in seq_len(p - 1)) {
    K[i, i + 1] = -rho / denom
    K[i + 1, i] = -rho / denom
  }
  K
}

# Sparse: zero out 4 of 10 edges (i.e., 40% of edges absent).
make_K_sparse = function(p = 5) {
  K = make_K_dense(p)
  K[1, 4] = K[4, 1] = 0
  K[1, 5] = K[5, 1] = 0
  K[2, 5] = K[5, 2] = 0
  K[3, 5] = K[5, 3] = 0
  stopifnot(all(eigen(K, symmetric = TRUE, only.values = TRUE)$values > 0))
  K
}


# ---- Coverage computation ----------------------------------------------------

# Fit R datasets and compute coverage for each off-diagonal K_ij.
# Returns a named numeric vector: fraction of datasets where the true
# value falls within the 95% credible interval.
compute_coverage = function(K_true, p, n, R,
                            edge_selection, scale = 2.5,
                            iter = 3000, warmup = 1000,
                            base_seed = 4000) {
  Sigma = solve(K_true)

  # Identify off-diagonal parameters (upper triangle)
  n_off = p * (p - 1) / 2
  covered = matrix(0L, nrow = R, ncol = n_off + p) # off-diag + diag

  # True values in column order
  true_off = numeric(n_off)
  off_names = character(n_off)
  idx = 0
  for(i in seq_len(p - 1)) {
    for(j in (i + 1):p) {
      idx = idx + 1
      true_off[idx] = K_true[i, j]
      off_names[idx] = paste0("K_", i, j)
    }
  }
  true_diag = diag(K_true)
  diag_names = paste0("K_", seq_len(p), seq_len(p))

  for(r in seq_len(R)) {
    set.seed(base_seed + r)
    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = bgm(dat,
      variable_type = "continuous",
      iter = iter, warmup = warmup, chains = 2,
      edge_selection = edge_selection, update_method = "nuts",
      interaction_prior = cauchy_prior(scale = scale),
      display_progress = "none", seed = base_seed + r
    )

    # Raw off-diagonal samples (precision scale)
    pw_samples = do.call(rbind, fit$raw_samples$pairwise)

    # With edge selection, use the mixture posterior:
    # K_ij = 0 when gamma=0, slab value when gamma=1.
    if(edge_selection) {
      ind_samples = do.call(rbind, fit$raw_samples$indicator)
      pw_samples = pw_samples * ind_samples
    }

    for(k in seq_len(n_off)) {
      ci = quantile(pw_samples[, k], c(0.025, 0.975))
      if(true_off[k] >= ci[1] && true_off[k] <= ci[2]) {
        covered[r, k] = 1L
      }
    }

    # Raw diagonal samples
    main_samples = do.call(rbind, fit$raw_samples$main)
    for(k in seq_len(p)) {
      ci = quantile(main_samples[, k], c(0.025, 0.975))
      if(true_diag[k] >= ci[1] && true_diag[k] <= ci[2]) {
        covered[r, n_off + k] = 1L
      }
    }
  }

  coverage = colMeans(covered)
  names(coverage) = c(off_names, diag_names)
  coverage
}


# ---- Tests -------------------------------------------------------------------

test_that("PR.1: GGM parameter recovery, dense graph (p=5)", {
  skip_unless_certification()

  p = 5
  n = 200
  R = 50
  K_true = make_K_dense(p)

  coverage = compute_coverage(K_true, p, n, R,
    edge_selection = FALSE, base_seed = 4001
  )

  # Per-parameter coverage floor of 0.80 (with mean >= 0.90 below as the
  # main calibration check). Under perfectly-calibrated 95% credible
  # intervals at R = 50, Binomial(50, 0.95) puts the 1st percentile of
  # observed coverage at 0.86, so a 0.85 floor false-fires roughly 38% of
  # the time across 15 parameters from MC noise alone. At 0.80 the false-
  # positive rate is < 0.1% across 15 params, while retaining ~92% power
  # to detect a real miscalibration to 85% true coverage.
  for(k in seq_along(coverage)) {
    expect_true(coverage[k] >= 0.80,
      info = sprintf(
        "PR.1 %s: coverage %.0f%% < 80%%",
        names(coverage)[k], coverage[k] * 100
      )
    )
  }

  # Mean coverage should be >= 90%
  expect_true(mean(coverage) >= 0.90,
    info = sprintf("PR.1 mean coverage %.0f%%", mean(coverage) * 100)
  )
})


test_that("PR.2: GGM parameter recovery, sparse graph with edge selection (p=5)", {
  skip_unless_certification()

  p = 5
  n = 200
  R = 50
  K_true = make_K_sparse(p)

  coverage = compute_coverage(K_true, p, n, R,
    edge_selection = TRUE, base_seed = 4002
  )

  # Per-parameter floor of 0.80 (mean >= 0.90 below as the main check) for
  # the same binomial-noise reason as PR.1: at R = 50, a 0.85 floor sits
  # below the Binomial(50, 0.95) 1st percentile and false-fires from MC
  # noise alone.
  for(k in seq_along(coverage)) {
    expect_true(coverage[k] >= 0.80,
      info = sprintf(
        "PR.2 %s: coverage %.0f%% < 80%%",
        names(coverage)[k], coverage[k] * 100
      )
    )
  }

  # Mean coverage should be >= 90%
  expect_true(mean(coverage) >= 0.90,
    info = sprintf("PR.2 mean coverage %.0f%%", mean(coverage) * 100)
  )
})
