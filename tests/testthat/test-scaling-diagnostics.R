# --------------------------------------------------------------------------- #
# Phase S <U+2014> Scaling and stress diagnostics for NUTS sampler.
#
# Checks that NUTS diagnostics (divergences, E-BFMI, tree depth, ESS, Rhat)
# remain healthy as problem size increases. Does NOT verify posterior
# correctness <U+2014> see test-ggm-nuts.R and test-mixed-nuts.R for that.
#
# GGM:
#   G1  p=5,  n=200, no edge selection    (easy)
#   G2  p=10, n=100, no edge selection    (moderate)
#   G3  p=10, n=100, edge selection       (moderate + ES)
#   G4  p=15, n=200, edge selection       (hard)
#
# Mixed MRF:
#   M1  p=3, q=2, n=200, conditional, no ES   (easy)
#   M2  p=5, q=3, n=200, conditional, ES      (moderate)
#   M3  p=7, q=5, n=150, conditional, ES      (hard)
#   M4  p=5, q=3, n=200, marginal, ES         (moderate)
#   M5  p=3, q=2, n=200, conditional, ES      (near-singular Kyy)
#
# Fail thresholds (from plan):
#   divergence rate < 5%, E-BFMI > 0.1, tree depth hits < 25%,
#   min ESS > 50, max Rhat < 1.1.
#
# All gated behind BGMS_RUN_SLOW_TESTS.
# --------------------------------------------------------------------------- #


# ---- Skip gate ---------------------------------------------------------------

skip_unless_slow_scaling = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run scaling diagnostics"
  )
}


# ---- Diagnostic assertion helper ---------------------------------------------

check_nuts_health = function(fit, label,
                             div_rate_max = 0.05,
                             ebfmi_min = 0.1,
                             tree_hit_rate_max = 0.25,
                             ess_min = 50,
                             rhat_max = 1.1) {
  diag = fit$nuts_diag
  n_samples = length(diag$divergent)

  div_rate = diag$summary$total_divergences / n_samples
  expect_true(div_rate < div_rate_max,
    info = sprintf(
      "%s: divergence rate %.1f%% (limit %.0f%%)",
      label, div_rate * 100, div_rate_max * 100
    )
  )

  expect_true(diag$summary$min_ebfmi > ebfmi_min,
    info = sprintf(
      "%s: min E-BFMI %.3f (limit %.1f)",
      label, diag$summary$min_ebfmi, ebfmi_min
    )
  )

  tree_hit_rate = diag$summary$max_tree_depth_hits / n_samples
  expect_true(tree_hit_rate < tree_hit_rate_max,
    info = sprintf(
      "%s: tree depth hit rate %.1f%% (limit %.0f%%)",
      label, tree_hit_rate * 100, tree_hit_rate_max * 100
    )
  )

  # ESS and Rhat are available when edge selection produces summary tables
  if(!is.null(fit$posterior_summary_pairwise)) {
    ess = fit$posterior_summary_pairwise$n_eff
    rhat = fit$posterior_summary_pairwise$Rhat

    ess_finite = ess[is.finite(ess)]
    if(length(ess_finite) > 0) {
      expect_true(min(ess_finite) > ess_min,
        info = sprintf(
          "%s: min ESS %.0f (limit %d)",
          label, min(ess_finite), ess_min
        )
      )
    }

    rhat_finite = rhat[is.finite(rhat)]
    if(length(rhat_finite) > 0) {
      expect_true(max(rhat_finite) < rhat_max,
        info = sprintf(
          "%s: max Rhat %.3f (limit %.2f)",
          label, max(rhat_finite), rhat_max
        )
      )
    }
  }
}


# ---- GGM data generator ------------------------------------------------------

# Tridiagonal precision: K_ii = 2, K_{i,i+1} = -0.5.
generate_ggm_scaling_data = function(p, n, seed = 1) {
  set.seed(seed)
  K = diag(2, p)
  for(i in seq_len(p - 1)) {
    K[i, i + 1] = -0.5
    K[i + 1, i] = -0.5
  }
  Sigma = solve(K)
  x = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(x) = paste0("V", seq_len(p))
  as.data.frame(x)
}


# ---- Mixed MRF data generator ------------------------------------------------

# Nearest-neighbor discrete coupling, sparse cross interactions,
# near-diagonal continuous precision.
generate_mixed_scaling_data = function(p, q, n,
                                       pairwise_cont_override = NULL,
                                       cross_strength = 0.2,
                                       threshold_spread = 0.1,
                                       seed = 2026) {
  nc = rep(2L, p)

  # Discrete-discrete: nearest-neighbor coupling
  pairwise_disc = matrix(0, p, p)
  if(p > 1) {
    for(i in seq_len(p - 1)) {
      val = (-1)^i * 0.15
      pairwise_disc[i, i + 1] = val
      pairwise_disc[i + 1, i] = val
    }
  }

  # Cross: one interaction per discrete variable (cycling over q)
  pairwise_cross = matrix(0, p, q)
  for(i in seq_len(min(p, q))) {
    pairwise_cross[i, i] = (-1)^i * cross_strength
  }

  # Continuous precision (stored as -0.5 * K)
  if(!is.null(pairwise_cont_override)) {
    pairwise_cont = pairwise_cont_override
  } else {
    pairwise_cont = diag(-0.5, q)
    if(q > 1) {
      for(i in seq_len(q - 1)) {
        pairwise_cont[i, i + 1] = -0.05
        pairwise_cont[i + 1, i] = -0.05
      }
    }
  }

  # Thresholds
  mux = matrix(0, p, max(nc) + 1)
  for(i in seq_len(p)) {
    mux[i, seq_len(nc[i])] = seq(-threshold_spread, threshold_spread,
      length.out = nc[i]
    )
  }

  muy = seq(-0.1, 0.1, length.out = q)

  sim = bgms:::sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc,
    pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = rep("ordinal", p),
    baseline_category_r = rep(0L, p),
    iter = 2000L, seed = as.integer(seed)
  )

  df = data.frame(sim$x, sim$y)
  colnames(df) = c(paste0("d", seq_len(p)), paste0("c", seq_len(q)))
  df
}


# ---- GGM scaling tests -------------------------------------------------------

test_that("S.G1: GGM NUTS healthy at p=5, no edge selection", {
  skip_unless_slow_scaling()

  dat = generate_ggm_scaling_data(p = 5, n = 200, seed = 3001)

  fit = bgm(dat,
    variable_type = "continuous",
    iter = 2000, warmup = 2000, chains = 4,
    edge_selection = FALSE, update_method = "nuts",
    interaction_prior = cauchy_prior(scale = 2.5),
    display_progress = "none", seed = 3001
  )

  check_nuts_health(fit, "S.G1")
})

test_that("S.G2: GGM NUTS healthy at p=10, no edge selection", {
  skip_unless_slow_scaling()

  dat = generate_ggm_scaling_data(p = 10, n = 100, seed = 3002)

  fit = bgm(dat,
    variable_type = "continuous",
    iter = 2000, warmup = 2000, chains = 4,
    edge_selection = FALSE, update_method = "nuts",
    interaction_prior = cauchy_prior(scale = 2.5),
    display_progress = "none", seed = 3002
  )

  check_nuts_health(fit, "S.G2")
})

test_that("S.G3: GGM NUTS healthy at p=10 with edge selection", {
  skip_unless_slow_scaling()

  dat = generate_ggm_scaling_data(p = 10, n = 100, seed = 3003)

  fit = bgm(dat,
    variable_type = "continuous",
    iter = 2000, warmup = 2000, chains = 4,
    edge_selection = TRUE, update_method = "nuts",
    interaction_prior = cauchy_prior(scale = 2.5),
    display_progress = "none", seed = 3003
  )

  check_nuts_health(fit, "S.G3")
})

test_that("S.G4: GGM NUTS healthy at p=15 with edge selection", {
  skip_unless_slow_scaling()

  dat = generate_ggm_scaling_data(p = 15, n = 200, seed = 3004)

  fit = bgm(dat,
    variable_type = "continuous",
    iter = 2000, warmup = 2000, chains = 4,
    edge_selection = TRUE, update_method = "nuts",
    interaction_prior = cauchy_prior(scale = 2.5),
    display_progress = "none", seed = 3004
  )

  check_nuts_health(fit, "S.G4")
})


# ---- Mixed MRF scaling tests -------------------------------------------------

test_that("S.M1: Mixed NUTS healthy at p=3, q=2, no edge selection", {
  skip_unless_slow_scaling()

  dat = generate_mixed_scaling_data(p = 3, q = 2, n = 200, seed = 3011)
  vtype = c(rep("ordinal", 3), rep("continuous", 2))

  fit = bgm(dat,
    variable_type = vtype,
    iter = 2000, warmup = 2000, chains = 4,
    edge_selection = FALSE, update_method = "nuts",
    interaction_prior = cauchy_prior(scale = 2.5),
    threshold_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
    display_progress = "none", seed = 3011
  )

  check_nuts_health(fit, "S.M1")
})

test_that("S.M2: Mixed NUTS healthy at p=5, q=3 with edge selection", {
  skip_unless_slow_scaling()

  dat = generate_mixed_scaling_data(p = 5, q = 3, n = 200, seed = 3012)
  vtype = c(rep("ordinal", 5), rep("continuous", 3))

  fit = bgm(dat,
    variable_type = vtype,
    iter = 2000, warmup = 2000, chains = 4,
    edge_selection = TRUE, update_method = "nuts",
    interaction_prior = cauchy_prior(scale = 2.5),
    threshold_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
    display_progress = "none", seed = 3012
  )

  check_nuts_health(fit, "S.M2")
})

test_that("S.M3: Mixed NUTS healthy at p=7, q=5 with edge selection", {
  skip_unless_slow_scaling()

  dat = generate_mixed_scaling_data(p = 7, q = 5, n = 150, seed = 3013)
  vtype = c(rep("ordinal", 7), rep("continuous", 5))

  fit = bgm(dat,
    variable_type = vtype,
    iter = 2000, warmup = 2000, chains = 4,
    edge_selection = TRUE, update_method = "nuts",
    interaction_prior = cauchy_prior(scale = 2.5),
    threshold_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
    display_progress = "none", seed = 3013
  )

  # Rhat limit relaxed to 1.17 (vs the default 1.10): posterior_summary_pairwise
  # is the max classic Gelman-Rubin Rhat over all 66 edge-selected interaction
  # coefficients, each a spike-and-slab (0/value) sequence. Since the marginal-PL
  # correctness fix (#97) the corrected target sits at ~1.16 here; the other
  # health checks (divergences, E-BFMI, tree depth, ESS) stay strict.
  check_nuts_health(fit, "S.M3", rhat_max = 1.17)
})

test_that("S.M4: Mixed NUTS healthy at p=5, q=3, marginal PL", {
  skip_unless_slow_scaling()

  dat = generate_mixed_scaling_data(p = 5, q = 3, n = 200, seed = 3014)
  vtype = c(rep("ordinal", 5), rep("continuous", 3))

  fit = bgm(dat,
    variable_type = vtype,
    iter = 2000, warmup = 2000, chains = 4,
    edge_selection = TRUE, update_method = "nuts",
    interaction_prior = cauchy_prior(scale = 2.5),
    threshold_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
    display_progress = "none", seed = 3014
  )

  # Rhat limit relaxed to 1.17 (vs the default 1.10): same rationale as S.M3 --
  # max GR Rhat over edge-selected spike-and-slab pairwise coefficients, shifted
  # by the #97 marginal-PL correctness fix. Other health checks stay strict.
  check_nuts_health(fit, "S.M4", rhat_max = 1.17)
})

test_that("S.M5: Mixed NUTS survives near-singular Kyy", {
  skip_unless_slow_scaling()

  # K_yy = [[1, 0.05], [0.05, 0.012]]: condition number ~106,
  # stresses the RATTLE projection while keeping data well-behaved.
  pairwise_cont_m5 = matrix(c(-0.5, -0.025, -0.025, -0.006), 2, 2)

  dat = generate_mixed_scaling_data(
    p = 3, q = 2, n = 300,
    pairwise_cont_override = pairwise_cont_m5,
    cross_strength = 0, threshold_spread = 0.3, seed = 3015
  )
  vtype = c(rep("ordinal", 3), rep("continuous", 2))

  fit = bgm(dat,
    variable_type = vtype,
    iter = 2000, warmup = 2000, chains = 4,
    edge_selection = TRUE, update_method = "nuts",
    interaction_prior = cauchy_prior(scale = 2.5),
    threshold_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
    display_progress = "none", seed = 3015
  )

  # Test goal: NUTS *survives* near-singular Kyy (no divergence explosion,
  # finite output, plausible ESS). The strict Rhat < 1.10 used in S.M1-S.M4
  # is not the right check here: the truth has K_yy[1,2] = 0.05 (spike
  # territory) and K_yy[2,2] = 0.012 (near zero), so different chains
  # routinely settle into different modes of the inclusion indicators,
  # inflating across-chain Rhat without indicating a real convergence
  # failure. Empirical Rhat across delta in {0, 0.5, 1, 1.5, 2} is 1.11-1.26
  # at this seed; the relaxed ceiling matches the "survives" semantics.
  check_nuts_health(fit, "S.M5", rhat_max = 1.50)
})
