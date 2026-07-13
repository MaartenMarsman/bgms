# Trust-gauge reference route (block_reference_logR): the block-local exact
# Monte-Carlo reference log R_e that the in-chain gauge compares against the
# deployed log J. These drive the engine reference in isolation.

test_that("the block reference is finite and lands near the deployed ratio", {
  skip_on_cran()
  p = 10
  delta = 0.5 * log(p)
  eta = 2
  zc = bgms:::zratio_constants(delta, eta = eta)
  set.seed(3)
  G = matrix(0L, p, p)
  ut = upper.tri(G)
  G[ut] = rbinom(sum(ut), 1, 0.6)
  G = G + t(G)
  diag(G) = 1L

  # First edge whose mediating block is non-trivial (m >= 2).
  edge = NULL
  for(i in seq_len(p - 1)) {
    for(j in (i + 1):p) {
      r = bgms:::zratio_test_reference(
        G, i, j, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
        delta, eta, 120L, 30L, 7L, FALSE
      )
      if(isTRUE(r$valid) && r$m >= 2) {
        edge = c(i, j)
        break
      }
    }
    if(!is.null(edge)) break
  }
  expect_false(is.null(edge))

  r = bgms:::zratio_test_reference(
    G, edge[1], edge[2], zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    delta, eta, 240L, 30L, 7L, FALSE
  )
  expect_true(r$ok)
  expect_true(is.finite(r$logR))
  expect_true(is.finite(r$mcse) && r$mcse > 0)
  # The deployed corrected ratio is a validated approximation of the exact
  # reference on an entangled edge: they agree to within a few MCSE.
  expect_lt(abs(r$logR - r$log_zratio), 0.05)
})

test_that("the reference MCSE shrinks with the draw count", {
  skip_on_cran()
  p = 10
  delta = 0.5 * log(p)
  eta = 2
  zc = bgms:::zratio_constants(delta, eta = eta)
  set.seed(3)
  G = matrix(0L, p, p)
  ut = upper.tri(G)
  G[ut] = rbinom(sum(ut), 1, 0.6)
  G = G + t(G)
  diag(G) = 1L

  # A high-degree endpoint pair has an entangled block.
  small = bgms:::zratio_test_reference(
    G, 9, 10, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    delta, eta, 60L, 30L, 11L, FALSE
  )
  large = bgms:::zratio_test_reference(
    G, 9, 10, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    delta, eta, 960L, 30L, 11L, FALSE
  )
  expect_true(small$valid && large$valid)
  # Two estimates from independent streams agree within combined MCSE.
  expect_lt(
    abs(small$logR - large$logR),
    5 * sqrt(small$mcse^2 + large$mcse^2)
  )
  expect_lt(large$mcse, small$mcse)
})

# Harm channel: standardized se (projected inclusion-probability distortion)
# on top of the flip-rate channel.

test_that("the harm channel computes the documented statistics", {
  gauge = list(
    flip_rate = 0.001, noise_floor = 1e-4, se_mean = 0.05, se_sd = 0.02,
    se_mcse = 0.002, n_ent = 100, n_ref = 25, n_capped = 0
  )
  chains = list(list(zratio = list(gauge = gauge)))

  # Beta-Bernoulli, exchangeable inclusion probabilities: the linearized
  # feedback gain is E * m / (theta (1 - theta) (a + b + E)).
  a = 9
  b = 1
  E = 190
  th = 0.9
  pip = rep(th, E)
  s = summarize_zratio_gauge(
    chains,
    verbose = FALSE,
    harm_inputs = list(pip = list(pip), a = a, b = b)
  )
  pc = s$per_chain
  expect_equal(pc$se_se, sqrt(0.002^2 + 0.02^2 / 25))
  m_bar = th * (1 - th)
  theta_hat = (a + E * th) / (a + b + E)
  gain = E * m_bar / (theta_hat * (1 - theta_hat) * (a + b + E))
  expect_equal(pc$amplification, 1 / (1 - min(gain, 0.98)))
  expect_equal(pc$harm_pred, abs(0.05) * m_bar * pc$amplification)
  # 0.05 nats, sensitivity 0.09, ~18x amplification: far over tolerance,
  # and resolved (|se_mean| = 0.05 > 2 * se_se).
  expect_true(pc$harm_flag)
  expect_true(s$flagged)

  # Fixed inclusion probability: no feedback, amplification 1.
  s_bern = summarize_zratio_gauge(
    chains,
    verbose = FALSE,
    harm_inputs = list(pip = list(pip), a = NULL, b = NULL)
  )
  expect_equal(s_bern$per_chain$amplification, 1)
  expect_equal(s_bern$per_chain$harm_pred, 0.05 * m_bar)

  # Unresolved se (large reference noise) must not flag regardless of size.
  gauge_noisy = gauge
  gauge_noisy$se_mcse = 0.2
  s_noisy = summarize_zratio_gauge(
    list(list(zratio = list(gauge = gauge_noisy))),
    verbose = FALSE,
    harm_inputs = list(pip = list(pip), a = a, b = b)
  )
  expect_false(s_noisy$per_chain$harm_flag)

  # No harm inputs: channel reports NA and never flags.
  s_off = summarize_zratio_gauge(chains, verbose = FALSE)
  expect_true(is.na(s_off$per_chain$harm_pred))
  expect_false(s_off$per_chain$harm_flag)

  # Unsupported edge priors disable the channel at the input builder.
  expect_null(zratio_harm_inputs(list(pip), "Stochastic-Block", a = 1, b = 1))
})

test_that("a known-biased evidence-free fit fires the harm channel", {
  skip_on_cran()
  # Bare additive kernel (calibration_window = 0) under a dense-leaning
  # Beta-Bernoulli prior with no data: the feedback-amplified regime where
  # the flip rate stays quiet but the projected distortion is first-order.
  f = sample_ggm_prior(
    p = 16L, n_samples = 1200L, n_warmup = 500L,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 6),
    spec = "hierarchical", edge_prior = beta_bernoulli_prior(9, 1),
    update_method = "gibbs", calibration_window = 0L,
    zratio_diagnostics = TRUE, seed = 7L, verbose = FALSE
  )
  pc = f$zratio_diagnostics$per_chain
  expect_true(is.finite(pc$se_mcse) && pc$se_mcse > 0)
  expect_true(is.finite(pc$harm_pred))
  expect_gt(pc$amplification, 5)
  expect_true(pc$harm_flag)
})

test_that("the harm channel weights errors by per-edge sensitivity", {
  # Per-pair audit stream present: harm_pred = |mean(m_e s_e)| * A with a
  # cluster-robust (per audited edge) resolution gate. Edge (0,1) is audited
  # twice to exercise the clustering.
  gauge = list(
    flip_rate = 0.001, noise_floor = 1e-4, se_mean = mean(c(0.05, 0.03, 0.05)),
    se_sd = 0.01, se_mcse = 0.002, n_ent = 100, n_ref = 3, n_capped = 0,
    pair_i = c(0L, 0L, 0L), pair_j = c(1L, 2L, 1L),
    pair_se = c(0.05, 0.03, 0.05), pair_mcse = c(0.002, 0.002, 0.002)
  )
  chains = list(list(zratio = list(gauge = gauge)))
  a = 9
  b = 1
  E = 190
  th = 0.9
  pip = rep(th, E)
  s = summarize_zratio_gauge(
    chains,
    verbose = FALSE,
    harm_inputs = list(pip = list(pip), a = a, b = b)
  )
  pc = s$per_chain

  m = th * (1 - th)
  theta_hat = (a + E * th) / (a + b + E)
  gain = E * m / (theta_hat * (1 - theta_hat) * (a + b + E))
  A = 1 / (1 - min(gain, 0.98))
  expect_equal(pc$kappa, m * A)
  x = m * c(0.05, 0.03, 0.05)
  expect_equal(pc$harm_pred, abs(mean(x)) * A)
  # Resolution gate: cluster-robust spread + reference noise.
  cr = ((sum(x[c(1, 3)]) - 2 * mean(x))^2 + (x[2] - mean(x))^2) / 9
  noise2 = sum((m * 0.002)^2 * 3) / 9
  expect_true(abs(mean(x)) > 2 * sqrt(cr + noise2))
  expect_true(pc$harm_flag)

  # Heterogeneous sensitivities: errors on pinned edges must not count.
  pip2 = rep(th, E)
  pip2[1] = 0.999 # edge (0,1) pinned: m_e ~ 0
  s2 = summarize_zratio_gauge(
    chains,
    verbose = FALSE,
    harm_inputs = list(pip = list(pip2), a = a, b = b)
  )
  m1 = 0.999 * (1 - 0.999)
  m2 = th * (1 - th)
  x2 = c(m1 * 0.05, m2 * 0.03, m1 * 0.05)
  expect_equal(
    s2$per_chain$harm_pred,
    abs(mean(x2)) * s2$per_chain$amplification
  )
  expect_lt(s2$per_chain$harm_pred, pc$harm_pred)
})
