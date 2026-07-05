# Tests for the Z-ratio alarm suite (hierarchical GGM prior specification):
# verdict thresholds, drift flags, the per-chain diagnostics block recorded
# by the sampler, and the audit channels on healthy and defective kernels.

test_that("zratio_tau maps regimes to the measured thresholds", {
  expect_equal(bgms:::zratio_tau(1), 0.01)
  expect_equal(bgms:::zratio_tau(1.9), 0.01)
  expect_equal(bgms:::zratio_tau(2), 0.02)
  expect_equal(bgms:::zratio_tau(3), 0.04)
  expect_equal(bgms:::zratio_tau(5), 0.04)
})

test_that("zratio_drift_flag separates flat from drifting streams", {
  n_pairs = 45
  set.seed(11)
  flat = 0.3 + rnorm(200, sd = 0.5 / sqrt(n_pairs) / 4)
  drifting = seq(0.5, 0.15, length.out = 200)
  expect_false(bgms:::zratio_drift_flag(flat, n_pairs))
  expect_true(bgms:::zratio_drift_flag(drifting, n_pairs))
  expect_true(is.na(bgms:::zratio_drift_flag(flat[1:20], n_pairs)))
})

test_that("zratio_indicator_graph rebuilds the adjacency", {
  p = 5
  G = matrix(0L, p, p)
  G[1, 3] = G[3, 1] = 1L
  G[2, 5] = G[5, 2] = 1L
  diag(G) = 1L
  column = integer(p * (p + 1) / 2)
  e = 1L
  for(i in seq_len(p)) {
    for(j in i:p) {
      column[e] = G[i, j]
      e = e + 1L
    }
  }
  expect_identical(bgms:::zratio_indicator_graph(column, p), G)
})

test_that("hierarchical chains attach the diagnostics block and pass quiet", {
  skip_on_cran()
  draws = sample_ggm_prior(
    p = 10, n_samples = 120, n_warmup = 200,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 2),
    spec = "hierarchical", update_method = "gibbs",
    edge_prior = beta_bernoulli_prior(2, 4),
    calibration_window = 80, seed = 7, verbose = FALSE
  )
  zd = draws$zratio_diagnostics
  expect_false(is.null(zd))
  pc = zd$per_chain
  expect_equal(nrow(pc), 1L)
  expect_true(pc$frozen)
  expect_gt(pc$n_oracle, 0)
  expect_equal(pc$n_anchors, pc$n_oracle)
  expect_gt(pc$n_audit, 0)
  # Healthy kernel at eta = 1: every audit channel sits under the verdict
  # threshold and the calibration stream is not drifting.
  expect_equal(zd$eta, 1)
  expect_equal(zd$tau, 0.01)
  expect_lt(pc$gate, zd$tau)
  expect_false(zd$verdict_flagged)
  expect_false(zd$calibration_incomplete)
  expect_false(pc$drift_density)
  expect_false(pc$drift_theta)
})

test_that("the audit flags a defective frozen kernel", {
  skip_on_cran()
  p = 10
  delta = 0.5 * log(p)
  zc = bgms:::zratio_constants(delta, sigma = 1, beta = 1)
  # Pack a fit whose correction is a constant 0.05 on every coupled-bridge
  # block, with a hull box wide enough that the clamp never engages. The
  # oracle discrepancy on healthy blocks is O(0.001), so every targeted
  # audit reads ~0.05 > tau = 0.01.
  addc = numeric(23)
  addc[1:6] = zc$addc[1:6]
  addc[7] = 0.05
  addc[13] = 1
  addc[14:23] = rep(c(0, 1e6), 5)
  set.seed(21)
  G = matrix(0L, p, p)
  ut = upper.tri(G)
  G[ut] = rbinom(sum(ut), 1, 0.5)
  G = G + t(G)
  diag(G) = 1L
  scan = bgms:::zratio_scan_graph(
    G, addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0
  )
  expect_gt(sum(scan[, 6] >= 2), 0) # coupled-bridge blocks present
  column = integer(p * (p + 1) / 2)
  e = 1L
  for(i in seq_len(p)) {
    for(j in i:p) {
      column[e] = G[i, j]
      e = e + 1L
    }
  }
  chain = list(
    indicator_samples = matrix(column, ncol = 4, nrow = length(column)),
    zratio = list(
      addc = addc,
      anchors_x = matrix(numeric(0), 0, 6),
      anchors_y = numeric(0),
      counters = c(
        n_hit = 0, n_miss = 0, n_pred = 0, n_add = 0, n_clamp = 0,
        n_oracle = 0, n_anchors = 0, cache_size = 0, frozen = 1
      ),
      warmup_density = numeric(0),
      warmup_theta = numeric(0)
    )
  )
  spec = list(
    tg = zc$tg, ihat = zc$ihat, ghat = zc$ghat, wt = zc$wt, psi0 = zc$psi0,
    delta = delta, sigma = 1, beta = 1
  )
  old = options(bgms.verbose = TRUE)
  on.exit(options(old), add = TRUE)
  expect_output(
    zd <- summarize_zratio_diagnostics(
      list(chain), spec,
      num_nodes = p, n_graphs = 2, seed = 3,
      verbose = TRUE
    ),
    "Audit verdict"
  )
  expect_true(zd$verdict_flagged)
  expect_gt(zd$per_chain$gate, zd$tau)
})
