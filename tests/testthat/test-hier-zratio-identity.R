# Graph-law identity gates for the hierarchical prior specification. Under
# spec = "hierarchical" the graph marginal is exactly the edge prior:
# Bernoulli(p) gives edge count ~ Binomial(E, p) and per-edge PIP = p;
# Beta-Bernoulli(a, b) gives theta ~ Beta(a, b) with clean conjugate
# updates (no C-correction on this path). The joint spec fails these
# identities (its graph marginal carries the Z(Gamma) tilt), which is the
# negative control. Fast identity smokes run under NOT_CRAN; the fuller
# cell battery is gated behind BGMS_RUN_SLOW_TESTS.

hier_prior_run = function(q, delta, sigma, p_inc, um, edge_prior = NULL,
                          spec = "hierarchical", n_samples = 4000L,
                          seed = 7L) {
  sample_ggm_prior(
    p = q, n_samples = n_samples, n_warmup = 1500L,
    interaction_prior = normal_prior(scale = sigma / 2),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    spec = spec, edge_inclusion_prob = p_inc,
    update_method = um, edge_prior = edge_prior,
    delta = delta, seed = seed, verbose = FALSE
  )
}

test_that("hierarchical graph marginal holds at gamma shapes", {
  skip_on_cran()
  # The generalized constants carry the diagonal Gamma shape through every
  # channel and the oracle sweep; the graph law must hold away from the
  # exponential (shape = 1) cell on both update methods.
  for(shape in c(0.5, 2)) {
    for(um in c("adaptive-metropolis", "gibbs")) {
      d = sample_ggm_prior(
        p = 6L, n_samples = 6000L, n_warmup = 1500L,
        interaction_prior = normal_prior(scale = 0.5),
        precision_scale_prior = gamma_prior(shape = shape, rate = 2),
        spec = "hierarchical", edge_inclusion_prob = 0.3,
        update_method = um, delta = 0.5 * log(6), seed = 7L,
        verbose = FALSE
      )
      expect_lt(
        abs(mean(d$edge_indicators) - 0.3), 0.02,
        label = paste0("shape = ", shape, ", ", um)
      )
    }
  }
})

test_that("gamma-shape constants build in the standardized cell", {
  skip_on_cran()
  # Same eta and shape, different frames: one constant set.
  d = 0.5 * log(6)
  a = bgms:::zratio_cell_constants(
    d,
    pairwise_scale = 0.5, scale_rate = 2, scale_shape = 2
  )
  b = bgms:::zratio_cell_constants(
    d,
    pairwise_scale = 0.25, scale_rate = 4, scale_shape = 2
  )
  expect_identical(a, b)
  expect_identical(a$alpha, 2)
})

test_that("hierarchical graph marginal matches Bernoulli(p); joint does not", {
  skip_on_cran()
  # rate = 2 (the eta = 1 default at scale 0.5) gives tau^2 = 2 beta sigma^2
  # = 2, where the joint spec's Z(Gamma) tilt separates cleanly from the
  # Bernoulli target (joint marginal ~0.22 at p = 0.3).
  run = function(um, spec) {
    sample_ggm_prior(
      p = 6L, n_samples = 4000L, n_warmup = 1500L,
      interaction_prior = normal_prior(scale = 0.5),
      spec = spec, edge_inclusion_prob = 0.3,
      update_method = um, seed = 7L, verbose = FALSE
    )
  }
  for(um in c("adaptive-metropolis", "gibbs")) {
    d = run(um, "hierarchical")
    expect_lt(abs(mean(d$edge_indicators) - 0.3), 0.02)
    dj = run(um, "joint")
    expect_gt(abs(mean(dj$edge_indicators) - 0.3), 0.05)
  }
})

test_that("Z-ratio constants build in the standardized cell", {
  skip_on_cran()
  # The between-graph ratio depends on (delta, eta) only, so every user
  # frame with the same eta = pairwise_scale * scale_rate must resolve to
  # one constant set, built at sigma = 1 (the frame the quadrature grids
  # are sized for). At the old bare-scale build the scale-2.5 frame drifted
  # the saddle by -0.10 in log J and the q = 6 identity by -0.008.
  d = 0.5 * log(6)
  a = bgms:::zratio_cell_constants(d, pairwise_scale = 0.5, scale_rate = 2)
  b = bgms:::zratio_cell_constants(d, pairwise_scale = 0.25, scale_rate = 4)
  e = bgms:::zratio_cell_constants(
    d,
    pairwise_scale = 2.5, scale_rate = 0.4, scale_eta = 1
  )
  expect_identical(a, b)
  expect_identical(a, e)
  expect_identical(a$eta, 1)
})

test_that("hierarchical graph marginal holds at a non-unit slab scale", {
  skip_on_cran()
  # pairwise scale 2.5 (the sample_ggm_prior default) with eta = 1 is the
  # same physical cell as scale 0.5 / rate 2; the graph law must hold there
  # identically.
  d = sample_ggm_prior(
    p = 6L, n_samples = 6000L, n_warmup = 1500L,
    interaction_prior = normal_prior(scale = 2.5),
    precision_scale_prior = gamma_prior(shape = 1, eta = 1),
    spec = "hierarchical", edge_inclusion_prob = 0.3,
    update_method = "adaptive-metropolis", delta = 0.5 * log(6),
    seed = 11L, verbose = FALSE
  )
  expect_lt(abs(mean(d$edge_indicators) - 0.3), 0.02)
})

test_that("hierarchical graph marginal holds for the Cauchy slab", {
  skip_on_cran()
  # The Cauchy slab needs its own (marginal-Cauchy) Z-ratio constants; with
  # the Normal tables this cell read 0.247 for a 0.30 edge prior.
  for(um in c("adaptive-metropolis", "gibbs")) {
    d = sample_ggm_prior(
      p = 6L, n_samples = 6000L, n_warmup = 1500L,
      interaction_prior = cauchy_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 1, rate = 2),
      spec = "hierarchical", edge_inclusion_prob = 0.3,
      update_method = um, delta = 0.5 * log(6), seed = 7L, verbose = FALSE
    )
    expect_lt(abs(mean(d$edge_indicators) - 0.3), 0.02, label = um)
  }
})

test_that("Cauchy graph law holds at dense high q (coupling regime)", {
  skip_on_cran()
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to run the dense high-q Cauchy identity"
  )
  # The single-edge memo flags the Gaussian-mixture slab as a coupling-
  # sensitive case (the additive moment approximation weakens on dense
  # blocks). At dense high q (q = 15, 20; p_inc = 0.5, 0.7) the deployed
  # system -- additive saddle + warm-up OLS correction on maxbd >= 2 blocks
  # (the default window engages at p >= 15) -- must still reproduce the edge
  # prior. The marginal is the release-relevant statistic; the per-block
  # alarm verdict is expected to flag here (it does for the Normal slab too)
  # and is not asserted.
  suppressMessages(library(parallel))
  cells = expand.grid(
    q = c(15L, 20L), p_inc = c(0.5, 0.7),
    um = c("adaptive-metropolis", "gibbs"), stringsAsFactors = FALSE
  )
  devs = unlist(mclapply(seq_len(nrow(cells)), function(r) {
    cl = cells[r, ]
    d = sample_ggm_prior(
      p = cl$q, n_samples = 6000L, n_warmup = 2000L,
      interaction_prior = cauchy_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 1, rate = 2),
      spec = "hierarchical", edge_inclusion_prob = cl$p_inc,
      update_method = cl$um, delta = 0.5 * log(cl$q),
      seed = 4000L + r, verbose = FALSE
    )
    mean(d$edge_indicators) - cl$p_inc
  }, mc.cores = min(8L, nrow(cells))))
  expect_lt(max(abs(devs)), 0.02)
})

test_that("hierarchical BB identity: theta ~ Beta(a, b), PIP = a/(a+b)", {
  skip_on_cran()
  d = hier_prior_run(
    6L, 0.5 * log(6), 1, 0.5, "adaptive-metropolis",
    edge_prior = beta_bernoulli_prior(2, 4), n_samples = 6000L
  )
  # theta mixes with IAT ~ 10-25; thin to near-independence for the KS and
  # size the mean tolerance at ~3 SE of the thinned mean.
  th_thin = d$theta[seq(1, length(d$theta), by = 60L)]
  expect_lt(abs(mean(d$theta) - 2 / 6), 0.035)
  expect_gt(suppressWarnings(ks.test(th_thin, pbeta, 2, 4)$p.value), 0.01)
  expect_lt(abs(mean(d$edge_indicators) - 1 / 3), 0.035)
})

test_that("hierarchical graph law across cells (slow battery)", {
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run the cell battery"
  )
  q = 10L
  E = q * (q - 1L) / 2L
  cells = list(
    c(0, 1, 0.2), c(0, 1, 0.5), c(0.5 * log(10), 1, 0.5),
    c(2, 2, 0.8), c(0.5, 0.5, 0.3)
  )
  for(um in c("adaptive-metropolis", "gibbs")) {
    for(cl in cells) {
      d = hier_prior_run(q, cl[1], cl[2], cl[3], um, n_samples = 6000L)
      ne = rowSums(d$edge_indicators)
      pip = colMeans(d$edge_indicators)
      lbl = sprintf("%s delta=%.2f sigma=%.1f p=%.1f", um, cl[1], cl[2], cl[3])
      expect_lt(
        abs(mean(ne) - E * cl[3]), 0.035 * E * cl[3] + 0.6,
        label = lbl
      )
      expect_lt(max(abs(pip - cl[3])), 0.05, label = lbl)
    }
  }
})
