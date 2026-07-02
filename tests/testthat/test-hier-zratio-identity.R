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

test_that("hierarchical spec requires the Normal slab and shape-1 diagonal", {
  expect_error(
    sample_ggm_prior(
      p = 4, n_samples = 5, n_warmup = 5, spec = "hierarchical",
      interaction_prior = cauchy_prior(scale = 1), verbose = FALSE
    ),
    "normal interaction"
  )
  expect_error(
    sample_ggm_prior(
      p = 4, n_samples = 5, n_warmup = 5, spec = "hierarchical",
      interaction_prior = normal_prior(scale = 1),
      precision_scale_prior = gamma_prior(shape = 2, rate = 1),
      verbose = FALSE
    ),
    "shape = 1"
  )
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
