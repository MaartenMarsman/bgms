# Tests for sample_graph_prior(): the hierarchical spec's ancestral graph
# law, hyperparameter conditioning, and the joint spec's Z(Gamma)-tilted
# negative control (mirroring test-hier-zratio-identity).

test_that("hierarchical Bernoulli graph law matches its prior", {
  g = sample_graph_prior(
    p = 6, n_samples = 4000,
    edge_prior = bernoulli_prior(0.3), seed = 11
  )
  expect_identical(dim(g$edge_indicators), c(4000L, 15L))
  expect_true(all(g$edge_indicators %in% 0:1))
  expect_lt(abs(mean(g$edge_indicators) - 0.3), 0.02)
  expect_null(g$theta)
  expect_null(g$allocations)
  expect_identical(g$spec, "hierarchical")
  expect_identical(g$edge_prior, "Bernoulli")
  expect_identical(g$pair_names[1:2], c("1-2", "1-3"))
})

test_that("hierarchical Beta-Bernoulli samples theta and coheres with it", {
  a = 2
  b = 4
  g = sample_graph_prior(
    p = 6, n_samples = 4000,
    edge_prior = beta_bernoulli_prior(a, b), seed = 5
  )
  expect_length(g$theta, 4000L)
  expect_true(all(g$theta > 0 & g$theta < 1))
  expect_lt(abs(mean(g$theta) - a / (a + b)), 0.02)
  expect_lt(abs(var(g$theta) - a * b / ((a + b)^2 * (a + b + 1))), 0.01)
  # The graph marginal integrates theta out to a/(a + b), and per draw the
  # edge mean tracks the sampled theta.
  expect_lt(abs(mean(g$edge_indicators) - a / (a + b)), 0.02)
  expect_gt(cor(rowMeans(g$edge_indicators), g$theta), 0.7)
})

test_that("hierarchical SBM separates within- from between-block edges", {
  g = sample_graph_prior(
    p = 8, n_samples = 1000,
    edge_prior = sbm_prior(
      alpha = 6, beta = 1, alpha_between = 1, beta_between = 6
    ),
    seed = 7
  )
  expect_identical(dim(g$allocations), c(1000L, 8L))
  expect_true(all(g$allocations >= 1L))

  pairs = which(upper.tri(matrix(0, 8, 8)), arr.ind = TRUE)
  pairs = pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
  within = matrix(
    g$allocations[, pairs[, 1L]] == g$allocations[, pairs[, 2L]],
    nrow = 1000L
  )
  ind = g$edge_indicators
  expect_gt(mean(ind[within]), mean(ind[!within]) + 0.3)
})

test_that("conditioning on theta fixes the graph law", {
  g = sample_graph_prior(
    p = 6, n_samples = 2000,
    edge_prior = beta_bernoulli_prior(1, 1), theta = 0.8, seed = 3
  )
  expect_lt(abs(mean(g$edge_indicators) - 0.8), 0.02)
  expect_null(g$theta)
})

test_that("conditioning on (allocations, block_probs) fixes the SBM law", {
  z = c(1L, 1L, 1L, 2L, 2L, 2L)
  bp = matrix(c(0.9, 0.05, 0.05, 0.9), 2, 2)
  g = sample_graph_prior(
    p = 6, n_samples = 2000,
    edge_prior = sbm_prior(), allocations = z, block_probs = bp, seed = 9
  )
  pairs = which(upper.tri(matrix(0, 6, 6)), arr.ind = TRUE)
  pairs = pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
  within = z[pairs[, 1L]] == z[pairs[, 2L]]
  expect_lt(abs(mean(g$edge_indicators[, within]) - 0.9), 0.03)
  expect_lt(abs(mean(g$edge_indicators[, !within]) - 0.05), 0.03)
  expect_null(g$allocations)
})

test_that("conditioning arguments are validated", {
  expect_error(
    sample_graph_prior(
      p = 5, n_samples = 10,
      edge_prior = sbm_prior(), theta = 0.5
    ),
    "allocations"
  )
  expect_error(
    sample_graph_prior(
      p = 5, n_samples = 10,
      edge_prior = sbm_prior(), allocations = rep(1L, 5)
    ),
    "block_probs"
  )
  expect_error(
    sample_graph_prior(
      p = 5, n_samples = 10,
      edge_prior = bernoulli_prior(0.5),
      allocations = rep(1L, 5),
      block_probs = matrix(0.5, 1, 1)
    ),
    "sbm_prior"
  )
  expect_error(
    sample_graph_prior(p = 5, n_samples = 10, theta = 1.5),
    "theta"
  )
  expect_error(
    sample_graph_prior(
      p = 5, n_samples = 10,
      edge_prior = sbm_prior(),
      allocations = c(1L, 1L, 2L, 2L, 3L),
      block_probs = matrix(0.5, 2, 2)
    ),
    "exceed"
  )
})

test_that("the same seed reproduces the draw and restores the RNG state", {
  set.seed(123)
  before = runif(1)
  set.seed(123)
  g1 = sample_graph_prior(p = 6, n_samples = 50, seed = 42)
  after = runif(1)
  g2 = sample_graph_prior(p = 6, n_samples = 50, seed = 42)
  expect_identical(g1$edge_indicators, g2$edge_indicators)
  expect_identical(before, after)
})

test_that("the joint spec produces the Z(Gamma)-tilted law, not the prior", {
  skip_on_cran()
  # Mirrors the negative control in test-hier-zratio-identity: at
  # tau^2 = 2 the joint marginal separates cleanly from Bernoulli(0.3)
  # (about 0.22), while the hierarchical law sits at 0.3.
  g = sample_graph_prior(
    p = 6, n_samples = 4000,
    edge_prior = bernoulli_prior(0.3), spec = "joint",
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 2),
    n_warmup = 1500, seed = 7, verbose = FALSE
  )
  expect_identical(dim(g$edge_indicators), c(4000L, 15L))
  expect_gt(abs(mean(g$edge_indicators) - 0.3), 0.05)

  h = sample_graph_prior(
    p = 6, n_samples = 4000,
    edge_prior = bernoulli_prior(0.3), seed = 7
  )
  expect_lt(abs(mean(h$edge_indicators) - 0.3), 0.02)
})
