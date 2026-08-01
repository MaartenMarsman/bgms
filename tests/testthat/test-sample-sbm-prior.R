# Tests for sample_sbm_prior(): ancestral draws from the MFM-SBM
# edge-prior hyperprior.

test_that("sample_sbm_prior returns coherent draws", {
  d = sample_sbm_prior(p = 8, n_samples = 500, seed = 4)
  expect_identical(dim(d$allocations), c(500L, 8L))
  expect_identical(dim(d$pair_probability), c(500L, 28L))
  expect_length(d$num_blocks, 500L)
  expect_true(all(d$allocations >= 1L))
  expect_true(all(d$pair_probability > 0 & d$pair_probability < 1))
  expect_true(all(d$num_blocks >= 1L & d$num_blocks <= 8L))
  # num_blocks counts the occupied blocks of each allocation.
  expect_identical(
    d$num_blocks,
    apply(d$allocations, 1, function(z) length(unique(z)))
  )
})

test_that("pair probabilities follow the block structure", {
  # Dense within (Beta(8, 1)), sparse between (Beta(1, 8)).
  d = sample_sbm_prior(
    p = 8, n_samples = 500,
    edge_prior = sbm_prior(
      alpha = 8, beta = 1, alpha_between = 1, beta_between = 8
    ),
    seed = 2
  )
  pairs = which(upper.tri(matrix(0, 8, 8)), arr.ind = TRUE)
  pairs = pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
  within = matrix(
    d$allocations[, pairs[, 1L]] == d$allocations[, pairs[, 2L]],
    nrow = 500L
  )
  expect_gt(mean(d$pair_probability[within]), 0.8)
  expect_lt(mean(d$pair_probability[!within]), 0.2)
})

test_that("lambda shifts the number of blocks", {
  few = sample_sbm_prior(
    p = 10, n_samples = 400,
    edge_prior = sbm_prior(lambda = 0.1), seed = 6
  )
  many = sample_sbm_prior(
    p = 10, n_samples = 400,
    edge_prior = sbm_prior(lambda = 8), seed = 6
  )
  expect_lt(mean(few$num_blocks), mean(many$num_blocks))
})

test_that("input validation and reproducibility", {
  expect_error(
    sample_sbm_prior(p = 8, n_samples = 10, edge_prior = bernoulli_prior()),
    "sbm_prior"
  )
  expect_error(sample_sbm_prior(p = 1, n_samples = 10), "'p'")
  d1 = sample_sbm_prior(p = 6, n_samples = 20, seed = 9)
  d2 = sample_sbm_prior(p = 6, n_samples = 20, seed = 9)
  expect_identical(d1, d2)
})
