# Regression tests for the 2026-07-03 audit high-confidence (H) fixes.

# ---------------------------------------------------------------------------
# H13 - mixed model with a single continuous variable
# ---------------------------------------------------------------------------

test_that("mixed model with one continuous variable does not crash (H13)", {
  skip_on_cran()

  set.seed(2)
  n = 60
  x = cbind(
    d1 = sample(0:2, n, replace = TRUE),
    d2 = sample(0:2, n, replace = TRUE),
    c1 = rnorm(n)
  )

  expect_no_error(
    bgm(
      x,
      variable_type = c("ordinal", "ordinal", "continuous"),
      edge_selection = FALSE, update_method = "adaptive-metropolis",
      iter = 30, warmup = 30, chains = 1, seed = 1, display_progress = "none"
    )
  )
})

# ---------------------------------------------------------------------------
# H15 - OMRF SBM allocation samples carry node names, not edge names
# ---------------------------------------------------------------------------

test_that("OMRF SBM allocation samples carry node names (H15)", {
  skip_on_cran()

  set.seed(3)
  p = 4
  x = matrix(sample(0:2, 60 * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("V", seq_len(p))

  fit = bgm(
    x,
    variable_type = "ordinal", edge_selection = TRUE,
    edge_prior = sbm_prior(),
    iter = 30, warmup = 30, chains = 1, seed = 1, display_progress = "none"
  )

  alloc_names = fit$raw_samples$parameter_names$allocations
  expect_equal(alloc_names, colnames(x))
})

# ---------------------------------------------------------------------------
# H16 - explicit integer baseline_category = 0L is honoured
# ---------------------------------------------------------------------------

test_that("explicit baseline_category = 0L works for Blume-Capel (H16)", {
  skip_on_cran()

  set.seed(4)
  p = 3
  x = matrix(sample(0:3, 60 * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("V", seq_len(p))

  # Passing integer 0L used to be read as "not provided" and errored with
  # "baseline_category is required for Blume-Capel variables".
  expect_no_error(
    bgm(
      x,
      variable_type = "blume-capel", baseline_category = 0L,
      edge_selection = FALSE, update_method = "adaptive-metropolis",
      iter = 20, warmup = 20, chains = 1, seed = 1, display_progress = "none"
    )
  )
})

# ---------------------------------------------------------------------------
# H9 - Blume-Capel prediction/simulation on the original category scale
# ---------------------------------------------------------------------------

test_that("Blume-Capel recode round-trips through the stored shift (H9)", {
  # Data coded 2..5 shifts to internal 0..3; the reverse map restores 2..5.
  x = matrix(c(2, 3, 4, 5), ncol = 1)
  shift = c(2)

  internal = recode_data_for_prediction(
    x,
    num_categories = 3L, is_ordinal = FALSE,
    category_levels = NULL, blume_capel_shift = shift
  )
  expect_equal(as.numeric(internal), c(0, 1, 2, 3))

  original = recode_simulated_to_original(
    internal,
    category_levels = list(NULL), blume_capel_shift = shift
  )
  expect_equal(as.numeric(original), c(2, 3, 4, 5))
})

test_that("simulate() returns Blume-Capel draws on the original scale (H9)", {
  skip_on_cran()

  set.seed(7)
  p = 3
  x = matrix(sample(1:5, 150 * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("B", seq_len(p))

  fit = bgm(
    x,
    variable_type = "blume-capel", baseline_category = 3,
    edge_selection = FALSE, update_method = "adaptive-metropolis",
    iter = 200, warmup = 200, chains = 1, seed = 1, display_progress = "none"
  )

  # The shift back to the original 1-based coding is stored per variable.
  expect_equal(extract_arguments(fit)$blume_capel_shift, rep(1, p))

  sim = simulate(fit, nsim = 200, seed = 2)
  # Original scale is 1..5; the internal 0-based scale would include 0.
  expect_gte(min(sim), 1)
  expect_lte(max(sim), 5)
})

# ---------------------------------------------------------------------------
# H10 - GGM prediction centers newdata on the training means
# ---------------------------------------------------------------------------

test_that("GGM predict uses stored training means, not newdata means (H10)", {
  skip_on_cran()

  set.seed(8)
  p = 4
  n = 120
  x = sweep(matrix(rnorm(n * p), n, p), 2, c(5, -3, 10, 0), "+")
  colnames(x) = paste0("V", seq_len(p))

  fit = bgm(
    x,
    variable_type = "continuous", update_method = "adaptive-metropolis",
    edge_selection = FALSE, iter = 300, warmup = 300, chains = 1, seed = 1,
    display_progress = "none"
  )

  # Training means are stored so prediction centers on the fitted scale.
  expect_equal(extract_arguments(fit)$column_means, colMeans(x),
    ignore_attr = TRUE
  )

  # With one row, centering by newdata's own mean gave back the input value
  # (a variable predicting itself). Centering on the training means does not.
  row1 = x[1, , drop = FALSE]
  pred = predict(fit, newdata = row1)
  pred_means = vapply(pred, function(m) m[, "mean"], numeric(1))
  expect_false(isTRUE(all.equal(unname(pred_means), unname(row1[1, ]))))
})

# ---------------------------------------------------------------------------
# H17 - Rhat splits each chain so it detects within-chain drift (split-Rhat)
# ---------------------------------------------------------------------------

test_that("split_chains halves each chain (H17)", {
  a = array(seq_len(10 * 2 * 3), dim = c(10, 2, 3))
  s = split_chains(a)
  expect_equal(dim(s), c(5L, 4L, 3L))
  # sub-chains 1-2 are the two halves of the original chain 1
  expect_equal(s[, 1, ], a[1:5, 1, ])
  expect_equal(s[, 2, ], a[6:10, 1, ])
  # an odd iteration count drops the middle draw
  odd = array(1:7, dim = c(7, 1, 1))
  expect_equal(as.vector(split_chains(odd)), c(1, 2, 3, 5, 6, 7))
})

test_that("split-Rhat flags a drifting chain that classic Rhat misses (H17)", {
  set.seed(17)
  # A single chain that drifts: classic Gelman-Rubin needs >1 chain and returns
  # NA, but split-Rhat compares the two halves and sees the trend.
  drift = array(cumsum(rnorm(2000, mean = 0.03)) + rnorm(2000), dim = c(2000, 1, 1))
  expect_true(is.na(.compute_rhat_cpp(drift)))
  expect_gt(.compute_rhat_cpp(split_chains(drift)), 1.1)

  # A stationary chain still passes.
  stat = array(rnorm(2000), dim = c(2000, 1, 1))
  expect_lt(.compute_rhat_cpp(split_chains(stat)), 1.05)
})

# ---------------------------------------------------------------------------
# H8 - bgmCompare difference prior family is independent of the interaction
# prior, defaulting to Cauchy
# ---------------------------------------------------------------------------

test_that("difference_family selects the difference prior independently (H8)", {
  skip_on_cran()

  set.seed(8)
  ng = 2
  npg = 40
  p = 3
  x = matrix(sample(0:2, ng * npg * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("V", seq_len(p))
  g = rep(seq_len(ng), each = npg)

  run = function(fam) {
    bgmCompare(
      x,
      group_indicator = g, difference_family = fam,
      difference_selection = FALSE, iter = 80, warmup = 80, chains = 1,
      seed = 3, display_progress = "none"
    )
  }

  cauchy = run("Cauchy")
  normal = run("Normal")

  # Same seed, different prior family -> different difference draws.
  expect_false(isTRUE(all.equal(
    cauchy$raw_samples$pairwise[[1]], normal$raw_samples$pairwise[[1]]
  )))

  # The default is Cauchy.
  default = run("Cauchy")
  expect_equal(
    default$raw_samples$pairwise[[1]], cauchy$raw_samples$pairwise[[1]]
  )

  # An unknown family is rejected.
  expect_error(
    bgmCompare(
      x,
      group_indicator = g, difference_family = "laplace",
      difference_selection = FALSE, iter = 10, warmup = 10, chains = 1,
      seed = 1, display_progress = "none"
    )
  )
})
