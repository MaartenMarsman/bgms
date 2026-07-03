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
