# Regression tests pinning fixed defects in imputation row bookkeeping,
# group-difference labels, summary scales, and mixed-model sampling paths.

# ---------------------------------------------------------------------------
# imputation row indices follow the group-sorted row order
# ---------------------------------------------------------------------------

test_that("imputation targets the right rows when groups are interleaved", {
  # 6 people, 2 ordinal variables, groups interleaved (1,2,1,2,1,2).
  x = data.frame(
    V1 = c(NA, 1, 2, 0, 1, 2),
    V2 = c(2, 1, 0, NA, 1, 0)
  )
  group_indicator = c(1, 2, 1, 2, 1, 2)

  spec = bgm_spec(
    x = x, model_type = "compare",
    group_indicator = group_indicator, na_action = "impute",
    edge_selection = FALSE, difference_selection = FALSE, seed = 1
  )

  mi = spec$missing$missing_index

  # Group 1 people (rows 1,3,5) move to sorted rows 1,2,3; group 2 people
  # (rows 2,4,6) move to sorted rows 4,5,6. The missing cell for person 1
  # is in sorted row 1 (0-based 0); the one for person 4 is in sorted row 5
  # (0-based 4). Column order: V1 first, then V2.
  expect_equal(mi[, 1], c(0, 4))
  expect_equal(mi[, 2], c(0, 1))
})

# ---------------------------------------------------------------------------
# difference labels attached to the right parameter for 3+ groups
# ---------------------------------------------------------------------------

test_that("pairwise difference summaries carry the right labels for 3 groups", {
  skip_on_cran()

  set.seed(11)
  n_groups = 3
  n_per_group = 20
  p = 3
  x = matrix(sample(0:2, n_groups * n_per_group * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("V", seq_len(p))
  group_indicator = rep(seq_len(n_groups), each = n_per_group)

  fit = bgmCompare(
    x = x, group_indicator = group_indicator,
    difference_selection = TRUE,
    iter = 50, warmup = 50, chains = 1, seed = 3,
    display_progress = "none"
  )

  actual = fit$posterior_summary_pairwise_differences$parameter

  # The rows run edge by edge (in i<j order), each edge repeated once per
  # group difference. The labels must line up with that order.
  expected = character()
  for(i in seq_len(p - 1L)) {
    for(j in (i + 1L):p) {
      for(h in seq_len(n_groups - 1L)) {
        expected = c(expected, sprintf("V%d-V%d (diff%d)", i, j, h))
      }
    }
  }

  expect_equal(actual, expected)
})

# ---------------------------------------------------------------------------
# GGM summary() pairwise means on the association scale
# ---------------------------------------------------------------------------

test_that("GGM summary() pairwise means match coef() (association scale)", {
  skip_on_cran()

  set.seed(5)
  p = 4
  n = 150
  x = matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) = paste0("V", seq_len(p))

  fit = bgm(
    x = x, variable_type = "continuous",
    update_method = "adaptive-metropolis", edge_selection = FALSE,
    iter = 400, warmup = 400, chains = 1, seed = 7,
    display_progress = "none"
  )

  summary_mean = fit$posterior_summary_pairwise$mean
  coef_pairwise = coef(fit)$pairwise
  coef_upper = coef_pairwise[upper.tri(coef_pairwise)]

  # summary() and coef() report the same association-scale pairwise means
  expect_equal(sort(summary_mean), sort(coef_upper), tolerance = 1e-8)

  # and not the raw precision scale (which is association * -2)
  expect_false(
    isTRUE(all.equal(sort(summary_mean), sort(-2 * coef_upper),
      tolerance = 1e-6
    ))
  )
})

# ---------------------------------------------------------------------------
# mixed model imputes discrete-only missing data
#
# A stale cached mean biases sampling without failing, so there is no cheap
# deterministic check; this exercises the path and checks the results are
# well formed.
# ---------------------------------------------------------------------------

test_that("mixed model runs with missing values in discrete columns only", {
  skip_on_cran()

  set.seed(303)
  n = 80
  x = cbind(
    d1 = sample(0:2, n, replace = TRUE),
    c1 = rnorm(n),
    d2 = sample(0:2, n, replace = TRUE),
    c2 = rnorm(n),
    d3 = sample(0:2, n, replace = TRUE)
  )
  x[sample(n, 6), "d1"] = NA
  x[sample(n, 6), "d3"] = NA

  fit = bgm(
    x = x,
    variable_type = c("ordinal", "continuous", "ordinal", "continuous", "ordinal"),
    na_action = "impute", edge_selection = FALSE,
    update_method = "adaptive-metropolis",
    iter = 100, warmup = 100, chains = 1, seed = 5,
    display_progress = "none"
  )

  expect_true(all(is.finite(coef(fit)$pairwise)))
  expect_true(all(is.finite(unlist(coef(fit)$main))))
})

# ---------------------------------------------------------------------------
# mixed NUTS gradient after imputation with edge selection
#
# Wrong gradients degrade mixing rather than erroring, so this exercises
# the NUTS + imputation + selection path and checks well-formed output.
# ---------------------------------------------------------------------------

test_that("mixed NUTS runs with imputation and edge selection", {
  skip_on_cran()

  set.seed(404)
  n = 80
  x = cbind(
    d1 = sample(0:2, n, replace = TRUE),
    c1 = rnorm(n),
    d2 = sample(0:2, n, replace = TRUE),
    c2 = rnorm(n),
    d3 = sample(0:2, n, replace = TRUE)
  )
  x[sample(n, 6), "d1"] = NA

  fit = bgm(
    x = x,
    variable_type = c("ordinal", "continuous", "ordinal", "continuous", "ordinal"),
    na_action = "impute", edge_selection = TRUE,
    update_method = "nuts",
    iter = 100, warmup = 100, chains = 1, seed = 6,
    display_progress = "none"
  )

  expect_true(all(is.finite(coef(fit)$pairwise)))
})

# ---------------------------------------------------------------------------
# fractional Dirichlet parameter can open new blocks
#
# The type fix is verified at compile time. Showing the behaviour (the sampler
# reaching three or more occupied blocks, which a value truncated to 0 forbids)
# needs data simulated with a clear three-block structure and a long enough run
# to be stable, which belongs in the SBC / block-recovery suite rather than a
# fast unit test.
# ---------------------------------------------------------------------------

test_that("fractional Dirichlet parameter can open new blocks", {
  skip("needs simulated three-block data + long run; covered by block-recovery suite")
})

# ---------------------------------------------------------------------------
# step-size search sign
#
# The fix only affects warmup cost, not the target distribution, so a fast pass
# / fail assertion is not available; warmup efficiency is measured by the
# benchmark suite.
# ---------------------------------------------------------------------------

test_that("step-size search does not collapse the step size", {
  skip("warmup-efficiency only; no effect on the posterior to assert on")
})
