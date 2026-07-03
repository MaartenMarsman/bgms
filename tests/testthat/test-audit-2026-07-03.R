# Regression tests for the 2026-07-03 codebase audit fixes (C1-C7, H1).

# ---------------------------------------------------------------------------
# C5 - GGM summary() pairwise means on the association scale
# ---------------------------------------------------------------------------

test_that("GGM summary() pairwise means match coef() (association scale) (C5)", {
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
      tolerance = 1e-6))
  )
})
