# ==============================================================================
# Core bgm() Function Tests
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Reproducibility, correlation with sufficient statistics
#
# These tests verify core bgm() functionality:
#   - Reproducibility: identical seeds produce identical MCMC chains
#   - Sanity: posterior means correlate with classical sufficient statistics
#
# INTEGRATION NOTE: Many sampler configurations (adaptive-metropolis,
# Blume-Capel, missing data imputation, standardization) are tested via the
# parameterized fixture approach in test-methods.R. See:
#   - helper-fixtures.R: Cached fit functions
#   - test-methods.R: get_bgms_fixtures() loops over all configurations
#
# This file focuses on tests that require special setup or unique assertions.
# ==============================================================================

# Tiers. The product surface stays local. The estimate-simulate-re-estimate
# cycle is parameter recovery -- it proves the settled numerics still recover
# what they planted -- so it runs nightly.

skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run the parameter-recovery cycle"
  )
}

test_that("bgm is reproducible", {
  # Use cached fixture as fit1, run one fresh fit as fit2 with same params
  fit1 = get_bgms_fit_ordinal()

  data("Wenchuan", package = "bgms")
  fit2 = bgm(
    Wenchuan[1:50, 1:4],
    iter = 50, warmup = 100, chains = 2,
    seed = 12345, # Same seed as fixture
    display_progress = "none"
  )

  testthat::expect_equal(fit1$raw_samples, fit2$raw_samples)
})

test_that("bgmCompare is reproducible", {
  # Use cached fixture as fit1, run one fresh fit as fit2 with same params
  fit1 = get_bgmcompare_fit_ordinal()

  data("Wenchuan", package = "bgms")
  x = Wenchuan[1:50, 1:4]
  group_ind = rep(1:2, each = 25)
  fit2 = without_support_warning(bgmCompare(
    x = x, group_indicator = group_ind,
    iter = 50, warmup = 100, chains = 2,
    seed = 54321, # Same seed as fixture
    display_progress = "none"
  ))

  combine_chains = function(fit) {
    pairs = do.call(rbind, fit$raw_samples$pairwise)
    mains = do.call(rbind, fit$raw_samples$main)
    cbind(mains, pairs)
  }

  testthat::expect_equal(combine_chains(fit1), combine_chains(fit2))
})

# ==============================================================================
# GGM Reproducibility and Structure Tests
# ==============================================================================

test_that("bgm GGM is reproducible", {
  fit1 = get_bgms_fit_ggm()

  set.seed(42)
  x = matrix(rnorm(200), nrow = 50, ncol = 4)
  colnames(x) = paste0("V", 1:4)

  fit2 = bgm(
    x = x, variable_type = "continuous",
    edge_selection = TRUE,
    iter = 50, warmup = 100, chains = 1,
    seed = 44442, # Same seed as fixture
    display_progress = "none"
  )

  testthat::expect_equal(fit1$raw_samples$main, fit2$raw_samples$main)
  testthat::expect_equal(fit1$raw_samples$pairwise, fit2$raw_samples$pairwise)
})

test_that("bgm GGM output has correct dimensions", {
  fit = get_bgms_fit_ggm()
  args = extract_arguments(fit)
  p = args$num_variables

  # GGM has no main effects; precision diagonal is in quadratic
  expect_null(fit$posterior_summary_main)
  expect_null(fit$posterior_mean_main)
  expect_equal(nrow(fit$posterior_summary_quadratic), p)

  # pairwise: p*(p-1)/2 off-diagonal elements
  n_edges = p * (p - 1) / 2
  expect_equal(nrow(fit$posterior_summary_pairwise), n_edges)
  expect_equal(nrow(fit$posterior_mean_pairwise), p)
  expect_equal(ncol(fit$posterior_mean_pairwise), p)

  # precision diagonal stored separately (positive for GGM)
  expect_true(all(fit$posterior_mean_residual_variance > 0))

  # pairwise: p*(p-1)/2 off-diagonal elements
  n_edges = p * (p - 1) / 2
  expect_equal(nrow(fit$posterior_summary_pairwise), n_edges)
  expect_equal(nrow(fit$posterior_mean_pairwise), p)
  expect_equal(ncol(fit$posterior_mean_pairwise), p)

  # indicators (edge selection = TRUE)
  expect_equal(nrow(fit$posterior_summary_indicator), n_edges)
  expect_equal(nrow(fit$posterior_mean_indicator), p)
  expect_equal(ncol(fit$posterior_mean_indicator), p)

  # raw samples
  expect_equal(ncol(fit$raw_samples$main[[1]]), p)
  expect_equal(ncol(fit$raw_samples$pairwise[[1]]), n_edges)
  expect_equal(nrow(fit$raw_samples$main[[1]]), args$iter)
})

test_that("bgm GGM without edge selection omits indicators", {
  fit = get_bgms_fit_ggm_no_es()

  expect_s3_class(fit, "bgms")
  expect_null(fit$posterior_summary_indicator)
  expect_null(fit$posterior_mean_indicator)
})

test_that("bgm GGM posterior precision diagonals are positive", {
  fit = get_bgms_fit_ggm_no_es()

  expect_true(all(fit$posterior_summary_quadratic$mean > 0))
})

# ==============================================================================
# Parameter Ordering Tests (p >= 4 required to detect row/column-major bugs)
# ==============================================================================
#
# For a symmetric p x p matrix, the off-diagonal elements can be vectorized in
# row-major or column-major upper-triangle order. These orderings are identical
# for p <= 3. At p = 4, the 3rd element differs: row-major gives (1,4) while
# column-major gives (2,3). Using p = 5 ensures robust detection.
#
# See helper-fixtures.R for check_summary_matrix_consistency() and
# check_extractor_matrix_consistency().
# ==============================================================================

test_that("bgm GGM output has correct parameter ordering", {
  skip_on_cran()

  # Precision matrix with distinctive values at swap positions (p=5):
  #   Row-major position 3 = (1,4), column-major position 3 = (2,3)
  #   Row-major position 7 = (2,5), column-major position 7 = (3,4)
  # Zero vs non-zero at these positions makes any swap detectable.
  p = 5
  omega = diag(p) * 2
  omega[1, 2] = omega[2, 1] = 0.6
  omega[1, 3] = omega[3, 1] = -0.4
  omega[1, 4] = omega[4, 1] = 0.0 # swaps with V2-V3 under wrong ordering
  omega[1, 5] = omega[5, 1] = 0.3
  omega[2, 3] = omega[3, 2] = 0.0 # swaps with V1-V4 under wrong ordering
  omega[2, 4] = omega[4, 2] = -0.5
  omega[2, 5] = omega[5, 2] = 0.0 # swaps with V3-V4 under wrong ordering
  omega[3, 4] = omega[4, 3] = 0.25
  omega[3, 5] = omega[5, 3] = 0.0
  omega[4, 5] = omega[5, 4] = -0.35

  n = 1000
  x = simulate_mrf(
    num_states = n, num_variables = p, pairwise = omega,
    variable_type = "continuous", seed = 42
  )
  colnames(x) = paste0("V", 1:p)

  fit = bgm(
    x,
    variable_type = "continuous",
    iter = 500, warmup = 500, chains = 1,
    edge_selection = FALSE, seed = 42,
    display_progress = "none"
  )

  # Summary names -> matrix positions (pairwise)
  # GGM: summary and matrix are both on the association scale.
  expect_true(
    all(check_summary_matrix_consistency(
      fit$posterior_summary_pairwise,
      fit$posterior_mean_pairwise
    )),
    info = "GGM pairwise summary names do not match matrix positions"
  )

  # Extractor column means -> matrix positions
  pw_means = colMeans(extract_pairwise_interactions(fit))
  expect_true(
    all(check_extractor_matrix_consistency(
      pw_means, fit$posterior_mean_pairwise
    )),
    info = paste(
      "GGM extract_pairwise_interactions()",
      "names do not match matrix positions"
    )
  )

  # Truth-based swap-position checks (GGM stores association-scale: A = -0.5 * precision):
  # V1-V4 (true precision = 0) should be near zero, not ~-0.125 (V3-V4's value)
  expect_true(
    abs(fit$posterior_mean_pairwise["V1", "V4"]) < 0.15,
    info = sprintf(
      "V1-V4 should be ~0 but is %.3f (possible swap with V2-V3)",
      fit$posterior_mean_pairwise["V1", "V4"]
    )
  )
  # V2-V3 (true precision = 0) should be near zero, not ~-0.3 (V1-V2's value)
  expect_true(
    abs(fit$posterior_mean_pairwise["V2", "V3"]) < 0.15,
    info = sprintf(
      "V2-V3 should be ~0 but is %.3f (possible swap with V1-V4)",
      fit$posterior_mean_pairwise["V2", "V3"]
    )
  )
  # V3-V4 (true precision = 0.25, true association = -0.125) should be negative
  expect_true(
    fit$posterior_mean_pairwise["V3", "V4"] < -0.05,
    info = sprintf(
      "V3-V4 should be ~-0.125 but is %.3f (possible swap with V2-V5)",
      fit$posterior_mean_pairwise["V3", "V4"]
    )
  )
  # V2-V4 (true precision = -0.5, true association = 0.25) should be positive
  expect_true(
    fit$posterior_mean_pairwise["V2", "V4"] > 0.15,
    info = sprintf(
      "V2-V4 should be ~0.25 but is %.3f",
      fit$posterior_mean_pairwise["V2", "V4"]
    )
  )
})

test_that("bgm OMRF output has correct parameter ordering", {
  # Ordering is structural: the chain only has to produce the vectors.
  skip_on_cran()

  data("Wenchuan", package = "bgms")
  x = na.omit(Wenchuan[, 1:5]) # p=5 to detect row/column-major bugs

  fit = bgm(
    x,
    iter = 60, warmup = 60, chains = 1,
    edge_selection = TRUE, seed = 42,
    display_progress = "none"
  )

  # Summary names -> matrix positions (pairwise)
  expect_true(
    all(check_summary_matrix_consistency(
      fit$posterior_summary_pairwise,
      fit$posterior_mean_pairwise
    )),
    info = "OMRF pairwise summary names do not match matrix positions"
  )

  # Extractor column means -> matrix positions
  pw_means = colMeans(extract_pairwise_interactions(fit))
  expect_true(
    all(check_extractor_matrix_consistency(
      pw_means, fit$posterior_mean_pairwise
    )),
    info = paste(
      "OMRF extract_pairwise_interactions()",
      "names do not match matrix positions"
    )
  )

  # Indicator summary names -> matrix positions
  expect_true(
    all(check_summary_matrix_consistency(
      fit$posterior_summary_indicator,
      fit$posterior_mean_indicator
    )),
    info = "OMRF indicator summary names do not match matrix positions"
  )
})


# ==============================================================================
# GGM Expanded Test Suite
# ==============================================================================
#
# Tests for GGM correctness, convergence, and edge detection.
# ==============================================================================


# --- Multi-chain convergence ---------------------------------------------

test_that("bgm GGM multi-chain produces valid Rhat", {
  skip_on_cran()

  p = 5
  omega = diag(p) * 2
  omega[1, 2] = omega[2, 1] = 0.5
  omega[2, 3] = omega[3, 2] = -0.4
  omega[3, 4] = omega[4, 3] = 0.3

  x = simulate_mrf(
    num_states = 200, num_variables = p, pairwise = omega,
    variable_type = "continuous", seed = 101
  )
  colnames(x) = paste0("V", 1:p)

  fit = bgm(
    x,
    variable_type = "continuous",
    edge_selection = FALSE,
    iter = 1000, warmup = 500, chains = 2,
    seed = 202, display_progress = "none"
  )

  # All Rhat values should be below 1.1 for converged chains
  rhat_quad = fit$posterior_summary_quadratic$Rhat
  rhat_pair = fit$posterior_summary_pairwise$Rhat

  expect_true(
    all(rhat_quad < 1.1),
    info = sprintf("Max quadratic Rhat = %.3f (expected < 1.1)", max(rhat_quad))
  )
  expect_true(
    all(rhat_pair < 1.1),
    info = sprintf("Max pairwise Rhat = %.3f (expected < 1.1)", max(rhat_pair))
  )
})


# --- Sufficient statistics / MLE convergence -----------------------------

test_that("bgm GGM posterior mean approaches MLE for large n", {
  # For large n without edge selection, the posterior mean should approach
  # the sample precision matrix, since the likelihood dominates the prior.
  p = 4
  omega_true = diag(p)
  omega_true[1, 2] = omega_true[2, 1] = 0.4
  omega_true[3, 4] = omega_true[4, 3] = -0.3

  n = 500
  x = simulate_mrf(
    num_states = n, num_variables = p, pairwise = omega_true,
    variable_type = "continuous", seed = 42
  )
  colnames(x) = paste0("V", 1:p)

  # MLE of the precision matrix = solve of sample covariance (centered)
  x_centered = scale(x, center = TRUE, scale = FALSE)
  S = crossprod(x_centered) / n
  mle_precision = solve(S)

  fit = bgm(
    x,
    variable_type = "continuous",
    edge_selection = FALSE,
    iter = 500, warmup = 500, chains = 1,
    seed = 43, display_progress = "none"
  )

  # Reconstruct posterior mean precision (precision = -2 * association)
  omega_hat = extract_precision(fit)

  # Posterior mean should correlate highly with MLE (likelihood dominates)
  cor_offdiag = cor(
    omega_hat[lower.tri(omega_hat)],
    mle_precision[lower.tri(mle_precision)]
  )
  cor_diag = cor(diag(omega_hat), diag(mle_precision))

  expect_true(
    cor_offdiag > 0.95,
    info = sprintf(
      "Off-diagonal cor with MLE = %.3f (expected > 0.95)",
      cor_offdiag
    )
  )
  expect_true(
    cor_diag > 0.95,
    info = sprintf("Diagonal cor with MLE = %.3f (expected > 0.95)", cor_diag)
  )
})


# --- Missing data handling -----------------------------------------------

test_that("bgm GGM with listwise deletion drops rows correctly", {
  set.seed(42)
  x = matrix(rnorm(200), nrow = 50, ncol = 4)
  colnames(x) = paste0("V", 1:4)

  # Introduce NAs
  x[5, 2] = NA
  x[10, 3] = NA
  x[20, 1] = NA

  fit = bgm(
    x,
    variable_type = "continuous",
    na_action = "listwise",
    edge_selection = FALSE,
    iter = 50, warmup = 100, chains = 1,
    seed = 44, display_progress = "none"
  )

  # 3 rows removed -> n = 47
  expect_equal(extract_arguments(fit)$num_cases, 47L)
  expect_s3_class(fit, "bgms")
})

test_that("GGM with na_action='impute' runs without error", {
  set.seed(1)
  n = 20
  p = 3
  x = matrix(rnorm(n * p), n, p)
  colnames(x) = paste0("V", 1:p)
  x[c(1, 5, 10)] = NA

  fit = bgm(x,
    iter = 50, warmup = 50, chains = 1,
    variable_type = "continuous",
    na_action = "impute", display_progress = "none"
  )

  expect_s3_class(fit, "bgms")
  expect_true(extract_arguments(fit)$na_impute)
  expect_equal(nrow(fit$raw_samples$pairwise[[1]]), 50)
})

test_that("GGM imputation preserves posterior accuracy", {
  skip_on_cran()
  set.seed(42)

  p = 6
  n = 300
  # Known sparse precision matrix
  Omega_true = diag(p)
  Omega_true[1, 2] = Omega_true[2, 1] = 0.4
  Omega_true[3, 4] = Omega_true[4, 3] = -0.3

  Sigma = solve(Omega_true)
  L = chol(Sigma)
  x_full = matrix(rnorm(n * p), n, p) %*% L
  colnames(x_full) = paste0("V", 1:p)

  # Fit on complete data
  fit_full = bgm(x_full,
    iter = 800, edge_selection = FALSE,
    variable_type = "continuous",
    chains = 1, display_progress = "none"
  )

  # Introduce 5% MCAR missing data
  x_miss = x_full
  miss_idx = sample(length(x_miss), size = round(0.05 * length(x_miss)))
  x_miss[miss_idx] = NA

  fit_miss = bgm(x_miss,
    iter = 800, edge_selection = FALSE,
    variable_type = "continuous",
    na_action = "impute", chains = 1,
    display_progress = "none"
  )

  # Posterior means should be correlated > 0.85
  cor_pairwise = cor(
    as.numeric(fit_full$posterior_mean_pairwise),
    as.numeric(fit_miss$posterior_mean_pairwise)
  )
  expect_gt(cor_pairwise, 0.85)
})

test_that("GGM imputation gives comparable results to listwise", {
  skip_on_cran()
  set.seed(123)

  p = 5
  n = 200
  x = matrix(rnorm(n * p), n, p)
  colnames(x) = paste0("V", 1:p)

  # Introduce 3% missing
  miss_idx = sample(length(x), size = round(0.03 * length(x)))
  x[miss_idx] = NA

  fit_listwise = bgm(x,
    iter = 500, na_action = "listwise",
    variable_type = "continuous",
    edge_selection = FALSE, chains = 1,
    display_progress = "none"
  )
  fit_impute = bgm(x,
    iter = 500, na_action = "impute",
    variable_type = "continuous",
    edge_selection = FALSE, chains = 1,
    display_progress = "none"
  )

  expect_s3_class(fit_listwise, "bgms")
  expect_s3_class(fit_impute, "bgms")

  cor_val = cor(
    as.numeric(fit_listwise$posterior_mean_pairwise),
    as.numeric(fit_impute$posterior_mean_pairwise)
  )
  expect_gt(cor_val, 0.80)
})

test_that("GGM imputation handles edge cases", {
  set.seed(1)
  n = 50
  p = 4
  x = matrix(rnorm(n * p), n, p)
  colnames(x) = paste0("V", 1:p)

  # Single missing value
  x[1, 1] = NA
  fit = bgm(x,
    iter = 50, warmup = 50, chains = 1,
    variable_type = "continuous",
    na_action = "impute", display_progress = "none"
  )
  expect_s3_class(fit, "bgms")

  # Multiple missing in same row
  x[2, 1:3] = NA
  fit2 = bgm(x,
    iter = 50, warmup = 50, chains = 1,
    variable_type = "continuous",
    na_action = "impute", display_progress = "none"
  )
  expect_s3_class(fit2, "bgms")
})

test_that("GGM impute: posterior samples remain finite and bounded", {
  set.seed(1)
  x = matrix(rnorm(500), 100, 5)
  colnames(x) = paste0("V", 1:5)
  x[sample(500, 25)] = NA
  fit = bgm(x,
    iter = 50, warmup = 50, chains = 1,
    variable_type = "continuous",
    na_action = "impute", edge_selection = FALSE,
    display_progress = "none"
  )
  samples = fit$raw_samples$pairwise[[1]]
  expect_true(all(is.finite(samples)))
  expect_true(all(abs(samples) < 100))
})

test_that("GGM with na_action='impute' but no NAs works transparently", {
  x = matrix(rnorm(100), 20, 5)
  colnames(x) = paste0("V", 1:5)
  fit = bgm(x,
    iter = 50, warmup = 50, chains = 1,
    variable_type = "continuous",
    na_action = "impute", display_progress = "none"
  )
  expect_s3_class(fit, "bgms")
})

test_that("GGM impute: entire-column-missing gives clear error", {
  x = matrix(rnorm(100), 20, 5)
  colnames(x) = paste0("V", 1:5)
  x[, 3] = NA
  expect_error(
    bgm(x,
      iter = 50, variable_type = "continuous",
      na_action = "impute", display_progress = "none"
    ),
    "no observed values"
  )
})

test_that("bgm GGM stores the training column means for prediction", {
  fit = get_bgms_fit_ggm()
  args = extract_arguments(fit)

  # predict() centers newdata on these means (the GGM analogue of the ordinal
  # recode map in category_levels), so they must be stored on the fit.
  expect_length(args$column_means, args$num_variables)
  expect_true(all(is.finite(args$column_means)))
})


# --- Larger p (Cholesky stability) ---------------------------------------

test_that("bgm GGM with p = 15 produces valid output", {
  skip_on_cran()

  p = 15
  # Diagonally dominant precision matrix (ensures positive definiteness)
  omega = diag(p) * 3
  set.seed(77)
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      if(runif(1) < 0.2) { # ~20% non-zero edges
        val = runif(1, -0.3, 0.3)
        omega[i, j] = omega[j, i] = val
      }
    }
  }

  x = simulate_mrf(
    num_states = 300, num_variables = p, pairwise = omega,
    variable_type = "continuous", seed = 88
  )
  colnames(x) = paste0("V", 1:p)

  fit = bgm(
    x,
    variable_type = "continuous",
    edge_selection = FALSE,
    iter = 500, warmup = 500, chains = 1,
    seed = 99, display_progress = "none"
  )

  # All precision diagonals should be positive
  expect_true(
    all(fit$posterior_summary_quadratic$mean > 0),
    info = "Some diagonal precision elements are non-positive"
  )

  # All values should be finite
  expect_true(all(is.finite(fit$posterior_summary_quadratic$mean)))
  expect_true(all(is.finite(fit$posterior_summary_pairwise$mean)))

  # Correct dimensions
  n_edges = p * (p - 1) / 2
  expect_equal(nrow(fit$posterior_summary_quadratic), p)
  expect_equal(nrow(fit$posterior_summary_pairwise), n_edges)
})


# --- Edge detection power ------------------------------------------------

test_that("bgm GGM edge selection discriminates true edges", {
  skip_on_cran()

  p = 6
  omega_true = diag(p) * 2
  # 5 true edges (sparse graph)
  omega_true[1, 2] = omega_true[2, 1] = 0.6
  omega_true[2, 3] = omega_true[3, 2] = -0.5
  omega_true[3, 4] = omega_true[4, 3] = 0.4
  omega_true[4, 5] = omega_true[5, 4] = -0.5
  omega_true[5, 6] = omega_true[6, 5] = 0.3

  # True adjacency
  adj_true = (omega_true != 0)
  diag(adj_true) = FALSE
  true_edges = adj_true[lower.tri(adj_true)]

  x = simulate_mrf(
    num_states = 500, num_variables = p, pairwise = omega_true,
    variable_type = "continuous", seed = 321
  )
  colnames(x) = paste0("V", 1:p)

  fit = bgm(
    x,
    variable_type = "continuous",
    edge_selection = TRUE,
    iter = 1000, warmup = 500, chains = 2,
    seed = 654, display_progress = "none"
  )

  # Posterior inclusion probabilities
  pip = fit$posterior_mean_indicator[lower.tri(fit$posterior_mean_indicator)]

  # True edges should have higher PIPs than non-edges on average
  mean_pip_true = mean(pip[true_edges])
  mean_pip_false = mean(pip[!true_edges])

  expect_true(
    mean_pip_true > mean_pip_false,
    info = sprintf(
      "Mean PIP for true edges (%.3f) should exceed non-edges (%.3f)",
      mean_pip_true, mean_pip_false
    )
  )

  # True edges should mostly have PIP > 0.5
  expect_true(
    mean(pip[true_edges] > 0.5) >= 0.6,
    info = sprintf(
      "Only %.0f%% of true edges have PIP > 0.5 (expected >= 60%%)",
      100 * mean(pip[true_edges] > 0.5)
    )
  )

  # Non-edges should mostly have PIP < 0.5
  expect_true(
    mean(pip[!true_edges] < 0.5) >= 0.6,
    info = sprintf(
      "Only %.0f%% of non-edges have PIP < 0.5 (expected >= 60%%)",
      100 * mean(pip[!true_edges] < 0.5)
    )
  )
})


# --- Conditional regression check ----------------------------------------

# ==============================================================================
# Mixed MRF End-to-End Tests
# ==============================================================================

test_that("bgm mixed MRF is reproducible", {
  fit1 = get_bgms_fit_mixed_mrf()

  set.seed(99)
  n = 80
  x = cbind(
    sample(0:2, n, replace = TRUE),
    rnorm(n),
    sample(0:2, n, replace = TRUE),
    rnorm(n),
    sample(0:2, n, replace = TRUE)
  )
  colnames(x) = c("d1", "c1", "d2", "c2", "d3")

  fit2 = bgm(
    x = x,
    variable_type = c(
      "ordinal", "continuous", "ordinal",
      "continuous", "ordinal"
    ),
    edge_selection = TRUE,
    iter = 50, warmup = 100, chains = 1,
    seed = 77771,
    display_progress = "none"
  )

  testthat::expect_equal(fit1$raw_samples$main, fit2$raw_samples$main)
  testthat::expect_equal(fit1$raw_samples$pairwise, fit2$raw_samples$pairwise)
})

test_that("bgm mixed MRF output has correct dimensions", {
  fit = get_bgms_fit_mixed_mrf()
  args = extract_arguments(fit)
  p_total = args$num_variables # 5
  p = 3L # discrete
  q = 2L # continuous

  # pairwise: p_total*(p_total-1)/2 edges
  n_edges = p_total * (p_total - 1) / 2
  expect_equal(nrow(fit$posterior_summary_pairwise), n_edges)
  expect_equal(nrow(fit$posterior_mean_pairwise), p_total)
  expect_equal(ncol(fit$posterior_mean_pairwise), p_total)

  # indicators (edge selection = TRUE)
  expect_equal(nrow(fit$posterior_summary_indicator), n_edges)
  expect_equal(nrow(fit$posterior_mean_indicator), p_total)
  expect_equal(ncol(fit$posterior_mean_indicator), p_total)

  # posterior_mean_main: list with discrete and continuous
  expect_true(is.list(fit$posterior_mean_main))
  expect_equal(nrow(fit$posterior_mean_main$discrete), p)
  expect_equal(nrow(fit$posterior_mean_main$continuous), q)
  expect_equal(ncol(fit$posterior_mean_main$continuous), 1) # mean only

  # raw samples
  expect_equal(ncol(fit$raw_samples$pairwise[[1]]), n_edges)
  expect_equal(nrow(fit$raw_samples$main[[1]]), args$iter)
})

test_that("bgm mixed MRF without edge selection omits indicators", {
  fit = get_bgms_fit_mixed_mrf_no_es()

  expect_s3_class(fit, "bgms")
  expect_null(fit$posterior_summary_indicator)
  expect_null(fit$posterior_mean_indicator)
})

test_that("bgm mixed MRF pairwise matrix has correct variable names", {
  fit = get_bgms_fit_mixed_mrf()

  # Interleaved order: d1, c1, d2, c2, d3
  expected_names = c("d1", "c1", "d2", "c2", "d3")
  expect_equal(rownames(fit$posterior_mean_pairwise), expected_names)
  expect_equal(colnames(fit$posterior_mean_pairwise), expected_names)
  expect_equal(rownames(fit$posterior_mean_indicator), expected_names)
  expect_equal(colnames(fit$posterior_mean_indicator), expected_names)
})

test_that("bgm mixed MRF pairwise matrix is symmetric", {
  fit = get_bgms_fit_mixed_mrf()
  expect_equal(fit$posterior_mean_pairwise, t(fit$posterior_mean_pairwise))
  expect_equal(fit$posterior_mean_indicator, t(fit$posterior_mean_indicator))
})

test_that("bgm mixed MRF summary-matrix consistency", {
  fit = get_bgms_fit_mixed_mrf()
  expect_true(
    all(check_summary_matrix_consistency(
      fit$posterior_summary_pairwise,
      fit$posterior_mean_pairwise
    )),
    info = "Mixed MRF pairwise summary names do not match matrix positions"
  )
  expect_true(
    all(check_summary_matrix_consistency(
      fit$posterior_summary_indicator,
      fit$posterior_mean_indicator
    )),
    info = "Mixed MRF indicator summary names do not match matrix positions"
  )
})

test_that("bgm mixed MRF residual variances are positive", {
  fit = get_bgms_fit_mixed_mrf_no_es()
  args = extract_arguments(fit)
  expect_true(all(fit$posterior_mean_residual_variance > 0))
})

test_that("bgm mixed MRF (no edge selection) runs", {
  fit = get_bgms_fit_mixed_mrf_marginal()
  expect_s3_class(fit, "bgms")
  expect_equal(nrow(fit$posterior_mean_pairwise), 5)
  expect_true(all(is.finite(fit$posterior_mean_pairwise)))
})

test_that("bgm mixed MRF with edge selection runs", {
  fit = get_bgms_fit_mixed_mrf_marginal_es()
  expect_s3_class(fit, "bgms")
  expect_equal(nrow(fit$posterior_mean_pairwise), 5)
  expect_equal(ncol(fit$posterior_mean_pairwise), 5)
  expect_true(all(is.finite(fit$posterior_mean_pairwise)))
  # Edge selection produces indicator matrix
  expect_false(is.null(fit$posterior_mean_indicator))
  expect_equal(nrow(fit$posterior_mean_indicator), 5)
  expect_true(all(fit$posterior_mean_indicator >= 0 &
    fit$posterior_mean_indicator <= 1))
})

test_that("bgm mixed MRF NUTS runs", {
  fit = get_bgms_fit_mixed_mrf_nuts()
  expect_s3_class(fit, "bgms")
  expect_equal(nrow(fit$posterior_mean_pairwise), 5)
  expect_equal(ncol(fit$posterior_mean_pairwise), 5)
  expect_true(all(is.finite(fit$posterior_mean_pairwise)))
  # Edge selection active
  expect_false(is.null(fit$posterior_mean_indicator))
})

test_that("bgm mixed MRF NUTS output dimensions", {
  fit = get_bgms_fit_mixed_mrf_nuts()
  args = extract_arguments(fit)
  p_total = args$num_variables
  n_edges = p_total * (p_total - 1) / 2

  expect_equal(nrow(fit$posterior_summary_pairwise), n_edges)
  expect_equal(nrow(fit$posterior_mean_pairwise), p_total)
  expect_equal(ncol(fit$posterior_mean_pairwise), p_total)
  expect_equal(nrow(fit$posterior_summary_indicator), n_edges)
  expect_equal(ncol(fit$raw_samples$pairwise[[1]]), n_edges)
  expect_equal(nrow(fit$raw_samples$main[[1]]), args$iter)
})

test_that("bgm mixed MRF NUTS is reproducible", {
  fit1 = get_bgms_fit_mixed_mrf_nuts()

  set.seed(99)
  n = 80
  x = cbind(
    sample(0:2, n, replace = TRUE),
    rnorm(n),
    sample(0:2, n, replace = TRUE),
    rnorm(n),
    sample(0:2, n, replace = TRUE)
  )
  colnames(x) = c("d1", "c1", "d2", "c2", "d3")

  fit2 = bgm(
    x = x,
    variable_type = c(
      "ordinal", "continuous", "ordinal",
      "continuous", "ordinal"
    ),
    edge_selection = TRUE,
    update_method = "nuts",
    iter = 50, warmup = 100, chains = 1,
    seed = 77775,
    display_progress = "none"
  )

  expect_equal(fit1$raw_samples$main, fit2$raw_samples$main)
  expect_equal(fit1$raw_samples$pairwise, fit2$raw_samples$pairwise)
})

test_that("bgm mixed MRF NUTS without edge selection runs", {
  fit = get_bgms_fit_mixed_mrf_nuts_no_es()
  expect_s3_class(fit, "bgms")
  expect_equal(nrow(fit$posterior_mean_pairwise), 5)
  expect_true(all(is.finite(fit$posterior_mean_pairwise)))
  expect_null(fit$posterior_summary_indicator)
  expect_null(fit$posterior_mean_indicator)
})

test_that("bgm mixed MRF Beta-Bernoulli prior runs", {
  fit = get_bgms_fit_mixed_mrf_beta_bernoulli()
  expect_s3_class(fit, "bgms")
  expect_equal(nrow(fit$posterior_mean_pairwise), 5)
  expect_true(all(is.finite(fit$posterior_mean_pairwise)))
  expect_false(is.null(fit$posterior_mean_indicator))
  expect_true(all(fit$posterior_mean_indicator >= 0 &
    fit$posterior_mean_indicator <= 1))
})

test_that("bgm mixed MRF Stochastic-Block prior runs", {
  fit = get_bgms_fit_mixed_mrf_sbm()
  expect_s3_class(fit, "bgms")
  expect_equal(nrow(fit$posterior_mean_pairwise), 5)
  expect_true(all(is.finite(fit$posterior_mean_pairwise)))
  expect_false(is.null(fit$posterior_mean_indicator))
  expect_true(all(fit$posterior_mean_indicator >= 0 &
    fit$posterior_mean_indicator <= 1))
})

test_that("bgm mixed MRF Blume-Capel + continuous runs", {
  fit = get_bgms_fit_mixed_mrf_bc()
  expect_s3_class(fit, "bgms")
  p_total = 4 # bc1, c1, bc2, c2
  expect_equal(nrow(fit$posterior_mean_pairwise), p_total)
  expect_equal(ncol(fit$posterior_mean_pairwise), p_total)
  expect_true(all(is.finite(fit$posterior_mean_pairwise)))
  expect_false(is.null(fit$posterior_mean_indicator))

  # Blume-Capel main effects: quadratic structure (linear + quadratic terms)
  expect_true(is.list(fit$posterior_mean_main))
  expect_true("discrete" %in% names(fit$posterior_mean_main))
  expect_true("continuous" %in% names(fit$posterior_mean_main))
})

test_that("bgm mixed MRF imputation runs", {
  fit = get_bgms_fit_mixed_mrf_impute()
  expect_s3_class(fit, "bgms")
  expect_equal(nrow(fit$posterior_mean_pairwise), 5)
  expect_true(all(is.finite(fit$posterior_mean_pairwise)))
  expect_false(is.null(fit$posterior_mean_indicator))
})

test_that("bgm mixed MRF multi-chain R-hat and ESS", {
  fit = get_bgms_fit_mixed_mrf_multichain()
  expect_s3_class(fit, "bgms")
  expect_equal(nrow(fit$posterior_mean_pairwise), 5)
  expect_true(all(is.finite(fit$posterior_mean_pairwise)))

  # Multi-chain produces multiple raw sample chains
  expect_equal(length(fit$raw_samples$pairwise), 2)
  expect_equal(length(fit$raw_samples$main), 2)

  # R-hat and ESS should be computable
  rhat = extract_rhat(fit)
  expect_true(is.numeric(rhat$pairwise))
  expect_true(all(is.na(rhat$pairwise) | rhat$pairwise > 0))
  ess = extract_ess(fit)
  expect_true(is.numeric(ess$pairwise))
  expect_true(all(is.na(ess$pairwise) | ess$pairwise > 0))
})

test_that("bgm mixed MRF output has correct parameter ordering", {
  skip_on_cran()

  # 5 interleaved variables: d1, c1, d2, c2, d3 (p=3 discrete, q=2 continuous)
  # Internal C++ order: d1, d2, d3, c1, c2
  #
  # For a 5x5 upper triangle, position 2 differs between orderings:
  #   Row-major position 2 = (1,4) = d1-c2
  #   Col-major position 2 = (2,3) = c1-d2
  # Strategic zeros make any swap detectable.
  p = 3L
  q = 2L
  n = 500L

  # Parameters in internal (dd/cc/dc block) order
  pairwise_disc = matrix(c(
    0, -0.2, 0.1,
    -0.2, 0, 0.0,
    0.1, 0.0, 0
  ), p, p, byrow = TRUE)

  pairwise_cross = matrix(c(
    0.3,  0.0, # d1-c1 = 0.3, d1-c2 = 0.0 (swap sentinel)
    0.5,  0.3, # d2-c1 = 0.5 (swap sentinel), d2-c2 = 0.3
    -0.3, 0.15 # d3-c1 = -0.3, d3-c2 = 0.15
  ), p, q, byrow = TRUE)

  pairwise_cont = diag(c(-0.75, -1.0))
  pairwise_cont[1, 2] = pairwise_cont[2, 1] = 0.0 # c1-c2 = 0 (swap sentinel)

  nc = c(2L, 2L, 2L)
  mux = matrix(0, p, max(nc) + 1)
  muy = rep(0, q)

  result = sample_mixed_mrf_gibbs(
    num_states = n, pairwise_disc_r = pairwise_disc,
    pairwise_cross_r = pairwise_cross, pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy, num_categories_r = nc,
    variable_type_r = rep("ordinal", p),
    baseline_category_r = rep(0L, p), iter = 500L, seed = 42L
  )

  # Reassemble in user (interleaved) order: d1, c1, d2, c2, d3
  x = data.frame(
    d1 = result$x[, 1],
    c1 = result$y[, 1],
    d2 = result$x[, 2],
    c2 = result$y[, 2],
    d3 = result$x[, 3]
  )

  fit = bgm(
    x,
    variable_type = c(
      "ordinal", "continuous", "ordinal",
      "continuous", "ordinal"
    ),
    iter = 400, warmup = 300, chains = 1,
    edge_selection = FALSE, seed = 42,
    display_progress = "none"
  )

  # Extractor column means -> matrix positions
  pw_means = colMeans(extract_pairwise_interactions(fit))
  expect_true(
    all(check_extractor_matrix_consistency(
      pw_means, fit$posterior_mean_pairwise
    )),
    info = paste(
      "Mixed MRF extract_pairwise_interactions()",
      "names do not match matrix positions"
    )
  )

  # Truth-based swap checks (user-order variable names):
  # d1-c2 (true = 0.0) should be near zero, not ~0.5 (d2-c1's value)
  expect_true(
    abs(fit$posterior_mean_pairwise["d1", "c2"]) < 0.2,
    info = sprintf(
      "d1-c2 should be ~0 but is %.3f (possible swap with c1-d2)",
      fit$posterior_mean_pairwise["d1", "c2"]
    )
  )
  # c1-d2 (true = 0.5) should be clearly positive, not ~0 (d1-c2's value)
  expect_true(
    fit$posterior_mean_pairwise["c1", "d2"] > 0.15,
    info = sprintf(
      "c1-d2 should be ~0.5 but is %.3f (possible swap with d1-c2)",
      fit$posterior_mean_pairwise["c1", "d2"]
    )
  )
  # c1-c2 (true = 0.0) should be near zero, not ~0.3 (d2-c2's value)
  expect_true(
    abs(fit$posterior_mean_pairwise["c1", "c2"]) < 0.2,
    info = sprintf(
      "c1-c2 should be ~0 but is %.3f (possible swap with d2-c2)",
      fit$posterior_mean_pairwise["c1", "c2"]
    )
  )
  # d2-d3 (true = 0.0) should be near zero
  expect_true(
    abs(fit$posterior_mean_pairwise["d2", "d3"]) < 0.2,
    info = sprintf(
      "d2-d3 should be ~0 but is %.3f",
      fit$posterior_mean_pairwise["d2", "d3"]
    )
  )
  # d1-d2 (true = -0.2) should be negative
  expect_true(
    fit$posterior_mean_pairwise["d1", "d2"] < -0.08,
    info = sprintf(
      "d1-d2 should be ~-0.2 but is %.3f",
      fit$posterior_mean_pairwise["d1", "d2"]
    )
  )
})


test_that("bgm GGM implied regression matches OLS for large n", {
  skip_on_cran()

  p = 5
  omega_true = diag(p) * 2
  omega_true[1, 2] = omega_true[2, 1] = 0.5
  omega_true[2, 3] = omega_true[3, 2] = -0.4
  omega_true[3, 4] = omega_true[4, 3] = 0.3
  omega_true[1, 5] = omega_true[5, 1] = 0.2

  n = 1000
  x = simulate_mrf(
    num_states = n, num_variables = p, pairwise = omega_true,
    variable_type = "continuous", seed = 111
  )
  colnames(x) = paste0("V", 1:p)

  fit = bgm(
    x,
    variable_type = "continuous",
    edge_selection = FALSE,
    iter = 800, warmup = 400, chains = 2,
    seed = 222, display_progress = "none"
  )

  # Reconstruct posterior mean precision matrix
  omega_hat = extract_precision(fit)

  # For each variable j, the implied regression coefficients are:
  #   beta_j = -omega_{j,-j} / omega_{jj}
  # This should match OLS coefficients from the data.
  x_centered = sweep(x, 2, colMeans(x))

  for(j in seq_len(p)) {
    rest = setdiff(seq_len(p), j)

    # Implied regression from precision matrix
    beta_implied = -omega_hat[rest, j] / omega_hat[j, j]

    # OLS regression on centered data
    ols_fit = lm(x_centered[, j] ~ x_centered[, rest] - 1)
    beta_ols = coef(ols_fit)

    # For large n, these should agree closely
    expect_true(
      cor(beta_implied, beta_ols) > 0.95,
      info = sprintf(
        "Variable %d: cor(beta_implied, beta_ols) = %.3f (expected > 0.95)",
        j, cor(beta_implied, beta_ols)
      )
    )
  }
})


# ==============================================================================
# Estimate-Simulate-Re-estimate Cycle Tests
# ==============================================================================
# Verify self-consistency: fit -> simulate -> re-fit -> compare.
# Posterior mean parameters from the re-fit should correlate with the original.

test_that("estimate-simulate-re-estimate cycle recovers parameters (OMRF)", {
  skip_unless_slow()
  skip_on_cran()

  data("Wenchuan", package = "bgms")
  fit1 = bgm(Wenchuan[1:100, 1:4],
    iter = 800, warmup = 400,
    edge_selection = FALSE, chains = 1, display_progress = "none"
  )
  sim = simulate(fit1, nsim = 200, method = "posterior-mean")
  fit2 = bgm(sim,
    iter = 800, warmup = 400,
    edge_selection = FALSE, chains = 1, display_progress = "none"
  )

  cor_pw = cor(
    as.numeric(fit1$posterior_mean_pairwise),
    as.numeric(fit2$posterior_mean_pairwise)
  )
  expect_gt(cor_pw, 0.7)
})

test_that("estimate-simulate-re-estimate cycle recovers parameters (GGM)", {
  skip_on_cran()

  set.seed(42)
  x = matrix(rnorm(400), nrow = 100, ncol = 4)
  colnames(x) = paste0("V", 1:4)

  fit1 = bgm(x,
    variable_type = "continuous",
    edge_selection = FALSE, iter = 800, warmup = 400,
    chains = 1, display_progress = "none"
  )
  sim = simulate(fit1, nsim = 200, method = "posterior-mean")
  fit2 = bgm(sim,
    variable_type = "continuous",
    edge_selection = FALSE, iter = 800, warmup = 400,
    chains = 1, display_progress = "none"
  )

  cor_pw = cor(
    as.numeric(fit1$posterior_mean_pairwise),
    as.numeric(fit2$posterior_mean_pairwise)
  )
  expect_gt(cor_pw, 0.7)
})

test_that("estimate-simulate-re-estimate cycle recovers parameters (mixed MRF)", {
  skip_on_cran()

  set.seed(99)
  n = 100
  x = data.frame(
    d1 = sample(0:2, n, replace = TRUE),
    c1 = rnorm(n),
    d2 = sample(0:2, n, replace = TRUE),
    c2 = rnorm(n)
  )
  vtypes = c("ordinal", "continuous", "ordinal", "continuous")

  fit1 = bgm(x,
    variable_type = vtypes,
    edge_selection = FALSE, iter = 800, warmup = 400,
    chains = 1, display_progress = "none"
  )
  sim = simulate(fit1, nsim = 200, method = "posterior-mean")
  fit2 = bgm(sim,
    variable_type = vtypes,
    edge_selection = FALSE, iter = 800, warmup = 400,
    chains = 1, display_progress = "none"
  )

  cor_pw = cor(
    as.numeric(fit1$posterior_mean_pairwise),
    as.numeric(fit2$posterior_mean_pairwise)
  )
  expect_gt(cor_pw, 0.7)
})
