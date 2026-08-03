# ==============================================================================
# Test Fixtures and Utilities for bgms Test Suite
#
# This file is automatically loaded before tests by testthat (files matching
# helper-*.R are sourced alphabetically before test files).
#
# Contents:
#   1. Pre-fitted model fixtures (loaded from inst/extdata)
#   2. Test data generators
#   3. Matrix validation helpers
#   4. Contract testing utilities
#   5. MCMC test helpers
#
# ==============================================================================
# TESTING PHILOSOPHY (see test-tolerance.R for the foundational approach)
# ==============================================================================
#
# All tests in this suite build on the "stochastic-robust" testing approach
# established in test-tolerance.R. Because bgms outputs are stochastic (MCMC),
# we avoid exact-value assertions and instead test:
#
#   1. RANGE INVARIANTS - Values within valid bounds
#      - Probabilities in [0, 1]
#      - Indicators are binary (0 or 1)
#      - Category predictions within valid range
#
#   2. SYMMETRY - Pairwise matrices should be symmetric
#      - Posterior mean pairwise effects
#      - Posterior inclusion probabilities
#
#   3. DIMENSION CONSISTENCY - Correct matrix sizes
#      - p x p for pairwise matrices
#      - p*(p-1)/2 edges for vectorized parameters
#      - n x p for simulated/predicted data
#
#   4. STRUCTURAL CONTRACTS - API stability for downstream packages
#      - Required fields in output objects
#      - Return types and structures
#
#   5. COARSE AGGREGATES (used sparingly) - Wide bounds on summary statistics
#      - Mean absolute interaction within [0, 0.8]
#      - Finite values where expected
#
# Helper functions below (is_symmetric, values_in_range, etc.) implement
# these testing patterns for reuse across all test files.
#
# ==============================================================================

# Ensure bgms package is loaded
library(bgms)

# Advisory output is quieted for the test run in setup.R, not here. Helper files
# are sourced by devtools::load_all() (setup files are not), so setting the
# option here would leak bgms.verbose = FALSE into interactive development
# sessions and silence fit-time messages and progress bars.

# ------------------------------------------------------------------------------
# 1. Session-Cached Model Fixtures
# ------------------------------------------------------------------------------
# These fixtures are computed once per test session using current code.
# This avoids stale RDS files while minimizing overhead.

.test_cache = new.env(parent = emptyenv())

#' @description Get cached bgms fit
#' (4 binary variables, edge selection, 2 chains)
get_bgms_fit = function() {
  if(is.null(.test_cache$bgms_fit)) {
    data("ADHD", package = "bgms")
    .test_cache$bgms_fit = bgm(
      ADHD[1:50, 2:5], # 4 binary symptom variables
      iter = 50, warmup = 100, chains = 2,
      seed = 12345,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit
}

#' @description Get cached bgms fit
#' (4 ordinal variables, edge selection, 2 chains)
get_bgms_fit_ordinal = function() {
  if(is.null(.test_cache$bgms_fit_ordinal)) {
    data("Wenchuan", package = "bgms")
    .test_cache$bgms_fit_ordinal = bgm(
      Wenchuan[1:50, 1:4], # 4 ordinal variables (0-4 scale)
      iter = 50, warmup = 100, chains = 2,
      seed = 12345,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_ordinal
}

#' @description Get cached bgmCompare fit
#' (4 binary variables, 2 groups, 2 chains)
get_bgmcompare_fit = function() {
  if(is.null(.test_cache$bgmcompare_fit)) {
    data("ADHD", package = "bgms")
    .test_cache$bgmcompare_fit = without_support_warning(bgmCompare(
        x = ADHD[, 2:5], # 4 binary symptom variables, full dataset
        group_indicator = ADHD[, "group"], # ADHD diagnosis group
        iter = 50, warmup = 100, chains = 2,
        seed = 54321,
        display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit
}

#' @description Get cached bgmCompare fit
#' using x,y interface (4 ordinal variables, 2 chains)
get_bgmcompare_fit_xy = function() {
  if(is.null(.test_cache$bgmcompare_fit_xy)) {
    data("Wenchuan", package = "bgms")
    x = Wenchuan[1:25, 1:4]
    y = Wenchuan[26:50, 1:4]
    .test_cache$bgmcompare_fit_xy = without_support_warning(bgmCompare(
        x = x, y = y,
        iter = 50, warmup = 100, chains = 2,
        seed = 1234,
        display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit_xy
}

#' @description Get cached bgmCompare fit
#' (4 ordinal variables, 2 groups, 2 chains)
get_bgmcompare_fit_ordinal = function() {
  if(is.null(.test_cache$bgmcompare_fit_ordinal)) {
    data("Wenchuan", package = "bgms")
    x = Wenchuan[1:50, 1:4] # 4 ordinal variables
    group_ind = rep(1:2, each = 25)
    .test_cache$bgmcompare_fit_ordinal = without_support_warning(bgmCompare(
        x = x, group_indicator = group_ind,
        iter = 50, warmup = 100, chains = 2,
        seed = 54321,
        display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit_ordinal
}

#' @description Get cached bgms fit with Blume-Capel variables (2 chains)
get_bgms_fit_blumecapel = function() {
  if(is.null(.test_cache$bgms_fit_blumecapel)) {
    data("Wenchuan", package = "bgms")
    .test_cache$bgms_fit_blumecapel = bgm(
      Wenchuan[1:50, 1:4], # 4 ordinal variables treated as Blume-Capel
      variable_type = "blume-capel",
      baseline_category = 2, # Middle category as baseline
      iter = 50, warmup = 100, chains = 2,
      seed = 11111,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_blumecapel
}

#' @description Get cached bgms fit with
#' single chain (for R-hat edge case testing)
get_bgms_fit_single_chain = function() {
  if(is.null(.test_cache$bgms_fit_single)) {
    data("ADHD", package = "bgms")
    .test_cache$bgms_fit_single = bgm(
      ADHD[1:50, 2:5],
      iter = 50, warmup = 100, chains = 1,
      seed = 99999,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_single
}

#' @description Get cached bgms fit using adaptive-metropolis sampler
get_bgms_fit_adaptive_metropolis = function() {
  if(is.null(.test_cache$bgms_fit_am)) {
    data("ADHD", package = "bgms")
    .test_cache$bgms_fit_am = bgm(
      ADHD[1:50, 2:5],
      update_method = "adaptive-metropolis",
      iter = 50, warmup = 100, chains = 2,
      seed = 77777,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_am
}

#' @description Get cached bgmCompare fit using adaptive-metropolis sampler
get_bgmcompare_fit_adaptive_metropolis = function() {
  if(is.null(.test_cache$bgmcompare_fit_am)) {
    data("ADHD", package = "bgms")
    .test_cache$bgmcompare_fit_am = without_support_warning(bgmCompare(
        x = ADHD[, 2:5],
        group_indicator = ADHD[, "group"],
        update_method = "adaptive-metropolis",
        iter = 50, warmup = 100, chains = 2,
        seed = 88888,
        display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit_am
}

#' @description Get cached bgmCompare fit with
#' main_difference_selection = TRUE +
#' Blume-Capel (1 chain)
#' Crosses Blume-Capel with difference_selection (Bernoulli prior)
get_bgmcompare_fit_main_selection = function() {
  if(is.null(.test_cache$bgmcompare_fit_main_sel)) {
    data("Boredom", package = "bgms")
    # Select 25 rows from each language group
    rows = c(1:25, 491:515)
    lang = Boredom[rows, "language"]
    .test_cache$bgmcompare_fit_main_sel = without_support_warning(bgmCompare(
        x = Boredom[rows, 2:5], # 4 ordinal variables (7 categories)
        group_indicator = lang,
        difference_selection = TRUE,
        main_difference_selection = TRUE,
        variable_type = "blume-capel",
        baseline_category = 3,
        iter = 25, warmup = 50, chains = 1,
        seed = 44444,
        display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit_main_sel
}

#' @description Get cached bgmCompare fit with
#' Beta-Bernoulli difference prior +
#' ordinal (1 chain)
#' Crosses Beta-Bernoulli prior with ordinal variables
get_bgmcompare_fit_beta_bernoulli = function() {
  if(is.null(.test_cache$bgmcompare_fit_bb)) {
    data("Wenchuan", package = "bgms")
    x = Wenchuan[1:25, 1:4]
    y = Wenchuan[26:50, 1:4]
    .test_cache$bgmcompare_fit_bb = without_support_warning(bgmCompare(
        x = x, y = y,
        difference_selection = TRUE,
        main_difference_selection = TRUE,
        difference_prior = beta_bernoulli_prior(alpha = 1, beta = 4),
        iter = 25, warmup = 50, chains = 1,
        seed = 55555,
        display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit_bb
}

#' @description Get cached bgms fit with Beta-Bernoulli edge prior (2 chains)
get_bgms_fit_beta_bernoulli = function() {
  if(is.null(.test_cache$bgms_fit_bb)) {
    data("ADHD", package = "bgms")
    .test_cache$bgms_fit_bb = bgm(
      ADHD[1:50, 2:5], # 4 binary symptom variables
      edge_prior = beta_bernoulli_prior(alpha = 1, beta = 4),
      iter = 50, warmup = 100, chains = 2,
      seed = 22222,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_bb
}

#' @description Get cached bgms fit with
#' Stochastic-Block Model edge prior (2 chains)
get_bgms_fit_sbm = function() {
  if(is.null(.test_cache$bgms_fit_sbm)) {
    data("ADHD", package = "bgms")
    .test_cache$bgms_fit_sbm = bgm(
      ADHD[1:50, 2:5], # 4 binary symptom variables
      edge_prior = sbm_prior(alpha = 1, beta = 1, dirichlet_alpha = 1),
      iter = 50, warmup = 100, chains = 2,
      seed = 33333,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_sbm
}

#' @description Get cached bgms fit with
#' adaptive-metropolis + Blume-Capel (1 chain)
get_bgms_fit_am_blumecapel = function() {
  if(is.null(.test_cache$bgms_fit_am_bc)) {
    data("Wenchuan", package = "bgms")
    .test_cache$bgms_fit_am_bc = bgm(
      Wenchuan[1:50, 1:4],
      update_method = "adaptive-metropolis",
      variable_type = "blume-capel",
      baseline_category = 1,
      iter = 25, warmup = 50, chains = 1,
      seed = 66666,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_am_bc
}

#' @description Get cached bgms fit with missing data imputation (1 chain)
get_bgms_fit_impute = function() {
  if(is.null(.test_cache$bgms_fit_impute)) {
    data("Wenchuan", package = "bgms")
    x = Wenchuan[1:50, 1:4]
    x[5, 2] = NA
    x[10, 3] = NA
    .test_cache$bgms_fit_impute = bgm(
      x,
      na_action = "impute",
      iter = 25, warmup = 50, chains = 1,
      seed = 77771,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_impute
}

#' @description Get cached bgmCompare fit with Blume-Capel variables (1 chain)
get_bgmcompare_fit_blumecapel = function() {
  if(is.null(.test_cache$bgmcompare_fit_bc)) {
    data("Boredom", package = "bgms")
    # Select 25 rows from each language group
    rows = c(1:25, 491:515)
    # The character column is a valid indicator; groups number by first appearance.
    lang = Boredom[rows, "language"]
    .test_cache$bgmcompare_fit_bc = without_support_warning(bgmCompare(
        x = Boredom[rows, 2:5], # 4 ordinal variables (7 categories)
        group_indicator = lang,
        variable_type = "blume-capel",
        baseline_category = 3,
        iter = 25, warmup = 50, chains = 1,
        seed = 99991,
        display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit_bc
}

#' @description Get cached bgmCompare fit
#' with adaptive-metropolis + Blume-Capel
#' (1 chain)
get_bgmcompare_fit_am_blumecapel = function() {
  if(is.null(.test_cache$bgmcompare_fit_am_bc)) {
    data("Boredom", package = "bgms")
    # Select 25 rows from each language group
    rows = c(1:25, 491:515)
    # The character column is a valid indicator; groups number by first appearance.
    lang = Boredom[rows, "language"]
    .test_cache$bgmcompare_fit_am_bc = without_support_warning(bgmCompare(
        x = Boredom[rows, 2:5], # 4 ordinal variables (7 categories)
        group_indicator = lang,
        update_method = "adaptive-metropolis",
        variable_type = "blume-capel",
        baseline_category = 3,
        iter = 25, warmup = 50, chains = 1,
        seed = 99992,
        display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit_am_bc
}

#' @description Get cached bgmCompare fit with missing data imputation (1 chain)
get_bgmcompare_fit_impute = function() {
  if(is.null(.test_cache$bgmcompare_fit_impute)) {
    data("Wenchuan", package = "bgms")
    x = Wenchuan[1:25, 1:4]
    y = Wenchuan[26:50, 1:4]
    x[5, 2] = NA
    y[10, 3] = NA
    .test_cache$bgmcompare_fit_impute = without_support_warning(bgmCompare(
        x = x, y = y,
        na_action = "impute",
        iter = 25, warmup = 50, chains = 1,
        seed = 11112,
        display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit_impute
}

#' @description Get cached bgmCompare fit
#' with Blume-Capel + missing data
#' imputation (1 chain)
get_bgmcompare_fit_blumecapel_impute = function() {
  if(is.null(.test_cache$bgmcompare_fit_bc_impute)) {
    data("Boredom", package = "bgms")
    # Select 25 rows from each language group
    rows = c(1:25, 491:515)
    x = Boredom[rows, 2:5] # 4 ordinal variables (7 categories)
    x[5, 2] = NA
    x[30, 3] = NA # Row in second group
    # The character column is a valid indicator; groups number by first appearance.
    lang = Boredom[rows, "language"]
    .test_cache$bgmcompare_fit_bc_impute = without_support_warning(bgmCompare(
        x = x,
        group_indicator = lang,
        variable_type = "blume-capel",
        baseline_category = 3,
        na_action = "impute",
        iter = 25, warmup = 50, chains = 1,
        seed = 11113,
        display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit_bc_impute
}

#' @description Get cached bgms fit for GGM
#' with edge selection
#' (4 continuous variables, 1 chain)
get_bgms_fit_ggm = function() {
  if(is.null(.test_cache$bgms_fit_ggm)) {
    set.seed(42)
    x = matrix(rnorm(200), nrow = 50, ncol = 4)
    colnames(x) = paste0("V", 1:4)
    .test_cache$bgms_fit_ggm = bgm(
      x = x,
      variable_type = "continuous",
      edge_selection = TRUE,
      iter = 50, warmup = 100, chains = 1,
      seed = 44442,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_ggm
}

#' @description Get cached bgms fit for GGM
#' without edge selection
#' (4 continuous variables, 1 chain)
get_bgms_fit_ggm_no_es = function() {
  if(is.null(.test_cache$bgms_fit_ggm_no_es)) {
    set.seed(42)
    x = matrix(rnorm(200), nrow = 50, ncol = 4)
    colnames(x) = paste0("V", 1:4)
    .test_cache$bgms_fit_ggm_no_es = bgm(
      x = x,
      variable_type = "continuous",
      edge_selection = FALSE,
      iter = 50, warmup = 100, chains = 1,
      seed = 44443,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_ggm_no_es
}

get_bgms_fit_mixed_mrf = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf)) {
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
    .test_cache$bgms_fit_mixed_mrf = bgm(
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
  }
  .test_cache$bgms_fit_mixed_mrf
}

get_bgms_fit_mixed_mrf_no_es = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf_no_es)) {
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
    .test_cache$bgms_fit_mixed_mrf_no_es = bgm(
      x = x,
      variable_type = c(
        "ordinal", "continuous", "ordinal",
        "continuous", "ordinal"
      ),
      edge_selection = FALSE,
      iter = 50, warmup = 100, chains = 1,
      seed = 77772,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_mixed_mrf_no_es
}

get_bgms_fit_mixed_mrf_marginal = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf_marginal)) {
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
    .test_cache$bgms_fit_mixed_mrf_marginal = bgm(
      x = x,
      variable_type = c(
        "ordinal", "continuous", "ordinal",
        "continuous", "ordinal"
      ),
      edge_selection = FALSE,
      iter = 50, warmup = 100, chains = 1,
      seed = 77773,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_mixed_mrf_marginal
}

get_bgms_fit_mixed_mrf_marginal_es = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf_marginal_es)) {
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
    .test_cache$bgms_fit_mixed_mrf_marginal_es = bgm(
      x = x,
      variable_type = c(
        "ordinal", "continuous", "ordinal",
        "continuous", "ordinal"
      ),
      edge_selection = TRUE,
      iter = 50, warmup = 100, chains = 1,
      seed = 77774,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_mixed_mrf_marginal_es
}

get_bgms_fit_mixed_mrf_nuts = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf_nuts)) {
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
    .test_cache$bgms_fit_mixed_mrf_nuts = bgm(
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
  }
  .test_cache$bgms_fit_mixed_mrf_nuts
}

get_bgms_fit_mixed_mrf_nuts_no_es = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf_nuts_no_es)) {
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
    .test_cache$bgms_fit_mixed_mrf_nuts_no_es = bgm(
      x = x,
      variable_type = c(
        "ordinal", "continuous", "ordinal",
        "continuous", "ordinal"
      ),
      edge_selection = FALSE,
      update_method = "nuts",
      iter = 50, warmup = 100, chains = 1,
      seed = 77776,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_mixed_mrf_nuts_no_es
}

get_bgms_fit_mixed_mrf_beta_bernoulli = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf_beta_bernoulli)) {
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
    .test_cache$bgms_fit_mixed_mrf_beta_bernoulli = bgm(
      x = x,
      variable_type = c(
        "ordinal", "continuous", "ordinal",
        "continuous", "ordinal"
      ),
      edge_selection = TRUE,
      edge_prior = beta_bernoulli_prior(alpha = 1, beta = 1),
      iter = 50, warmup = 100, chains = 1,
      seed = 77777,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_mixed_mrf_beta_bernoulli
}

get_bgms_fit_mixed_mrf_sbm = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf_sbm)) {
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
    # Two continuous variables: the block-model correction warns that its
    # slope curve is not resolvable at a single tilted pair and keeps the
    # plain conjugate updates (asserted in test-mixed-correction.R).
    .test_cache$bgms_fit_mixed_mrf_sbm = suppressWarnings(bgm(
      x = x,
      variable_type = c(
        "ordinal", "continuous", "ordinal",
        "continuous", "ordinal"
      ),
      edge_selection = TRUE,
      edge_prior = sbm_prior(),
      iter = 50, warmup = 100, chains = 1,
      seed = 77778,
      display_progress = "none"
    ))
  }
  .test_cache$bgms_fit_mixed_mrf_sbm
}

get_bgms_fit_mixed_mrf_bc = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf_bc)) {
    set.seed(99)
    n = 80
    x = cbind(
      sample(0:4, n, replace = TRUE),
      rnorm(n),
      sample(0:4, n, replace = TRUE),
      rnorm(n)
    )
    colnames(x) = c("bc1", "c1", "bc2", "c2")
    .test_cache$bgms_fit_mixed_mrf_bc = bgm(
      x = x,
      variable_type = c(
        "blume-capel", "continuous",
        "blume-capel", "continuous"
      ),
      baseline_category = 2L,
      edge_selection = TRUE,
      iter = 50, warmup = 100, chains = 1,
      seed = 77779,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_mixed_mrf_bc
}

get_bgms_fit_mixed_mrf_impute = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf_impute)) {
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
    # Insert NAs in both discrete and continuous columns
    x[1, 1] = NA
    x[2, 2] = NA
    .test_cache$bgms_fit_mixed_mrf_impute = bgm(
      x = x,
      variable_type = c(
        "ordinal", "continuous", "ordinal",
        "continuous", "ordinal"
      ),
      edge_selection = TRUE,
      na_action = "impute",
      iter = 50, warmup = 100, chains = 1,
      seed = 77780,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_mixed_mrf_impute
}

get_bgms_fit_mixed_mrf_multichain = function() {
  if(is.null(.test_cache$bgms_fit_mixed_mrf_multichain)) {
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
    .test_cache$bgms_fit_mixed_mrf_multichain = bgm(
      x = x,
      variable_type = c(
        "ordinal", "continuous", "ordinal",
        "continuous", "ordinal"
      ),
      edge_selection = TRUE,
      iter = 50, warmup = 100, chains = 2,
      seed = 77781,
      display_progress = "none"
    )
  }
  .test_cache$bgms_fit_mixed_mrf_multichain
}

# The user-facing check surface -- verdicts(), extract_centrality(),
# calibration_check(), plot() -- reads a settled posterior rather than driving
# the sampler, so one fit per cell serves every test that inspects it. The
# chains are long enough for the fragility flag to see both of its states, which
# is what these fits are sized for.

#' @description Get cached bgms fit
#' (6 ordinal Wenchuan variables, edge selection, 2 chains)
get_bgms_fit_wenchuan6 = function() {
  if(is.null(.test_cache$bgms_fit_wenchuan6)) {
    data("Wenchuan", package = "bgms")
    .test_cache$bgms_fit_wenchuan6 = bgm(
      Wenchuan[, 1:6],
      chains = 2, iter = 400, warmup = 400, cores = 2, seed = 1,
      display_progress = "none", verbose = FALSE
    )
  }
  .test_cache$bgms_fit_wenchuan6
}

#' @description Get cached bgms fit
#' (6 ordinal Wenchuan variables, NO edge selection, 2 chains).
#' The panel of such a fit reads its evidence off the prior and posterior
#' ordinates at zero, because there is no indicator to Rao-Blackwellize.
get_bgms_fit_wenchuan6_noselection = function() {
  if(is.null(.test_cache$bgms_fit_wenchuan6_nosel)) {
    data("Wenchuan", package = "bgms")
    .test_cache$bgms_fit_wenchuan6_nosel = bgm(
      Wenchuan[, 1:6], edge_selection = FALSE,
      chains = 2, iter = 400, warmup = 400, cores = 2, seed = 1,
      display_progress = "none", verbose = FALSE
    )
  }
  .test_cache$bgms_fit_wenchuan6_nosel
}

#' @description Get cached bgms fit
#' (5 ordinal Wenchuan variables, edge selection, 2 chains)
get_bgms_fit_wenchuan5 = function() {
  if(is.null(.test_cache$bgms_fit_wenchuan5)) {
    data("Wenchuan", package = "bgms")
    .test_cache$bgms_fit_wenchuan5 = bgm(
      Wenchuan[, 1:5],
      chains = 2, iter = 300, warmup = 300, cores = 2, seed = 7,
      display_progress = "none", verbose = FALSE
    )
  }
  .test_cache$bgms_fit_wenchuan5
}

#' @description Get cached bgmCompare fit
#' (5 ordinal Wenchuan variables, 2 groups, difference selection, 2 chains)
get_bgmcompare_fit_wenchuan5 = function() {
  if(is.null(.test_cache$bgmcompare_fit_wenchuan5)) {
    data("Wenchuan", package = "bgms")
    .test_cache$bgmcompare_fit_wenchuan5 = without_support_warning(bgmCompare(
        x = Wenchuan[1:120, 1:5], group_indicator = rep(1:2, each = 60),
        iter = 300, warmup = 300, chains = 2, cores = 2, seed = 13,
        difference_selection = TRUE, display_progress = "none"
    ))
  }
  .test_cache$bgmcompare_fit_wenchuan5
}

# ------------------------------------------------------------------------------
# 2. Prediction Data Helpers
# ------------------------------------------------------------------------------

#' Get prediction data matching the binary bgms fixture
get_prediction_data_binary = function(n = 10) {
  data("ADHD", package = "bgms")
  ADHD[51:(50 + n), 2:5] # Use different rows than training, 4 variables
}

#' Get prediction data matching the ordinal bgms fixture
get_prediction_data_ordinal = function(n = 10) {
  data("Wenchuan", package = "bgms")
  Wenchuan[51:(50 + n), 1:4] # Use different rows than training, 4 variables
}

#' Get prediction data matching the binary bgmCompare fixture
get_prediction_data_bgmcompare_binary = function(n = 10) {
  data("ADHD", package = "bgms")
  ADHD[sample(nrow(ADHD), n), 2:5] # Random sample, 4 variables
}

#' Get prediction data matching the ordinal bgmCompare fixture
get_prediction_data_bgmcompare_ordinal = function(n = 10) {
  data("Wenchuan", package = "bgms")
  Wenchuan[sample(nrow(Wenchuan), n), 1:4] # Random sample, 4 variables
}

#' Get prediction data matching the Blume-Capel bgmCompare fixture (Boredom)
get_prediction_data_bgmcompare_blumecapel = function(n = 10) {
  data("Boredom", package = "bgms")
  Boredom[26:35, 2:5] # Use different rows than training, 4 ordinal variables
}

#' Get prediction data matching the GGM bgms fixture (continuous)
get_prediction_data_ggm = function(n = 10) {
  set.seed(99)
  x = matrix(rnorm(n * 4), nrow = n, ncol = 4)
  colnames(x) = paste0("V", 1:4)
  x
}

#' Get prediction data matching the mixed MRF bgms fixture
#' Columns: d1 (ordinal 0-2), c1 (continuous), d2 (ordinal 0-2),
#'          c2 (continuous), d3 (ordinal 0-2)
get_prediction_data_mixed = function(n = 10) {
  set.seed(199)
  x = cbind(
    sample(0:2, n, replace = TRUE),
    rnorm(n),
    sample(0:2, n, replace = TRUE),
    rnorm(n),
    sample(0:2, n, replace = TRUE)
  )
  colnames(x) = c("d1", "c1", "d2", "c2", "d3")
  x
}

#' Get prediction data matching the mixed MRF Blume-Capel bgms fixture
#' Columns: bc1 (ordinal 0-4), c1 (continuous), bc2 (ordinal 0-4),
#'          c2 (continuous)
get_prediction_data_mixed_bc = function(n = 10) {
  set.seed(299)
  x = cbind(
    sample(0:4, n, replace = TRUE),
    rnorm(n),
    sample(0:4, n, replace = TRUE),
    rnorm(n)
  )
  colnames(x) = c("bc1", "c1", "bc2", "c2")
  x
}

# ------------------------------------------------------------------------------
# 3. Test Data Generators
# ------------------------------------------------------------------------------

#' Generate small test dataset for quick MCMC runs
#' @param n Number of observations
#' @param p Number of variables
#' @param seed Random seed
generate_test_data = function(n = 30, p = 4, seed = 42) {
  set.seed(seed)
  # Binary/ordinal data with values 0, 1, 2
  data = matrix(sample(0:2, n * p, replace = TRUE), nrow = n, ncol = p)
  colnames(data) = paste0("V", seq_len(p))
  as.data.frame(data)
}

#' Generate grouped test data for bgmCompare
#' @param n_per_group Observations per group
#' @param p Number of variables
#' @param n_groups Number of groups
#' @param seed Random seed
generate_grouped_test_data = function(n_per_group = 20, p = 4, n_groups = 2,
                                      seed = 42) {
  set.seed(seed)
  total_n = n_per_group * n_groups
  data = matrix(sample(0:2, total_n * p, replace = TRUE),
    nrow = total_n, ncol = p
  )
  colnames(data) = paste0("V", seq_len(p))
  list(
    x = as.data.frame(data),
    group_indicator = rep(seq_len(n_groups), each = n_per_group)
  )
}


# ------------------------------------------------------------------------------
# 3. Matrix Validation Helpers
# ------------------------------------------------------------------------------

#' Check if matrix is symmetric within tolerance
is_symmetric = function(M, tol = 1e-10) {
  if(!is.matrix(M)) {
    return(FALSE)
  }
  if(nrow(M) != ncol(M)) {
    return(FALSE)
  }
  max(abs(M - t(M)), na.rm = TRUE) <= tol
}

#' Check if all values in matrix are within bounds
values_in_range = function(M, lower = -Inf, upper = Inf) {
  vals = as.vector(M)
  vals = vals[!is.na(vals)]
  all(vals >= lower & vals <= upper)
}

#' Get upper triangle values (for pairwise parameters)
upper_vals = function(M) {
  M[upper.tri(M)]
}

#' Check that named summary entries match
#' matrix positions (ordering consistency)
#'
#' For each row of summary_df (named "Vi-Vj"), verify that summary_df$mean[k]
#' equals matrix_val[Vi, Vj]. Returns a logical vector (TRUE = match).
#' Requires p >= 4 to detect row-major vs column-major ordering bugs.
check_summary_matrix_consistency = function(summary_df, matrix_val) {
  matches = logical(nrow(summary_df))
  for(k in seq_len(nrow(summary_df))) {
    parts = strsplit(rownames(summary_df)[k], "-")[[1]]
    matches[k] = abs(
      summary_df$mean[k] -
        matrix_val[parts[1], parts[2]]
    ) < 1e-10
  }
  matches
}

#' Check that extractor column means match
#' matrix positions (ordering consistency)
#'
#' For each named element of extracted_means (named "Vi-Vj"), verify that
#' the value matches matrix_val[Vi, Vj].
#' Returns a logical vector (TRUE = match).
check_extractor_matrix_consistency = function(extracted_means, matrix_val) {
  matches = logical(length(extracted_means))
  for(k in seq_along(extracted_means)) {
    parts = strsplit(names(extracted_means)[k], "-")[[1]]
    matches[k] = abs(extracted_means[k] - matrix_val[parts[1], parts[2]]) < 1e-6
  }
  matches
}


# ------------------------------------------------------------------------------
# 4. Contract Testing Utilities
# ------------------------------------------------------------------------------
# These helpers verify that extractor functions return objects with expected
# structure, enabling contract testing for downstream packages like easybgm.

#' Verify extractor output structure
#' @param obj Output from an extractor function
#' @param type Expected type: "matrix", "data.frame", "list", "numeric", etc.
#' @param expected_dim Expected dimensions (for matrix/data.frame)
#' @param expected_names Expected column/row names or list names
expect_extractor_structure = function(obj, type, expected_dim = NULL,
                                      expected_names = NULL) {
  # Type check
  expect_true(
    inherits(obj, type),
    info = sprintf(
      "Expected class %s, got %s",
      type, paste(class(obj), collapse = ", ")
    )
  )

  # Dimension check
  if(!is.null(expected_dim)) {
    if(is.matrix(obj) || is.data.frame(obj)) {
      expect_equal(dim(obj), expected_dim,
        info = sprintf(
          "Expected dim %s, got %s",
          paste(expected_dim, collapse = "x"),
          paste(dim(obj), collapse = "x")
        )
      )
    }
  }

  # Names check
  if(!is.null(expected_names)) {
    if(is.matrix(obj)) {
      expect_true(
        all(expected_names %in% colnames(obj)) ||
          all(expected_names %in% rownames(obj)),
        info = "Expected names not found in matrix row/colnames"
      )
    } else if(is.list(obj)) {
      expect_true(
        all(expected_names %in% names(obj)),
        info = sprintf(
          "Expected list names %s, got %s",
          paste(expected_names, collapse = ", "),
          paste(names(obj), collapse = ", ")
        )
      )
    }
  }
}

#' Check that function errors with expected message pattern
expect_error_pattern = function(expr, pattern) {
  expect_error(expr, regexp = pattern)
}


# ------------------------------------------------------------------------------
# 5. MCMC Test Helpers
# ------------------------------------------------------------------------------

#' Get appropriate number of cores for testing
test_cores = function() {
  on_ci = isTRUE(as.logical(Sys.getenv("CI", "false")))
  if(on_ci) 2L else min(2L, parallel::detectCores())
}

#' Quick MCMC settings for testing (minimal iterations)
quick_mcmc_args = function() {
  list(
    iter = 100,
    warmup = 100,
    chains = 1,
    display_progress = "none"
  )
}

#' Moderate MCMC settings for more thorough testing
moderate_mcmc_args = function() {
  list(
    iter = 500,
    warmup = 500,
    chains = 2,
    display_progress = "none"
  )
}


# ==============================================================================
# 6. Consolidated Fixture Spec Lists (single source of truth)
# ==============================================================================
# These are used by test-methods.R, test-simulate-predict-regression.R,
# and test-extractor-functions.R. Each entry is a named list describing one
# fixture: label, get_fit, var_type, and boolean flags (is_continuous, is_mixed).
# Entries used by simulate/predict also carry get_prediction_data.

# ------------------------------------------------------------------
# get_bgms_fixtures
# ------------------------------------------------------------------
# All bgms fit variants for parameterized testing.
# Entries carry get_prediction_data for simulate/predict loops.
#
# Returns: list of fixture spec lists.
# ------------------------------------------------------------------
get_bgms_fixtures = function() {
  list(
    list(
      label = "binary",
      get_fit = get_bgms_fit,
      get_prediction_data = get_prediction_data_binary,
      var_type = "binary",
      is_continuous = FALSE
    ),
    list(
      label = "ordinal",
      get_fit = get_bgms_fit_ordinal,
      get_prediction_data = get_prediction_data_ordinal,
      var_type = "ordinal",
      is_continuous = FALSE
    ),
    list(
      label = "single-chain",
      get_fit = get_bgms_fit_single_chain,
      get_prediction_data = get_prediction_data_binary,
      var_type = "binary",
      is_continuous = FALSE
    ),
    list(
      label = "blume-capel",
      get_fit = get_bgms_fit_blumecapel,
      get_prediction_data = get_prediction_data_ordinal,
      var_type = "blume-capel",
      is_continuous = FALSE
    ),
    list(
      label = "adaptive-metropolis",
      get_fit = get_bgms_fit_adaptive_metropolis,
      get_prediction_data = get_prediction_data_binary,
      var_type = "binary",
      is_continuous = FALSE
    ),
    list(
      label = "am-blumecapel",
      get_fit = get_bgms_fit_am_blumecapel,
      get_prediction_data = get_prediction_data_ordinal,
      var_type = "blume-capel",
      is_continuous = FALSE
    ),
    list(
      label = "impute",
      get_fit = get_bgms_fit_impute,
      get_prediction_data = get_prediction_data_ordinal,
      var_type = "ordinal",
      is_continuous = FALSE
    ),
    list(
      label = "beta-bernoulli",
      get_fit = get_bgms_fit_beta_bernoulli,
      get_prediction_data = get_prediction_data_binary,
      var_type = "binary",
      is_continuous = FALSE
    ),
    list(
      label = "sbm",
      get_fit = get_bgms_fit_sbm,
      get_prediction_data = get_prediction_data_binary,
      var_type = "binary",
      is_continuous = FALSE
    ),
    list(
      label = "ggm",
      get_fit = get_bgms_fit_ggm,
      get_prediction_data = get_prediction_data_ggm,
      var_type = "continuous",
      is_continuous = TRUE
    ),
    list(
      label = "ggm-no-es",
      get_fit = get_bgms_fit_ggm_no_es,
      get_prediction_data = get_prediction_data_ggm,
      var_type = "continuous",
      is_continuous = TRUE
    ),
    list(
      label = "mixed-mrf",
      get_fit = get_bgms_fit_mixed_mrf,
      get_prediction_data = get_prediction_data_mixed,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    ),
    list(
      label = "mixed-mrf-no-es",
      get_fit = get_bgms_fit_mixed_mrf_no_es,
      get_prediction_data = get_prediction_data_mixed,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    ),
    list(
      label = "mixed-mrf-marginal",
      get_fit = get_bgms_fit_mixed_mrf_marginal,
      get_prediction_data = get_prediction_data_mixed,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    ),
    list(
      label = "mixed-mrf-marginal-es",
      get_fit = get_bgms_fit_mixed_mrf_marginal_es,
      get_prediction_data = get_prediction_data_mixed,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    ),
    list(
      label = "mixed-mrf-nuts",
      get_fit = get_bgms_fit_mixed_mrf_nuts,
      get_prediction_data = get_prediction_data_mixed,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    ),
    list(
      label = "mixed-mrf-nuts-no-es",
      get_fit = get_bgms_fit_mixed_mrf_nuts_no_es,
      get_prediction_data = get_prediction_data_mixed,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    ),
    list(
      label = "mixed-mrf-beta-bernoulli",
      get_fit = get_bgms_fit_mixed_mrf_beta_bernoulli,
      get_prediction_data = get_prediction_data_mixed,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    ),
    list(
      label = "mixed-mrf-sbm",
      get_fit = get_bgms_fit_mixed_mrf_sbm,
      get_prediction_data = get_prediction_data_mixed,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    ),
    list(
      label = "mixed-mrf-bc",
      get_fit = get_bgms_fit_mixed_mrf_bc,
      get_prediction_data = get_prediction_data_mixed_bc,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    ),
    list(
      label = "mixed-mrf-impute",
      get_fit = get_bgms_fit_mixed_mrf_impute,
      get_prediction_data = get_prediction_data_mixed,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    ),
    list(
      label = "mixed-mrf-multichain",
      get_fit = get_bgms_fit_mixed_mrf_multichain,
      get_prediction_data = get_prediction_data_mixed,
      var_type = "mixed",
      is_continuous = FALSE,
      is_mixed = TRUE
    )
  )
}

# ------------------------------------------------------------------
# get_bgmcompare_fixtures
# ------------------------------------------------------------------
# All bgmCompare fit variants for parameterized testing.
#
# Returns: list of fixture spec lists.
# ------------------------------------------------------------------
get_bgmcompare_fixtures = function() {
  list(
    list(
      label = "binary",
      get_fit = get_bgmcompare_fit,
      get_prediction_data = get_prediction_data_bgmcompare_binary,
      var_type = "binary"
    ),
    list(
      label = "ordinal",
      get_fit = get_bgmcompare_fit_ordinal,
      get_prediction_data = get_prediction_data_bgmcompare_ordinal,
      var_type = "ordinal"
    ),
    list(
      label = "adaptive-metropolis",
      get_fit = get_bgmcompare_fit_adaptive_metropolis,
      get_prediction_data = get_prediction_data_bgmcompare_binary,
      var_type = "binary"
    ),
    list(
      label = "blume-capel",
      get_fit = get_bgmcompare_fit_blumecapel,
      get_prediction_data = get_prediction_data_bgmcompare_blumecapel,
      var_type = "blume-capel"
    ),
    list(
      label = "am-blume-capel",
      get_fit = get_bgmcompare_fit_am_blumecapel,
      get_prediction_data = get_prediction_data_bgmcompare_blumecapel,
      var_type = "blume-capel"
    ),
    list(
      label = "impute",
      get_fit = get_bgmcompare_fit_impute,
      get_prediction_data = get_prediction_data_bgmcompare_ordinal,
      var_type = "ordinal"
    ),
    list(
      label = "blume-capel-impute",
      get_fit = get_bgmcompare_fit_blumecapel_impute,
      get_prediction_data = get_prediction_data_bgmcompare_blumecapel,
      var_type = "blume-capel"
    ),
    list(
      label = "beta-bernoulli",
      get_fit = get_bgmcompare_fit_beta_bernoulli,
      get_prediction_data = get_prediction_data_bgmcompare_ordinal,
      var_type = "ordinal"
    ),
    list(
      label = "xy",
      get_fit = get_bgmcompare_fit_xy,
      get_prediction_data = get_prediction_data_bgmcompare_ordinal,
      var_type = "ordinal"
    ),
    list(
      label = "main-selection",
      get_fit = get_bgmcompare_fit_main_selection,
      get_prediction_data = get_prediction_data_bgmcompare_blumecapel,
      var_type = "blume-capel"
    )
  )
}

# ------------------------------------------------------------------
# get_extractor_fixtures
# ------------------------------------------------------------------
# Representative subset for extractor-function contract tests.
# Covers all model families (OMRF, GGM, mixed) and both object classes.
#
# Returns: list of fixture spec lists.
# ------------------------------------------------------------------
get_extractor_fixtures = function() {
  list(
    list(
      label = "bgms_binary",
      get_fit = get_bgms_fit,
      type = "bgms",
      var_type = "binary"
    ),
    list(
      label = "bgms_ordinal",
      get_fit = get_bgms_fit_ordinal,
      type = "bgms",
      var_type = "ordinal"
    ),
    list(
      label = "bgms_blumecapel",
      get_fit = get_bgms_fit_blumecapel,
      type = "bgms",
      var_type = "blume-capel"
    ),
    list(
      label = "bgms_ggm",
      get_fit = get_bgms_fit_ggm,
      type = "bgms",
      var_type = "continuous",
      is_continuous = TRUE
    ),
    list(
      label = "bgms_mixed",
      get_fit = get_bgms_fit_mixed_mrf,
      type = "bgms",
      var_type = "mixed",
      is_mixed = TRUE
    ),
    list(
      label = "bgmCompare_binary",
      get_fit = get_bgmcompare_fit,
      type = "bgmCompare",
      var_type = "binary"
    ),
    list(
      label = "bgmCompare_ordinal",
      get_fit = get_bgmcompare_fit_ordinal,
      type = "bgmCompare",
      var_type = "ordinal"
    ),
    list(
      label = "bgmCompare_blumecapel",
      get_fit = get_bgmcompare_fit_blumecapel,
      type = "bgmCompare",
      var_type = "blume-capel"
    ),
    list(
      label = "bgmCompare_main_sel",
      get_fit = get_bgmcompare_fit_main_selection,
      type = "bgmCompare",
      var_type = "blume-capel"
    )
  )
}

# ------------------------------------------------------------------------------
# Group-support warning
# ------------------------------------------------------------------------------
#
# bgmCompare() warns, by design, whenever a retained category has no
# observations in some group (F-075/F-110). Small fixtures -- 25 rows per group
# of a five-point scale -- trigger it routinely, and for a test about
# reproducibility or object structure it is incidental noise. This muffles that
# ONE condition class and nothing else, so a genuine new warning still fails the
# run. Tests that are about the support behaviour assert on it directly instead;
# see test-collapse-categories.R.
without_support_warning = function(expr) {
  withCallingHandlers(
    expr,
    bgms_group_support_warning = function(w) invokeRestart("muffleWarning")
  )
}
