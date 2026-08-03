# ==============================================================================
# Contract Tests for Extractor Functions (Parameterized)
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Range invariants, symmetry checks, dimension consistency
#
# This file uses parameterized testing (specs + loop) to reduce code repetition.
# Each extractor is tested across multiple fixture types with shared assertions.
#
# IMPORTANT: Changes to extractor function output structure may break easybgm!
# ==============================================================================

# ------------------------------------------------------------------------------
# Fixture Specifications <U+2014> defined in helper-fixtures.R
# get_extractor_fixtures()
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# extract_arguments() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_arguments returns complete argument list for all fit types", {
  fixtures = get_extractor_fixtures()

  for(spec in fixtures) {
    ctx = sprintf("[%s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    # Basic structure
    expect_true(is.list(args), info = paste(ctx, "should be list"))
    expect_true(length(args) > 0, info = ctx)

    # Essential fields present
    expect_true("num_variables" %in% names(args), info = paste(ctx, "missing num_variables"))
    expect_true("num_cases" %in% names(args), info = paste(ctx, "missing num_cases"))
    expect_true("data_columnnames" %in% names(args), info = paste(ctx, "missing data_columnnames"))

    # Values are sensible
    expect_true(args$num_variables >= 1, info = ctx)
    expect_true(args$num_cases >= 1, info = ctx)

    # Type-specific fields
    if(spec$type == "bgms") {
      expect_true(is.logical(args$edge_selection), info = paste(ctx, "edge_selection should be logical"))
    } else {
      expect_true(args$num_groups >= 2, info = paste(ctx, "bgmCompare should have >= 2 groups"))
    }
  }
})

test_that("extract_arguments errors on non-bgms objects", {
  expect_error(extract_arguments(list()), class = "error")
  expect_error(extract_arguments(data.frame()), class = "error")
})

test_that("extract_arguments exposes main_effect_indices for bgmCompare fits", {
  # Regression: the main-effect row layout lived only in the internal cache,
  # so callers that need to map main-effect parameter rows back to variables
  # got NULL from extract_arguments() and errored downstream.
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)

  mei = args$main_effect_indices
  expect_false(is.null(mei), info = "bgmCompare should expose main_effect_indices")
  expect_true(is.matrix(mei))
  expect_equal(nrow(mei), args$num_variables)
  expect_equal(ncol(mei), 2L)

  # Zero-based, contiguous, non-overlapping blocks covering every parameter row.
  expect_equal(mei[1, 1], 0L)
  expect_true(all(mei[, 2] >= mei[, 1]))
  if(nrow(mei) > 1) {
    expect_equal(mei[-1, 1], mei[-nrow(mei), 2] + 1L)
  }

  # Each block is as wide as that variable's parameter count: one row per
  # category for ordinal variables, two for Blume-Capel.
  widths = mei[, 2] - mei[, 1] + 1L
  expected = ifelse(args$is_ordinal_variable, args$num_categories, 2L)
  expect_equal(as.integer(widths), as.integer(expected))

  # Together the blocks tile the baseline block of the main-effect samples,
  # which is followed by the difference columns.
  num_baseline = mei[nrow(mei), 2] + 1L
  expect_equal(num_baseline, sum(as.integer(expected)))
  expect_lte(num_baseline, ncol(fit$raw_samples$main[[1]]))
})


# ------------------------------------------------------------------------------
# extract_pairwise_interactions() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_pairwise_interactions returns valid matrix for all fit types", {
  fixtures = get_extractor_fixtures()

  for(spec in fixtures) {
    ctx = sprintf("[%s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    pairwise = extract_pairwise_interactions(fit)

    # Structure checks
    expect_true(is.matrix(pairwise), info = paste(ctx, "should be matrix"))

    p = args$num_variables
    expected_cols = p * (p - 1) / 2
    expect_equal(ncol(pairwise), expected_cols,
      info = paste(ctx, "wrong number of edge columns")
    )

    # Values finite
    expect_true(all(is.finite(pairwise)), info = paste(ctx, "should have finite values"))

    # Has column names
    expect_true(!is.null(colnames(pairwise)), info = paste(ctx, "should have column names"))
  }
})


# ------------------------------------------------------------------------------
# extract_main_effects() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_main_effects returns valid output for all fit types", {
  fixtures = get_extractor_fixtures()

  for(spec in fixtures) {
    ctx = sprintf("[%s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    main = extract_main_effects(fit)

    if(isTRUE(args$is_continuous)) {
      # GGM: no main effects; returns NULL silently
      main_null = extract_main_effects(fit)
      expect_null(main_null, info = paste(ctx, "GGM should return NULL"))
    } else if(isTRUE(args$is_mixed)) {
      # Mixed MRF returns a list
      expect_true(is.list(main), info = paste(ctx, "should be list for mixed"))
      expect_true(is.matrix(main$discrete), info = paste(ctx, "$discrete should be matrix"))
      expect_true(is.matrix(main$continuous), info = paste(ctx, "$continuous should be matrix"))
    } else {
      # OMRF / Blume-Capel return matrix
      expect_true(is.matrix(main), info = paste(ctx, "should be matrix"))
      vals = main[!is.na(main)]
      expect_true(all(is.finite(vals)), info = paste(ctx, "non-NA values should be finite"))
    }
  }
})

test_that("extract_category_thresholds emits deprecation warning", {
  fit = get_bgms_fit()
  expect_warning(
    extract_category_thresholds(fit),
    "extract_main_effects"
  )
})


# ------------------------------------------------------------------------------
# extract_indicators() and extract_posterior_inclusion_probabilities() Tests
# ------------------------------------------------------------------------------
# These only apply to fits with edge_selection = TRUE

test_that("extract_indicators returns binary matrix for edge-selection fits", {
  # Only test fixtures with edge selection
  fixtures = list(
    list(label = "bgms_binary", get_fit = get_bgms_fit)
  )

  for(spec in fixtures) {
    ctx = sprintf("[%s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    if(!isTRUE(args$edge_selection)) {
      next
    }

    indicators = extract_indicators(fit)

    # Structure
    expect_true(is.matrix(indicators), info = ctx)

    p = args$num_variables
    expected_cols = p * (p - 1) / 2
    expect_equal(ncol(indicators), expected_cols, info = paste(ctx, "wrong indicator columns"))

    # Binary values
    expect_true(all(indicators %in% c(0, 1)),
      info = paste(ctx, "indicators should be 0 or 1")
    )
  }
})

test_that("extract_posterior_inclusion_probabilities returns symmetric PIP matrix", {
  # Only test fixtures with edge selection
  fixtures = list(
    list(label = "bgms_binary", get_fit = get_bgms_fit)
  )

  for(spec in fixtures) {
    ctx = sprintf("[%s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    if(!isTRUE(args$edge_selection)) {
      next
    }

    pip = extract_posterior_inclusion_probabilities(fit)
    p = args$num_variables

    # Structure
    expect_true(is.matrix(pip), info = ctx)
    expect_equal(dim(pip), c(p, p), info = paste(ctx, "should be p x p"))

    # Symmetry
    expect_true(is_symmetric(pip), info = paste(ctx, "should be symmetric"))

    # Range [0, 1]
    expect_true(values_in_range(pip, 0, 1), info = paste(ctx, "PIPs should be in [0,1]"))

    # Diagonal is zero (no self-loops)
    expect_true(all(diag(pip) == 0), info = paste(ctx, "diagonal should be 0"))

    # Has variable names
    expect_equal(colnames(pip), args$data_columnnames, info = ctx)
  }
})

test_that("extract_indicators errors when edge_selection = FALSE", {
  data = generate_test_data(n = 20, p = 3)
  args = c(list(x = data, edge_selection = FALSE), quick_mcmc_args())
  fit = do.call(bgm, args)

  expect_error(extract_indicators(fit), regexp = "edge_selection")
})


# ------------------------------------------------------------------------------
# extract_rhat() and extract_ess() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_rhat returns valid diagnostics for all fit types", {
  fixtures = get_extractor_fixtures()

  for(spec in fixtures) {
    ctx = sprintf("[%s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    rhat = extract_rhat(fit)

    expect_true(is.list(rhat), info = paste(ctx, "should be list"))

    if(spec$type == "bgms") {
      expect_true("pairwise" %in% names(rhat), info = paste(ctx, "missing pairwise"))
      expect_true(is.numeric(rhat$pairwise), info = ctx)
      expect_true(all(is.na(rhat$pairwise) | rhat$pairwise > 0),
        info = paste(ctx, "R-hat should be positive")
      )
    } else {
      expect_true("pairwise_baseline" %in% names(rhat), info = paste(ctx, "missing pairwise_baseline"))
      expect_true(is.numeric(rhat$pairwise_baseline), info = ctx)
      expect_true(all(is.na(rhat$pairwise_baseline) | rhat$pairwise_baseline > 0), info = ctx)
    }
  }
})

test_that("extract_ess returns valid diagnostics for all fit types", {
  fixtures = get_extractor_fixtures()

  for(spec in fixtures) {
    ctx = sprintf("[%s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    ess = extract_ess(fit)

    expect_true(is.list(ess), info = paste(ctx, "should be list"))

    if(spec$type == "bgms") {
      expect_true("pairwise" %in% names(ess), info = paste(ctx, "missing pairwise"))
      expect_true(is.numeric(ess$pairwise), info = ctx)
      expect_true(all(is.na(ess$pairwise) | ess$pairwise > 0),
        info = paste(ctx, "ESS should be positive")
      )
    } else {
      expect_true("pairwise_baseline" %in% names(ess), info = paste(ctx, "missing pairwise_baseline"))
      expect_true(is.numeric(ess$pairwise_baseline), info = ctx)
    }
  }
})

test_that("extract_ess indicators default to the RB n_eff and honour estimator", {
  fixtures = get_extractor_fixtures()
  checked = 0L

  for(spec in fixtures) {
    ctx = sprintf("[%s]", spec$label)
    fit = spec$get_fit()
    ind = summary(fit)$indicator

    if(is.null(ind) || !("n_eff" %in% names(ind))) next
    checked = checked + 1L

    rb = extract_ess(fit)$indicator

    # The default is the RB n_eff column, matching summary() and the R-hat that
    # extract_rhat() reports.
    expect_equal(unname(rb), ind$n_eff, info = paste(ctx, "default is RB n_eff"))
    expect_equal(unname(extract_ess(fit, estimator = "rb")$indicator), ind$n_eff,
      info = paste(ctx, "explicit rb is RB n_eff")
    )

    # NA masking is inherited from the summary table.
    expect_identical(is.na(unname(rb)), is.na(ind$n_eff),
      info = paste(ctx, "RB NA positions")
    )

    # The transition ESS is deprecated and no longer a summary column; the
    # estimator still resolves, warns, and recomputes it from the raw draws.
    expect_false("n_eff_mixt" %in% names(ind), info = paste(ctx, "no transition column"))
    expect_warning(
      mixt <- extract_ess(fit, estimator = "mixt")$indicator,
      class = "lifecycle_warning_deprecated"
    )
    raw = bgms:::get_fit_cache(fit)$raw
    ref = bgms:::.compute_indicator_ess_cpp(
      bgms:::combine_chains(raw, "indicator_samples")
    )[, "n_eff_mixt"]
    expect_equal(unname(mixt), unname(ref), info = paste(ctx, "mixt is the transition ESS"))

    expect_error(extract_ess(fit, estimator = "nonsense"))
  }

  expect_gt(checked, 0L)
})

test_that("extract_ess indicators fall back to the transition ESS without RB draws", {
  # Summary tables from fits predating the RB regime (bgms < 0.2.0.0) carry no
  # n_eff column: a default call falls back, an explicit "rb" errors, and an
  # explicit "mixt" reads the legacy column.
  rb_table = data.frame(n_eff = c(100, NA))
  legacy_table = data.frame(n_eff_mixt = c(80, NA))
  default_arg = c("rb", "mixt")

  expect_equal(indicator_ess_column(NULL, rb_table, default_arg, TRUE), c(100, NA))
  expect_equal(indicator_ess_column(NULL, legacy_table, default_arg, TRUE), c(80, NA))
  expect_equal(indicator_ess_column(NULL, legacy_table, "mixt", FALSE), c(80, NA))
  expect_error(indicator_ess_column(NULL, legacy_table, "rb", FALSE), "Rao-Blackwellized")
})

test_that("extract_ess indicator names and other elements ignore estimator", {
  fit = get_bgms_fit()
  rb = extract_ess(fit)
  expect_warning(mixt <- extract_ess(fit, estimator = "mixt"),
    class = "lifecycle_warning_deprecated"
  )

  expect_identical(names(rb), names(mixt))
  expect_identical(names(rb$indicator), names(mixt$indicator))
  expect_equal(rb$pairwise, mixt$pairwise)
})

test_that("extract_rhat and extract_ess error on non-bgms objects", {
  expect_error(extract_rhat(list()), class = "error")
  expect_error(extract_rhat(data.frame()), class = "error")
  expect_error(extract_ess(list()), class = "error")
  expect_error(extract_ess(data.frame()), class = "error")
})


# ------------------------------------------------------------------------------
# extract_indicator_priors() Tests
# ------------------------------------------------------------------------------

test_that("extract_indicator_priors returns prior specification", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)

  if(!isTRUE(args$edge_selection)) {
    skip("Fit object does not have edge_selection = TRUE")
  }

  priors = extract_indicator_priors(fit)

  expect_type(priors, "list")
  expect_true("type" %in% names(priors))

  valid_types = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block")
  expect_true(priors$type %in% valid_types)

  # Type-specific checks
  if(priors$type == "Bernoulli") {
    expect_true("prior_inclusion_probability" %in% names(priors))
    pip = priors$prior_inclusion_probability
    expect_true(all(pip >= 0 & pip <= 1))
  }

  if(priors$type == "Beta-Bernoulli") {
    expect_true(all(c("alpha", "beta") %in% names(priors)))
    expect_true(priors$alpha > 0 && priors$beta > 0)
  }
})

test_that("extract_indicator_priors errors when no selection performed", {
  data = generate_test_data(n = 20, p = 3)
  args = c(list(x = data, edge_selection = FALSE), quick_mcmc_args())
  fit = do.call(bgm, args)

  expect_error(extract_indicator_priors(fit), regexp = "selection")
})


# ------------------------------------------------------------------------------
# bgmCompare-specific Tests
# ------------------------------------------------------------------------------

test_that("extract_group_params returns group-level parameters", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)

  group_params = extract_group_params(fit)

  expect_type(group_params, "list")
  expect_true("main_effects_groups" %in% names(group_params))
  expect_true("pairwise_effects_groups" %in% names(group_params))

  # Dimensions match number of groups
  n_groups = args$num_groups
  expect_equal(ncol(group_params$main_effects_groups), n_groups)
  expect_equal(ncol(group_params$pairwise_effects_groups), n_groups)

  # Values finite
  expect_true(all(is.finite(group_params$main_effects_groups)))
  expect_true(all(is.finite(group_params$pairwise_effects_groups)))
})


# ------------------------------------------------------------------------------
# Cross-Function Consistency Tests
# ------------------------------------------------------------------------------

test_that("extractor outputs are dimensionally consistent", {
  fixtures = list(
    list(label = "bgms_binary", get_fit = get_bgms_fit)
  )

  for(spec in fixtures) {
    ctx = sprintf("[%s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    if(!isTRUE(args$edge_selection)) {
      next
    }

    p = args$num_variables
    n_edges = p * (p - 1) / 2

    # All should agree on number of variables/edges
    pip = extract_posterior_inclusion_probabilities(fit)
    expect_equal(nrow(pip), p, info = paste(ctx, "PIP rows"))

    indicators = extract_indicators(fit)
    expect_equal(ncol(indicators), n_edges, info = paste(ctx, "indicator cols"))

    pairwise = extract_pairwise_interactions(fit)
    expect_equal(ncol(pairwise), n_edges, info = paste(ctx, "pairwise cols"))

    thresholds = suppressWarnings(extract_category_thresholds(fit))
    expect_equal(nrow(thresholds), p, info = paste(ctx, "threshold rows"))
  }
})


# ------------------------------------------------------------------------------
# Contract Tests for easybgm Integration
# ------------------------------------------------------------------------------

test_that("bgms fit contains all fields accessed by easybgm", {
  fixtures = list(
    list(label = "bgms", get_fit = get_bgms_fit, type = "bgms"),
    list(label = "bgmCompare", get_fit = get_bgmcompare_fit, type = "bgmCompare")
  )

  for(spec in fixtures) {
    ctx = sprintf("[%s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    if(spec$type == "bgms") {
      expect_true("posterior_summary_pairwise" %in% names(fit), info = ctx)
      expect_true(is.data.frame(fit$posterior_summary_pairwise), info = ctx)
      expect_true("Rhat" %in% names(fit$posterior_summary_pairwise), info = ctx)
      expect_true("n_eff" %in% names(fit$posterior_summary_pairwise), info = ctx)

      if(isTRUE(args$edge_selection)) {
        expect_true("posterior_summary_indicator" %in% names(fit), info = ctx)
        # Inclusion summary now reports the Rao-Blackwellized estimate with
        # continuous ESS/Rhat (n_eff, not the binary n_eff_mixt).
        expect_true("n_eff" %in% names(fit$posterior_summary_indicator), info = ctx)
        expect_true("Rhat" %in% names(fit$posterior_summary_indicator), info = ctx)
      }
    } else {
      expect_true("posterior_summary_pairwise_baseline" %in% names(fit), info = ctx)
      expect_true(is.data.frame(fit$posterior_summary_pairwise_baseline), info = ctx)
      expect_true("Rhat" %in% names(fit$posterior_summary_pairwise_baseline), info = ctx)
      expect_true("n_eff" %in% names(fit$posterior_summary_pairwise_baseline), info = ctx)
    }
  }
})


# ------------------------------------------------------------------------------
# extract_indicators.bgmCompare Tests
# ------------------------------------------------------------------------------

test_that("extract_indicators.bgmCompare returns indicator matrix for difference selection fits", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)

  if(!isTRUE(args$difference_selection)) {
    skip("Fit object does not have difference_selection = TRUE")
  }

  indicators = extract_indicators(fit)

  # Structure
  expect_true(is.matrix(indicators), info = "should be matrix")

  # Binary values
  expect_true(all(indicators %in% c(0, 1)),
    info = "indicators should be 0 or 1"
  )

  # Has column names
  expect_true(!is.null(colnames(indicators)), info = "should have column names")
})

test_that("extract_indicators.bgmCompare errors when difference_selection = FALSE", {
  data = generate_grouped_test_data(n_per_group = 15, p = 3)
  args = c(
    list(x = data$x, group_indicator = data$group_indicator, difference_selection = FALSE),
    quick_mcmc_args()
  )
  fit = do.call(bgmCompare, args)

  expect_error(extract_indicators(fit), regexp = "difference_selection")
})


# ------------------------------------------------------------------------------
# extract_posterior_inclusion_probabilities.bgmCompare Tests
# ------------------------------------------------------------------------------

test_that("extract_posterior_inclusion_probabilities.bgmCompare returns symmetric PIP matrix", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)

  if(!isTRUE(args$difference_selection)) {
    skip("Fit object does not have difference_selection = TRUE")
  }

  pip = extract_posterior_inclusion_probabilities(fit)
  p = args$num_variables

  # Structure
  expect_true(is.matrix(pip), info = "should be matrix")
  expect_equal(dim(pip), c(p, p), info = "should be p x p")

  # Symmetry
  expect_true(is_symmetric(pip), info = "should be symmetric")

  # Range [0, 1]
  expect_true(values_in_range(pip, 0, 1), info = "PIPs should be in [0,1]")

  # Has variable names
  expect_equal(colnames(pip), args$data_columnnames, info = "should have column names")
})

test_that("extract_posterior_inclusion_probabilities.bgmCompare errors when difference_selection = FALSE", {
  data = generate_grouped_test_data(n_per_group = 15, p = 3)
  args = c(
    list(x = data$x, group_indicator = data$group_indicator, difference_selection = FALSE),
    quick_mcmc_args()
  )
  fit = do.call(bgmCompare, args)

  expect_error(extract_posterior_inclusion_probabilities(fit), regexp = "difference_selection")
})


# ------------------------------------------------------------------------------
# extract_indicator_priors.bgmCompare Tests
# ------------------------------------------------------------------------------

test_that("extract_indicator_priors.bgmCompare returns prior specification", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)

  if(!isTRUE(args$difference_selection)) {
    skip("Fit object does not have difference_selection = TRUE")
  }

  priors = extract_indicator_priors(fit)

  # Returns the difference_prior from arguments
  expect_true(!is.null(priors), info = "should return prior specification")
})

test_that("extract_indicator_priors.bgmCompare errors when difference_selection = FALSE", {
  data = generate_grouped_test_data(n_per_group = 15, p = 3)
  args = c(
    list(x = data$x, group_indicator = data$group_indicator, difference_selection = FALSE),
    quick_mcmc_args()
  )
  fit = do.call(bgmCompare, args)

  expect_error(extract_indicator_priors(fit), regexp = "selection")
})


# ------------------------------------------------------------------------------
# main_difference_selection Tests
# ------------------------------------------------------------------------------

test_that("bgmCompare with main_difference_selection = TRUE produces valid output", {
  fit = get_bgmcompare_fit_main_selection()
  args = extract_arguments(fit)

  # Verify main_difference_selection is TRUE in arguments

  expect_true(isTRUE(args$main_difference_selection),
    info = "main_difference_selection should be TRUE"
  )
  expect_true(isTRUE(args$difference_selection),
    info = "difference_selection should be TRUE"
  )
})

test_that("extract_indicators works with main_difference_selection = TRUE", {
  fit = get_bgmcompare_fit_main_selection()
  args = extract_arguments(fit)

  indicators = extract_indicators(fit)

  # Structure
  expect_true(is.matrix(indicators), info = "should be matrix")

  # Binary values
  expect_true(all(indicators %in% c(0, 1)),
    info = "indicators should be 0 or 1"
  )

  # With main_difference_selection = TRUE, there should be more indicator columns
  # than just pairwise (includes main effect indicators)
  p = args$num_variables
  n_pairwise = p * (p - 1) / 2

  # Indicator dimensions should include main effects + pairwise
  # Exact count depends on number of categories, but should be > n_pairwise
  expect_true(ncol(indicators) >= n_pairwise,
    info = "indicators should include at least pairwise effects"
  )
})

test_that("extract_posterior_inclusion_probabilities works with main_difference_selection = TRUE", {
  fit = get_bgmcompare_fit_main_selection()
  args = extract_arguments(fit)

  pip = extract_posterior_inclusion_probabilities(fit)
  p = args$num_variables

  # Structure - should be p x p matrix
  expect_true(is.matrix(pip), info = "should be matrix")
  expect_equal(dim(pip), c(p, p), info = "should be p x p")

  # Symmetry
  expect_true(is_symmetric(pip), info = "should be symmetric")

  # Range [0, 1]
  expect_true(values_in_range(pip, 0, 1), info = "PIPs should be in [0,1]")
})

test_that("extract_group_params works with main_difference_selection = TRUE", {
  fit = get_bgmcompare_fit_main_selection()
  args = extract_arguments(fit)

  group_params = extract_group_params(fit)

  expect_type(group_params, "list")
  expect_true("main_effects_groups" %in% names(group_params))
  expect_true("pairwise_effects_groups" %in% names(group_params))

  # Dimensions match number of groups
  n_groups = args$num_groups
  expect_equal(ncol(group_params$main_effects_groups), n_groups)
  expect_equal(ncol(group_params$pairwise_effects_groups), n_groups)

  # Values finite
  expect_true(all(is.finite(group_params$main_effects_groups)))
  expect_true(all(is.finite(group_params$pairwise_effects_groups)))
})


# ------------------------------------------------------------------------------
# extract_sbm.bgms Tests (Stochastic Block Model)
# ------------------------------------------------------------------------------

test_that("extract_sbm.bgms returns SBM summaries for Stochastic-Block prior", {
  fit = get_bgms_fit_sbm()
  args = extract_arguments(fit)

  sbm = extract_sbm(fit)

  # Structure
  expect_type(sbm, "list")

  # Required fields
  expect_true("posterior_num_blocks" %in% names(sbm),
    info = "should have posterior_num_blocks"
  )
  expect_true("posterior_mean_allocations" %in% names(sbm),
    info = "should have posterior_mean_allocations"
  )
  expect_true("posterior_mode_allocations" %in% names(sbm),
    info = "should have posterior_mode_allocations"
  )
  expect_true("posterior_mean_coclustering_matrix" %in% names(sbm),
    info = "should have posterior_mean_coclustering_matrix"
  )

  # Coclustering matrix should be symmetric
  ccm = sbm$posterior_mean_coclustering_matrix
  expect_true(is.matrix(ccm), info = "coclustering matrix should be matrix")
  expect_true(is_symmetric(ccm), info = "coclustering matrix should be symmetric")

  # Values in [0, 1] for coclustering probabilities
  expect_true(values_in_range(ccm, 0, 1),
    info = "coclustering probabilities should be in [0,1]"
  )
})

test_that("extract_sbm.bgms errors for non-SBM prior", {
  fit = get_bgms_fit() # Uses default Bernoulli prior

  expect_error(extract_sbm(fit), regexp = "Stochastic-Block")
})


# ------------------------------------------------------------------------------
# extract_indicator_priors with Beta-Bernoulli Prior Tests
# ------------------------------------------------------------------------------

test_that("extract_indicator_priors returns Beta-Bernoulli parameters", {
  fit = get_bgms_fit_beta_bernoulli()
  args = extract_arguments(fit)

  priors = extract_indicator_priors(fit)

  # Type check
  expect_type(priors, "list")
  expect_equal(priors$type, "Beta-Bernoulli")

  # Required parameters
  expect_true("alpha" %in% names(priors), info = "should have alpha parameter")
  expect_true("beta" %in% names(priors), info = "should have beta parameter")

  # Positive values

  expect_true(priors$alpha > 0, info = "alpha should be positive")
  expect_true(priors$beta > 0, info = "beta should be positive")
})

test_that("extract_indicator_priors returns Stochastic-Block parameters", {
  fit = get_bgms_fit_sbm()

  priors = extract_indicator_priors(fit)

  # Type check
  expect_type(priors, "list")
  expect_equal(priors$type, "Stochastic-Block")

  # Required parameters
  expect_true("beta_bernoulli_alpha" %in% names(priors),
    info = "should have beta_bernoulli_alpha"
  )
  expect_true("beta_bernoulli_beta" %in% names(priors),
    info = "should have beta_bernoulli_beta"
  )
  expect_true("dirichlet_alpha" %in% names(priors),
    info = "should have dirichlet_alpha"
  )
})

# ==============================================================================
# Legacy Format Support Tests
# ==============================================================================
# These tests verify backward compatibility with fit objects from older bgms versions.
# Legacy fixtures are stored in tests/testthat/fixtures/legacy/ (NOT shipped with package).
#
# To generate fixtures, run: Rscript tests/fixtures/generate_legacy_fixtures.R
#
# Tests skip on CRAN since fixtures aren't available in installed package.
#
# PATTERN: Unified fixture specs for both bgm and bgmCompare, mirroring get_extractor_fixtures()
# ==============================================================================
# Legacy Format Compatibility Tests
# ==============================================================================
#
# These tests verify backward compatibility with fit objects from older bgms
# versions. They require legacy fixture files (*.rds) that are:
#   - Generated by tests/fixtures/generate_legacy_fixtures.R
#   - Stored in tests/testthat/fixtures/legacy/
#   - NOT shipped to CRAN (excluded via .Rbuildignore)
#   - Skipped on CRAN via skip_on_cran() in get_legacy_dir()
#
# Format evolution:
#   bgm:
#     - pre-0.1.4: $gamma (defunct), $interactions, $thresholds
#     - 0.1.4-0.1.5: $indicator at top level (deprecated)
#     - 0.1.6+: $raw_samples$indicator (current)
#   bgmCompare:
#     - 0.1.4-0.1.5: $pairwise_difference_indicator, $interactions, $thresholds
#     - 0.1.6+: $raw_samples$indicator, $raw_samples$pairwise, $raw_samples$main
# ==============================================================================

# ------------------------------------------------------------------------------
# Legacy Fixture Infrastructure
# ------------------------------------------------------------------------------

# Get the legacy fixtures directory path
# NOTE: skip_on_cran() here ensures ALL legacy tests are skipped on CRAN
get_legacy_dir = function() {
  skip_on_cran() # Legacy fixtures not shipped to CRAN

  # Try relative path first (for devtools::test())
  legacy_dir = file.path("fixtures", "legacy")
  if(!dir.exists(legacy_dir)) {
    # Try from package root (for testthat::test_file())
    legacy_dir = file.path("tests", "testthat", "fixtures", "legacy")
  }
  if(!dir.exists(legacy_dir)) {
    return(NULL)
  }
  legacy_dir
}

# Load a legacy fixture by filename
load_legacy_fixture = function(filename) {
  legacy_dir = get_legacy_dir()
  if(is.null(legacy_dir)) {
    skip("Legacy fixtures directory not found - run tests/fixtures/generate_legacy_fixtures.R")
  }

  path = file.path(legacy_dir, paste0(filename, ".rds"))
  if(!file.exists(path)) {
    skip(paste("Legacy fixture not found:", filename, "- run tests/fixtures/generate_legacy_fixtures.R"))
  }
  readRDS(path)
}

# Categorize version by format era (works for both bgm and bgmCompare)
categorize_version = function(version, type = "bgm") {
  v = numeric_version(version)

  if(type == "bgm") {
    if(v < "0.1.4") {
      return("pre-0.1.4") # Defunct: $gamma field <U+2192> error
    } else if(v < "0.1.6") {
      return("0.1.4-0.1.5") # Deprecated: $indicator at top level <U+2192> warning
    } else {
      return("0.1.6+") # Current: $raw_samples$indicator <U+2192> no warning
    }
  } else {
    # bgmCompare (introduced in 0.1.4)
    if(v < "0.1.6") {
      return("0.1.4-0.1.5") # Deprecated: top-level fields <U+2192> warning
    } else {
      return("0.1.6+") # Current: $raw_samples$* <U+2192> no warning
    }
  }
}

# Build legacy fixture specs from available files
# Returns list of specs like get_extractor_fixtures(), with:
#   label, version, type (bgm/bgmCompare), era, get_fit
get_legacy_fixture_specs = function() {
  legacy_dir = get_legacy_dir()
  if(is.null(legacy_dir)) {
    return(list())
  }

  specs = list()

  # bgm fixtures: fit_v*.rds
  # Use local() to properly capture variables in closures
  bgm_files = list.files(legacy_dir, pattern = "^fit_v.*\\.rds$")
  for(file in bgm_files) {
    specs[[length(specs) + 1]] = local({
      version = gsub("^fit_v(.*)\\.rds$", "\\1", file)
      fn = gsub("\\.rds$", "", file)
      list(
        label = paste0("bgm_v", version),
        version = version,
        type = "bgm",
        era = categorize_version(version, "bgm"),
        get_fit = function() load_legacy_fixture(fn)
      )
    })
  }

  # bgmCompare fixtures: bgmcompare_v*.rds
  bgmcompare_files = list.files(legacy_dir, pattern = "^bgmcompare_v.*\\.rds$")
  for(file in bgmcompare_files) {
    specs[[length(specs) + 1]] = local({
      version = gsub("^bgmcompare_v(.*)\\.rds$", "\\1", file)
      fn = gsub("\\.rds$", "", file)
      list(
        label = paste0("bgmCompare_v", version),
        version = version,
        type = "bgmCompare",
        era = categorize_version(version, "bgmCompare"),
        get_fit = function() load_legacy_fixture(fn)
      )
    })
  }

  specs
}

# Helper to filter specs by type and/or era
filter_legacy_specs = function(specs, type = NULL, era = NULL) {
  Filter(function(spec) {
    type_match = is.null(type) || spec$type == type
    era_match = is.null(era) || spec$era == era
    type_match && era_match
  }, specs)
}

# ------------------------------------------------------------------------------
# Legacy Lifecycle Tests (Parameterized)
# ------------------------------------------------------------------------------

test_that("pre-0.1.4 bgm formats throw defunct errors for indicator extraction", {
  specs = filter_legacy_specs(get_legacy_fixture_specs(), type = "bgm", era = "pre-0.1.4")
  skip_if(length(specs) == 0, "No pre-0.1.4 bgm fixtures available")

  for(spec in specs) {
    fit = spec$get_fit()

    expect_error(extract_indicators(fit), "defunct",
      info = paste(spec$label, "extract_indicators should error (defunct)")
    )
    expect_error(extract_posterior_inclusion_probabilities(fit), "defunct",
      info = paste(spec$label, "extract_pip should error (defunct)")
    )
  }
})

test_that("pre-0.1.4 bgm formats emit deprecation warnings for pairwise/thresholds", {
  specs = filter_legacy_specs(get_legacy_fixture_specs(), type = "bgm", era = "pre-0.1.4")
  skip_if(length(specs) == 0, "No pre-0.1.4 bgm fixtures available")

  for(spec in specs) {
    fit = spec$get_fit()

    expect_warning(extract_pairwise_interactions(fit), "deprecated",
      info = paste(spec$label, "extract_pairwise should warn")
    )
    threshold_warnings = capture_warnings(extract_category_thresholds(fit))
    expect_true(
      any(grepl("extract_main_effects", threshold_warnings)),
      info = paste(spec$label, "extract_thresholds should warn about rename")
    )
    expect_true(
      any(grepl("deprecated", threshold_warnings)),
      info = paste(spec$label, "extract_thresholds should warn about legacy format")
    )
  }
})

test_that("0.1.4-0.1.5 formats emit deprecation warnings", {
  specs = filter_legacy_specs(get_legacy_fixture_specs(), era = "0.1.4-0.1.5")
  skip_if(length(specs) == 0, "No 0.1.4-0.1.5 fixtures available")

  for(spec in specs) {
    fit = spec$get_fit()

    expect_warning(extract_indicators(fit), "deprecated",
      info = paste(spec$label, "extract_indicators should warn")
    )
    expect_warning(extract_posterior_inclusion_probabilities(fit), "deprecated",
      info = paste(spec$label, "extract_pip should warn")
    )
    expect_warning(extract_pairwise_interactions(fit), "deprecated",
      info = paste(spec$label, "extract_pairwise should warn")
    )
    threshold_warnings = capture_warnings(extract_category_thresholds(fit))
    expect_true(
      any(grepl("extract_main_effects", threshold_warnings)),
      info = paste(spec$label, "extract_thresholds should warn about rename")
    )
    expect_true(
      any(grepl("deprecated", threshold_warnings)),
      info = paste(spec$label, "extract_thresholds should warn about legacy format")
    )

    # bgmCompare also has extract_group_params
    if(spec$type == "bgmCompare") {
      expect_warning(extract_group_params(fit), "deprecated",
        info = paste(spec$label, "extract_group_params should warn")
      )
    }
  }
})

test_that("0.1.6+ formats work without deprecation warnings", {
  specs = filter_legacy_specs(get_legacy_fixture_specs(), era = "0.1.6+")
  skip_if(length(specs) == 0, "No 0.1.6+ fixtures available")

  for(spec in specs) {
    fit = spec$get_fit()

    # expect_no_warning doesn't support info= parameter, so use labeled tests
    expect_no_warning(extract_indicators(fit))
    expect_no_warning(extract_posterior_inclusion_probabilities(fit))
    expect_no_warning(extract_pairwise_interactions(fit))
    expect_warning(extract_category_thresholds(fit), "extract_main_effects")
    expect_no_warning(extract_main_effects(fit))
  }
})

# ------------------------------------------------------------------------------
# Legacy Functional Tests (Parameterized)
# ------------------------------------------------------------------------------

test_that("extract_indicators works with deprecated formats (0.1.4-0.1.5)", {
  specs = filter_legacy_specs(get_legacy_fixture_specs(), era = "0.1.4-0.1.5")
  skip_if(length(specs) == 0, "No 0.1.4-0.1.5 fixtures available")

  for(spec in specs) {
    fit = spec$get_fit()
    result = suppressWarnings(extract_indicators(fit))

    expect_true(is.matrix(result), info = paste(spec$label, "should return matrix"))
    expect_true(nrow(result) > 0, info = paste(spec$label, "should have rows"))
    expect_true(ncol(result) > 0, info = paste(spec$label, "should have columns"))
    expect_true(all(result %in% c(0, 1)), info = paste(spec$label, "should have binary values"))
  }
})

test_that("extract_posterior_inclusion_probabilities works with deprecated formats", {
  specs = filter_legacy_specs(get_legacy_fixture_specs(), era = "0.1.4-0.1.5")
  skip_if(length(specs) == 0, "No 0.1.4-0.1.5 fixtures available")

  for(spec in specs) {
    fit = spec$get_fit()
    result = suppressWarnings(extract_posterior_inclusion_probabilities(fit))

    expect_true(is.matrix(result), info = paste(spec$label, "should return matrix"))
    expect_true(isSymmetric(result), info = paste(spec$label, "should be symmetric"))
    expect_true(all(result >= 0 & result <= 1), info = paste(spec$label, "should be in [0,1]"))
  }
})

test_that("extract_pairwise_interactions works with pre-0.1.6 formats", {
  specs = filter_legacy_specs(get_legacy_fixture_specs(), era = NULL)
  specs = Filter(function(s) s$era != "0.1.6+", specs)
  skip_if(length(specs) == 0, "No pre-0.1.6 fixtures available")

  for(spec in specs) {
    fit = spec$get_fit()
    result = suppressWarnings(extract_pairwise_interactions(fit))

    expect_true(is.matrix(result), info = paste(spec$label, "should return matrix"))
    expect_true(nrow(result) > 0, info = paste(spec$label, "should have rows"))
  }
})

test_that("extract_category_thresholds works with pre-0.1.6 formats", {
  specs = filter_legacy_specs(get_legacy_fixture_specs(), era = NULL)
  specs = Filter(function(s) s$era != "0.1.6+", specs)
  skip_if(length(specs) == 0, "No pre-0.1.6 fixtures available")

  for(spec in specs) {
    fit = spec$get_fit()
    result = suppressWarnings(extract_category_thresholds(fit))

    expect_true(is.matrix(result), info = paste(spec$label, "should return matrix"))
    expect_true(nrow(result) > 0, info = paste(spec$label, "should have rows"))
  }
})

test_that("extract_group_params works with deprecated bgmCompare formats", {
  specs = filter_legacy_specs(get_legacy_fixture_specs(), type = "bgmCompare", era = "0.1.4-0.1.5")
  skip_if(length(specs) == 0, "No 0.1.4-0.1.5 bgmCompare fixtures available")

  for(spec in specs) {
    fit = spec$get_fit()
    result = suppressWarnings(extract_group_params(fit))

    expect_type(result, "list")
    expect_true("main_effects_groups" %in% names(result),
      info = paste(spec$label, "should have main_effects_groups")
    )
    expect_true("pairwise_effects_groups" %in% names(result),
      info = paste(spec$label, "should have pairwise_effects_groups")
    )
    expect_equal(ncol(result$main_effects_groups), 2,
      info = paste(spec$label, "should have 2 groups")
    )
    expect_equal(ncol(result$pairwise_effects_groups), 2,
      info = paste(spec$label, "should have 2 groups")
    )
  }
})

test_that("extract_arguments works with all legacy versions", {
  specs = get_legacy_fixture_specs()
  skip_if(length(specs) == 0, "No legacy fixtures available")

  for(spec in specs) {
    fit = spec$get_fit()
    args = extract_arguments(fit)

    expect_type(args, "list")
    expect_true(
      "no_variables" %in% names(args) || "num_variables" %in% names(args),
      info = paste(spec$label, "should have variable count")
    )
    expect_true("data_columnnames" %in% names(args),
      info = paste(spec$label, "should have column names")
    )

    # bgmCompare specific
    if(spec$type == "bgmCompare") {
      expect_true("difference_selection" %in% names(args),
        info = paste(spec$label, "should have difference_selection")
      )
    }
  }
})


# ==============================================================================
# Value-Checking Tests for Association-Scale Extractors
# ==============================================================================
#
# These tests build synthetic bgms objects with known parameter values and
# verify that extractors produce correct output.
# ==============================================================================

# ------------------------------------------------------------------
# Helper: build a minimal synthetic bgms object
# ------------------------------------------------------------------
make_synthetic_bgms = function(associations, residual_variance = NULL,
                               is_continuous = FALSE, is_mixed = FALSE,
                               discrete_indices = NULL,
                               continuous_indices = NULL,
                               data_columnnames_discrete = NULL,
                               data_columnnames_continuous = NULL) {
  p = nrow(associations)
  obj = list(
    posterior_mean_pairwise = associations,
    posterior_mean_residual_variance = residual_variance,
    arguments = list(
      num_variables = p,
      num_cases = 100,
      data_columnnames = colnames(associations),
      is_continuous = is_continuous,
      is_mixed = is_mixed,
      discrete_indices = discrete_indices,
      continuous_indices = continuous_indices,
      data_columnnames_discrete = data_columnnames_discrete,
      data_columnnames_continuous = data_columnnames_continuous
    )
  )
  class(obj) = "bgms"
  obj
}


# ------------------------------------------------------------------------------
# extract_log_odds() Value Tests
# ------------------------------------------------------------------------------

test_that("extract_log_odds returns twice the associations for OMRF", {
  associations = matrix(c(0, 0.3, -0.2, 0.3, 0, 0.5, -0.2, 0.5, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(paste0("X", 1:3), paste0("X", 1:3))
  )
  fit = make_synthetic_bgms(associations)
  result = extract_log_odds(fit)
  expect_equal(result, 2 * associations)
})

test_that("extract_log_odds returns NULL for GGM", {
  associations = matrix(c(0, -0.5, -0.5, 0),
    nrow = 2,
    dimnames = list(c("Y1", "Y2"), c("Y1", "Y2"))
  )
  fit = make_synthetic_bgms(associations,
    residual_variance = c(1.2, 0.8),
    is_continuous = TRUE
  )
  result = extract_log_odds(fit)
  expect_null(result)
})

test_that("extract_log_odds extracts discrete block for mixed MRF", {
  # 2 discrete (d1, d2) + 1 continuous (c1)
  associations = matrix(0, 3, 3,
    dimnames = list(c("d1", "c1", "d2"), c("d1", "c1", "d2"))
  )
  associations[1, 3] = 0.4
  associations[3, 1] = 0.4 # pairwise_disc between d1 and d2
  associations[1, 2] = 0.1
  associations[2, 1] = 0.1 # cross
  associations[2, 3] = -0.2
  associations[3, 2] = -0.2 # cross

  fit = make_synthetic_bgms(associations,
    residual_variance = c(0.5),
    is_mixed = TRUE,
    discrete_indices = c(1, 3),
    continuous_indices = 2,
    data_columnnames_discrete = c("d1", "d2"),
    data_columnnames_continuous = "c1"
  )

  result = extract_log_odds(fit)
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(dimnames(result), list(c("d1", "d2"), c("d1", "d2")))
  # Value: 2 * pairwise_disc block
  expect_equal(result["d1", "d2"], 2 * 0.4)
  expect_equal(result["d2", "d1"], 2 * 0.4)
  expect_equal(result["d1", "d1"], 0)
  expect_equal(result["d2", "d2"], 0)
})


# ------------------------------------------------------------------------------
# extract_precision() Value Tests
# ------------------------------------------------------------------------------

test_that("extract_precision reconstructs precision = -2 * association for GGM", {
  # Known association matrix (A = -0.5 * precision)
  associations = matrix(c(0, -0.3, 0.1, -0.3, 0, 0.2, 0.1, 0.2, 0),
    nrow = 3,
    dimnames = list(paste0("Y", 1:3), paste0("Y", 1:3))
  )
  rv = c(0.5, 0.8, 0.6) # residual variance = 1/precision_ii
  names(rv) = paste0("Y", 1:3)

  fit = make_synthetic_bgms(associations, residual_variance = rv, is_continuous = TRUE)
  Theta = extract_precision(fit)

  # Off-diagonal: Theta_ij = -2 * A_ij (precision = -2 * association)
  expect_equal(Theta[1, 2], -2 * associations[1, 2])
  expect_equal(Theta[1, 3], -2 * associations[1, 3])
  expect_equal(Theta[2, 3], -2 * associations[2, 3])

  # Diagonal: precision_ii = 1/rv_i
  expect_equal(unname(diag(Theta)), unname(1 / rv))

  # Symmetric
  expect_equal(Theta, t(Theta))
})

test_that("extract_precision returns NULL for OMRF", {
  associations = matrix(c(0, 0.3, 0.3, 0),
    nrow = 2,
    dimnames = list(c("X1", "X2"), c("X1", "X2"))
  )
  fit = make_synthetic_bgms(associations)
  result = extract_precision(fit)
  expect_null(result)
})

test_that("extract_precision extracts continuous block for mixed MRF", {
  # 2 discrete + 2 continuous
  associations = matrix(0, 4, 4,
    dimnames = list(c("d1", "c1", "d2", "c2"), c("d1", "c1", "d2", "c2"))
  )
  associations[2, 4] = -0.25
  associations[4, 2] = -0.25 # pairwise_cont between c1, c2

  rv = c(0.5, 0.4)
  names(rv) = c("c1", "c2")

  fit = make_synthetic_bgms(associations,
    residual_variance = rv,
    is_mixed = TRUE,
    discrete_indices = c(1, 3),
    continuous_indices = c(2, 4),
    data_columnnames_discrete = c("d1", "d2"),
    data_columnnames_continuous = c("c1", "c2")
  )

  Theta = extract_precision(fit)

  expect_equal(nrow(Theta), 2)
  expect_equal(dimnames(Theta), list(c("c1", "c2"), c("c1", "c2")))

  # Off-diagonal: -2 * associations
  expect_equal(Theta["c1", "c2"], -2 * (-0.25))
  expect_equal(Theta["c2", "c1"], -2 * (-0.25))

  # Diagonal: 1/rv
  expect_equal(Theta["c1", "c1"], 1 / 0.5)
  expect_equal(Theta["c2", "c2"], 1 / 0.4)
})


# ------------------------------------------------------------------------------
# extract_partial_correlations() Value Tests
# ------------------------------------------------------------------------------

test_that("extract_partial_correlations derives from precision for GGM", {
  associations = matrix(c(0, -0.3, 0.1, -0.3, 0, 0.2, 0.1, 0.2, 0),
    nrow = 3,
    dimnames = list(paste0("Y", 1:3), paste0("Y", 1:3))
  )
  rv = c(0.5, 0.8, 0.6)
  names(rv) = paste0("Y", 1:3)

  fit = make_synthetic_bgms(associations, residual_variance = rv, is_continuous = TRUE)
  pcor = extract_partial_correlations(fit)
  Theta = extract_precision(fit)

  # rho_ij = -precision_ij / sqrt(precision_ii * precision_jj)
  for(i in 1:3) {
    for(j in 1:3) {
      if(i == j) {
        expect_equal(pcor[i, j], 1)
      } else {
        expected = -Theta[i, j] / sqrt(Theta[i, i] * Theta[j, j])
        expect_equal(pcor[i, j], expected,
          info = sprintf("pcor[%d,%d]", i, j)
        )
      }
    }
  }

  # Symmetric
  expect_equal(pcor, t(pcor))
})

test_that("extract_partial_correlations returns NULL for OMRF", {
  associations = matrix(c(0, 0.3, 0.3, 0),
    nrow = 2,
    dimnames = list(c("X1", "X2"), c("X1", "X2"))
  )
  fit = make_synthetic_bgms(associations)
  result = extract_partial_correlations(fit)
  expect_null(result)
})

test_that("extract_partial_correlations bounded in [-1, 1] for mixed MRF", {
  associations = matrix(0, 4, 4,
    dimnames = list(c("d1", "c1", "d2", "c2"), c("d1", "c1", "d2", "c2"))
  )
  associations[2, 4] = -0.25
  associations[4, 2] = -0.25

  rv = c(0.5, 0.4)
  names(rv) = c("c1", "c2")

  fit = make_synthetic_bgms(associations,
    residual_variance = rv,
    is_mixed = TRUE,
    discrete_indices = c(1, 3),
    continuous_indices = c(2, 4),
    data_columnnames_discrete = c("d1", "d2"),
    data_columnnames_continuous = c("c1", "c2")
  )

  pcor = extract_partial_correlations(fit)

  expect_equal(nrow(pcor), 2)
  expect_true(all(pcor >= -1 & pcor <= 1))
  expect_equal(diag(pcor), c(c1 = 1, c2 = 1))
})


# ------------------------------------------------------------------------------
# Mutual Consistency Tests
# ------------------------------------------------------------------------------

test_that("GGM extractors are mutually consistent: association -> precision -> pcor", {
  associations = matrix(c(0, -0.4, 0.15, -0.4, 0, -0.1, 0.15, -0.1, 0),
    nrow = 3,
    dimnames = list(paste0("Y", 1:3), paste0("Y", 1:3))
  )
  rv = c(0.5, 1.0, 0.8)
  names(rv) = paste0("Y", 1:3)

  fit = make_synthetic_bgms(associations, residual_variance = rv, is_continuous = TRUE)

  # association -> precision -> pcor chain
  Theta = extract_precision(fit)
  pcor = extract_partial_correlations(fit)
  log_odds = extract_log_odds(fit)

  # precision = -2 * association
  expect_equal(Theta[lower.tri(Theta)], -2 * associations[lower.tri(associations)])

  # pcor from Theta
  d = sqrt(diag(Theta))
  expected_pcor = -Theta / outer(d, d)
  diag(expected_pcor) = 1
  expect_equal(pcor, expected_pcor)

  # GGM has no log-odds
  expect_null(log_odds)
})

test_that("OMRF extractors: log_odds present, precision/pcor NULL", {
  associations = matrix(c(0, 0.3, 0.3, 0),
    nrow = 2,
    dimnames = list(c("X1", "X2"), c("X1", "X2"))
  )
  fit = make_synthetic_bgms(associations)

  expect_equal(extract_log_odds(fit), 2 * associations)
  expect_null(extract_precision(fit))
  expect_null(extract_partial_correlations(fit))
})

test_that("Mixed MRF: log_odds and precision from disjoint blocks", {
  associations = matrix(0, 4, 4,
    dimnames = list(c("d1", "c1", "d2", "c2"), c("d1", "c1", "d2", "c2"))
  )
  associations[1, 3] = 0.4
  associations[3, 1] = 0.4 # pairwise_disc
  associations[2, 4] = -0.3
  associations[4, 2] = -0.3 # pairwise_cont
  associations[1, 2] = 0.1
  associations[2, 1] = 0.1 # cross

  rv = c(0.5, 0.6)
  names(rv) = c("c1", "c2")

  fit = make_synthetic_bgms(associations,
    residual_variance = rv,
    is_mixed = TRUE,
    discrete_indices = c(1, 3),
    continuous_indices = c(2, 4),
    data_columnnames_discrete = c("d1", "d2"),
    data_columnnames_continuous = c("c1", "c2")
  )

  lo = extract_log_odds(fit)
  Theta = extract_precision(fit)
  pcor = extract_partial_correlations(fit)

  # Discrete block: log_odds = 2 * associations
  expect_equal(lo["d1", "d2"], 2 * 0.4)

  # Continuous block: precision = -2 * association
  expect_equal(Theta["c1", "c2"], -2 * (-0.3))

  # Dimensions are independent
  expect_equal(dim(lo), c(2, 2))
  expect_equal(dim(Theta), c(2, 2))
  expect_equal(dim(pcor), c(2, 2))
})


# ==============================================================================
# Convention / Relationship Tests on Real Fit Objects
# ==============================================================================
#
# Verify association-scale invariants hold on real bgm() output (not synthetic objects).
# ==============================================================================

test_that("GGM: associations diagonal is zero", {
  fit = get_bgms_fit_ggm()
  associations = fit$posterior_mean_pairwise
  expect_true(all(diag(associations) == 0))
})

test_that("GGM: residual variance is positive", {
  fit = get_bgms_fit_ggm()
  rv = fit$posterior_mean_residual_variance
  expect_true(all(rv > 0))
})

test_that("GGM: precision diagonal = 1/residual_variance", {
  fit = get_bgms_fit_ggm()
  Theta = extract_precision(fit)
  rv = fit$posterior_mean_residual_variance
  expect_equal(unname(diag(Theta)), unname(1 / rv))
})

test_that("GGM: precision is symmetric", {
  fit = get_bgms_fit_ggm()
  Theta = extract_precision(fit)
  expect_equal(Theta, t(Theta))
})

test_that("GGM: partial correlations in [-1, 1] with unit diagonal", {
  fit = get_bgms_fit_ggm()
  pcor = extract_partial_correlations(fit)
  expect_true(all(pcor >= -1 & pcor <= 1))
  expect_equal(unname(diag(pcor)), rep(1, nrow(pcor)))
})

test_that("OMRF: associations matrix is symmetric with zero diagonal", {
  fit = get_bgms_fit()
  associations = fit$posterior_mean_pairwise
  expect_equal(associations, t(associations))
  expect_true(all(diag(associations) == 0))
})

test_that("OMRF: log_odds = 2 * associations", {
  fit = get_bgms_fit()
  lo = extract_log_odds(fit)
  associations = fit$posterior_mean_pairwise
  expect_equal(lo, 2 * associations)
})

test_that("Mixed MRF: associations diagonal is zero", {
  fit = get_bgms_fit_mixed_mrf()
  associations = fit$posterior_mean_pairwise
  expect_true(all(diag(associations) == 0))
})

test_that("Mixed MRF: residual variance is positive", {
  fit = get_bgms_fit_mixed_mrf()
  rv = fit$posterior_mean_residual_variance
  expect_true(all(rv > 0))
})

test_that("Mixed MRF: precision diagonal = 1/residual_variance", {
  fit = get_bgms_fit_mixed_mrf()
  Theta = extract_precision(fit)
  rv = fit$posterior_mean_residual_variance
  expect_equal(unname(diag(Theta)), unname(1 / rv))
})

test_that("Mixed MRF: discrete block log_odds = 2 * associations", {
  fit = get_bgms_fit_mixed_mrf()
  args = extract_arguments(fit)
  associations = fit$posterior_mean_pairwise
  disc_idx = args$discrete_indices
  lo = extract_log_odds(fit)
  expected = 2 * associations[disc_idx, disc_idx]
  dimnames(expected) = dimnames(lo)
  expect_equal(lo, expected)
})

test_that("Mixed MRF: continuous precision = -2 * associations", {
  fit = get_bgms_fit_mixed_mrf()
  args = extract_arguments(fit)
  associations = fit$posterior_mean_pairwise
  cont_idx = args$continuous_indices
  Theta = extract_precision(fit)
  rv = fit$posterior_mean_residual_variance
  expected_offdiag = -2 * associations[cont_idx, cont_idx]
  diag(expected_offdiag) = unname(1 / rv)
  dimnames(expected_offdiag) = dimnames(Theta)
  expect_equal(Theta, expected_offdiag)
})


# ---------------------------------------------------------------------------
# SBM number-of-blocks summary: p(K | t) convention
# ---------------------------------------------------------------------------

test_that("p(K | t) matches the shifted-Poisson MFM generative prior", {
  skip_on_cran()

  q = 5L
  lambda = 1
  dirichlet_alpha = 1
  log_Vn = compute_Vn_mfm_sbm(q, dirichlet_alpha, q + 10L, lambda)

  # Simulate the sampler's generative prior: K - 1 ~ Poisson(lambda),
  # symmetric Dirichlet(alpha) weights, iid allocations; t = #occupied.
  set.seed(21)
  n_sims = 3e5
  K = rpois(n_sims, lambda) + 1L
  t_obs = vapply(K, function(k) {
    w = rgamma(k, dirichlet_alpha)
    length(unique(sample.int(k, q, replace = TRUE, prob = w)))
  }, integer(1))

  for(t0 in 1:3) {
    # The summary reports the conditional restricted to K <= q.
    sel = t_obs == t0 & K <= q
    empirical = tabulate(K[sel], nbins = q) / sum(sel)
    analytic = compute_p_k_given_t(t0, log_Vn, dirichlet_alpha, q, lambda)
    expect_lt(0.5 * sum(abs(empirical - analytic)), 0.02)
  }
})
