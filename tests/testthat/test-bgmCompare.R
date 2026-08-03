# ==============================================================================
# Tests for bgmCompare() - Multi-group MRF Comparison
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Reproducibility, range invariants, dimension consistency
#
# These tests parallel the structure of test-bgm.R for consistency.
# Tests both reproducibility and basic output structure.
#
# INTEGRATION NOTE: Many configurations (Blume-Capel, missing data imputation,
# standardization) are tested via the parameterized fixture approach in
# test-methods.R. See:
#   - helper-fixtures.R: Cached fit functions (get_bgmcompare_fit_blumecapel, etc.)
#   - test-methods.R: get_bgmcompare_fixtures() loops over all configurations
#
# This file focuses on tests that require special setup or unique assertions.
# ==============================================================================

# ------------------------------------------------------------------------------
# Reproducibility Tests (using fixtures to save one model fit)
# ------------------------------------------------------------------------------

test_that("bgmCompare is reproducible with seed (x, y interface)", {
  # Use cached fixture as fit1, run one fresh fit as fit2 with same params
  fit1 = get_bgmcompare_fit_xy()

  data("Wenchuan", package = "bgms")
  x = Wenchuan[1:25, 1:4]
  y = Wenchuan[26:50, 1:4]

  fit2 = without_support_warning(bgmCompare(x = x, y = y, iter = 50, warmup = 100, chains = 2, seed = 1234, display_progress = "none"))

  combine_chains = function(fit) {
    pairs = do.call(rbind, fit$raw_samples$pairwise)
    mains = do.call(rbind, fit$raw_samples$main)
    cbind(mains, pairs)
  }

  expect_equal(combine_chains(fit1), combine_chains(fit2))
})


test_that("bgmCompare accepts its default update_method", {
  data("Wenchuan", package = "bgms")
  x = Wenchuan[1:20, 1:3]
  y = Wenchuan[21:40, 1:3]
  expect_error(
    without_support_warning(bgmCompare(
      x = x, y = y, iter = 20, warmup = 20, chains = 1, cores = 1,
      seed = 1, display_progress = "none", verbose = FALSE
    )),
    NA
  )
})


# ------------------------------------------------------------------------------
# Output Structure Tests (using saved fit)
# ------------------------------------------------------------------------------

test_that("bgmCompare output has expected structure", {
  fit = get_bgmcompare_fit()

  expect_s3_class(fit, "bgmCompare")

  # Should have key components
  expect_true("arguments" %in% names(fit))
  expect_true("raw_samples" %in% names(fit))

  # Raw samples should have required components
  expect_true("pairwise" %in% names(fit$raw_samples))
  expect_true("main" %in% names(fit$raw_samples))
})

test_that("bgmCompare stores correct number of groups", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)

  expect_true("num_groups" %in% names(args))
  expect_true(args$num_groups >= 2)
})

test_that("bgmCompare posterior summaries have expected format", {
  fit = get_bgmcompare_fit()

  # Should have baseline summaries
  expect_true(!is.null(fit$posterior_summary_pairwise_baseline))
  expect_true(!is.null(fit$posterior_mean_pairwise_baseline))
})


# ------------------------------------------------------------------------------
# Tolerance/Sanity Tests (Stochastic-robust)
# ------------------------------------------------------------------------------

test_that("bgmCompare outputs are numerically sane", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)
  p = args$num_variables

  # Check baseline pairwise
  M = fit$posterior_mean_pairwise_baseline

  expect_true(is.matrix(M))
  expect_equal(dim(M), c(p, p))

  # Symmetry check
  asym = max(abs(M - t(M)), na.rm = TRUE)
  expect_true(asym <= 1e-8, info = sprintf("Asymmetry too large: %g", asym))

  # Values should be finite
  expect_true(all(is.finite(M)))

  # Check group params
  group_params = extract_group_params(fit)

  expect_true(all(is.finite(group_params$main_effects_groups)))
  expect_true(all(is.finite(group_params$pairwise_effects_groups)))
})


# ------------------------------------------------------------------------------
# Fresh Fit Tests
# ------------------------------------------------------------------------------

test_that("bgmCompare without selection produces valid estimates", {
  data = generate_grouped_test_data(n_per_group = 20, p = 3, n_groups = 2, seed = 42)

  fit = bgmCompare(
    x = data$x,
    group_indicator = data$group_indicator,
    difference_selection = FALSE,
    iter = 50,
    warmup = 100,
    chains = 1,
    display_progress = "none"
  )

  expect_s3_class(fit, "bgmCompare")

  # Should have posterior means
  expect_true(!is.null(fit$posterior_mean_pairwise_baseline))
  expect_true(!is.null(fit$posterior_mean_main_baseline))
})

test_that("bgmCompare with selection produces valid indicators", {
  data = generate_grouped_test_data(n_per_group = 20, p = 3, n_groups = 2, seed = 123)

  # A single chain exercises the single-chain path of summarize_mixture_effect
  fit = bgmCompare(
    x = data$x,
    group_indicator = data$group_indicator,
    difference_selection = TRUE,
    iter = 50,
    warmup = 100,
    chains = 1,
    display_progress = "none"
  )

  expect_s3_class(fit, "bgmCompare")

  # Should have indicator samples
  expect_true(!is.null(fit$raw_samples$indicator))

  # Indicators should be binary
  ind_samples = do.call(rbind, fit$raw_samples$indicator)
  expect_true(all(ind_samples %in% c(0, 1)))
})


# ------------------------------------------------------------------------------
# Method Variations Tests
# ------------------------------------------------------------------------------

test_that("bgmCompare works with different update methods", {
  data = generate_grouped_test_data(n_per_group = 20, p = 3, n_groups = 2, seed = 99)

  methods_to_test = c("adaptive-metropolis")
  # Note: Could add "hmc", "nuts" if testing more thoroughly

  for(method in methods_to_test) {
    fit = tryCatch(
      without_support_warning(bgmCompare(
        x = data$x,
        group_indicator = data$group_indicator,
        update_method = method,
        iter = 25,
        warmup = 50,
        chains = 1,
        display_progress = "none"
      )),
      error = function(e) e
    )

    if(!inherits(fit, "error")) {
      expect_s3_class(fit, "bgmCompare")
    }
  }
})


# ------------------------------------------------------------------------------
# More Than Two Groups
# ------------------------------------------------------------------------------

test_that("bgmCompare handles more than 2 groups", {
  data = generate_grouped_test_data(
    n_per_group = 15, p = 3, n_groups = 3, seed = 456
  )

  # Use difference_selection = FALSE to avoid summary computation issues
  # with very short chains
  fit = bgmCompare(
    x = data$x,
    group_indicator = data$group_indicator,
    difference_selection = FALSE,
    iter = 25,
    warmup = 50,
    chains = 1,
    display_progress = "none"
  )

  expect_s3_class(fit, "bgmCompare")

  args = extract_arguments(fit)
  expect_equal(args$num_groups, 3)

  # Group-specific effects should have 3 columns
  group_params = extract_group_params(fit)
  expect_equal(ncol(group_params$main_effects_groups), 3)
  expect_equal(ncol(group_params$pairwise_effects_groups), 3)
})


# ==============================================================================
# Parameter Ordering Test (p >= 4 required to detect row/column-major bugs)
# ==============================================================================
#
# See test-bgm.R header comment and helper-fixtures.R for background on
# row-major vs column-major ordering bugs.
# ==============================================================================

test_that("bgmCompare output has correct parameter ordering", {
  skip_on_cran()

  data("Wenchuan", package = "bgms")
  x = na.omit(Wenchuan[, 1:5]) # p=5 to detect row/column-major bugs
  group_ind = rep(1:2, length.out = nrow(x))

  fit = bgmCompare(
    x = x, group_indicator = group_ind,
    difference_selection = TRUE,
    iter = 1000, warmup = 500, chains = 1,
    seed = 42,
    display_progress = "none"
  )

  # Summary mean vector -> matrix lower triangle (same row-major order)
  M = fit$posterior_mean_pairwise_baseline
  expect_true(
    all(abs(fit$posterior_summary_pairwise_baseline$mean - M[lower.tri(M)]) < 1e-10),
    info = "bgmCompare pairwise baseline summary means do not match matrix lower triangle"
  )

  # Extractor column means -> matrix positions (uses named "Vi-Vj" columns)
  pw_means = colMeans(extract_pairwise_interactions(fit))
  expect_true(
    all(check_extractor_matrix_consistency(
      pw_means, fit$posterior_mean_pairwise_baseline
    )),
    info = "bgmCompare extract_pairwise_interactions() names do not match matrix positions"
  )
})


# ------------------------------------------------------------------------------
# Stochastic-Block difference prior
# ------------------------------------------------------------------------------

test_that("bgmCompare accepts sbm_prior() and surfaces allocations", {
  data("Wenchuan", package = "bgms")
  x = na.omit(Wenchuan[, 1:4])
  group_ind = rep(1:2, length.out = nrow(x))

  fit = bgmCompare(
    x = x, group_indicator = group_ind,
    difference_prior = sbm_prior(),
    difference_selection = TRUE,
    iter = 100, warmup = 100, chains = 2,
    seed = 1, display_progress = "none"
  )

  expect_s3_class(fit, "bgmCompare")

  # Posterior allocation summaries are populated for SBM.
  expect_false(is.null(fit$posterior_mean_allocations))
  expect_length(fit$posterior_mean_allocations, ncol(x))
  expect_length(fit$posterior_mode_allocations, ncol(x))

  # Coclustering matrix is square with the right dim and 1's on the diagonal.
  cm = fit$posterior_mean_coclustering_matrix
  expect_equal(dim(cm), c(ncol(x), ncol(x)))
  expect_equal(unname(diag(cm)), rep(1, ncol(x)))
  expect_true(all(cm >= 0 & cm <= 1))

  # Per-iteration allocations land in raw_samples with shape iter x p.
  expect_equal(
    dim(fit$raw_samples$allocations[[1]]),
    c(nrow(fit$raw_samples$indicator[[1]]), ncol(x))
  )
})

test_that("bgmCompare without sbm_prior() does not produce allocation fields", {
  data("Wenchuan", package = "bgms")
  x = na.omit(Wenchuan[, 1:4])
  group_ind = rep(1:2, length.out = nrow(x))

  fit = bgmCompare(
    x = x, group_indicator = group_ind,
    difference_prior = beta_bernoulli_prior(alpha = 1, beta = 1),
    difference_selection = TRUE,
    iter = 100, warmup = 100, chains = 1,
    seed = 1, display_progress = "none"
  )

  expect_null(fit$posterior_mean_allocations)
  expect_null(fit$posterior_mode_allocations)
  expect_null(fit$posterior_mean_coclustering_matrix)
  expect_null(fit$raw_samples$allocations)
})

# ------------------------------------------------------------------------------
# Pairwise parameterization
# ------------------------------------------------------------------------------

test_that("bgmCompare pairwise effects are on the association scale", {
  # simulate_mrf() is the package's reference for the association scale: a
  # rest score carries 2 * omega * x. Both groups are drawn from the same
  # omega, so each group's estimate must recover omega, not 2 * omega.
  p = 3
  omega = matrix(0, p, p)
  omega[upper.tri(omega)] = c(0.5, 0.0, 0.45)
  omega = omega + t(omega)
  main = matrix(c(0, -0.5), nrow = p, ncol = 2, byrow = TRUE)

  draw = function(seed) {
    simulate_mrf(
      400, p, num_categories = 2, pairwise = omega, main = main,
      variable_type = "ordinal", iter = 50, seed = seed
    )
  }

  fit = bgmCompare(
    rbind(draw(11), draw(12)), group = rep(1:2, each = 400),
    iter = 600, warmup = 300, chains = 1, seed = 1234,
    difference_selection = FALSE, display_progress = "none"
  )

  groups = extract_group_params(fit)$pairwise_effects_groups
  target = omega[t(utils::combn(p, 2))]

  # The scale is read as the slope of the estimate on the target, which is
  # bounded from BOTH sides: a doubled parameterization gives slope 2 and a
  # halved one slope 0.5, and neither is inside the band. An rmse bound
  # against 2 * target alone tests only the doubling direction, and at this
  # data size it cannot see the halving one at all -- target and 0.5 * target
  # are 0.194 apart in rmse while the run's own error reaches 0.19.
  #
  # The band comes from the estimator's spread, not from taste: over 12
  # alternative data/fit seeds at this construction the slope stayed within
  # [0.851, 1.369] on both groups, so [0.6, 1.6] holds it with room while
  # still refusing 0.5 and 2.
  #
  # Both groups are drawn from the same omega, so both must recover it. Only
  # group 1 was ever checked, which left a difference parameterization that
  # mis-signs or mis-scales the second group's reconstruction unexamined.
  slope = function(estimate) sum(estimate * target) / sum(target^2)
  rmse = function(estimate, x) sqrt(mean((estimate - x)^2))

  for(g in 1:2) {
    estimate = groups[, g]
    ctx = paste("group", g)
    expect_lt(rmse(estimate, target), 0.2)
    expect_gt(slope(estimate), 0.6)
    expect_lt(slope(estimate), 1.6)
    # And the correct scale has to fit better than the doubled one, which is
    # the direction this data size does resolve.
    expect_lt(rmse(estimate, target), 0.8 * rmse(estimate, 2 * target),
      label = ctx)
  }
})

test_that("bgmCompare recovers a planted group difference at its planted size", {
  skip_on_cran()
  # The guard above fits two groups drawn from the SAME omega, so every
  # difference in it is zero and a difference parameterization off by a factor
  # would still pass. Nothing else in the every-run tier looks at the size of
  # a nonzero difference. This plants one and reads it back.
  #
  # The two group matrices swap their two nonzero pairs, so the planted
  # difference is large (+-0.5) while both groups stay in a well-identified
  # coupling range -- planting a large difference by inflating one group
  # instead pushes that group toward deterministic data, where the posterior
  # is wide and the read-back is noisier than the difference being measured.
  p = 3
  pairs = t(utils::combn(p, 2))
  symmetric = function(values) {
    m = matrix(0, p, p)
    m[upper.tri(m)] = values
    m + t(m)
  }
  omega_1 = symmetric(c(0.60, 0.00, 0.10))
  omega_2 = symmetric(c(0.10, 0.00, 0.60))
  main = matrix(c(0, -0.5), nrow = p, ncol = 2, byrow = TRUE)

  draw = function(omega, seed) {
    simulate_mrf(
      1200, p, num_categories = 2, pairwise = omega, main = main,
      variable_type = "ordinal", iter = 50, seed = seed
    )
  }

  fit = bgmCompare(
    rbind(draw(omega_1, 101), draw(omega_2, 201)),
    group = rep(1:2, each = 1200),
    iter = 400, warmup = 250, chains = 1, seed = 1,
    difference_selection = FALSE, display_progress = "none"
  )

  groups = extract_group_params(fit)$pairwise_effects_groups
  recovered = groups[, 2] - groups[, 1]
  planted = omega_2[pairs] - omega_1[pairs]

  # Tolerances from the estimator's own spread: over 12 alternative
  # data/fit seeds at this construction the rmse of the recovered difference
  # reached 0.089 and no single pair was off by more than 0.125, while the
  # slope stayed within [0.813, 1.167]. The bounds below sit above those and
  # below what a mis-scaled parameterization produces: at half or double the
  # difference the slope is 0.5 or 2, and the rmse against the planted
  # difference is 0.204 or 0.408.
  expect_lt(sqrt(mean((recovered - planted)^2)), 0.15)
  expect_lt(max(abs(recovered - planted)), 0.20)

  slope = sum(recovered * planted) / sum(planted^2)
  expect_gt(slope, 0.7)
  expect_lt(slope, 1.4)

  # The signs are the qualitative half of the same claim: the pair the second
  # group loses and the pair it gains must come back with opposite signs, and
  # the pair with no planted difference must not acquire one.
  expect_lt(recovered[1], 0)
  expect_gt(recovered[3], 0)
  expect_lt(abs(recovered[2]), 0.20)
})


# Weekly certification (T2): this is the only test in the suite that runs
# bgmCompare() at its shipped defaults on a full shipped dataset, and it costs
# ~14 min on the 2-core CI runner -- by itself a quarter of the nightly budget.
# A full-defaults end-to-end fit is a product-surface check, which the tier
# contract puts in T1, but at that size it is not a smoke. The cheap end of the
# same surface stays local: the label-propagation test below fits the same data
# at iter = 50.
test_that("the shipped data's own language column works as the group indicator", {
  skip_on_cran()
  skip_unless_certification()
  data("Boredom", package = "bgms")
  fit = bgmCompare(Boredom[, -1], group_indicator = Boredom$language)
  expect_s3_class(fit, "bgmCompare")
  expect_equal(
    tabulate(extract_arguments(fit)$group),
    tabulate(match(Boredom$language, unique(Boredom$language)))
  )
})


test_that("predict.bgmCompare centers Blume-Capel variables at the fit's baseline", {
  fit = get_bgmcompare_fit_blumecapel()
  arguments = extract_arguments(fit)
  # The fit's own baseline (category 3, shifted to the 0-based scale) must be
  # stored; without it predict() silently centered every Blume-Capel term at 0.
  expect_equal(arguments$baseline_category, rep(2L, 4L))

  data("Boredom", package = "bgms")
  newdata = Boredom[c(1:4, 494:497), 2:5]
  probs = predict(fit, newdata = newdata, group = 1)

  # Manual reference, the sampler's own convention: category c contributes
  # exp(lin*(c-ref) + quad*(c-ref)^2 + (c-ref)*rest), with the rest score
  # summing 2 * (x_v - ref_v) * pairwise[v, j] over the other variables.
  shift = arguments$blume_capel_shift
  ref = arguments$baseline_category
  x0 = sweep(data.matrix(newdata), 2, shift)
  gp = extract_group_params(fit)
  p = arguments$num_variables
  pw = matrix(0, p, p)
  pw[lower.tri(pw)] = gp$pairwise_effects_groups[, 1]
  pw = pw + t(pw)
  main = matrix(gp$main_effects_groups[, 1], ncol = 2, byrow = TRUE)
  for(j in seq_len(p)) {
    rest = as.numeric((sweep(x0[, -j, drop = FALSE], 2, ref[-j])) %*% (2 * pw[-j, j]))
    cats = 0:arguments$num_categories[j] - ref[j]
    expected = t(vapply(rest, function(r) {
      e = exp(main[j, 1] * cats + main[j, 2] * cats^2 + cats * r)
      e / sum(e)
    }, numeric(length(cats))))
    expect_equal(unname(probs[[j]]), unname(expected), tolerance = 1e-10)
  }
})


# ---- group labels on human displays (F-072) ---------------------------------
# Groups are numbered by first appearance in the indicator and every extractor
# keys on that number; the original labels ride along so the displays a person
# reads can name the group as well as number it.

test_that("a compare fit stores the group indicator's own labels", {
  data("Boredom", package = "bgms")
  # The shipped data is fr-first, so first-appearance numbering makes fr group 1
  # even though "en" sorts first -- exactly the confusion the labels remove.
  expect_identical(Boredom$language[1], "fr")
  fit = bgmCompare(
    x = Boredom[, 2:5], group_indicator = Boredom$language,
    iter = 50, warmup = 100, chains = 2, seed = 8, display_progress = "none"
  )
  arguments = extract_arguments(fit)
  expect_identical(arguments$group_labels, c("fr", "en"))
  expect_equal(tabulate(arguments$group), c(490L, 496L))

  # print() and summary() name the groups; the numbers stay the key.
  expect_output(print(fit), "groups: 1 = fr \\(n = 490\\), 2 = en \\(n = 496\\)")
  expect_output(
    print(summary(fit)),
    "groups: 1 = fr \\(n = 490\\), 2 = en \\(n = 496\\)"
  )

  # Centrality labels carry them too.
  expect_match(
    attr(extract_centrality(fit, group = 2), "label"), "group 2 (en)",
    fixed = TRUE
  )
  expect_match(
    attr(extract_centrality(fit, group = c(1, 2)), "label"),
    "(group 1 (fr) - group 2 (en))",
    fixed = TRUE
  )

  # The extractor contract is numeric and must not move: easybgm and JASP read
  # these column names.
  expect_identical(
    colnames(extract_group_params(fit)$pairwise_effects_groups),
    c("group1", "group2")
  )
})

test_that("the x/y path labels the groups x and y", {
  fit = get_bgmcompare_fit_xy()
  expect_identical(extract_arguments(fit)$group_labels, c("x", "y"))
  # The counts are the fit's own, after listwise deletion -- which is the point
  # of reporting them next to the labels.
  n = tabulate(extract_arguments(fit)$group)
  expect_output(
    print(fit),
    sprintf("groups: 1 = x \\(n = %d\\), 2 = y \\(n = %d\\)", n[1], n[2])
  )
  expect_match(
    attr(extract_centrality(fit, group = 1), "label"), "group 1 (x)",
    fixed = TRUE
  )
})

test_that("a fit without the stored labels degrades to bare group numbers", {
  # Fits made before the field existed have no group_labels; every display path
  # goes through these two helpers, so this is the whole degradation contract.
  legacy = extract_arguments(get_bgmcompare_fit())
  legacy$group_labels = NULL

  expect_null(compare_group_labels(legacy))
  expect_null(group_mapping_line(legacy))
  expect_identical(group_tag(NULL, 2L), "group 2")
  expect_identical(group_tag(c("fr", "en"), 2L), "group 2 (en)")

  # A label vector that cannot name every group is refused wholesale rather
  # than half-applied.
  expect_null(compare_group_labels(list(group_labels = "only-one"), 2L))
  expect_null(compare_group_labels(list(group_labels = c("a", NA))))
})
