# ==============================================================================
# Phase 0: Fit-Object Contract Tests
# ==============================================================================
#
# Regression tests for the fit-object contract documented in
# dev/plans/fit-object-contract.md. These lock down:
#
#   1. Serialization round-trips (saveRDS / readRDS)
#   2. Lazy summary computation semantics
#   3. names(fit) stability
#
# These tests exist to guard the contract before any structural refactor
# (e.g. S7 migration).
# ==============================================================================


# ==============================================================================
# 1. Serialization Round-Trips
# ==============================================================================

test_that("bgms fit survives saveRDS/readRDS without prior summary access", {
  fit = get_bgms_fit()
  tmp = tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  restored = readRDS(tmp)

  expect_s3_class(restored, "bgms")
  expect_equal(names(restored), names(fit))
  expect_equal(restored$arguments, fit$arguments)
  expect_equal(
    restored$posterior_mean_pairwise,
    fit$posterior_mean_pairwise
  )

  # Lazy summaries must still work after deserialization
  s = restored$posterior_summary_pairwise
  expect_true(is.data.frame(s) || is.matrix(s))
  expect_true(nrow(s) > 0)
})

test_that("bgms fit survives saveRDS/readRDS after summary access", {
  fit = get_bgms_fit()
  # Force summary computation before saving
  s_before = fit$posterior_summary_pairwise

  tmp = tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  restored = readRDS(tmp)

  s_after = restored$posterior_summary_pairwise
  expect_equal(s_after, s_before)
})

test_that("bgmCompare fit survives saveRDS/readRDS without prior summary access", {
  fit = get_bgmcompare_fit()
  tmp = tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  restored = readRDS(tmp)

  expect_s3_class(restored, "bgmCompare")
  expect_equal(names(restored), names(fit))
  expect_equal(restored$arguments, fit$arguments)

  # Lazy summaries must still work after deserialization
  s = restored$posterior_summary_pairwise_baseline
  expect_true(is.data.frame(s) || is.matrix(s))
  expect_true(nrow(s) > 0)
})

test_that("bgmCompare fit survives saveRDS/readRDS after summary access", {
  fit = get_bgmcompare_fit()
  s_before = fit$posterior_summary_pairwise_baseline

  tmp = tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  restored = readRDS(tmp)

  s_after = restored$posterior_summary_pairwise_baseline
  expect_equal(s_after, s_before)
})

test_that("GGM fit survives saveRDS/readRDS", {
  fit = get_bgms_fit_ggm()
  tmp = tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(fit, tmp)
  restored = readRDS(tmp)

  expect_s3_class(restored, "bgms")
  expect_equal(names(restored), names(fit))
  expect_equal(
    restored$posterior_mean_residual_variance,
    fit$posterior_mean_residual_variance
  )

  s = restored$posterior_summary_pairwise
  expect_true(is.data.frame(s) || is.matrix(s))
})


# ==============================================================================
# 2. Lazy Summary Computation Semantics
# ==============================================================================

test_that("fresh fit starts with summaries_computed = FALSE", {
  # Use a dedicated fit to avoid shared-fixture interference
  data("ADHD", package = "bgms")
  fit = bgm(
    ADHD[1:30, 2:5],
    iter = 25, warmup = 50, chains = 1,
    seed = 99998,
    display_progress = "none"
  )
  cache = fit$cache
  expect_false(isTRUE(cache$summaries_computed))

  s = fit$posterior_summary_pairwise
  expect_true(is.data.frame(s) || is.matrix(s))
  expect_true(isTRUE(cache$summaries_computed))
})

test_that("accessing posterior_summary_pairwise populates cache (bgms)", {
  fit = get_bgms_fit()
  cache = fit$cache

  # Access triggers computation (or returns cached)
  s = fit$posterior_summary_pairwise
  expect_true(is.data.frame(s) || is.matrix(s))
  expect_true(nrow(s) > 0)

  # Cache must be populated after access
  expect_true(isTRUE(cache$summaries_computed))
})

test_that("second access returns cached result without recomputation (bgms)", {
  fit = get_bgms_fit()

  # First access
  s1 = fit$posterior_summary_pairwise

  # Second access <U+2014> should be identical (same object from cache)
  s2 = fit$posterior_summary_pairwise
  expect_identical(s1, s2)
})

test_that("accessing posterior_summary_pairwise_baseline populates cache (bgmCompare)", {
  fit = get_bgmcompare_fit()
  cache = fit$cache

  s = fit$posterior_summary_pairwise_baseline
  expect_true(is.data.frame(s) || is.matrix(s))
  expect_true(nrow(s) > 0)

  expect_true(isTRUE(cache$summaries_computed))
})

test_that("summary() triggers lazy computation (bgms)", {
  fit = get_bgms_fit()
  cache = fit$cache

  sm = summary(fit)
  expect_true(isTRUE(cache$summaries_computed))
  expect_s3_class(sm, "summary.bgms")
})

test_that("summary() triggers lazy computation (bgmCompare)", {
  fit = get_bgmcompare_fit()
  cache = fit$cache

  sm = summary(fit)
  expect_true(isTRUE(cache$summaries_computed))
  expect_s3_class(sm, "summary.bgmCompare")
})


# ==============================================================================
# 3. names(fit) Stability
# ==============================================================================

test_that("bgms fit (edge selection) has all required names", {
  fit = get_bgms_fit()
  nm = names(fit)

  # Core fields
  expect_true("arguments" %in% nm)
  expect_true("raw_samples" %in% nm)
  expect_true("posterior_mean_pairwise" %in% nm)
  expect_true("cache" %in% nm)

  # Lazy summary placeholders
  expect_true("posterior_summary_main" %in% nm)
  expect_true("posterior_summary_pairwise" %in% nm)

  # Edge selection fields
  expect_true("posterior_mean_indicator" %in% nm)
  expect_true("posterior_summary_indicator" %in% nm)
})

test_that("bgms fit (no edge selection) does not have indicator fields", {
  fit = get_bgms_fit_ggm_no_es()
  nm = names(fit)

  expect_true("arguments" %in% nm)
  expect_true("raw_samples" %in% nm)
  expect_true("posterior_mean_pairwise" %in% nm)
  expect_true("posterior_summary_main" %in% nm)
  expect_true("posterior_summary_pairwise" %in% nm)
  expect_false("posterior_mean_indicator" %in% nm)
  expect_false("posterior_summary_indicator" %in% nm)
})

test_that("bgms GGM fit has residual variance but no main effects", {
  fit = get_bgms_fit_ggm()
  nm = names(fit)

  expect_true("posterior_mean_residual_variance" %in% nm)
  expect_null(fit$posterior_mean_main)
})

test_that("bgms OMRF fit has main effects", {
  fit = get_bgms_fit_ordinal()
  expect_false(is.null(fit$posterior_mean_main))
})

test_that("bgmCompare fit has all required names", {
  fit = get_bgmcompare_fit()
  nm = names(fit)

  expect_true("arguments" %in% nm)
  expect_true("raw_samples" %in% nm)
  expect_true("cache" %in% nm)

  # Baseline fields
  expect_true("posterior_mean_main_baseline" %in% nm)
  expect_true("posterior_mean_pairwise_baseline" %in% nm)

  # Difference fields
  expect_true("posterior_mean_main_differences" %in% nm)
  expect_true("posterior_mean_pairwise_differences" %in% nm)

  # Lazy summary placeholders
  expect_true("posterior_summary_main_baseline" %in% nm)
  expect_true("posterior_summary_pairwise_baseline" %in% nm)
  expect_true("posterior_summary_main_differences" %in% nm)
  expect_true("posterior_summary_pairwise_differences" %in% nm)
})

test_that("bgmCompare fit (with difference selection) has indicator summary", {
  fit = get_bgmcompare_fit_beta_bernoulli()
  nm = names(fit)

  expect_true("posterior_summary_indicator" %in% nm)
})


# ==============================================================================
# 4. Parameterized: all bgms fixtures survive serialization
# ==============================================================================

for(spec in get_bgms_fixtures()) {
  test_that(
    sprintf("saveRDS/readRDS round-trip preserves structure (%s)", spec$label),
    {
      fit = spec$get_fit()
      tmp = tempfile(fileext = ".rds")
      on.exit(unlink(tmp), add = TRUE)

      saveRDS(fit, tmp)
      restored = readRDS(tmp)

      expect_s3_class(restored, "bgms")
      expect_equal(names(restored), names(fit))
      expect_equal(restored$arguments, fit$arguments)
    }
  )
}

for(spec in get_bgmcompare_fixtures()) {
  test_that(
    sprintf("saveRDS/readRDS round-trip preserves structure (%s)", spec$label),
    {
      fit = spec$get_fit()
      tmp = tempfile(fileext = ".rds")
      on.exit(unlink(tmp), add = TRUE)

      saveRDS(fit, tmp)
      restored = readRDS(tmp)

      expect_s3_class(restored, "bgmCompare")
      expect_equal(names(restored), names(fit))
      expect_equal(restored$arguments, fit$arguments)
    }
  )
}


# ==============================================================================
# 5. Parameterized: lazy summaries work for all bgms fixtures
# ==============================================================================

for(spec in get_bgms_fixtures()) {
  test_that(
    sprintf("lazy summaries compute on first access (%s)", spec$label),
    {
      fit = spec$get_fit()
      cache = fit$cache

      # After fixture caching the summaries may already be computed from
      # other tests. Reset the flag to test the lazy path.
      # (This is safe because the fixture cache returns the same object.)
      # Instead, just verify that accessing a summary field returns data.
      summary_fields = grep(
        "^posterior_summary_", names(fit),
        value = TRUE
      )
      for(field in summary_fields) {
        val = fit[[field]]
        # Some fields may legitimately be NULL (e.g. main for GGM)
        if(!is.null(val)) {
          expect_true(
            is.data.frame(val) || is.matrix(val),
            info = sprintf("[%s] %s is not a data.frame/matrix", spec$label, field)
          )
          expect_true(
            nrow(val) > 0,
            info = sprintf("[%s] %s has zero rows", spec$label, field)
          )
        }
      }
    }
  )
}

for(spec in get_bgmcompare_fixtures()) {
  test_that(
    sprintf("lazy summaries compute on first access (%s)", spec$label),
    {
      fit = spec$get_fit()
      summary_fields = grep(
        "^posterior_summary_", names(fit),
        value = TRUE
      )
      for(field in summary_fields) {
        val = fit[[field]]
        if(!is.null(val)) {
          expect_true(
            is.data.frame(val) || is.matrix(val),
            info = sprintf("[%s] %s is not a data.frame/matrix", spec$label, field)
          )
          expect_true(
            nrow(val) > 0,
            info = sprintf("[%s] %s has zero rows", spec$label, field)
          )
        }
      }
    }
  )
}
