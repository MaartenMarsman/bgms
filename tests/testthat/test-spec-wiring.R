# ==============================================================================
# Unit tests for bgm_spec wiring into bgm() and bgmCompare()
# Phase B.7-B.8 of the R scaffolding refactor.
#
# Verifies that the .bgm_spec object attached by bgm()/bgmCompare() produces
# build_arguments() output matching the output$arguments list.
# ==============================================================================

# ==============================================================================
# Helper: compare build_arguments(spec) vs output$arguments
# ==============================================================================
#
# Compares all common fields except 'version' (always identical but created
# at different moments) and fields that legitimately differ between paths.
# ==============================================================================
assert_arguments_match = function(fit, skip = "version") {
  spec = fit$.bgm_spec
  if(is.null(spec)) {
    fail("No .bgm_spec attached to fit object")
    return(invisible())
  }
  spec_args = build_arguments(spec)
  out_args = fit$arguments

  common = intersect(names(spec_args), names(out_args))
  common = setdiff(common, skip)

  for(f in common) {
    info_msg = sprintf("field '%s'", f)
    expect_equal(spec_args[[f]], out_args[[f]], info = info_msg)
  }
  invisible(TRUE)
}


# ==============================================================================
# 1.  bgm() <U+2014> spec attached
# ==============================================================================

test_that("bgm() attaches .bgm_spec to output (binary/ordinal)", {
  fit = get_bgms_fit()
  expect_false(is.null(fit$.bgm_spec))
  expect_s3_class(fit$.bgm_spec, "bgm_spec")
})

test_that("bgm() attaches .bgm_spec to output (ordinal)", {
  fit = get_bgms_fit_ordinal()
  expect_false(is.null(fit$.bgm_spec))
  expect_s3_class(fit$.bgm_spec, "bgm_spec")
})


# ==============================================================================
# 2.  bgm() <U+2014> build_arguments matches output$arguments
# ==============================================================================

test_that("bgm() binary: build_arguments matches output$arguments", {
  fit = get_bgms_fit()
  assert_arguments_match(fit)
})

test_that("bgm() ordinal: build_arguments matches output$arguments", {
  fit = get_bgms_fit_ordinal()
  assert_arguments_match(fit)
})

test_that("bgm() blume-capel: build_arguments matches output$arguments", {
  fit = get_bgms_fit_blumecapel()
  assert_arguments_match(fit)
})

test_that("bgm() Beta-Bernoulli: build_arguments matches output$arguments", {
  fit = get_bgms_fit_beta_bernoulli()
  assert_arguments_match(fit)
})

test_that("bgm() SBM: build_arguments matches output$arguments", {
  fit = get_bgms_fit_sbm()
  assert_arguments_match(fit)
})

# ==============================================================================
# 3.  bgmCompare() <U+2014> spec attached
# ==============================================================================

test_that("bgmCompare() attaches .bgm_spec (group_indicator)", {
  fit = get_bgmcompare_fit()
  expect_false(is.null(fit$.bgm_spec))
  expect_s3_class(fit$.bgm_spec, "bgm_spec")
})

test_that("bgmCompare() attaches .bgm_spec (x/y interface)", {
  fit = get_bgmcompare_fit_xy()
  expect_false(is.null(fit$.bgm_spec))
  expect_s3_class(fit$.bgm_spec, "bgm_spec")
})


# ==============================================================================
# 4.  bgmCompare() <U+2014> build_arguments matches output$arguments
# ==============================================================================

test_that("bgmCompare() group_indicator: build_arguments matches", {
  fit = get_bgmcompare_fit()
  assert_arguments_match(fit)
})

test_that("bgmCompare() x/y: build_arguments matches", {
  fit = get_bgmcompare_fit_xy()
  assert_arguments_match(fit)
})

test_that("bgmCompare() ordinal: build_arguments matches", {
  fit = get_bgmcompare_fit_ordinal()
  assert_arguments_match(fit)
})

test_that("bgmCompare() Beta-Bernoulli: build_arguments matches", {
  fit = get_bgmcompare_fit_beta_bernoulli()
  assert_arguments_match(fit)
})


# ==============================================================================
# 5.  Spec model_type correctness
# ==============================================================================

test_that("bgm() binary <U+2192> model_type = 'omrf'", {
  fit = get_bgms_fit()
  expect_equal(fit$.bgm_spec$model_type, "omrf")
})

test_that("bgm() ordinal <U+2192> model_type = 'omrf'", {
  fit = get_bgms_fit_ordinal()
  expect_equal(fit$.bgm_spec$model_type, "omrf")
})

test_that("bgmCompare() <U+2192> model_type = 'compare'", {
  fit = get_bgmcompare_fit()
  expect_equal(fit$.bgm_spec$model_type, "compare")
})
