# ==============================================================================
# easybgm S3 compatibility mode
# ==============================================================================
#
# bgm() and bgmCompare() have two user-visible return shapes. The default is
# the S7 object; the second appears when easybgm < 0.5.0 is loaded, because
# those versions overwrite class(fit), which an S7 object does not survive.
# needs_easybgm_s3_compat() (R/build_output.R) decides between them, and the
# builders return the plain S3 list instead of calling the converter
# s3_list_to_bgms() / s3_list_to_bgmCompare() (R/class_s7.R).
#
# Two things are tested here:
#   1. the gate -- easybgm absent, easybgm >= 0.5.0, easybgm < 0.5.0;
#   2. field parity between the two shapes, in both directions, so that a
#      field added to one and forgotten in the other fails loudly.
#
# The compat shape can only be produced by a fit BUILT while the gate is TRUE,
# so the session fixtures in helper-fixtures.R (all S7) cannot supply it. The
# two fits below therefore mirror the cheapest fixtures (get_bgms_fit /
# get_bgmcompare_fit) at chains = 1 and are cached for the file.
# ==============================================================================


# ------------------------------------------------------------------------------
# File-local compat fits
# ------------------------------------------------------------------------------

.compat_cache = new.env(parent = emptyenv())

#' @description Get a cached bgm() fit built in easybgm S3 compatibility mode.
get_compat_bgm_list = function() {
  if(is.null(.compat_cache$bgm)) {
    local_mocked_bindings(needs_easybgm_s3_compat = function() TRUE)
    data("ADHD", package = "bgms")
    .compat_cache$bgm = suppressWarnings(bgm(
      ADHD[1:50, 2:5], # 4 binary symptom variables
      iter = 50, warmup = 100, chains = 1,
      seed = 12345,
      display_progress = "none"
    ))
  }
  .compat_cache$bgm
}

#' @description Get a cached bgmCompare() fit built in compatibility mode.
get_compat_bgmcompare_list = function() {
  if(is.null(.compat_cache$compare)) {
    local_mocked_bindings(needs_easybgm_s3_compat = function() TRUE)
    data("ADHD", package = "bgms")
    .compat_cache$compare = suppressWarnings(bgmCompare(
      x = ADHD[, 2:5],
      group_indicator = ADHD[, "group"],
      iter = 50, warmup = 100, chains = 1,
      seed = 54321,
      display_progress = "none"
    ))
  }
  .compat_cache$compare
}

#' @description Names of the properties a converter has to set explicitly.
#' Getter-backed properties are computed from the cache and are excluded, as
#' is .field_names, which the converter derives from names(results).
settable_properties = function(s7_class) {
  props = s7_class@properties
  computed = vapply(props, function(p) !is.null(p$getter), logical(1))
  setdiff(names(props)[!computed], ".field_names")
}

#' @description Assert that one field is carried across unchanged: names,
#' classes, dims, and value.
#'
#' The comparison is at accessor level, which is the level the contract lives
#' at. Both shapes store the posterior_summary_* fields as NULL and fill them
#' from the cache on first read -- the S7 object through its property getters,
#' the S3 list through `[[.bgms` (R/methods_bgms.R) -- so `[[` is what has to
#' agree, not the stored slot.
expect_field_carried = function(s3_value, s7_value, label) {
  expect_identical(class(s7_value), class(s3_value), label = label)
  expect_identical(dim(s7_value), dim(s3_value), label = label)
  expect_identical(names(s7_value), names(s3_value), label = label)
  expect_identical(dimnames(s7_value), dimnames(s3_value), label = label)
  expect_identical(s7_value, s3_value, label = label)
}


# ==============================================================================
# 1. The gate
# ==============================================================================

test_that("needs_easybgm_s3_compat() is FALSE when easybgm is not loaded", {
  # The gate reads loadedNamespaces(), not the installed set: easybgm may well
  # be installed on the machine running the suite, and that is not the trigger.
  # The real set is captured BEFORE the mock is installed -- reading it from
  # inside the replacement would call the replacement.
  loaded = loadedNamespaces()
  local_mocked_bindings(
    loadedNamespaces = function() setdiff(loaded, "easybgm"),
    .package = "base"
  )
  expect_false(expect_no_warning(bgms:::needs_easybgm_s3_compat()))
})

test_that("needs_easybgm_s3_compat() is FALSE for easybgm >= 0.5.0", {
  # 0.5.0 moved easybgm to the extractor functions and it no longer overwrites
  # class(fit), so the S7 object travels unharmed and no compat shape is built.
  loaded = loadedNamespaces()
  local_mocked_bindings(
    loadedNamespaces = function() union(loaded, "easybgm"),
    .package = "base"
  )
  local_mocked_bindings(
    packageVersion = function(...) package_version("0.5.0"),
    .package = "utils"
  )
  expect_false(expect_no_warning(bgms:::needs_easybgm_s3_compat()))

  local_mocked_bindings(
    packageVersion = function(...) package_version("1.2.0"),
    .package = "utils"
  )
  expect_false(expect_no_warning(bgms:::needs_easybgm_s3_compat()))
})

test_that("needs_easybgm_s3_compat() warns and is TRUE for easybgm < 0.5.0", {
  loaded = loadedNamespaces()
  local_mocked_bindings(
    loadedNamespaces = function() union(loaded, "easybgm"),
    .package = "base"
  )
  local_mocked_bindings(
    packageVersion = function(...) package_version("0.4.0"),
    .package = "utils"
  )
  # The warning names the version it found and the version to move to; both
  # are what makes it actionable, so both are asserted.
  expect_warning(
    result <- bgms:::needs_easybgm_s3_compat(),
    regexp = "easybgm 0\\.4\\.0 is not compatible.*0\\.5\\.0 or later"
  )
  expect_true(result)
})


# ==============================================================================
# 2. The compatibility return shape
# ==============================================================================

test_that("bgm() returns a plain classed list in compatibility mode", {
  fit = get_compat_bgm_list()
  # A base list, NOT an S7 object: this is exactly what makes `class(fit) <- .`
  # in old easybgm survivable.
  expect_identical(typeof(fit), "list")
  expect_identical(class(fit), "bgms")
  expect_false(inherits(fit, "S7_object"))

  # The read paths still work on it, and give what the S7 shape gives: the
  # same fit through the converter has to be indistinguishable to a caller.
  expect_identical(class(summary(fit)), "summary.bgms")
  s7 = bgms:::s3_list_to_bgms(fit)
  expect_identical(
    extract_pairwise_interactions(fit), extract_pairwise_interactions(s7)
  )
  expect_identical(extract_indicators(fit), extract_indicators(s7))
  expect_identical(extract_arguments(fit), extract_arguments(s7))
})

test_that("bgmCompare() returns a plain classed list in compatibility mode", {
  fit = get_compat_bgmcompare_list()
  expect_identical(typeof(fit), "list")
  expect_identical(class(fit), "bgmCompare")
  expect_false(inherits(fit, "S7_object"))
  expect_identical(class(summary(fit)), "summary.bgmCompare")

  s7 = bgms:::s3_list_to_bgmCompare(fit)
  expect_identical(extract_main_effects(fit), extract_main_effects(s7))
  expect_identical(extract_arguments(fit), extract_arguments(s7))
})


# ==============================================================================
# 3. Field parity: the S3 list vs the S7 object it converts to
# ==============================================================================

test_that("s3_list_to_bgms() carries every field of a real compat fit", {
  s3 = get_compat_bgm_list()
  s7 = bgms:::s3_list_to_bgms(s3)

  expect_true(inherits(s7, "S7_object"))
  # .field_names is the record of what the builder produced; it is what the
  # print and summary methods key off, so it has to be the list's own names.
  expect_identical(S7::prop(s7, ".field_names"), names(s3))

  known = names(bgms:::bgms_class@properties)
  for(nm in names(s3)) {
    # Direction 1: a field the builder produced that the S7 class does not
    # declare would be silently dropped by the converter.
    expect_true(nm %in% known, label = paste0("bgms_class declares '", nm, "'"))
    expect_field_carried(s3[[nm]], S7::prop(s7, nm), label = nm)
  }
})

test_that("s3_list_to_bgmCompare() carries every field of a real compat fit", {
  s3 = get_compat_bgmcompare_list()
  s7 = bgms:::s3_list_to_bgmCompare(s3)

  expect_true(inherits(s7, "S7_object"))
  expect_identical(S7::prop(s7, ".field_names"), names(s3))

  known = names(bgms:::bgmCompare_class@properties)
  for(nm in names(s3)) {
    expect_true(
      nm %in% known,
      label = paste0("bgmCompare_class declares '", nm, "'")
    )
    expect_field_carried(s3[[nm]], S7::prop(s7, nm), label = nm)
  }
})


# ==============================================================================
# 4. Field parity: every declared property reaches the converter
# ==============================================================================
#
# The walk above only sees the fields a particular fit happens to populate --
# SBM allocations, residual variances and the deprecated easybgm aliases are
# all absent from a plain binary edge-selection fit. A sentinel list with one
# entry per settable property closes that hole: it fails the moment a property
# is added to the S7 class and forgotten in the converter, whatever the
# fixture looks like.
# ------------------------------------------------------------------------------

test_that("s3_list_to_bgms() sets every settable property of bgms_class", {
  fields = settable_properties(bgms:::bgms_class)
  sentinels = lapply(fields, function(nm) structure(nm, class = "compat_sentinel"))
  names(sentinels) = fields

  s7 = bgms:::s3_list_to_bgms(sentinels)

  for(nm in fields) {
    expect_identical(
      S7::prop(s7, nm), sentinels[[nm]],
      label = paste0("s3_list_to_bgms() carries '", nm, "'")
    )
  }
})

test_that("s3_list_to_bgmCompare() sets every settable property", {
  fields = settable_properties(bgms:::bgmCompare_class)
  sentinels = lapply(fields, function(nm) structure(nm, class = "compat_sentinel"))
  names(sentinels) = fields

  s7 = bgms:::s3_list_to_bgmCompare(sentinels)

  for(nm in fields) {
    expect_identical(
      S7::prop(s7, nm), sentinels[[nm]],
      label = paste0("s3_list_to_bgmCompare() carries '", nm, "'")
    )
  }
})
