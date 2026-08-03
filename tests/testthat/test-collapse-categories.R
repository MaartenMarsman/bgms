# ==============================================================================
# Tests for collapse_categories_across_groups()
#
# Semantics (F-075): an ordinal variable is recoded onto the contiguous UNION
# of the category values observed in any group. A category observed by at
# least one group is always retained; only a true gap -- a value no group
# observes -- is dropped, and the survivors are renumbered contiguously.
# Retaining a category that some group never observes makes that cell a
# structural zero, which is reported by a warning. Blume-Capel variables are
# exempt.
#
# Before the F-075 fix the function kept only categories observed in EVERY
# group, so cases 2, 3, 4, 5 and 6 below asserted merged categories. Their
# expectations changed with the semantics, not because they were wrong about
# the old code.
# ==============================================================================


# ==============================================================================
# 1. Nothing to do --- all categories in all groups
# ==============================================================================

test_that("no recoding when all categories appear in all groups", {
  x = matrix(c(
    0, 1, 2, 0, 1, 2,
    0, 1, 2, 0, 1, 2
  ), nrow = 6, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(2, 2)
  bc = c(0, 0)

  expect_silent(
    result <- collapse_categories_across_groups(
      x, group, is_ordinal,
      num_categories, bc
    )
  )
  expect_equal(result$x, x)
  expect_equal(result$num_categories, c(2, 2))
  expect_equal(result$baseline_category, bc)
  # Support table is reported for every ordinal variable.
  expect_equal(dim(result$category_support[[1]]), c(3L, 2L))
  expect_equal(unname(result$category_support[[1]][, 1]), c(1L, 1L, 1L))
})


# ==============================================================================
# 2. A category missing from one group is RETAINED, and warned about
# ==============================================================================

test_that("category missing from one group is retained, not merged", {
  # Variable 1: group 1 has {0,1,2}, group 2 has {0,2}. Category 1 is observed
  # by group 1, so it survives; group 2's cell for it is a structural zero.
  x = matrix(c(
    0, 1, 2, 0, 2, 0,
    0, 1, 2, 0, 1, 2
  ), nrow = 6, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(2, 2)
  bc = c(0, 0)

  expect_warning(
    result <- collapse_categories_across_groups(
      x, group, is_ordinal,
      num_categories, bc
    ),
    "category 1, group 2"
  )
  expect_equal(result$x[, 1], c(0, 1, 2, 0, 2, 0))
  expect_equal(result$num_categories[1], 2)
  # Variable 2 unchanged and unremarked (all categories in both groups)
  expect_equal(result$x[, 2], c(0, 1, 2, 0, 1, 2))
  expect_equal(result$num_categories[2], 2)

  sup = result$category_support[[1]]
  expect_equal(unname(sup[, 1]), c(1L, 1L, 1L)) # group 1: 0, 1, 2
  expect_equal(unname(sup[, 2]), c(2L, 0L, 1L)) # group 2: no category 1
})

test_that("several categories missing from one group are all retained", {
  # Variable: group 1 has {0,1,2,3}, group 2 has {0,3}. Categories 1 and 2 are
  # observed by group 1 and survive.
  x = matrix(c(0, 1, 2, 3, 0, 3, 0, 1, 0, 1, 0, 1), nrow = 6, ncol = 2)
  group = c(1, 1, 1, 1, 2, 2)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(3, 1)
  bc = c(0, 0)

  expect_warning(
    result <- collapse_categories_across_groups(
      x, group, is_ordinal,
      num_categories, bc
    ),
    "category 1, group 2"
  )
  expect_equal(result$x[, 1], c(0, 1, 2, 3, 0, 3))
  expect_equal(result$num_categories[1], 3)
  expect_equal(unname(result$category_support[[1]][, 2]), c(1L, 0L, 0L, 1L))
})


# ==============================================================================
# 3. Blume-Capel variables pass through unchanged, and are never warned about
# ==============================================================================

test_that("BC variables are not affected by the cross-group pass", {
  # Variable 1: BC, values 0,1,2,3 (group 2 observes only 0 and 1)
  # Variable 2: ordinal, values 0,1,2
  x = matrix(c(
    0, 1, 2, 3, 0, 1,
    0, 1, 2, 0, 2, 0
  ), nrow = 6, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2)
  is_ordinal = c(FALSE, TRUE)
  num_categories = c(3, 2)
  bc = c(1, 0)

  expect_warning(
    result <- collapse_categories_across_groups(
      x, group, is_ordinal,
      num_categories, bc
    ),
    "variable 2" # only the ordinal variable is reported
  )
  # BC variable unchanged: values, count, baseline, and no support table
  expect_equal(result$x[, 1], c(0, 1, 2, 3, 0, 1))
  expect_equal(result$num_categories[1], 3)
  expect_equal(result$baseline_category[1], 1)
  expect_null(result$category_support[[1]])
  # Ordinal variable: group 2 lacks category 1, which is retained regardless
  expect_equal(result$x[, 2], c(0, 1, 2, 0, 2, 0))
  expect_equal(result$num_categories[2], 2)
})

test_that("a BC variable alone raises nothing", {
  x = matrix(c(
    0, 1, 2, 3, 0, 1,
    0, 1, 2, 3, 0, 1
  ), nrow = 6, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2)

  expect_silent(
    result <- collapse_categories_across_groups(
      x, group,
      is_ordinal = c(FALSE, FALSE),
      num_categories = c(3, 3),
      baseline_category = c(1, 1)
    )
  )
  expect_equal(result$x, x)
  expect_equal(result$num_categories, c(3, 3))
})


# ==============================================================================
# 4. More than two groups: the union runs over all K
# ==============================================================================

test_that("three groups take the union, not the intersection", {
  # Variable 1: group 1 has {0,1,2}, group 2 has {0,1}, group 3 has {0,2}.
  # Only category 0 is in all three groups; under the union all three survive.
  x = matrix(c(
    0, 1, 2, 0, 1, 0, 2, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 1, 0
  ), nrow = 9, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(2, 1)
  bc = c(0, 0)

  expect_warning(
    result <- collapse_categories_across_groups(
      x, group, is_ordinal,
      num_categories, bc
    ),
    "category 2, group 2"
  )
  expect_equal(result$x, x)
  expect_equal(result$num_categories, c(2, 1))
  sup = result$category_support[[1]]
  expect_equal(dim(sup), c(3L, 3L))
  expect_equal(unname(sup[, 2]), c(2L, 1L, 0L)) # group 2 lacks category 2
  expect_equal(unname(sup[, 3]), c(2L, 0L, 1L)) # group 3 lacks category 1
})

test_that("three groups with all categories present", {
  x = matrix(c(
    0, 1, 2, 0, 1, 2, 0, 1, 2,
    0, 1, 0, 0, 1, 0, 0, 1, 0
  ), nrow = 9, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(2, 1)
  bc = c(0, 0)

  expect_silent(
    result <- collapse_categories_across_groups(
      x, group, is_ordinal,
      num_categories, bc
    )
  )
  expect_equal(result$x, x)
  expect_equal(result$num_categories, c(2, 1))
})


# ==============================================================================
# 5. Disjoint supports: still the union, no error
# ==============================================================================

test_that("groups with disjoint supports keep every observed category", {
  # Group 1 has {0,1}, group 2 has {2,3}. The intersection is empty; under the
  # union all four categories survive and every cell of one group is zero.
  x = matrix(c(
    0, 1, 2, 3,
    0, 1, 0, 1
  ), nrow = 4, ncol = 2)
  group = c(1, 1, 2, 2)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(3, 1)
  bc = c(0, 0)

  expect_warning(
    result <- collapse_categories_across_groups(
      x, group, is_ordinal,
      num_categories, bc
    ),
    "category 2, group 1"
  )
  expect_equal(result$x[, 1], c(0, 1, 2, 3))
  expect_equal(result$num_categories[1], 3)
  expect_equal(unname(result$category_support[[1]][, 1]), c(1L, 1L, 0L, 0L))
  expect_equal(unname(result$category_support[[1]][, 2]), c(0L, 0L, 1L, 1L))
})

test_that("a variable with a single observed value is still an error", {
  x = matrix(c(
    2, 2, 2, 2,
    0, 1, 0, 1
  ), nrow = 4, ncol = 2)
  expect_error(
    collapse_categories_across_groups(
      x,
      group = c(1, 1, 2, 2), is_ordinal = c(TRUE, TRUE),
      num_categories = c(0, 1), baseline_category = c(0, 0)
    ),
    "Only one value was observed for variable 1"
  )
})


# ==============================================================================
# 6. True gaps: a value NO group observes is dropped, with a message
# ==============================================================================

test_that("a value observed in no group is dropped and the rest renumbered", {
  # Values 0, 1, 3, 4 are observed; 2 is a true gap. The four survivors are
  # renumbered 0, 1, 2, 3.
  withr::local_options(bgms.verbose = TRUE)
  x = cbind(
    v1 = c(rep(0, 5), rep(1, 5), rep(3, 5), rep(4, 5)),
    v2 = rep(0:1, 10)
  )
  group = rep(1:2, each = 10)

  expect_message(
    suppressWarnings(
      result <- collapse_categories_across_groups(
        x, group,
        is_ordinal = c(TRUE, TRUE),
        num_categories = c(4, 1), baseline_category = c(0, 0)
      )
    ),
    "not used by any group"
  )
  expect_equal(sort(unique(result$x[, 1])), 0:3)
  expect_equal(result$num_categories[1], 3)
  # The observed values keep their order and none of them are merged.
  expect_equal(result$x[, 1], rep(0:3, each = 5))
})

test_that("a gap common to all groups is dropped without touching support", {
  # Both groups observe {0, 2}; 1 is a true gap. No structural zeros here, so
  # the message fires and the warning does not.
  withr::local_options(bgms.verbose = TRUE)
  x = cbind(
    v1 = c(0, 0, 2, 2, 0, 0, 2, 2),
    v2 = c(0, 1, 0, 1, 0, 1, 0, 1)
  )
  group = rep(1:2, each = 4)

  expect_message(
    result <- collapse_categories_across_groups(
      x, group,
      is_ordinal = c(TRUE, TRUE),
      num_categories = c(2, 1), baseline_category = c(0, 0)
    ),
    "No observed category was merged"
  )
  expect_equal(result$x[, 1], c(0, 0, 1, 1, 0, 0, 1, 1))
  expect_equal(result$num_categories[1], 1)
  expect_false(any(result$category_support[[1]] == 0L))
})


# ==============================================================================
# 7. The report-16 reproducer, as a unit test of the new semantics
#
# Four categories, 40 observations each; group 1 spans {0,1,2} and group 2
# spans {1,2,3}. Under the pre-fix intersection rule this variable became
# BINARY (num_categories 1) with the recode 0->0 1->0 2->1 3->1, merging 80
# well-observed group-1 observations into one category. Under the union rule
# all four levels survive untouched.
# ==============================================================================

test_that("report-16 reproducer: four levels survive, nothing is merged", {
  g1 = cbind(v1 = rep(0:2, each = 40), v2 = rep(0:3, times = 30))
  g2 = cbind(v1 = rep(1:3, each = 40), v2 = rep(0:3, times = 30))
  x = rbind(g1, g2)
  group = rep(1:2, each = 120)

  expect_warning(
    result <- collapse_categories_across_groups(
      x = x, group = group,
      is_ordinal = c(TRUE, TRUE),
      num_categories = c(3, 3),
      baseline_category = c(0L, 0L)
    ),
    "not used by every group"
  )

  # Four levels, not two.
  expect_equal(result$num_categories[1], 3)
  # No observation moved: the recode is the identity.
  expect_identical(result$x, x)
  # Both groups keep all 40 observations in each category they use.
  sup = result$category_support[[1]]
  expect_equal(unname(sup[, 1]), c(40L, 40L, 40L, 0L))
  expect_equal(unname(sup[, 2]), c(0L, 40L, 40L, 40L))
  # v2 is observed everywhere and stays out of it.
  expect_equal(result$num_categories[2], 3)
  expect_false(any(result$category_support[[2]] == 0L))
})

test_that("the structural-zero warning names variable, category, and group", {
  g1 = cbind(v1 = rep(0:2, each = 40), v2 = rep(0:3, times = 30))
  g2 = cbind(v1 = rep(1:3, each = 40), v2 = rep(0:3, times = 30))
  x = rbind(g1, g2)
  group = rep(1:2, each = 120)

  w = tryCatch(
    collapse_categories_across_groups(
      x = x, group = group, is_ordinal = c(TRUE, TRUE),
      num_categories = c(3, 3), baseline_category = c(0L, 0L)
    ),
    warning = conditionMessage
  )
  expect_match(w, "variable 'v1', category 3, group 1", fixed = TRUE)
  expect_match(w, "variable 'v1', category 0, group 2", fixed = TRUE)
  expect_false(grepl("v2", w, fixed = TRUE))
  # and says what it means for the user, in words
  expect_match(w, "set by the prior, not")
})


# ==============================================================================
# 8. Integration: the compare data pipeline end to end
# ==============================================================================

test_that("full pipeline retains categories one group never observes", {
  # Variable 1: group 1 has {1,2,3}, group 2 has {1,3}. After the ordinal
  # recode both live on {0,1,2}; category 1 is observed by group 1 only and
  # survives.
  x = matrix(c(
    1, 2, 3, 1, 2, 3, 1, 3, 1, 3,
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 10, ncol = 2)
  group = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2)
  variable_bool = c(TRUE, TRUE)
  bc = c(0, 0)

  md = validate_missing_data(x,
    na_action = "listwise", is_continuous = FALSE,
    group = group
  )
  ord = reformat_ordinal_data(md$x,
    is_ordinal = variable_bool,
    baseline_category = bc
  )
  expect_warning(
    result <- collapse_categories_across_groups(
      x = ord$x, group = md$group, is_ordinal = variable_bool,
      num_categories = ord$num_categories,
      baseline_category = ord$baseline_category
    ),
    "category 1, group 2"
  )

  # Variable 1: 1->0, 2->1, 3->2, with nothing merged
  expect_equal(result$x[1:6, 1], c(0, 1, 2, 0, 1, 2))
  expect_equal(result$x[7:10, 1], c(0, 2, 0, 2))
  expect_equal(result$num_categories[1], 2)

  # Variable 2: 0,1 both in all groups -> unchanged
  expect_equal(result$x[, 2], c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1))
  expect_equal(result$num_categories[2], 1)
})

test_that("full pipeline preserves BC variables", {
  x = matrix(c(
    0, 1, 2, 3, 0, 1, 2, 3,
    0, 1, 2, 0, 0, 1, 2, 0
  ), nrow = 8, ncol = 2)
  group = c(1, 1, 1, 1, 2, 2, 2, 2)
  variable_bool = c(FALSE, TRUE)
  bc = c(1, 0)

  md = validate_missing_data(x,
    na_action = "listwise", is_continuous = FALSE,
    group = group
  )
  ord = reformat_ordinal_data(md$x,
    is_ordinal = variable_bool,
    baseline_category = bc
  )
  result = collapse_categories_across_groups(
    x = ord$x, group = md$group, is_ordinal = variable_bool,
    num_categories = ord$num_categories,
    baseline_category = ord$baseline_category
  )

  # BC variable: no cross-group pass, just the standard ordinal recode
  expect_equal(result$x[, 1], c(0, 1, 2, 3, 0, 1, 2, 3))
  expect_equal(result$num_categories[1], 3)
  expect_equal(result$baseline_category[1], 1)

  # Ordinal variable: 0,1,2 all in both groups -> unchanged
  expect_equal(result$x[, 2], c(0, 1, 2, 0, 0, 1, 2, 0))
  expect_equal(result$num_categories[2], 2)
})


# ==============================================================================
# 9. bgmCompare() surfaces the support table on the fitted object
# ==============================================================================

test_that("bgmCompare() reports a value no group used", {
  # reformat_ordinal_data() runs before the cross-group pass and has already
  # closed the gap by then, so this message has to come from the compare spec
  # builder, which still has the supplied values. Value 2 is never used here.
  withr::local_options(bgms.verbose = TRUE)
  x = cbind(
    a = c(rep(0, 20), rep(1, 20), rep(3, 20)),
    b = rep(0:2, 20), c = rep(0:1, 30)
  )
  group = rep(1:2, each = 30)

  expect_message(
    suppressWarnings(
      spec <- bgm_spec(x = x, group_indicator = group, model = "compare")
    ),
    "not used by any group"
  )
  expect_equal(spec$data$num_categories, c(2L, 2L, 1L))
  # 0 -> 0, 1 -> 1, 3 -> 2: nothing merged, the gap simply closed.
  expect_equal(unname(spec$data$category_levels[[1]]), c(0, 1, 2))
  expect_equal(names(spec$data$category_levels[[1]]), c("0", "1", "3"))
})

test_that("a scale that merely starts above zero is not called a gap", {
  # Values 1..5 are contiguous; recoding them to 0..4 is an offset, not a
  # dropped category, and is not worth telling anyone about.
  withr::local_options(bgms.verbose = TRUE)
  x = cbind(a = rep(1:5, 24), b = rep(1:3, 40), c = rep(1:2, 60))
  group = rep(1:2, each = 60)

  expect_no_message(
    spec <- bgm_spec(x = x, group_indicator = group, model = "compare")
  )
  expect_equal(spec$data$num_categories, c(4L, 2L, 1L))
})

test_that("bgmCompare() records the per-group support and warns once", {
  set.seed(19)
  n = 60
  mk = function(lo, hi) {
    cbind(
      a = sample(lo:hi, n, TRUE), b = sample(0:2, n, TRUE),
      c = sample(0:2, n, TRUE)
    )
  }
  x = rbind(mk(0, 1), mk(1, 2))
  group = rep(1:2, each = n)

  expect_warning(
    fit <- bgmCompare(
      x, group_indicator = group, iter = 60, warmup = 60,
      chains = 1, cores = 1, seed = 3, display_progress = "none",
      verbose = FALSE
    ),
    "not used by every group"
  )

  args = extract_arguments(fit)
  expect_equal(args$num_categories[1], 2)
  sup = args$category_support[[1]]
  expect_equal(dim(sup), c(3L, 2L))
  expect_equal(unname(sup[3, 1]), 0L) # group 1 never uses category 2
  expect_equal(unname(sup[1, 2]), 0L) # group 2 never uses category 0
  expect_equal(sum(sup), 2L * n)
})


# ==============================================================================
# 10. Estimate survival (T1) --- the assertion this file never had
#
# The mechanics cases above check what the recode does; none of them checks
# that the resulting fit is any good, which is why a function that turned a
# five-level variable into a binary one passed its own test suite (report 16,
# finding 16-3).
#
# Two groups drawn from the SAME network, then given deliberately different
# response supports: group 1 never uses the top option, group 2 never the
# bottom one. The union is all four options; the intersection is two, so
# before the F-075 fix each variable was reduced to a BINARY one. The
# assertion is that bgmCompare()'s group-specific pairwise estimates agree
# with a matched pair of bgm() fits on the same rows.
#
# Tolerance is derived, not chosen: the two estimates are posterior means of
# the same quantity from the same data, so their difference is compared
# against sqrt(sd_compare^2 + sd_bgm^2), the scale of one posterior standard
# deviation of the difference under independence. Three of those is a loose
# bound for two fits that should agree closely; it is far inside what the
# pre-fix collapse produced.
# ==============================================================================

test_that("group-differing support: pairwise estimates match matched bgm() fits", {
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to run the estimate-survival check (T1)"
  )

  p = 5
  K = 3
  n = 400
  set.seed(75)
  omega = matrix(0, p, p)
  omega[lower.tri(omega)] = c(
    0.3258, 0.1583, 0.1089, 0.0631, 0.2457,
    0.0605, 0.0554, 0.0497, 0.0970, 0.3859
  )
  omega = omega + t(omega)
  main = rbind(
    c(0.4919, -1.8439, -4.7586),
    c(-0.7551, -4.1471, -7.7251),
    c(-0.4009, -3.2720, -6.6058),
    c(0.0876, -2.1691, -5.1081),
    c(-0.7857, -3.7272, -7.2905)
  )

  x1 = simulate_mrf(n, p,
    num_categories = K, pairwise = omega, main = main,
    variable_type = "ordinal", iter = 200, seed = 41
  )
  x2 = simulate_mrf(n, p,
    num_categories = K, pairwise = omega, main = main,
    variable_type = "ordinal", iter = 200, seed = 42
  )
  x1[x1 == 3] = 2 # group 1 never uses the top option
  x2[x2 == 0] = 1 # group 2 never uses the bottom option
  expect_setequal(unique(as.vector(x1)), 0:2)
  expect_setequal(unique(as.vector(x2)), 1:3)

  expect_warning(
    fit <- bgmCompare(
      rbind(x1, x2),
      group_indicator = rep(1:2, each = n),
      interaction_prior = cauchy_prior(1),
      iter = 800, warmup = 800, chains = 2, cores = 2, seed = 91,
      display_progress = "none", verbose = FALSE
    ),
    "not used by every group"
  )

  args = extract_arguments(fit)
  # All four response options survive: 3 free thresholds, not 1.
  expect_equal(args$num_categories, rep(3L, p))
  expect_equal(unname(args$category_support[[1]][4, 1]), 0L)
  expect_equal(unname(args$category_support[[1]][1, 2]), 0L)

  fit1 = bgm(x1,
    edge_selection = FALSE, interaction_prior = cauchy_prior(1),
    iter = 800, warmup = 800, chains = 2, cores = 2, seed = 92,
    display_progress = "none", verbose = FALSE
  )
  fit2 = bgm(x2,
    edge_selection = FALSE, interaction_prior = cauchy_prior(1),
    iter = 800, warmup = 800, chains = 2, cores = 2, seed = 93,
    display_progress = "none", verbose = FALSE
  )

  draws1 = extract_pairwise_interactions(fit1)
  draws2 = extract_pairwise_interactions(fit2)

  # Group-specific compare draws: baseline + projection weight * difference.
  raw = bgms:::get_raw_samples(fit)$pairwise
  np = ncol(raw[[1]]) / args$num_groups
  base = do.call(rbind, lapply(raw, function(m) m[, seq_len(np), drop = FALSE]))
  diff = do.call(rbind, lapply(raw, function(m) m[, np + seq_len(np), drop = FALSE]))
  proj = args$projection
  cmp1 = base + proj[1, 1] * diff
  cmp2 = base + proj[2, 1] * diff

  z = function(cmp, bgm_draws) {
    (colMeans(cmp) - colMeans(bgm_draws)) /
      sqrt(apply(cmp, 2, var) + apply(bgm_draws, 2, var))
  }
  expect_lt(max(abs(z(cmp1, draws1))), 3)
  expect_lt(max(abs(z(cmp2, draws2))), 3)
  # Not one lucky edge: the typical disagreement is well inside one posterior
  # standard deviation too.
  expect_lt(mean(abs(z(cmp1, draws1))), 1)
  expect_lt(mean(abs(z(cmp2, draws2))), 1)

  # And they track each other edge by edge, not merely on average.
  expect_gt(cor(colMeans(cmp1), colMeans(draws1)), 0.9)
  expect_gt(cor(colMeans(cmp2), colMeans(draws2)), 0.9)
})
