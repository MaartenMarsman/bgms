# ==============================================================================
# Marking prior-only rows in the bgmCompare summary (F-117)
#
# bgmCompare() keeps the union of the ordinal categories the groups observe, so
# a retained category can have no observations at all in some group. That
# group's threshold for the category is then set by the prior, not by the data,
# and the reported threshold difference is large and very uncertain without
# being evidence of anything. The fit warns about it once; the printed summary
# marks the affected rows so that the reader of a saved fit, who never saw the
# warning, is told as well.
#
# A row is marked on either of two conditions: its own cell is empty, or the
# variable's REFERENCE cell (category 0) is empty in some group. The second
# condition marks every threshold of that variable, because every threshold is
# identified relative to category 0 and a group that never used it has nothing
# fixing the level of its threshold vector.
#
# The snapshot drops the MCMC numbers (see `labels_only()`): what is under test
# is the printed FORMAT -- which rows carry the mark, where the mark sits, and
# the footnote text -- and posterior means are not stable enough across
# platforms to belong in a snapshot.
# ==============================================================================

# Drop the trailing run of numeric cells from every table row, keeping the row
# name, the mark gutter and the parameter label. Headers, the truncation line
# and the notes all end in a word or a ")", so they survive whole.
labels_only = function(lines) {
  sub("(\\s+-?[0-9.]+)+\\s*$", "", lines)
}

# The `n` table rows printed under one block heading (the heading is followed
# by the column header, then the rows), trimmed of the alignment padding.
block_rows = function(out, heading, n) {
  start = which(out == heading)
  stopifnot(length(start) == 1L)
  trimws(out[seq.int(start + 2L, start + 1L + n)])
}

# Two groups, three ordinal variables. Variable A is four-level, and the groups
# split its support: group 1 uses {0, 1, 2}, group 2 uses {1, 2, 3}. The union
# keeps all four, so A has an empty cell in each group. B and C are observed on
# {0, 1, 2} by both groups and must stay unmarked.
adverse_compare_fit = function() {
  set.seed(117)
  n = 80
  draw = function(lo, hi) {
    cbind(
      A = sample(lo:hi, n, TRUE),
      B = sample(0:2, n, TRUE),
      C = sample(0:2, n, TRUE)
    )
  }
  x1 = draw(0, 2)
  x2 = draw(1, 3)
  x2[, "B"] = sample(0:2, n, TRUE)
  x2[, "C"] = sample(0:2, n, TRUE)

  bgmCompare(
    rbind(x1, x2),
    group_indicator = rep(1:2, each = n),
    iter = 200, warmup = 200, chains = 1, cores = 1, seed = 117,
    display_progress = "none", verbose = FALSE
  )
}

shared_support_compare_fit = function() {
  set.seed(118)
  n = 80
  draw = function() {
    cbind(
      A = sample(0:3, n, TRUE),
      B = sample(0:2, n, TRUE),
      C = sample(0:2, n, TRUE)
    )
  }
  bgmCompare(
    rbind(draw(), draw()),
    group_indicator = rep(1:2, each = n),
    iter = 200, warmup = 200, chains = 1, cores = 1, seed = 118,
    display_progress = "none", verbose = FALSE
  )
}


FOOTNOTE = paste0(
  "* a group lacks observations in this category or in the reference ",
  "category; the estimate reflects the prior, not the data"
)


test_that("an empty group-by-category cell marks its difference rows", {
  # The classed condition is named rather than muffled wholesale, so an
  # unrelated warning from the fit would still fail this test.
  expect_warning(
    fit <- adverse_compare_fit(),
    class = "bgms_group_support_warning"
  )

  support = extract_arguments(fit)$category_support
  # The fixture is what it claims to be: A loses category 3 in group 1 and
  # category 0 in group 2; B and C are complete.
  expect_equal(dim(support[[1]]), c(4L, 2L))
  expect_equal(unname(support[[1]][4L, 1L]), 0L)
  expect_equal(unname(support[[1]][1L, 2L]), 0L)
  expect_false(any(support[[2]] == 0L))
  expect_false(any(support[[3]] == 0L))

  out = capture.output(print(summary(fit)))
  rows = block_rows(out, "Group differences (main effects):", 6L)

  # All three of A's thresholds are marked, not just threshold 3. Threshold 3
  # has its own empty cell in group 1; thresholds 1 and 2 are marked because
  # group 2 never uses the reference category, which unfixes the level of its
  # whole threshold vector. B and C are backed by observations throughout.
  expect_identical(
    grepl("^\\*", rows),
    c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  )
  expect_true(all(startsWith(
    sub("^\\* ", "", rows[1:3]), paste0("A (diff1; ", 1:3, ")")
  )))

  expect_true(any(out == FOOTNOTE))

  # The mark is a display device: the data frame behind it is untouched.
  expect_false(any(grepl("*", summary(fit)$main_diff$parameter, fixed = TRUE)))

  expect_snapshot(print(summary(fit)), transform = labels_only)
})


test_that("a fit whose groups share their support prints unmarked", {
  expect_no_warning(fit <- shared_support_compare_fit())

  support = extract_arguments(fit)$category_support
  expect_false(any(vapply(support, function(s) any(s == 0L), logical(1))))

  out = capture.output(print(summary(fit)))
  rows = block_rows(out, "Group differences (main effects):", 6L)

  expect_false(any(grepl("^\\*", rows)))
  expect_false(any(startsWith(out, "* a group lacks")))
  # No mark anywhere means no gutter either: the labels print as they are.
  expect_true(all(startsWith(
    rows, head(summary(fit)$main_diff$parameter, 6L)
  )))
})


test_that("the mark follows the row layout of both difference summaries", {
  # The two summarizers lay the difference rows out differently: with selection
  # the contrasts of one threshold are adjacent (variable-major), without it
  # the whole main-effects matrix repeats once per contrast (contrast-major).
  # Two groups cannot tell the two orderings apart -- there is only one
  # contrast -- so this uses three, and checks the mark against the labels the
  # marking code itself never looks at.
  set.seed(219)
  n = 70
  draw = function(lo, hi) {
    cbind(
      A = sample(lo:hi, n, TRUE),
      B = sample(0:2, n, TRUE),
      C = sample(0:2, n, TRUE)
    )
  }
  x = rbind(draw(0, 2), draw(1, 3), draw(0, 3))
  group = rep(1:3, each = n)

  for(selection in c(TRUE, FALSE)) {
    expect_warning(
      fit <- bgmCompare(
        x, group_indicator = group, difference_selection = selection,
        iter = 200, warmup = 200, chains = 1, cores = 1, seed = 9,
        display_progress = "none", verbose = FALSE
      ),
      class = "bgms_group_support_warning"
    )

    labels = summary(fit)$main_diff$parameter
    flags = compare_prior_only_main_diff(
      extract_arguments(fit), length(labels)
    )
    # Group 2 never uses A's reference category, so every threshold of A is
    # marked, in every contrast. B and C are untouched.
    expect_identical(flags, grepl("^A \\(diff\\d+; \\d+\\)$", labels))

    rows = block_rows(
      capture.output(print(summary(fit))),
      "Group differences (main effects):", 6L
    )
    expect_identical(grepl("^\\*", rows), head(flags, 6L))
  }
})


test_that("Blume-Capel rows are never marked", {
  # Blume-Capel variables are exempt from the union recode -- their two
  # parameters are functions of the category SCORE, so an unobserved score is
  # still a meaningful point on the scale -- and they carry no support matrix.
  # The exemption is by construction, which is exactly why it wants pinning:
  # the mapping walks all variables and has to hand back two unmarked rows for
  # a variable whose `category_support` entry is NULL.
  set.seed(311)
  n = 90
  draw = function(lo, hi) {
    cbind(
      A = sample(lo:hi, n, TRUE),   # ordinal, split support -> marked
      D = sample(lo:hi, n, TRUE),   # same data, Blume-Capel -> exempt
      B = sample(0:2, n, TRUE)
    )
  }
  x = rbind(draw(0, 2), draw(1, 3))
  x[(n + 1L):(2L * n), "B"] = sample(0:2, n, TRUE)

  expect_warning(
    fit <- bgmCompare(
      x, group_indicator = rep(1:2, each = n),
      variable_type = c("ordinal", "blume-capel", "ordinal"),
      baseline_category = 0L,
      iter = 200, warmup = 200, chains = 1, cores = 1, seed = 311,
      display_progress = "none", verbose = FALSE
    ),
    class = "bgms_group_support_warning"
  )

  args = extract_arguments(fit)
  expect_false(args$is_ordinal_variable[2L])
  # D's support is never computed: the recode skips Blume-Capel variables.
  expect_null(args$category_support[[2L]])
  # A's is, and it is adverse in both directions.
  expect_equal(unname(args$category_support[[1L]][4L, 1L]), 0L)
  expect_equal(unname(args$category_support[[1L]][1L, 2L]), 0L)

  labels = summary(fit)$main_diff$parameter
  flags = compare_prior_only_main_diff(args, length(labels))

  # D contributes exactly two rows per contrast, `linear` and `quadratic`, and
  # neither is marked -- while A, on the same data, is marked throughout.
  expect_identical(flags, grepl("^A \\(", labels))
  expect_true(any(grepl("^D \\(diff1; linear\\)$", labels)))
  expect_true(any(grepl("^D \\(diff1; quadratic\\)$", labels)))
  expect_false(any(flags[grepl("^D \\(", labels)]))
})


test_that("a fit that carries no category_support prints exactly as before", {
  expect_warning(
    fit <- adverse_compare_fit(),
    class = "bgms_group_support_warning"
  )

  # A fit saved before the field existed: the arguments simply lack it.
  old = summary(fit)
  old$arguments$category_support = NULL
  expect_null(old$arguments$category_support)

  expect_no_error(out <- capture.output(print(old)))
  expect_no_warning(capture.output(print(old)))

  rows = block_rows(out, "Group differences (main effects):", 6L)
  expect_false(any(grepl("^\\*", rows)))
  expect_false(any(startsWith(out, "* a group lacks")))
  expect_true(all(startsWith(rows, head(old$main_diff$parameter, 6L))))
})
