test_that("two_state_se_logit reproduces the Jeffreys-smoothed reference", {
  # Independent reference: the pooled two-state summary of the calibration
  # study, transitions counted within chains, Jeffreys-smoothed rates.
  reference = function(draws) {
    n_iter = nrow(draws)
    n01 = n10 = n0 = n1 = 0
    for(ch in seq_len(ncol(draws))) {
      x = draws[, ch]
      from = x[-n_iter]
      to = x[-1L]
      n0 = n0 + sum(from == 0)
      n1 = n1 + sum(from == 1)
      n01 = n01 + sum(from == 0 & to == 1)
      n10 = n10 + sum(from == 1 & to == 0)
    }
    a = (n01 + 0.5) / (n0 + 1)
    b = (n10 + 0.5) / (n1 + 1)
    p = a / (a + b)
    ess = n_iter * ncol(draws) * (a + b) / (2 - a - b)
    sqrt(1 / (ess * p * (1 - p)))
  }

  set.seed(11)
  n_iter = 500
  mixing = cbind(rbinom(n_iter, 1, 0.4), rbinom(n_iter, 1, 0.4))
  sticky = cbind(c(rep(0, 300), rep(1, 200)), c(rep(1, 250), rep(0, 250)))
  stuck = cbind(rep(1, n_iter), rep(1, n_iter))

  chains = list(
    cbind(mixing[, 1], sticky[, 1], stuck[, 1]),
    cbind(mixing[, 2], sticky[, 2], stuck[, 2])
  )
  se = two_state_se_logit(chains)

  for(k in 1:3) {
    expect_equal(se[k], reference(cbind(chains[[1]][, k], chains[[2]][, k])),
      tolerance = 1e-12, ignore_attr = TRUE
    )
  }

  # The smoothing is what makes the standard error exist at zero flips, and it
  # saturates at sqrt(3) there rather than diverging or vanishing.
  expect_true(is.finite(se[3]))
  expect_equal(unname(se[3]), sqrt(3), tolerance = 1e-3)
})

test_that("boundary_distance measures both boundaries in standard errors", {
  lthr = log(10)
  # Evidence sitting exactly on the presence boundary is zero standard errors
  # away from it, whatever the standard error.
  expect_equal(boundary_distance(lthr, 0.5, lthr), 0)
  expect_equal(boundary_distance(-lthr, 0.5, lthr), 0)
  # The Bayes factor and the standard errors share the natural log scale, so a
  # gap of two standard errors reads as two.
  expect_equal(boundary_distance(lthr + 1, 0.5, lthr), 2)
  # A saturated Bayes factor is infinitely far from either boundary.
  expect_equal(boundary_distance(Inf, 0.5, lthr), Inf)
  # No standard error, no distance.
  expect_true(is.na(boundary_distance(0.5, NA_real_, lthr)))
})

test_that("the reported evidence is the natural log Bayes factor", {
  # The classification boundaries sit at +/- log(t), not at +/- log10(t), and
  # the two evidence columns are exp/log inverses of each other.
  parameter = c("a-b", "a-c", "b-c", "a-d")
  log_bf = c(log(50), log(0.01), 0, -1.474)
  set.seed(4)
  draws = list(matrix(rbinom(4L * 60L, 1L, 0.5), nrow = 60L))

  v = build_verdicts(
    parameter = parameter, log_bf = log_bf, pip = c(0.9, 0.02, 0.5, 0.2),
    mcse = rep(0.01, 4L), draws = draws, evidence_threshold = 10,
    flag_validated = TRUE
  )

  expect_true("log_bf" %in% names(v))
  expect_false("log10_bf" %in% names(v))
  expect_equal(v$log_bf, log_bf)
  expect_equal(v$bf, exp(log_bf))
  # -1.474 is below the log10 boundary but inside the natural-log one.
  expect_equal(
    as.character(v$verdict),
    c("presence", "absence", "undecided", "undecided")
  )

  # The printed header names the boundaries in the unit the table reports.
  out = paste(utils::capture.output(print(v)), collapse = "\n")
  expect_match(out, "presence: log BF > 2.30; absence: log BF < -2.30")
  expect_match(out, "log_bf")
  expect_false(grepl("log10", out, fixed = TRUE))
})

test_that("format_log_bf caps the magnitude it prints", {
  expect_equal(format_log_bf(301.4), "= 301.4")
  expect_equal(format_log_bf(-2.3456), "= -2.3")
  expect_equal(format_log_bf(1e5), "> 10,000")
  expect_equal(format_log_bf(-1e5), "< -10,000")
  expect_equal(format_log_bf(Inf), "> 10,000")
  expect_equal(format_log_bf(-Inf), "< -10,000")
  expect_equal(format_log_bf(NA_real_), "NA")

  # A value that rounds away carries no sign the run established, so it must
  # not print as "-0.0" (F-105). print.bgms_verdicts reaches this directly.
  expect_equal(format_log_bf(-0.004), "= 0.0")
  expect_equal(format_log_bf(-0.04), "= 0.0")
  expect_equal(format_log_bf(0.004), "= 0.0")
  expect_equal(format_log_bf(-0.06), "= -0.1") # still rounds to a real value
})

test_that("compare_indicator_index lays out main then pairwise per variable", {
  idx = compare_indicator_index(3L)
  expect_equal(nrow(idx), 6L)
  expect_equal(idx[, 1], c(1, 1, 1, 2, 2, 3))
  expect_equal(idx[, 2], c(1, 2, 3, 2, 3, 3))
})

test_that("verdicts() reports one row per edge with a three-way verdict", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()
  v = verdicts(fit)

  expect_s3_class(v, "bgms_verdicts")
  expect_equal(nrow(v), 15L)
  expect_equal(v$parameter, rownames(fit$posterior_summary_indicator))
  expect_equal(levels(v$verdict), c("presence", "undecided", "absence"))
  expect_false(anyNA(v$verdict))
  expect_identical(attr(v, "evidence_threshold"), 10)
  expect_true(attr(v, "flag_validated"))

  # The verdict is the reading of the Bayes factor at the threshold, and the
  # Bayes factor is the one extract_inclusion_bf() reports.
  bf = extract_inclusion_bf(fit)
  edges = strsplit(v$parameter, "-", fixed = TRUE)
  reported = vapply(edges, function(e) bf[e[1], e[2]], numeric(1))
  expect_equal(v$bf, reported, tolerance = 1e-8)
  expect_true(all(v$verdict[reported > 10] == "presence"))
  expect_true(all(v$verdict[reported < 0.1] == "absence"))
  expect_true(all(v$verdict[reported >= 0.1 & reported <= 10] == "undecided"))

  # The inclusion probability is the Rao-Blackwellized one the fit reports.
  expect_equal(v$pip, fit$posterior_summary_indicator$mean, tolerance = 1e-10)
})

test_that("the fragility flag is the union of the two standard errors", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()
  v = verdicts(fit)

  expected = (!is.na(v$distance_two_state) & v$distance_two_state < 2) |
    (!is.na(v$distance_rb) & v$distance_rb < 2)
  expect_identical(v$fragile, expected)

  # An edge either standard error calls close is flagged; agreement is not
  # required. Both directions of the union are exercised by construction.
  expect_true(all(v$fragile[!is.na(v$distance_rb) & v$distance_rb < 2]))
  expect_true(all(v$fragile[!is.na(v$distance_two_state) & v$distance_two_state < 2]))

  # The two-state standard error is always available; it is the one that keeps
  # zero-flip edges from abstaining when the RB draws are constant.
  expect_false(anyNA(v$se_two_state))
  expect_true(any(is.na(v$se_rb)))
  expect_false(anyNA(v$fragile))
})

test_that("verdicts() moves the boundaries with evidence_threshold", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()
  loose = verdicts(fit, evidence_threshold = 3)
  strict = verdicts(fit, evidence_threshold = 100)

  # A stricter threshold can only move verdicts toward undecided.
  expect_gte(sum(strict$verdict == "undecided"), sum(loose$verdict == "undecided"))
  expect_identical(attr(strict, "evidence_threshold"), 100)

  # The standard errors are properties of the chain, not of the threshold.
  expect_equal(loose$se_two_state, strict$se_two_state)
  expect_equal(loose$se_rb, strict$se_rb)

  expect_error(verdicts(fit, evidence_threshold = 1), "greater than 1")
  expect_error(verdicts(fit, evidence_threshold = c(10, 30)), "single number")
})

test_that("verdicts() errors without selection and covers bgmCompare", {
  # The no-selection guard reads the spec, before any draw.
  skip_on_cran()
  no_selection = bgm(Wenchuan[, 1:4],
    chains = 1, iter = 20, warmup = 20, seed = 2,
    edge_selection = FALSE, display_progress = "none", verbose = FALSE
  )
  expect_error(verdicts(no_selection), "edge selection")

  x = Wenchuan[1:80, 1:4]
  fit = bgmCompare(
    x = x, group_indicator = rep(1:2, each = 40),
    iter = 300, warmup = 300, chains = 2, cores = 2, seed = 13,
    difference_selection = TRUE, display_progress = "none"
  )
  v = verdicts(fit)
  expect_s3_class(v, "bgms_verdicts")
  # The compare producer calls this element `indicators`. verdicts.bgmCompare
  # used to reach it as `$indicator`, which only resolved through `$` partial
  # matching -- correct today, NULL the day a sibling name shares the prefix.
  # The labels are unchanged by spelling it exactly.
  names_compare = get_raw_samples(fit)$parameter_names
  expect_true("indicators" %in% names(names_compare))
  expect_false("indicator" %in% names(names_compare))
  expect_equal(v$parameter, names_compare[["indicators"]])
  expect_equal(v$parameter, names_compare$indicator) # the pre-fix partial match
  expect_equal(nrow(v), 10L)

  # The pairwise differences were selected and carry verdicts; the main-effect
  # differences are not under selection (main_difference_selection defaults to
  # FALSE), so the table keeps their rows as NA while the print leaves them
  # out of the counts and says why.
  is_main = grepl("(main)", v$parameter, fixed = TRUE)
  expect_true(all(is.na(v$verdict[is_main])))
  expect_false(anyNA(v$verdict[!is_main]))
  expect_false(any(v$fragile[is_main]))
  expect_equal(attr(v, "unselected_main"), is_main)

  # The fragility flag's operating point was measured on single-network edge
  # indicators only; the difference-indicator print must say so rather than
  # borrow the single-network numbers.
  expect_false(attr(v, "flag_validated"))
  out = paste(utils::capture.output(print(v)), collapse = "\n")
  expect_match(out, "not under selection")
  expect_match(out, "6 indicators")
  expect_false(grepl("never updated and carry no verdict", out))
  expect_match(out, "not validated for difference indicators")
  expect_match(out, "scale-contingent")
})

test_that("print.bgms_verdicts tallies verdicts and warns once when fragile", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()
  v = verdicts(fit)

  out = paste(utils::capture.output(print(v)), collapse = "\n")
  expect_match(out, "presence \\d+ \\| undecided \\d+ \\| absence \\d+")
  expect_match(out, "15 indicators")
  # Edge indicators are what the flag was calibrated on, so no caveat here.
  expect_false(grepl("not validated", out))
  if(any(v$fragile)) {
    expect_match(out, "Monte-Carlo fragile")
    expect_match(out, "Consider a")
  }

  # No fragile edges, no advice line.
  v_none = v
  v_none$fragile = rep(FALSE, nrow(v))
  quiet = paste(utils::capture.output(print(v_none)), collapse = "\n")
  expect_false(grepl("Monte-Carlo fragile", quiet))

  # Selecting columns keeps the class but not the table; printing what is left
  # must degrade to the plain data frame rather than fail on a missing column.
  subset_columns = v[, c("parameter", "log_bf", "verdict")]
  expect_s3_class(subset_columns, "bgms_verdicts")
  expect_silent(plain <- utils::capture.output(print(subset_columns)))
  expect_false(any(grepl("Edge verdicts at", plain)))
  expect_true(any(grepl("intrusion", plain)))

  # Selecting rows keeps the full table, and the header with it.
  subset_rows = v[v$verdict == "absence", ]
  expect_match(
    paste(utils::capture.output(print(subset_rows)), collapse = "\n"),
    "Edge verdicts at"
  )

  # subset() supplies `j`, so it keeps every display column but drops the
  # evidence_threshold attribute -- the one combination a column-only guard
  # misses, and it used to fail in log(threshold) (F-091).
  expect_silent(dropped <- utils::capture.output(print(subset(v, fragile))))
  expect_false(any(grepl("Edge verdicts at", dropped)))
})

test_that("the print's tally pluralizes and states the absence threshold cleanly", {
  # Hand-built so the counts and the threshold are exact on every platform.
  one = structure(
    data.frame(
      parameter = "a-b",
      pip = 0.99,
      log_bf = 4.6,
      verdict = factor("presence", levels = c("presence", "undecided", "absence")),
      fragile = FALSE,
      stringsAsFactors = FALSE
    ),
    class = c("bgms_verdicts", "data.frame"),
    evidence_threshold = 30
  )

  out = paste(utils::capture.output(print(one)), collapse = "\n")
  # One row is one indicator, and 1/30 is stated at a readable precision
  # rather than as %g's 0.0333333.
  expect_match(out, "\\(1 indicator\\)")
  expect_false(grepl("1 indicators", out, fixed = TRUE))
  expect_match(out, "Bayes factor of 30 \\(and 0\\.0333 for absence\\)")

  two = one[c(1, 1), ]
  attr(two, "evidence_threshold") = 10
  out_two = paste(utils::capture.output(print(two)), collapse = "\n")
  expect_match(out_two, "\\(2 indicators\\)")
  expect_match(out_two, "\\(and 0\\.1 for absence\\)")
})

test_that("indicator_pair_index follows the fit's own indicator layout", {
  skip_on_cran()
  # Interleaved types, so the block layout and the row-major upper triangle are
  # genuinely different permutations.
  set.seed(21)
  n = 120
  x = cbind(
    sample(0:2, n, TRUE), rnorm(n), sample(0:2, n, TRUE), rnorm(n), rnorm(n)
  )
  colnames(x) = c("d1", "c1", "d2", "c2", "c3")
  fit = bgm(x,
    variable_type = c(
      "ordinal", "continuous", "ordinal", "continuous",
      "continuous"
    ),
    chains = 2, iter = 300, warmup = 300, cores = 2, seed = 6,
    display_progress = "none", verbose = FALSE
  )

  names_out = colnames(x)
  idx = indicator_pair_index(fit, 5L)
  labels = paste(names_out[idx[, 1]], names_out[idx[, 2]], sep = "-")
  flipped = paste(names_out[idx[, 2]], names_out[idx[, 1]], sep = "-")
  raw_names = get_raw_samples(fit)$parameter_names$indicator
  expect_equal(nrow(idx), 10L)
  expect_true(all(labels == raw_names | flipped == raw_names))

  # Every reported number belongs to the edge the row names, which the
  # row-major upper triangle got wrong for a mixed fit.
  v = verdicts(fit)
  expect_equal(v$parameter, raw_names)
  bf = extract_inclusion_bf(fit, log = TRUE)
  pip = extract_posterior_inclusion_probabilities(fit)
  ends = strsplit(v$parameter, "-", fixed = TRUE)
  expect_equal(v$log_bf, vapply(ends, function(e) bf[e[1], e[2]], numeric(1)))
  expect_equal(v$pip, vapply(ends, function(e) pip[e[1], e[2]], numeric(1)))
})

test_that("a single-type fit keeps the row-major upper-triangle layout", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan6()
  idx = indicator_pair_index(fit, 6L)
  expected = which(upper.tri(matrix(0, 6L, 6L)), arr.ind = TRUE)
  expected = expected[order(expected[, "row"], expected[, "col"]), , drop = FALSE]
  expect_equal(unname(idx), unname(expected))
})


# ------------------------------------------------------------------------------
# Printed report for bgmCompare difference verdicts (F-060)
# ------------------------------------------------------------------------------

# A hand-built verdicts table prints deterministically on every platform,
# which is what a snapshot of the layout needs; fitted numbers would not.
make_compare_verdicts = function(main_selected) {
  parameter = c(
    "A (main)", "A-B (pairwise)", "A-C (pairwise)",
    "B (main)", "B-C (pairwise)", "C (main)"
  )
  is_main = grepl("(main)", parameter, fixed = TRUE)
  pip = c(NaN, 0.95, 0.10, NaN, 0.52, NaN)
  log_bf = c(NA, 2.94, -2.20, NA, 0.08, NA)
  if(main_selected) {
    pip[is_main] = c(0.05, 0.90, 0.50)
    log_bf[is_main] = c(-2.94, 2.20, 0)
  }
  out = data.frame(
    parameter = parameter,
    pip = pip,
    bf = exp(log_bf),
    log_bf = log_bf,
    verdict = factor(
      verdict_from_lbf(log_bf, log(10)),
      levels = c("presence", "undecided", "absence")
    ),
    se_two_state = 0.1, se_rb = 0.1,
    distance_two_state = 5, distance_rb = 5,
    fragile = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  structure(
    out,
    class = c("bgms_verdicts", "data.frame"),
    evidence_threshold = 10,
    flag_validated = FALSE,
    difference_fit = TRUE,
    unselected_main = if(main_selected) NULL else is_main
  )
}

test_that("the compare print leaves unselected main differences out of the table", {
  expect_snapshot(print(make_compare_verdicts(main_selected = FALSE)))
})

test_that("the compare print tabulates selected main differences as rows", {
  expect_snapshot(print(make_compare_verdicts(main_selected = TRUE)))
})

test_that("the fragility footer suggests a longer run without commanding one", {
  out = capture.output(print(make_compare_verdicts(main_selected = FALSE)))
  expect_true(any(grepl("Consider a", out)))
  expect_false(any(grepl("Run longer", out)))
})
