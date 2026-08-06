test_that("fitted_observed_data returns the fitted data on the input scale", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan5()
  observed = fitted_observed_data(fit)

  complete = Wenchuan[stats::complete.cases(Wenchuan[, 1:5]), 1:5]
  expect_equal(dim(observed), dim(complete))
  expect_equal(colnames(observed), colnames(complete))
  expect_equal(unname(as.matrix(observed)), unname(as.matrix(complete)))

  # It is the default newdata of calibration_check(), so it has to be on the
  # scale predict() expects; the fit itself stores the data zero-based.
  expect_false(identical(unname(get_fit_spec(fit)$data$x), unname(observed)))
  expect_silent(stats::predict(fit,
    newdata = observed, type = "probabilities", method = "posterior-mean"
  ))
})

test_that("calibration_check fits isotonic curves inside a nested-event band", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan5()
  check = calibration_check(fit, nrep = 50, seed = 2)

  expect_s3_class(check, "bgms_calibration")
  expect_equal(nrow(check$summary), 5L)
  expect_equal(sort(unique(check$curves$variable)), sort(colnames(Wenchuan)[1:5]))
  expect_equal(nrow(check$curves), 5L * 101L)

  expect_true(all(check$curves$lower <= check$curves$upper))
  expect_true(all(check$curves$curve >= 0 & check$curves$curve <= 1))
  # Isotonic regression is monotone by construction, per variable.
  for(v in unique(check$curves$variable)) {
    curve = check$curves$curve[check$curves$variable == v]
    expect_false(is.unsorted(curve))
  }

  # Worst variable first, and the summary is a function of the curves.
  expect_false(is.unsorted(rev(check$summary$max_dev)))
  for(v in check$summary$variable) {
    df = check$curves[check$curves$variable == v, ]
    row = check$summary[check$summary$variable == v, ]
    expect_equal(row$max_dev, max(abs(df$curve - df$grid)))
    expect_equal(row$mean_dev, mean(abs(df$curve - df$grid)))
    expect_equal(row$share_outside_band, mean(df$curve < df$lower | df$curve > df$upper))
  }
})

test_that("calibration_check is reproducible and accepts newdata", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan5()
  first = calibration_check(fit, nrep = 40, seed = 9)
  second = calibration_check(fit, nrep = 40, seed = 9)
  expect_equal(first$curves, second$curves)

  # The isotonic fit itself does not depend on the band's resampling.
  other_seed = calibration_check(fit, nrep = 40, seed = 10)
  expect_equal(first$curves$curve, other_seed$curves$curve)
  expect_false(isTRUE(all.equal(first$curves$lower, other_seed$curves$lower)))

  held_out = Wenchuan[stats::complete.cases(Wenchuan[, 1:5]), 1:5][1:100, ]
  out_of_sample = calibration_check(fit, newdata = held_out, nrep = 20, seed = 9)
  expect_equal(nrow(out_of_sample$summary), 5L)

  expect_error(calibration_check(fit, nrep = 0), "nrep")
  expect_error(calibration_check(fit, probs = c(0.5, 0.5)), "two increasing")
})

test_that("uniform_ecdf_band brackets the diagonal and depends on n alone", {
  set.seed(3)
  grid = seq(0, 1, length.out = 21)
  band = uniform_ecdf_band(200, 300, c(0.025, 0.975), grid)

  expect_equal(dim(band), c(2L, length(grid)))
  expect_true(all(band[1, ] <= band[2, ]))
  # A calibrated variable's curve is the diagonal, which the band contains.
  expect_true(all(band[1, ] <= grid + 1e-9 & band[2, ] >= grid - 1e-9))
  # The endpoints are degenerate: every uniform is at or below 1.
  expect_equal(unname(band[, 1]), c(0, 0))
  expect_equal(unname(band[, length(grid)]), c(1, 1))

  # It narrows with n, since the empirical distribution function does.
  set.seed(3)
  wide = uniform_ecdf_band(50, 300, c(0.025, 0.975), grid)
  mid = which.min(abs(grid - 0.5))
  expect_lt(diff(band[, mid]), diff(wide[, mid]))
})

test_that("continuous variables get a PIT panel on the same geometry", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  set.seed(5)
  precision = diag(4)
  precision[1, 2] = precision[2, 1] = -0.3
  y = MASS::mvrnorm(300, rep(0, 4), solve(precision))
  colnames(y) = paste0("v", 1:4)
  ggm = bgm(y,
    variable_type = "continuous", chains = 2, iter = 400, warmup = 400,
    cores = 2, seed = 6, display_progress = "none", verbose = FALSE
  )
  check = calibration_check(ggm, nrep = 40, ndraws = 100, seed = 4)

  expect_s3_class(check, "bgms_calibration")
  expect_equal(nrow(check$summary), 4L)
  expect_true(all(check$summary$kind == "pit"))
  # Same unit square and the same summary columns as the isotonic panel.
  expect_true(all(check$curves$curve >= 0 & check$curves$curve <= 1))
  expect_true(all(check$curves$lower <= check$curves$upper))
  # An empirical distribution function is monotone.
  for(v in unique(check$curves$variable)) {
    expect_false(is.unsorted(check$curves$curve[check$curves$variable == v]))
  }
  # The band does not depend on the variable: under the transform the null is
  # uniform whatever the conditional density was.
  bands = split(check$curves$lower, check$curves$variable)
  expect_equal(length(unique(bands)), 1L)

  # Data whose residual spread is 1.6x what the model predicts leaves the band;
  # correctly scaled held-out data does not. Without that contrast the panel
  # would be decoration.
  set.seed(77)
  n = 300
  root = chol(solve(precision))
  correct = matrix(rnorm(n * 4), n, 4) %*% root
  inflated = correct * 1.6
  dev_of = function(d) {
    max(calibration_check(ggm, newdata = d, nrep = 40, ndraws = 100, seed = 4)$
      summary$share_outside_band)
  }
  expect_lt(dev_of(correct), 0.4)
  expect_gt(dev_of(inflated), 0.4)
})

test_that("a mixed fit produces both panel kinds in one check", {
  skip_on_cran()
  set.seed(9)
  n = 250
  a = sample(0:2, n, TRUE)
  data = data.frame(
    a = a, b = sample(0:2, n, TRUE),
    y1 = rnorm(n) + 0.5 * a, y2 = rnorm(n)
  )
  fit = bgm(data,
    variable_type = c("ordinal", "ordinal", "continuous", "continuous"),
    chains = 2, iter = 400, warmup = 400, cores = 2, seed = 5,
    display_progress = "none", verbose = FALSE
  )
  check = calibration_check(fit, nrep = 30, ndraws = 80, seed = 1)

  expect_equal(nrow(check$summary), 4L)
  expect_equal(
    check$summary$kind[match(c("a", "b", "y1", "y2"), check$summary$variable)],
    c("pav", "pav", "pit", "pit")
  )
  # One figure, one summary table, whatever the mix.
  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(check))

  out = paste(utils::capture.output(print(check)), collapse = "\n")
  expect_match(out, "'pit' panels are continuous variables")
  expect_match(out, "'pav' panels are discrete variables")
})

test_that("plot.bgms_calibration draws small multiples and checks variables", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan5()
  check = calibration_check(fit, nrep = 20, seed = 2)

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(check))
  expect_invisible(plot(check, variables = c("intrusion", "dreams")))
  expect_error(plot(check, variables = "nonesuch"), "not in the calibration check")

  out = paste(utils::capture.output(print(check)), collapse = "\n")
  expect_match(out, "consistency band from 20 resamples")
  expect_match(out, "share_outside_band")
})

test_that("calibration_check handles Blume-Capel and mixed Blume-Capel fits", {
  skip_on_cran()
  x = Wenchuan[stats::complete.cases(Wenchuan[, 1:5]), 1:5]

  # Blume-Capel columns carry an additive shift instead of a recode map, so
  # the decode and the observed-category lookup both have to read the shift.
  bc = bgm(x,
    variable_type = "blume-capel", baseline_category = 1,
    chains = 2, iter = 300, warmup = 300, cores = 2, seed = 3,
    display_progress = "none", verbose = FALSE
  )
  expect_true(all(vapply(
    extract_arguments(bc)$category_levels, is.null, logical(1)
  )))
  observed = fitted_observed_data(bc)
  expect_equal(unname(observed), unname(as.matrix(x)))

  check = calibration_check(bc, nrep = 20, seed = 1)
  expect_s3_class(check, "bgms_calibration")
  expect_equal(nrow(check$summary), 5L)
  expect_true(all(check$summary$kind == "pav"))
  expect_false(anyNA(check$curves$curve))

  mixed = bgm(x,
    variable_type = c(
      "continuous", "continuous", "ordinal", "blume-capel",
      "blume-capel"
    ),
    baseline_category = 1,
    chains = 2, iter = 300, warmup = 300, cores = 2, seed = 3,
    display_progress = "none", verbose = FALSE
  )
  expect_equal(unname(fitted_observed_data(mixed)), unname(as.matrix(x)))
  mixed_check = calibration_check(mixed, nrep = 20, ndraws = 60, seed = 1)
  expect_equal(nrow(mixed_check$summary), 5L)
  expect_equal(
    mixed_check$summary$kind[
      match(colnames(x), mixed_check$summary$variable)
    ],
    c("pit", "pit", "pav", "pav", "pav")
  )
  expect_false(anyNA(mixed_check$curves$curve))
})

test_that("calibration_check reads non-contiguous ordinal category scores", {
  skip_on_cran()
  x = Wenchuan[stats::complete.cases(Wenchuan[, 1:5]), 1:5] * 2L + 1L

  fit = bgm(x,
    variable_type = "ordinal",
    chains = 2, iter = 300, warmup = 300, cores = 2, seed = 3,
    display_progress = "none", verbose = FALSE
  )
  expect_equal(
    extract_arguments(fit)$category_levels[[1]], sort(unique(x[, 1]))
  )
  expect_equal(unname(fitted_observed_data(fit)), unname(as.matrix(x)))

  check = calibration_check(fit, nrep = 20, seed = 1)
  expect_equal(nrow(check$summary), 5L)
  expect_false(anyNA(check$curves$curve))
})

test_that("plot.bgms_calibration pages a wide fit and reports the paging", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan5()
  check = calibration_check(fit, nrep = 20, seed = 2)

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Five variables at three panels a page is two pages, and the message says so.
  withr::local_options(bgms.verbose = TRUE)
  expect_message(plot(check, max_panels = 3L), "Showing page 1 of 2")
  expect_message(plot(check, max_panels = 3L, page = 2L), "Showing page 2 of 2")
  expect_error(plot(check, max_panels = 3L, page = 3L), "make 2 pages")
  # One page, nothing to page through, nothing to say.
  expect_silent(plot(check))
  expect_invisible(plot(check, variables = "intrusion", max_panels = 1L))
})

test_that("the printed calibration summary states the share's units", {
  skip_on_cran()
  fit = get_bgms_fit_wenchuan5()
  check = calibration_check(fit, nrep = 20, seed = 2)

  expect_true(all(check$summary$share_outside_band >= 0 &
    check$summary$share_outside_band <= 1))
  out = paste(utils::capture.output(print(check)), collapse = "\n")
  expect_match(out, "proportion of the curve that does not")
  expect_match(out, "0 to 1 scale")
})


# ------------------------------------------------------------------------------
# bgmCompare fits (F-021)
# ------------------------------------------------------------------------------

test_that("discrete_category_index reads a compare fit's collapsed lookup", {
  # bgmCompare stores a named lookup from the original value to the final
  # 0-based category; a cross-group collapse makes it many-to-one. Reading it
  # with match() against the lookup's values (the bgm layout) would attach
  # observations to the wrong predicted column.
  lookup = c("1" = 0L, "3" = 1L, "7" = 1L, "9" = 2L)
  expect_equal(
    discrete_category_index(c(1, 3, 7, 9, 3), lookup, NA_real_, "v"),
    c(1L, 2L, 2L, 3L, 2L)
  )
  # The bgm layout (sorted original values, unnamed) is unchanged.
  expect_equal(
    discrete_category_index(c(2, 5, 9), c(2, 5, 9), NA_real_, "v"),
    c(1L, 2L, 3L)
  )
})

test_that("fitted_observed_data round-trips a compare fit, non-contiguous scores included", {
  set.seed(8)
  x = matrix(sample(c(2, 5, 9), 150, replace = TRUE), 50, 3)
  g = rep(1:2, length.out = 50) # interleaved on purpose
  fit = bgmCompare(
    x = x, group_indicator = g,
    iter = 25, warmup = 150, chains = 1, seed = 31, display_progress = "none"
  )
  observed = fitted_observed_data(fit)
  # The fit stores its cases sorted by group (stable), on the recoded scale;
  # the decode must return the original scores in that order.
  expect_equal(unname(observed), unname(x[order(g), , drop = FALSE]))
})

test_that("fitted_observed_data puts a compare fit's Blume-Capel baseline back", {
  data("Boredom", package = "bgms")
  rows = c(1:25, 491:515)
  fit = get_bgmcompare_fit_blumecapel()
  observed = fitted_observed_data(fit)
  original = data.matrix(Boredom[rows, 2:5])
  lang = Boredom[rows, "language"]
  expect_equal(
    unname(observed),
    unname(original[order(match(lang, unique(lang))), , drop = FALSE])
  )
})

test_that("calibration_check runs per group on a compare fit", {
  fit = get_bgmcompare_fit_wenchuan5()
  check = calibration_check(fit, nrep = 20, seed = 2)

  expect_s3_class(check, "bgms_calibration")
  expect_equal(check$groups, 1:2)
  expect_true("group" %in% names(check$summary))
  # One curve per variable per group.
  expect_equal(nrow(check$summary), 5L * 2L)
  expect_true(all(table(check$summary$group) == 5L))
  expect_true(all(check$summary$kind == "pav"))
  expect_true(all(check$curves$curve >= 0 & check$curves$curve <= 1))

  # A single group on request.
  one = calibration_check(fit, nrep = 20, seed = 2, group = 2)
  expect_equal(one$groups, 2L)
  expect_equal(nrow(one$summary), 5L)
  expect_error(calibration_check(fit, group = 3), "between 1 and 2")

  out = paste(utils::capture.output(print(check)), collapse = "\n")
  expect_match(out, "group 1, group 2")
  expect_match(out, "own predictions")

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(check))
  withr::local_options(bgms.verbose = TRUE)
  expect_message(plot(check, max_panels = 4L), "Showing page 1 of 3")
  expect_message(plot(check, max_panels = 4L), "panels, worst departure first")
})

test_that("a compare check accepts the fitted rows as newdata, in input order", {
  set.seed(8)
  x = matrix(sample(c(2, 5, 9), 150, replace = TRUE), 50, 3)
  g = rep(1:2, length.out = 50) # interleaved: input order != internal order
  fit = bgmCompare(
    x = x, group_indicator = g,
    iter = 25, warmup = 150, chains = 1, seed = 31, display_progress = "none"
  )
  insample = calibration_check(fit, nrep = 15, seed = 4)
  supplied = calibration_check(fit, nrep = 15, seed = 4, newdata = x)
  expect_equal(supplied$curves, insample$curves)
  expect_equal(supplied$summary, insample$summary)

  expect_error(
    calibration_check(fit, newdata = x[-1, , drop = FALSE]),
    "one row per case"
  )
})

test_that("calibration_check covers a Blume-Capel compare fit", {
  fit = get_bgmcompare_fit_blumecapel()
  check = calibration_check(fit, nrep = 15, seed = 5)
  expect_s3_class(check, "bgms_calibration")
  expect_equal(nrow(check$summary), 4L * 2L)
  expect_true(all(is.finite(check$summary$max_dev)))
})


test_that("an observed value the fit has no category for is named, not absorbed", {
  # Both lookup branches used to fail quietly: the recode map returns NA, which
  # reaches isoreg() as a missing outcome, and the Blume-Capel shift returns an
  # index outside the category range, which the threshold comparison reads as
  # "below all" or "above all".

  # Recode-map branch (regular ordinal, and the compare fit's named lookup).
  expect_error(
    discrete_category_index(c(1, 2, 9), c(1, 2, 3), NA_real_, "worry"),
    "'worry'"
  )
  expect_error(
    discrete_category_index(c(1, 2, 9), c(1, 2, 3), NA_real_, "worry"),
    "no category for"
  )
  expect_error(
    discrete_category_index(c(1, 3, 4), c("1" = 0L, "3" = 1L), NA_real_, "v"),
    "no category for"
  )
  # In support, nothing changes.
  expect_equal(
    discrete_category_index(c(1, 2, 3), c(1, 2, 3), NA_real_, "v", 3L),
    c(1L, 2L, 3L)
  )

  # Blume-Capel branch: the index is the shifted score, and only the category
  # count can say it left the range.
  expect_equal(
    discrete_category_index(c(1, 2, 5), NULL, 1, "v", 5L),
    c(1L, 2L, 5L)
  )
  expect_error(
    discrete_category_index(c(1, 2, 9), NULL, 1, "bc_var", 5L),
    "'bc_var'"
  )
  expect_error(
    discrete_category_index(c(0, 2), NULL, 1, "bc_var", 5L),
    "no category for"
  )
  # The offending value itself is named, not just the variable.
  expect_error(
    discrete_category_index(c(1, 2, 9), NULL, 1, "bc_var", 5L),
    "9"
  )
})


test_that("calibration_check rejects out-of-support newdata in both branches", {
  skip_on_cran()
  x = Wenchuan[stats::complete.cases(Wenchuan[, 1:4]), 1:4]

  ordinal_fit = bgm(x,
    chains = 1, iter = 200, warmup = 200, seed = 5,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  bad = as.matrix(x)
  bad[1, 1] = max(x[, 1]) + 3
  # predict() already warns that the cell is unobserved; the check must go
  # further and refuse, rather than push NA into isoreg().
  expect_error(
    suppressWarnings(
      calibration_check(ordinal_fit, newdata = bad, nrep = 5, seed = 1)
    ),
    "no category for"
  )
  expect_s3_class(
    calibration_check(ordinal_fit, newdata = as.matrix(x), nrep = 5, seed = 1),
    "bgms_calibration"
  )

  bc_fit = bgm(x,
    variable_type = "blume-capel", baseline_category = 1,
    chains = 1, iter = 200, warmup = 200, seed = 5,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  expect_error(
    suppressWarnings(
      calibration_check(bc_fit, newdata = bad, nrep = 5, seed = 1)
    ),
    "no category for"
  )
  expect_s3_class(
    calibration_check(bc_fit, newdata = as.matrix(x), nrep = 5, seed = 1),
    "bgms_calibration"
  )
})
