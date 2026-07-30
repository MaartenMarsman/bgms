calibration_fit = function() {
  bgm(Wenchuan[, 1:5],
    chains = 2, iter = 300, warmup = 300, seed = 7,
    display_progress = "none", verbose = FALSE
  )
}

test_that("fitted_observed_data returns the fitted data on the input scale", {
  skip_on_cran()
  fit = calibration_fit()
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
  fit = calibration_fit()
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
  fit = calibration_fit()
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
    seed = 6, display_progress = "none", verbose = FALSE
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
    chains = 2, iter = 400, warmup = 400, seed = 5,
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
  fit = calibration_fit()
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
