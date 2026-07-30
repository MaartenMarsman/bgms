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
  expect_true(all(check$curves$pav >= 0 & check$curves$pav <= 1))
  # Isotonic regression is monotone by construction, per variable.
  for(v in unique(check$curves$variable)) {
    curve = check$curves$pav[check$curves$variable == v]
    expect_false(is.unsorted(curve))
  }

  # Worst variable first, and the summary is a function of the curves.
  expect_false(is.unsorted(rev(check$summary$max_dev)))
  for(v in check$summary$variable) {
    df = check$curves[check$curves$variable == v, ]
    row = check$summary[check$summary$variable == v, ]
    expect_equal(row$max_dev, max(abs(df$pav - df$grid)))
    expect_equal(row$mean_dev, mean(abs(df$pav - df$grid)))
    expect_equal(row$share_outside_band, mean(df$pav < df$lower | df$pav > df$upper))
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
  expect_equal(first$curves$pav, other_seed$curves$pav)
  expect_false(isTRUE(all.equal(first$curves$lower, other_seed$curves$lower)))

  held_out = Wenchuan[stats::complete.cases(Wenchuan[, 1:5]), 1:5][1:100, ]
  out_of_sample = calibration_check(fit, newdata = held_out, nrep = 20, seed = 9)
  expect_equal(nrow(out_of_sample$summary), 5L)

  expect_error(calibration_check(fit, nrep = 0), "nrep")
  expect_error(calibration_check(fit, probs = c(0.5, 0.5)), "two increasing")
})

test_that("calibration_check errors on continuous variables and says why", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  set.seed(5)
  precision = diag(4)
  precision[1, 2] = precision[2, 1] = -0.3
  y = MASS::mvrnorm(120, rep(0, 4), solve(precision))
  colnames(y) = paste0("v", 1:4)
  ggm = bgm(y,
    variable_type = "continuous", chains = 2, iter = 200, warmup = 200,
    seed = 6, display_progress = "none", verbose = FALSE
  )

  # The isotonic curve reads an observed category against its predicted
  # probability; a continuous prediction is a density and has neither.
  expect_error(calibration_check(ggm), "densities rather than")
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
