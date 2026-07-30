ppc_fit = function() {
  bgm(Wenchuan[, 1:5],
    chains = 2, iter = 300, warmup = 300, seed = 7,
    display_progress = "none", verbose = FALSE
  )
}

test_that("ppc_observed_data returns the fitted data on simulate()'s scale", {
  skip_on_cran()
  fit = ppc_fit()
  observed = ppc_observed_data(fit)

  complete = Wenchuan[stats::complete.cases(Wenchuan[, 1:5]), 1:5]
  expect_equal(dim(observed), dim(complete))
  expect_equal(colnames(observed), colnames(complete))
  expect_equal(unname(as.matrix(observed)), unname(as.matrix(complete)))

  # Same scale as the replicates: a mismatch here would make every observed
  # statistic incomparable with its predictive distribution.
  replicate = stats::simulate(fit,
    nsim = 20, method = "posterior-sample", ndraws = 1,
    seed = 1, display_progress = "none"
  )[[1]]
  expect_true(all(replicate %in% unique(as.vector(observed))))
})

test_that("ppc statistics line up between observed data and replicates", {
  skip_on_cran()
  fit = ppc_fit()
  observed = ppc_observed_data(fit)
  levels_list = get_fit_spec(fit)$data$category_levels

  margins = ppc_statistic_fun("margins", observed, levels_list)
  expect_equal(nrow(margins$labels), sum(lengths(levels_list)))
  # The proportions of one variable sum to one.
  by_variable = tapply(margins$fun(observed), margins$labels$variable, sum)
  expect_true(all(abs(by_variable - 1) < 1e-12))

  pairwise = ppc_statistic_fun("pairwise", observed, levels_list)
  expect_equal(nrow(pairwise$labels), 10L)
  reference = stats::cor(observed, method = "spearman")
  expect_equal(
    pairwise$fun(observed),
    reference[cbind(
      match(pairwise$labels$variable1, colnames(observed)),
      match(pairwise$labels$variable2, colnames(observed))
    )]
  )

  sumscore = ppc_statistic_fun("sumscore", observed, levels_list)
  values = sumscore$fun(observed)
  expect_equal(sum(values), 1, tolerance = 1e-12)
  expect_equal(sumscore$labels$score[which.max(values)],
    as.numeric(names(sort(table(rowSums(observed)), decreasing = TRUE))[1])
  )
})

test_that("posterior_predictive_check reports coverage per statistic", {
  skip_on_cran()
  fit = ppc_fit()
  check = posterior_predictive_check(fit,
    ndraws = 60, seed = 3, cores = 2, display_progress = "none"
  )

  expect_s3_class(check, "bgms_ppc")
  expect_named(check$statistics, c("margins", "pairwise", "sumscore"))
  expect_equal(check$ndraws, 60L)
  expect_equal(check$nsim, nrow(ppc_observed_data(fit)))

  for(s in names(check$statistics)) {
    df = check$statistics[[s]]
    expect_true(all(c("observed", "predicted", "lower", "upper", "covered") %in% names(df)))
    expect_true(all(df$lower <= df$upper))
    expect_identical(df$covered, df$observed >= df$lower & df$observed <= df$upper)
    expect_equal(check$coverage$covered[check$coverage$statistic == s], sum(df$covered))
    expect_equal(check$coverage$total[check$coverage$statistic == s], nrow(df))
  }

  # Expected coverage is the interval level times the element count, which is
  # the number the printed table compares the realised coverage against.
  expect_equal(check$coverage$expected, 0.95 * check$coverage$total)

  # A subset of statistics is honoured, and only the requested ones are run.
  one = posterior_predictive_check(fit,
    statistic = "pairwise", ndraws = 20, seed = 3, cores = 2,
    display_progress = "none"
  )
  expect_named(one$statistics, "pairwise")
  expect_equal(nrow(one$coverage), 1L)
})

test_that("posterior_predictive_check validates its arguments and model type", {
  skip_on_cran()
  fit = ppc_fit()
  expect_error(
    posterior_predictive_check(fit, statistic = "entropy", display_progress = "none"),
    "should be one of"
  )
  expect_error(
    posterior_predictive_check(fit, probs = 0.95, display_progress = "none"),
    "two increasing probabilities"
  )
  expect_error(
    posterior_predictive_check(fit, probs = c(0.9, 0.1), display_progress = "none"),
    "two increasing probabilities"
  )

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
  expect_error(
    posterior_predictive_check(ggm, display_progress = "none"),
    "discrete variables only"
  )
})

test_that("plot.bgms_ppc draws one panel per statistic", {
  skip_on_cran()
  fit = ppc_fit()
  check = posterior_predictive_check(fit,
    ndraws = 20, seed = 3, cores = 2, display_progress = "none"
  )

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(check))
  expect_invisible(plot(check, max_labels = 0L))

  out = paste(utils::capture.output(print(check)), collapse = "\n")
  expect_match(out, "20 replicated datasets")
  expect_match(out, "margins")
  expect_match(out, "sumscore")
})

test_that("calibration_check fits isotonic curves inside a nested-event band", {
  skip_on_cran()
  fit = ppc_fit()
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
  fit = ppc_fit()
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

test_that("plot.bgms_calibration draws small multiples and checks variables", {
  skip_on_cran()
  fit = ppc_fit()
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
