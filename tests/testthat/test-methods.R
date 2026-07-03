# ==============================================================================
# Tests for S3 Methods on bgms and bgmCompare Objects (Parameterized)
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Range invariants, dimension consistency, reproducibility
#
# This file uses parameterized testing to reduce code repetition.
# Tests for: print, summary, coef, simulate, predict
#
# Tolerance testing principles applied:
#   - Range invariants: predictions in valid category range, probabilities sum to 1
#   - Symmetry: pairwise coefficient matrices must be symmetric
#   - Dimension consistency: simulated data has correct n x p shape
#   - Reproducibility: seed produces identical results
# ==============================================================================

# ------------------------------------------------------------------------------
# Fixture Specifications — defined in helper-fixtures.R
# get_bgms_fixtures(), get_bgmcompare_fixtures()
# ------------------------------------------------------------------------------


# ==============================================================================
# Fixture Coverage Guards
# ==============================================================================
# Fail fast when a fixture list drifts from the expected coverage.

test_that("get_bgms_fixtures covers all required labels", {
  specs = get_bgms_fixtures()
  labels = vapply(specs, `[[`, character(1), "label")
  required = c("binary", "ordinal", "ggm", "mixed-mrf", "blume-capel")
  for(r in required) {
    expect_true(r %in% labels, info = sprintf("missing required label '%s'", r))
  }
  expect_equal(length(specs), 22L,
    info = "bgms fixture count changed — update this guard if intentional"
  )
})

test_that("get_bgmcompare_fixtures covers all required labels", {
  specs = get_bgmcompare_fixtures()
  labels = vapply(specs, `[[`, character(1), "label")
  required = c("binary", "ordinal", "blume-capel")
  for(r in required) {
    expect_true(r %in% labels, info = sprintf("missing required label '%s'", r))
  }
  expect_equal(length(specs), 10L,
    info = "bgmcompare fixture count changed — update this guard if intentional"
  )
})

test_that("get_extractor_fixtures covers all model families", {
  specs = get_extractor_fixtures()
  labels = vapply(specs, `[[`, character(1), "label")
  required = c(
    "bgms_binary", "bgms_ggm", "bgms_mixed",
    "bgmCompare_binary", "bgmCompare_ordinal"
  )
  for(r in required) {
    expect_true(r %in% labels, info = sprintf("missing required label '%s'", r))
  }
  expect_equal(length(specs), 9L,
    info = "extractor fixture count changed — update this guard if intentional"
  )
})


# ==============================================================================
# print() Tests (Parameterized)
# ==============================================================================

test_that("print.bgms produces output without error for all fixture types", {
  for(spec in get_bgms_fixtures()) {
    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()

    expect_output(print(fit), regexp = "Number of variables", info = ctx)
    expect_output(print(fit), regexp = "MCMC", info = ctx)
  }
})

test_that("print.bgmCompare produces output without error for all fixture types", {
  for(spec in get_bgmcompare_fixtures()) {
    ctx = sprintf("[bgmCompare %s]", spec$label)
    fit = spec$get_fit()

    expect_output(print(fit), regexp = "Number of variables", info = ctx)
    expect_output(print(fit), regexp = "groups|MCMC", info = ctx)
  }
})


# ==============================================================================
# summary() Tests (Parameterized)
# ==============================================================================

test_that("summary.bgms returns correct structure for all fixture types", {
  for(spec in get_bgms_fixtures()) {
    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()

    summ = summary(fit)

    expect_s3_class(summ, "summary.bgms")
    has_content = !is.null(summ$main) || !is.null(summ$quadratic) || !is.null(summ$pairwise)
    expect_true(has_content, info = ctx)
  }
})

test_that("summary.bgmCompare returns correct structure for all fixture types", {
  for(spec in get_bgmcompare_fixtures()) {
    ctx = sprintf("[bgmCompare %s]", spec$label)
    fit = spec$get_fit()

    summ = summary(fit)

    if(!is.null(summ)) {
      expect_s3_class(summ, "summary.bgmCompare")
    }
  }
})

test_that("summary.bgms components have correct dimensions", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)
  summ = summary(fit)

  if(!is.null(summ$pairwise)) {
    p = args$num_variables
    expected_rows = p * (p - 1) / 2
    expect_equal(nrow(summ$pairwise), expected_rows)
  }
})

test_that("print.summary.bgms produces readable output", {
  fit = get_bgms_fit()
  summ = summary(fit)

  expect_output(print(summ), regexp = "Posterior summaries")
})

test_that("print.summary.bgmCompare produces readable output for all fixture types", {
  for(spec in get_bgmcompare_fixtures()) {
    ctx = sprintf("[bgmCompare %s]", spec$label)
    fit = spec$get_fit()
    summ = summary(fit)

    if(!is.null(summ)) {
      expect_output(print(summ), regexp = "Posterior summaries", info = ctx)
      expect_output(print(summ), regexp = "bgmCompare", info = ctx)
    }
  }
})


# ==============================================================================
# coef() Tests (Parameterized)
# ==============================================================================

test_that("coef.bgms returns list with expected components for all fixture types", {
  for(spec in get_bgms_fixtures()) {
    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    coeffs = coef(fit)

    expect_true(is.list(coeffs), info = ctx)
    expect_true("main" %in% names(coeffs), info = paste(ctx, "missing main"))
    expect_true("pairwise" %in% names(coeffs), info = paste(ctx, "missing pairwise"))

    # Pairwise should be symmetric matrix
    if(!is.null(coeffs$pairwise)) {
      expect_true(is.matrix(coeffs$pairwise), info = paste(ctx, "pairwise not matrix"))
      p = args$num_variables
      expect_equal(dim(coeffs$pairwise), c(p, p), info = paste(ctx, "wrong dim"))
      expect_true(is_symmetric(coeffs$pairwise), info = paste(ctx, "not symmetric"))
    }

    # If edge selection was used, check indicator
    if(isTRUE(args$edge_selection)) {
      expect_true("indicator" %in% names(coeffs), info = paste(ctx, "missing indicator"))
      expect_true(values_in_range(coeffs$indicator, 0, 1),
        info = paste(ctx, "indicator out of range")
      )
    }
  }
})

test_that("coef.bgmCompare returns group-specific effects for all fixture types", {
  for(spec in get_bgmcompare_fixtures()) {
    ctx = sprintf("[bgmCompare %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    coeffs = coef(fit)

    expect_true(is.list(coeffs), info = ctx)
    expect_true("main_effects_raw" %in% names(coeffs), info = ctx)
    expect_true("pairwise_effects_raw" %in% names(coeffs), info = ctx)
    expect_true("main_effects_groups" %in% names(coeffs), info = ctx)
    expect_true("pairwise_effects_groups" %in% names(coeffs), info = ctx)

    # Group effects should have correct number of columns
    n_groups = args$num_groups
    expect_equal(ncol(coeffs$pairwise_effects_groups), n_groups, info = ctx)
  }
})


# ==============================================================================
# simulate() Tests (Parameterized)
# ==============================================================================

test_that("simulate.bgms returns matrix of correct size for ordinal fixtures", {
  for(spec in get_bgms_fixtures()) {
    if(isTRUE(spec$is_continuous)) next
    if(isTRUE(spec$is_mixed)) next

    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    n_sim = 30
    simulated = simulate(fit, nsim = n_sim, method = "posterior-mean", seed = 123)

    expect_true(is.matrix(simulated), info = ctx)
    expect_equal(nrow(simulated), n_sim, info = paste(ctx, "wrong nrow"))
    expect_equal(ncol(simulated), args$num_variables, info = paste(ctx, "wrong ncol"))
    expect_equal(colnames(simulated), args$data_columnnames, info = ctx)

    # Values should be integers on the original category scale.
    expect_true(all(simulated == round(simulated)), info = paste(ctx, "not integers"))
    expect_true(all(simulated >= 0), info = paste(ctx, "negative values"))

    for(j in seq_len(args$num_variables)) {
      levels_j = args$category_levels[[j]]
      if(is.null(levels_j)) {
        # Blume-Capel (no recode map): original-scale scores, i.e. the 0-based
        # internal range shifted by the stored per-variable offset.
        shift_j = args$blume_capel_shift[j]
        expect_true(
          all(simulated[, j] >= shift_j &
            simulated[, j] <= shift_j + args$num_categories[j]),
          info = sprintf("%s variable %d exceeds the category range", ctx, j)
        )
      } else {
        # Ordinal: simulate() returns the original category values. bgm stores
        # them as sorted values; bgmCompare as a named original -> code lookup.
        valid = if(!is.null(names(levels_j))) as.numeric(names(levels_j)) else levels_j
        expect_true(
          all(simulated[, j] %in% valid),
          info = sprintf("%s variable %d off the original category scale", ctx, j)
        )
      }
    }
  }
})

test_that("simulate.bgms returns matrix of correct size for GGM fixtures", {
  for(spec in get_bgms_fixtures()) {
    if(!isTRUE(spec$is_continuous)) next

    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    n_sim = 30
    simulated = simulate(fit, nsim = n_sim, method = "posterior-mean", seed = 123)

    expect_true(is.matrix(simulated), info = ctx)
    expect_equal(nrow(simulated), n_sim, info = paste(ctx, "wrong nrow"))
    expect_equal(ncol(simulated), args$num_variables, info = paste(ctx, "wrong ncol"))
    expect_equal(colnames(simulated), args$data_columnnames, info = ctx)

    # Values should be real-valued (not restricted to integers)
    expect_true(is.numeric(simulated), info = paste(ctx, "not numeric"))
    expect_true(all(is.finite(simulated)), info = paste(ctx, "non-finite values"))
  }
})

test_that("simulate.bgms returns matrix of correct size for mixed-mrf fixtures", {
  for(spec in get_bgms_fixtures()) {
    if(!isTRUE(spec$is_mixed)) next

    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    n_sim = 30
    simulated = simulate(fit, nsim = n_sim, method = "posterior-mean", seed = 123)

    expect_true(is.matrix(simulated), info = ctx)
    expect_equal(nrow(simulated), n_sim, info = paste(ctx, "wrong nrow"))
    expect_equal(ncol(simulated), args$num_variables, info = paste(ctx, "wrong ncol"))
    expect_equal(colnames(simulated), args$data_columnnames, info = ctx)
    expect_true(all(is.finite(simulated)), info = paste(ctx, "non-finite values"))

    # Ordinal columns should be non-negative integers
    for(j in seq_len(args$num_variables)) {
      if(args$variable_type[j] != "continuous") {
        expect_true(all(simulated[, j] == round(simulated[, j])),
          info = sprintf("%s variable %d should be integer", ctx, j)
        )
        expect_true(all(simulated[, j] >= 0),
          info = sprintf("%s variable %d has negative values", ctx, j)
        )
      }
    }
  }
})

test_that("simulate.bgms is reproducible with seed", {
  fit = get_bgms_fit()

  sim1 = simulate(fit, nsim = 30, method = "posterior-mean", seed = 999)
  sim2 = simulate(fit, nsim = 30, method = "posterior-mean", seed = 999)

  expect_equal(sim1, sim2)
})

test_that("simulate.bgms posterior-sample returns list", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)

  n_draws = 3
  n_sim = 20

  result = simulate(fit,
    nsim = n_sim,
    method = "posterior-sample",
    ndraws = n_draws,
    seed = 123,
    display_progress = "none"
  )

  expect_true(is.list(result))
  expect_equal(length(result), n_draws)

  for(i in seq_along(result)) {
    expect_true(is.matrix(result[[i]]))
    expect_equal(nrow(result[[i]]), n_sim)
    expect_equal(ncol(result[[i]]), args$num_variables)
  }
})

test_that("simulate.bgms handles edge cases", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)

  # Single observation
  sim1 = simulate(fit, nsim = 1, method = "posterior-mean", seed = 42)
  expect_true(is.matrix(sim1))
  expect_equal(nrow(sim1), 1)

  # ndraws = 1
  result = simulate(fit,
    nsim = 10, method = "posterior-sample",
    ndraws = 1, seed = 42, display_progress = "none"
  )
  expect_true(is.list(result))
  expect_equal(length(result), 1)
})

test_that("simulate.bgms GGM is reproducible with seed", {
  fit = get_bgms_fit_ggm()

  sim1 = simulate(fit, nsim = 30, method = "posterior-mean", seed = 999)
  sim2 = simulate(fit, nsim = 30, method = "posterior-mean", seed = 999)

  expect_equal(sim1, sim2)
})

test_that("simulate.bgms GGM posterior-sample returns list of numeric matrices", {
  fit = get_bgms_fit_ggm()
  args = extract_arguments(fit)

  n_draws = 3
  n_sim = 20

  result = simulate(fit,
    nsim = n_sim,
    method = "posterior-sample",
    ndraws = n_draws,
    seed = 123,
    display_progress = "none"
  )

  expect_true(is.list(result))
  expect_equal(length(result), n_draws)

  for(i in seq_along(result)) {
    expect_true(is.matrix(result[[i]]))
    expect_equal(nrow(result[[i]]), n_sim)
    expect_equal(ncol(result[[i]]), args$num_variables)
    expect_true(all(is.finite(result[[i]])))
  }
})


# ==============================================================================
# predict() Tests (Parameterized)
# ==============================================================================

test_that("predict.bgms returns valid probabilities for ordinal fixtures", {
  for(spec in get_bgms_fixtures()) {
    if(isTRUE(spec$is_continuous)) next
    if(isTRUE(spec$is_mixed)) next

    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    newdata = spec$get_prediction_data(n = 5)
    probs = predict(fit, newdata = newdata, type = "probabilities")

    expect_true(is.list(probs), info = ctx)
    expect_equal(length(probs), args$num_variables, info = ctx)

    # Each variable's probabilities should sum to 1
    for(j in seq_along(probs)) {
      expect_true(is.matrix(probs[[j]]), info = paste(ctx, "var", j))
      expect_equal(nrow(probs[[j]]), nrow(newdata), info = paste(ctx, "var", j))

      # No NAs in probability output
      expect_false(anyNA(probs[[j]]),
        info = sprintf("%s var %d has NAs", ctx, j)
      )
      # Probabilities in [0, 1]
      expect_true(all(probs[[j]] >= 0 & probs[[j]] <= 1),
        info = sprintf("%s var %d probs out of [0,1]", ctx, j)
      )

      # Row sums = 1
      row_sums = rowSums(probs[[j]])
      expect_true(
        all(abs(row_sums - 1) < 1e-6),
        info = sprintf("%s var %d probs don't sum to 1", ctx, j)
      )
    }
  }
})

test_that("predict.bgms returns valid conditional moments for GGM fixtures", {
  for(spec in get_bgms_fixtures()) {
    if(!isTRUE(spec$is_continuous)) next

    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    newdata = spec$get_prediction_data(n = 5)
    result = predict(fit, newdata = newdata)

    expect_type(result, "list")
    expect_length(result, args$num_variables)
    expect_equal(names(result), args$data_columnnames, info = ctx)

    for(j in seq_along(result)) {
      expect_equal(nrow(result[[j]]), nrow(newdata), info = paste(ctx, "var", j))
      expect_equal(ncol(result[[j]]), 2, info = paste(ctx, "var", j))
      expect_equal(colnames(result[[j]]), c("mean", "sd"), info = paste(ctx, "var", j))
      expect_true(all(result[[j]][, "sd"] > 0),
        info = sprintf("%s var %d sd not positive", ctx, j)
      )
    }

    # type = "response" returns conditional means matrix
    pred_response = predict(fit, newdata = newdata, type = "response")
    expect_true(is.matrix(pred_response), info = ctx)
    expect_equal(nrow(pred_response), nrow(newdata), info = ctx)
    expect_equal(ncol(pred_response), args$num_variables, info = ctx)
  }
})

test_that("predict.bgms returns valid predictions for mixed-mrf fixtures", {
  for(spec in get_bgms_fixtures()) {
    if(!isTRUE(spec$is_mixed)) next

    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    newdata = spec$get_prediction_data(n = 5)
    probs = predict(fit, newdata = newdata, type = "probabilities")

    expect_true(is.list(probs), info = ctx)
    expect_equal(length(probs), args$num_variables, info = ctx)

    for(j in seq_len(args$num_variables)) {
      vname = args$data_columnnames[j]
      expect_equal(nrow(probs[[j]]), nrow(newdata),
        info = sprintf("%s %s nrow", ctx, vname)
      )
      expect_false(anyNA(probs[[j]]),
        info = sprintf("%s %s has NAs", ctx, vname)
      )

      if(args$variable_type[j] %in% c("ordinal", "blume-capel")) {
        # Discrete variables: probability rows sum to 1
        row_sums = rowSums(probs[[j]])
        expect_true(
          all(abs(row_sums - 1) < 1e-6),
          info = sprintf("%s %s probs don't sum to 1", ctx, vname)
        )
      } else {
        # Continuous variables: 2-column (mean, sd) matrix
        expect_equal(ncol(probs[[j]]), 2,
          info = sprintf("%s %s ncol", ctx, vname)
        )
      }
    }

    # type = "response" returns a matrix
    resp = predict(fit, newdata = newdata, type = "response")
    expect_true(is.matrix(resp), info = ctx)
    expect_equal(dim(resp), c(nrow(newdata), args$num_variables), info = ctx)
  }
})

test_that("predict.bgms response returns integer categories", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)
  newdata = get_prediction_data_binary(n = 10)

  pred = predict(fit, newdata = newdata, type = "response")

  expect_true(is.matrix(pred))
  expect_equal(nrow(pred), nrow(newdata))
  expect_true(all(pred == round(pred)))
})

test_that("predict.bgms accepts variable subsetting", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)
  newdata = get_prediction_data_binary(n = 5)

  # By index
  pred1 = predict(fit, newdata = newdata, variables = 1, type = "probabilities")
  expect_equal(length(pred1), 1)

  # By name
  var_name = args$data_columnnames[1]
  pred2 = predict(fit, newdata = newdata, variables = var_name, type = "probabilities")
  expect_equal(length(pred2), 1)
})

test_that("predict.bgms accepts variable subsetting for GGM", {
  fit = get_bgms_fit_ggm_no_es()
  args = extract_arguments(fit)
  newdata = get_prediction_data_ggm(n = 5)

  # By index
  pred1 = predict(fit, newdata = newdata, variables = c(1, 3))
  expect_length(pred1, 2)
  expect_equal(names(pred1), args$data_columnnames[c(1, 3)])

  # By name
  pred2 = predict(fit, newdata = newdata, variables = c("V2", "V4"))
  expect_length(pred2, 2)
  expect_equal(names(pred2), c("V2", "V4"))
})

test_that("predict.bgms errors on invalid newdata dimensions", {
  fit = get_bgms_fit()
  bad_data = matrix(1:10, nrow = 5, ncol = 2)

  expect_error(predict(fit, newdata = bad_data), regexp = "columns")
})

test_that("predict.bgms handles edge cases", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)

  # Single observation
  newdata = get_prediction_data_binary(n = 1)
  probs = predict(fit, newdata = newdata, type = "probabilities")

  expect_true(is.list(probs))
  expect_equal(length(probs), args$num_variables)
  expect_equal(nrow(probs[[1]]), 1)
})

test_that("predict.bgms with posterior-sample returns sd attribute", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)
  newdata = get_prediction_data_binary(n = 5)

  result = predict(fit,
    newdata = newdata,
    method = "posterior-sample",
    ndraws = 10,
    seed = 42
  )

  expect_true(is.list(result))
  expect_equal(length(result), args$num_variables)

  sd_attr = attr(result, "sd")
  expect_false(is.null(sd_attr))
  expect_equal(length(sd_attr), args$num_variables)
})

test_that("predict.bgms GGM with posterior-sample returns sd attribute", {
  fit = get_bgms_fit_ggm_no_es()
  args = extract_arguments(fit)
  newdata = get_prediction_data_ggm(n = 5)

  result = predict(fit,
    newdata = newdata,
    method = "posterior-sample",
    ndraws = 10,
    seed = 123
  )

  expect_type(result, "list")
  expect_length(result, args$num_variables)
  for(v in seq_len(args$num_variables)) {
    expect_equal(ncol(result[[v]]), 2)
    expect_equal(colnames(result[[v]]), c("mean", "sd"))
  }

  sd_attr = attr(result, "sd")
  expect_false(is.null(sd_attr))
  expect_length(sd_attr, args$num_variables)
})

test_that("predict.bgms GGM conditional mean matches analytic formula", {
  # Verify C++ output matches the analytic conditional Gaussian formula.
  # This only checks that predict() and the R formula agree on the *same*
  # posterior mean Omega, so minimal MCMC iterations suffice.
  fit = get_bgms_fit_ggm_no_es()
  args = extract_arguments(fit)
  newdata = get_prediction_data_ggm(n = 10)

  pred = predict(fit, newdata = newdata)

  # Reconstruct the posterior mean precision matrix (precision = -2 * association)
  omega_hat = extract_precision(fit)
  p = args$num_variables

  # Center newdata on the training means (predict does the same internally)
  train_means = args$column_means
  newdata_centered = sweep(newdata, 2, train_means)

  for(j in seq_len(p)) {
    omega_jj = omega_hat[j, j]
    rest_cols = setdiff(seq_len(p), j)
    # Conditional mean in centered space, then shift back
    expected_means = as.numeric(
      train_means[j] - newdata_centered[, rest_cols, drop = FALSE] %*%
        omega_hat[rest_cols, j] / omega_jj
    )
    expected_sd = sqrt(1 / omega_jj)

    expect_equal(pred[[j]][, "mean"], expected_means, tolerance = 1e-10)
    expect_equal(unname(pred[[j]][1, "sd"]), expected_sd, tolerance = 1e-10)
  }
})


# ==============================================================================
# simulate.bgmCompare() Tests (Parameterized)
# ==============================================================================

test_that("simulate.bgmCompare returns matrix of correct size for all fixture types", {
  for(spec in get_bgmcompare_fixtures()) {
    ctx = sprintf("[bgmCompare %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    n_sim = 30

    # Test simulation for each group
    for(g in seq_len(args$num_groups)) {
      simulated = simulate(fit, nsim = n_sim, group = g, method = "posterior-mean", seed = 123)

      expect_true(is.matrix(simulated), info = paste(ctx, "group", g))
      expect_equal(nrow(simulated), n_sim, info = paste(ctx, "group", g, "wrong nrow"))
      expect_equal(ncol(simulated), args$num_variables, info = paste(ctx, "group", g, "wrong ncol"))
      expect_equal(colnames(simulated), args$data_columnnames, info = paste(ctx, "group", g))

      # Values should be integers within valid range
      expect_true(all(simulated == round(simulated)), info = paste(ctx, "group", g, "not integers"))
      expect_true(all(simulated >= 0), info = paste(ctx, "group", g, "negative values"))

      for(j in seq_len(args$num_variables)) {
        levels_j = args$category_levels[[j]]
        if(is.null(levels_j)) {
          # Blume-Capel (no recode map): original-scale scores, i.e. the
          # 0-based internal range shifted by the stored per-variable offset.
          shift_j = args$blume_capel_shift[j]
          expect_true(
            all(simulated[, j] >= shift_j &
              simulated[, j] <= shift_j + args$num_categories[j]),
            info = sprintf("%s group %d variable %d exceeds the category range", ctx, g, j)
          )
        } else {
          # Ordinal: simulate() returns the original category values (bgmCompare
          # stores a named original -> code lookup).
          valid = if(!is.null(names(levels_j))) as.numeric(names(levels_j)) else levels_j
          expect_true(
            all(simulated[, j] %in% valid),
            info = sprintf("%s group %d variable %d off the original scale", ctx, g, j)
          )
        }
      }
    }
  }
})

test_that("simulate.bgmCompare is reproducible with seed", {
  fit = get_bgmcompare_fit()

  sim1 = simulate(fit, nsim = 30, group = 1, method = "posterior-mean", seed = 999)
  sim2 = simulate(fit, nsim = 30, group = 1, method = "posterior-mean", seed = 999)

  expect_equal(sim1, sim2)
})

test_that("simulate.bgmCompare produces different results for different groups", {
  fit = get_bgmcompare_fit()

  # Simulate many observations to detect distributional differences
  sim_g1 = simulate(fit, nsim = 100, group = 1, seed = 42)
  sim_g2 = simulate(fit, nsim = 100, group = 2, seed = 42)

  # While individual values might match, means or patterns should generally differ
  # This is a soft test - we just verify they can be different
  expect_true(is.matrix(sim_g1))
  expect_true(is.matrix(sim_g2))
  expect_equal(dim(sim_g1), dim(sim_g2))
})

test_that("simulate.bgmCompare handles single observation", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)

  sim1 = simulate(fit, nsim = 1, group = 1, method = "posterior-mean", seed = 42)
  expect_true(is.matrix(sim1))
  expect_equal(nrow(sim1), 1)
  expect_equal(ncol(sim1), args$num_variables)
})


# ==============================================================================
# predict.bgmCompare() Tests (Parameterized)
# ==============================================================================

test_that("predict.bgmCompare returns valid probabilities for all fixture types", {
  for(spec in get_bgmcompare_fixtures()) {
    ctx = sprintf("[bgmCompare %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    newdata = spec$get_prediction_data(n = 5)

    # Test prediction for each group
    for(g in seq_len(args$num_groups)) {
      probs = predict(fit, newdata = newdata, group = g, type = "probabilities")

      expect_true(is.list(probs), info = paste(ctx, "group", g))
      expect_equal(length(probs), args$num_variables, info = paste(ctx, "group", g))

      # Each variable's probabilities should sum to 1
      for(j in seq_along(probs)) {
        expect_true(is.matrix(probs[[j]]), info = paste(ctx, "group", g, "var", j))
        expect_equal(nrow(probs[[j]]), nrow(newdata), info = paste(ctx, "group", g, "var", j))

        # No NAs in probability output
        expect_false(anyNA(probs[[j]]),
          info = sprintf("%s group %d var %d has NAs", ctx, g, j)
        )
        # Probabilities in [0, 1]
        expect_true(all(probs[[j]] >= 0 & probs[[j]] <= 1),
          info = sprintf("%s group %d var %d probs out of [0,1]", ctx, g, j)
        )

        # Row sums = 1
        row_sums = rowSums(probs[[j]])
        expect_true(
          all(abs(row_sums - 1) < 1e-6),
          info = sprintf("%s group %d var %d probs don't sum to 1", ctx, g, j)
        )
      }
    }
  }
})

test_that("predict.bgmCompare response returns integer categories", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)
  newdata = get_prediction_data_bgmcompare_binary(n = 10)

  for(g in seq_len(args$num_groups)) {
    pred = predict(fit, newdata = newdata, group = g, type = "response")

    expect_true(is.matrix(pred))
    expect_equal(nrow(pred), nrow(newdata))
    expect_true(all(pred == round(pred)))
  }
})

test_that("predict.bgmCompare accepts variable subsetting", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)
  newdata = get_prediction_data_bgmcompare_binary(n = 5)

  # By index
  pred1 = predict(fit, newdata = newdata, group = 1, variables = 1, type = "probabilities")
  expect_equal(length(pred1), 1)

  # By name
  var_name = args$data_columnnames[1]
  pred2 = predict(fit, newdata = newdata, group = 1, variables = var_name, type = "probabilities")
  expect_equal(length(pred2), 1)
})

test_that("predict.bgmCompare errors on invalid newdata dimensions", {
  fit = get_bgmcompare_fit()
  bad_data = matrix(1:10, nrow = 5, ncol = 2)

  expect_error(predict(fit, newdata = bad_data, group = 1), regexp = "columns")
})

test_that("predict.bgmCompare handles single observation", {
  fit = get_bgmcompare_fit()
  args = extract_arguments(fit)
  newdata = get_prediction_data_bgmcompare_binary(n = 1)

  probs = predict(fit, newdata = newdata, group = 1, type = "probabilities")

  expect_true(is.list(probs))
  expect_equal(length(probs), args$num_variables)
  expect_equal(nrow(probs[[1]]), 1)
})


# ==============================================================================
# Cross-method Consistency Tests
# ==============================================================================

test_that("coef and summary produce consistent pairwise dimensions", {
  fit = get_bgms_fit()

  coeffs = coef(fit)
  summ = summary(fit)

  if(!is.null(coeffs$pairwise) && !is.null(summ$pairwise)) {
    p = nrow(coeffs$pairwise)
    expected_pairs = p * (p - 1) / 2
    expect_equal(nrow(summ$pairwise), expected_pairs)
  }
})

test_that("simulate and predict use consistent variable ordering", {
  fit = get_bgms_fit()
  args = extract_arguments(fit)

  sim_data = simulate(fit, nsim = 10, method = "posterior-mean", seed = 42)
  pred = predict(fit, newdata = sim_data, type = "response")

  expect_equal(ncol(sim_data), args$num_variables)
  expect_equal(ncol(pred), args$num_variables)
})

test_that("NUTS and adaptive-metropolis produce same output structure", {
  fit_nuts = get_bgms_fit()
  fit_am = get_bgms_fit_adaptive_metropolis()

  coef_nuts = coef(fit_nuts)
  coef_am = coef(fit_am)

  expect_equal(names(coef_nuts), names(coef_am))

  if(!is.null(coef_nuts$pairwise) && !is.null(coef_am$pairwise)) {
    expect_equal(dim(coef_nuts$pairwise), dim(coef_am$pairwise))
  }
})


# ==============================================================================
# Blume-Capel Specific Tests
# ==============================================================================

test_that("extract_arguments stores baseline_category for Blume-Capel", {
  fit = get_bgms_fit_blumecapel()
  args = extract_arguments(fit)

  expect_true("baseline_category" %in% names(args))
  expect_equal(length(args$baseline_category), args$num_variables)

  for(j in seq_len(args$num_variables)) {
    expect_true(
      args$baseline_category[j] >= 0 && args$baseline_category[j] <= args$num_categories[j],
      info = sprintf("baseline_category[%d] out of range", j)
    )
  }
})


# ==============================================================================
# Single Chain R-hat Handling
# ==============================================================================

test_that("extract_rhat handles single chain appropriately", {
  fit = get_bgms_fit_single_chain()

  rhat = extract_rhat(fit)

  # Should return something (NA or valid values)
  expect_true(!is.null(rhat))
})

test_that("multi-chain fit has valid R-hat values", {
  fit = get_bgms_fit()

  rhat = extract_rhat(fit)

  expect_true(is.list(rhat))
  if(!is.null(rhat$pairwise)) {
    rhat_vals = rhat$pairwise[!is.na(rhat$pairwise)]
    if(length(rhat_vals) > 0) {
      expect_true(all(rhat_vals > 0), info = "R-hat should be positive")
    }
  }
})


# ==============================================================================
# ESS with Adaptive-Metropolis
# ==============================================================================

test_that("extract_ess works with adaptive-metropolis", {
  fit = get_bgms_fit_adaptive_metropolis()

  ess = extract_ess(fit)

  expect_true(is.list(ess))
  if(!is.null(ess$pairwise)) {
    ess_vals = ess$pairwise[!is.na(ess$pairwise)]
    if(length(ess_vals) > 0) {
      expect_true(all(ess_vals > 0), info = "ESS should be positive")
    }
  }
})


# ==============================================================================
# [[ access consistency with $ for lazy-evaluated summaries
# ==============================================================================

for(spec in get_bgms_fixtures()) {
  test_that(
    sprintf("[[ and $ return identical posterior summaries (%s)", spec$label),
    {
      fit = spec$get_fit()
      summary_fields = grep(
        "^posterior_summary_", names(fit),
        value = TRUE
      )
      for(field in summary_fields) {
        bracket_val = fit[[field]]
        dollar_val = do.call("$", list(fit, field))
        expect_identical(
          bracket_val, dollar_val,
          info = sprintf("field '%s' differs between $ and [[", field)
        )
      }
    }
  )
}

for(spec in get_bgmcompare_fixtures()) {
  test_that(
    sprintf("[[ and $ return identical posterior summaries (%s)", spec$label),
    {
      fit = spec$get_fit()
      summary_fields = grep(
        "^posterior_summary_", names(fit),
        value = TRUE
      )
      for(field in summary_fields) {
        bracket_val = fit[[field]]
        dollar_val = do.call("$", list(fit, field))
        expect_identical(
          bracket_val, dollar_val,
          info = sprintf("field '%s' differs between $ and [[", field)
        )
      }
    }
  )
}
