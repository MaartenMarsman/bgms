# ==============================================================================
# What the R layer around the prediction kernels owes its caller.
#
# The C++ kernels take the recoded observation matrix at face value and know
# nothing about the fit's category scale, so four contracts live entirely in R:
#
#   - a row whose conditioning set is incomplete comes back as NA, not as a
#     confident answer computed from R's NA_integer_ sentinel;
#   - responses and probability labels are on the ORIGINAL category scale --
#     the one simulate() returns and predict() reads back;
#   - one row of newdata is a shape like any other;
#   - the posterior-mean precision matrix is a posterior mean throughout,
#     diagonal included;
#   - newdata columns are matched to the fit deliberately, not by width alone.
# ==============================================================================

.contract_cache = new.env(parent = emptyenv())

# Three recode maps, none of them the identity: an ordinal variable on 1..3, a
# non-contiguous ordinal on {1, 5, 9}, and a Blume-Capel variable on 2..5.
get_original_scale_fit = function() {
  if(is.null(.contract_cache$fit)) {
    set.seed(1)
    n = 200
    x = cbind(
      a = sample(c(1, 2, 3), n, replace = TRUE),
      b = sample(c(1, 5, 9), n, replace = TRUE),
      bc = sample(2:5, n, replace = TRUE)
    )
    .contract_cache$fit = bgm(
      x = x,
      variable_type = c("ordinal", "ordinal", "blume-capel"),
      baseline_category = c(1L, 1L, 3L),
      iter = 200, warmup = 200, chains = 1, cores = 1,
      update_method = "adaptive-metropolis",
      display_progress = "none", seed = 7
    )
  }
  .contract_cache$fit
}

collect_warnings = function(expr) {
  seen = character()
  withCallingHandlers(
    force(expr),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  seen
}


# ------------------------------------------------------------------------------
# Original category scale (simulate <-> predict round trip)
# ------------------------------------------------------------------------------

test_that("simulate() and predict(type = 'response') meet on the original scale", {
  skip_on_cran()
  fit = get_original_scale_fit()

  sim = simulate(fit, nsim = 30, method = "posterior-mean", seed = 3)
  expect_true(all(sim[, "a"] %in% c(1, 2, 3)))
  expect_true(all(sim[, "b"] %in% c(1, 5, 9)))
  expect_true(all(sim[, "bc"] %in% 2:5))

  resp = predict(fit, newdata = sim, type = "response")
  expect_equal(colnames(resp), c("a", "b", "bc"))
  expect_true(all(resp[, "a"] %in% c(1, 2, 3)))
  expect_true(all(resp[, "b"] %in% c(1, 5, 9)))
  expect_true(all(resp[, "bc"] %in% 2:5))

  # The round trip closes: the responses are themselves valid newdata, with no
  # unseen-category warning and no NA anywhere.
  expect_no_warning(again <- predict(fit, newdata = resp, type = "response"))
  expect_false(anyNA(again))
})

test_that("probability columns are labelled with the original categories", {
  skip_on_cran()
  fit = get_original_scale_fit()
  sim = simulate(fit, nsim = 10, method = "posterior-mean", seed = 4)

  probs = predict(fit, newdata = sim, type = "probabilities")
  expect_equal(colnames(probs$a), paste0("cat_", c(1, 2, 3)))
  expect_equal(colnames(probs$b), paste0("cat_", c(1, 5, 9)))
  expect_equal(colnames(probs$bc), paste0("cat_", 2:5))

  # The reported response is the value its own label names.
  resp = predict(fit, newdata = sim, type = "response")
  for(v in names(probs)) {
    labels = as.numeric(sub("^cat_", "", colnames(probs[[v]])))
    expect_equal(
      unname(labels[max.col(probs[[v]], ties.method = "first")]),
      unname(resp[, v])
    )
  }
})

test_that("posterior-sample prediction carries the original labels too", {
  skip_on_cran()
  fit = get_original_scale_fit()
  sim = simulate(fit, nsim = 8, method = "posterior-mean", seed = 11)

  probs = predict(
    fit,
    newdata = sim, type = "probabilities",
    method = "posterior-sample", ndraws = 5, seed = 2
  )
  expect_equal(colnames(probs$b), paste0("cat_", c(1, 5, 9)))
  expect_equal(colnames(attr(probs, "sd")$b), paste0("cat_", c(1, 5, 9)))
})


# ------------------------------------------------------------------------------
# Incomplete conditioning
# ------------------------------------------------------------------------------

test_that("an unseen category or a plain NA blanks the row's other variables", {
  skip_on_cran()
  fit = get_original_scale_fit()
  sim = simulate(fit, nsim = 6, method = "posterior-mean", seed = 5)

  newdata = sim
  newdata[1, "a"] = 99 # a value never observed in training
  newdata[2, "b"] = NA # a plain NA

  probs = suppressWarnings(
    predict(fit, newdata = newdata, type = "probabilities")
  )
  # The other variables in a corrupted row are undefined ...
  expect_true(all(is.na(probs$b[1, ])))
  expect_true(all(is.na(probs$bc[1, ])))
  expect_true(all(is.na(probs$a[2, ])))
  expect_true(all(is.na(probs$bc[2, ])))
  # ... while the variable whose own value is missing is not: its conditional
  # never reads that value. Clean rows are untouched.
  expect_false(anyNA(probs$a[1, ]))
  expect_false(anyNA(probs$b[2, ]))
  expect_false(anyNA(probs$a[3:6, ]))
  expect_false(anyNA(probs$bc[3:6, ]))

  resp = suppressWarnings(predict(fit, newdata = newdata, type = "response"))
  expect_true(all(is.na(resp[1, c("b", "bc")])))
  expect_true(all(is.na(resp[2, c("a", "bc")])))
  expect_false(anyNA(resp[3:6, ]))
})

test_that("incomplete conditioning warns once, with the affected row count", {
  skip_on_cran()
  fit = get_original_scale_fit()
  sim = simulate(fit, nsim = 6, method = "posterior-mean", seed = 5)

  newdata = sim
  newdata[1, "a"] = 99
  newdata[2, "b"] = NA

  seen = collect_warnings(predict(fit, newdata = newdata))
  hits = grep("conditioning variable is missing", seen, fixed = TRUE)
  expect_length(hits, 1L)
  expect_match(seen[hits], "2 row\\(s\\)")

  # A plain NA on its own used to pass in silence.
  only_na = sim
  only_na[4, "bc"] = NA
  seen_na = collect_warnings(predict(fit, newdata = only_na))
  expect_length(grep("conditioning variable is missing", seen_na, fixed = TRUE), 1L)
  expect_match(seen_na, "1 row\\(s\\)", all = FALSE)
})

test_that("a mixed fit blanks continuous predictions on a missing discrete cell", {
  skip_on_cran()
  fit = get_bgms_fit_mixed_mrf()
  newdata = get_prediction_data_mixed(n = 5)
  newdata[1, "d1"] = NA

  probs = suppressWarnings(
    predict(fit, newdata = newdata, type = "probabilities")
  )
  # The discrete rest score is what the sentinel corrupts, and it feeds the
  # continuous conditional as well as the discrete ones.
  expect_true(all(is.na(probs$c1[1, ])))
  expect_true(all(is.na(probs$d2[1, ])))
  expect_false(anyNA(probs$d1[1, ]))
  expect_false(anyNA(probs$c1[2:5, ]))

  resp = suppressWarnings(predict(fit, newdata = newdata, type = "response"))
  expect_true(all(is.na(resp[1, c("c1", "d2", "c2", "d3")])))
  expect_false(anyNA(resp[2:5, ]))
})


# ------------------------------------------------------------------------------
# One row of newdata
# ------------------------------------------------------------------------------

test_that("a single row of newdata predicts for every model family", {
  skip_on_cran()

  single_row_ok = function(fit, newdata, ...) {
    arguments = extract_arguments(fit)
    resp = predict(fit, newdata = newdata, type = "response", ...)
    expect_true(is.matrix(resp))
    expect_equal(dim(resp), c(1L, arguments$num_variables))
    expect_equal(colnames(resp), arguments$data_columnnames)
    expect_false(anyNA(resp))
  }

  single_row_ok(get_bgms_fit(), get_prediction_data_binary(n = 1))
  single_row_ok(get_bgms_fit_ggm(), get_prediction_data_ggm(n = 1))
  single_row_ok(get_bgms_fit_mixed_mrf(), get_prediction_data_mixed(n = 1))
  single_row_ok(
    get_bgmcompare_fit(), get_prediction_data_bgmcompare_binary(n = 1),
    group = 1
  )
})


# ------------------------------------------------------------------------------
# The posterior-mean precision diagonal
# ------------------------------------------------------------------------------

test_that("the GGM precision diagonal is the pooled mean of the raw draws", {
  skip_on_cran()
  fit = get_bgms_fit_ggm()

  raw_diagonal = do.call(rbind, get_raw_samples(fit)$main)
  omega = bgms:::reconstruct_precision(
    get_posterior_mean(fit, "pairwise"),
    bgms:::posterior_mean_precision_diagonal(fit)
  )
  expect_equal(
    unname(diag(omega)), unname(colMeans(raw_diagonal)),
    tolerance = 1e-12
  )

  # Not the reciprocal of the residual variance, which is a harmonic mean of
  # the same draws and so a different number.
  expect_false(isTRUE(all.equal(
    unname(diag(omega)),
    unname(1 / get_posterior_mean(fit, "residual_variance")),
    tolerance = 1e-8
  )))

  # Omega is now exactly the mean of the per-draw precision matrices, and a
  # mean of positive-definite matrices is positive definite.
  raw_pairwise = do.call(rbind, get_raw_samples(fit)$pairwise)
  p = extract_arguments(fit)$num_variables
  pooled = matrix(0, p, p)
  for(i in seq_len(nrow(raw_pairwise))) {
    pooled = pooled + bgms:::build_precision_from_draw(
      raw_pairwise[i, ], raw_diagonal[i, ], p
    )
  }
  pooled = pooled / nrow(raw_pairwise)
  expect_equal(unname(omega), unname(pooled), tolerance = 1e-12)
  expect_true(all(eigen(omega, symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("the mixed continuous diagonal comes from its own draws", {
  skip_on_cran()
  fit = get_bgms_fit_mixed_mrf()
  arguments = extract_arguments(fit)

  params = bgms:::build_mixed_params_mean(fit, arguments)
  draws = bgms:::mixed_cont_diagonal_draws(fit, arguments)
  expect_equal(
    unname(diag(params$pairwise_cont)), unname(colMeans(draws)),
    tolerance = 1e-12
  )

  precision = -2 * params$pairwise_cont
  expect_true(
    all(eigen(precision, symmetric = TRUE, only.values = TRUE)$values > 0)
  )
})


# ------------------------------------------------------------------------------
# newdata columns
# ------------------------------------------------------------------------------

test_that("named newdata must carry the fit's variables, in its order", {
  fit = get_bgms_fit()
  nodes = extract_arguments(fit)$data_columnnames
  newdata = as.matrix(get_prediction_data_binary(n = 4))

  renamed = newdata
  colnames(renamed) = paste0("x", seq_along(nodes))
  expect_error(
    predict(fit, newdata = renamed),
    "column names do not match the fitted model"
  )

  reordered = newdata[, rev(seq_along(nodes)), drop = FALSE]
  expect_error(predict(fit, newdata = reordered), "in a different order")

  expect_no_error(predict(fit, newdata = newdata))
})

test_that("unnamed newdata is accepted by position, and says so", {
  fit = get_bgms_fit()
  newdata = as.matrix(get_prediction_data_binary(n = 4))
  dimnames(newdata) = NULL

  expect_warning(
    predict(fit, newdata = newdata),
    "matched to the fitted model by position"
  )
})

test_that("bgmCompare applies the same newdata column rules", {
  fit = get_bgmcompare_fit()
  nodes = extract_arguments(fit)$data_columnnames
  newdata = as.matrix(get_prediction_data_bgmcompare_binary(n = 4))

  renamed = newdata
  colnames(renamed) = paste0("x", seq_along(nodes))
  expect_error(
    predict(fit, newdata = renamed, group = 1),
    "column names do not match the fitted model"
  )

  dimnames(newdata) = NULL
  expect_warning(
    predict(fit, newdata = newdata, group = 1),
    "matched to the fitted model by position"
  )
})


# ------------------------------------------------------------------------------
# simulate() entry validation
# ------------------------------------------------------------------------------

test_that("simulate.bgms validates nsim and iter for every path", {
  fit = get_bgms_fit()
  expect_error(simulate(fit, nsim = 0), "nsim")
  expect_error(simulate(fit, nsim = 10.5), "nsim")
  expect_error(simulate(fit, nsim = 10, iter = -1), "iter")
})

test_that("simulate.bgms survives a machine that cannot report its cores", {
  # parallel::detectCores() is documented to return NA when it cannot tell, and
  # that is the argument's default.
  fit = get_bgms_fit()
  expect_no_error(simulate(fit, nsim = 5, iter = 10, cores = NA_integer_))
})
