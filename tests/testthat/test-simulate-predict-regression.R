# ==============================================================================
# Phase C.6: simulate / predict regression tests
#
# Verifies that simulate() and predict() work correctly with fit objects
# produced by the refactored bgm() / bgmCompare() pipeline.
#
# Test groups:
#   1. $arguments contract <U+2014> every field simulate/predict extract is present
#   2. Fit-object structure <U+2014> posterior_mean_*, raw_samples present
#   4. Functional roundtrip <U+2014> simulate <U+2192> predict for each model type
#   5. Posterior-sample method regression
#   6. Field type and value invariants
# ==============================================================================


# ------------------------------------------------------------------------------
# Fixture Specifications <U+2014> defined in helper-fixtures.R
# get_bgms_fixtures(), get_bgmcompare_fixtures()
# ------------------------------------------------------------------------------


# ==============================================================================
# 1. $arguments contract tests
# ==============================================================================
# simulate.bgms / predict.bgms read these fields from extract_arguments():
#   num_variables, num_categories, variable_type, data_columnnames,
#   baseline_category, is_continuous (GGM only, via isTRUE guard)
#
# simulate.bgmCompare / predict.bgmCompare read:
#   num_groups, num_variables, num_categories, is_ordinal_variable,
#   data_columnnames, projection, baseline_category (NULL-safe)
# ==============================================================================

# Fields required by ALL bgms simulate/predict paths:
BGMS_COMMON_FIELDS = c(
  "num_variables", "variable_type", "data_columnnames"
)

# Additional fields required only for OMRF (ordinal / blume-capel):
BGMS_OMRF_FIELDS = c(
  "num_categories", "baseline_category"
)

COMPARE_SIM_PRED_FIELDS = c(
  "num_groups", "num_variables", "num_categories",
  "is_ordinal_variable", "data_columnnames", "projection"
)

test_that("bgms $arguments contains all fields needed by simulate/predict", {
  for(spec in get_bgms_fixtures()) {
    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    for(field in BGMS_COMMON_FIELDS) {
      expect_true(
        field %in% names(args),
        info = sprintf("%s: missing arguments$%s", ctx, field)
      )
    }

    if(isTRUE(spec$is_continuous)) {
      # GGM fits must carry is_continuous = TRUE
      expect_true(
        isTRUE(args$is_continuous),
        info = sprintf("%s: is_continuous should be TRUE for GGM", ctx)
      )
    } else if(isTRUE(spec$is_mixed)) {
      # Mixed MRF: OMRF fields plus mixed-specific fields
      for(field in BGMS_OMRF_FIELDS) {
        expect_true(
          field %in% names(args),
          info = sprintf("%s: missing arguments$%s", ctx, field)
        )
      }
      for(field in c(
        "is_mixed", "discrete_indices", "continuous_indices",
        "num_discrete", "num_continuous", "is_ordinal",
        "data_columnnames_discrete", "data_columnnames_continuous"
      )) {
        expect_true(
          field %in% names(args),
          info = sprintf("%s: missing mixed arguments$%s", ctx, field)
        )
      }
      expect_true(isTRUE(args$is_mixed),
        info = sprintf("%s: is_mixed should be TRUE", ctx)
      )
    } else {
      # OMRF fits must also carry num_categories and baseline_category
      for(field in BGMS_OMRF_FIELDS) {
        expect_true(
          field %in% names(args),
          info = sprintf("%s: missing arguments$%s", ctx, field)
        )
      }
    }
  }
})

test_that(paste(
  "bgmCompare $arguments contains all fields",
  "needed by simulate/predict"
), {
  for(spec in get_bgmcompare_fixtures()) {
    ctx = sprintf("[bgmCompare %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    for(field in COMPARE_SIM_PRED_FIELDS) {
      expect_true(
        field %in% names(args),
        info = sprintf("%s: missing arguments$%s", ctx, field)
      )
    }
  }
})


# ==============================================================================
# 2. Fit-object structure tests
# ==============================================================================
# simulate/predict also read directly from the fit object:
#   - posterior_mean_pairwise, posterior_mean_main  (posterior-mean method)
#   - raw_samples$pairwise, raw_samples$main       (posterior-sample method)
# ==============================================================================

test_that("bgms fit objects have posterior_mean fields for simulate/predict", {
  for(spec in get_bgms_fixtures()) {
    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)
    p = args$num_variables

    expect_false(is.null(fit$posterior_mean_pairwise),
      info = paste(ctx, "missing posterior_mean_pairwise")
    )
    expect_true(is.matrix(fit$posterior_mean_pairwise),
      info = paste(ctx, "posterior_mean_pairwise not a matrix")
    )
    expect_equal(nrow(fit$posterior_mean_pairwise), p,
      info = paste(ctx, "posterior_mean_pairwise wrong nrow")
    )
    expect_equal(ncol(fit$posterior_mean_pairwise), p,
      info = paste(ctx, "posterior_mean_pairwise wrong ncol")
    )

    if(isTRUE(args$is_continuous)) {
      # GGM: no main effects; precision diagonal stored separately
      expect_null(fit$posterior_mean_main,
        info = paste(ctx, "GGM posterior_mean_main should be NULL")
      )
      expect_true(all(fit$posterior_mean_residual_variance > 0),
        info = paste(ctx, "GGM residual variance should be positive (1/precision_ii)")
      )
    } else if(isTRUE(spec$is_mixed)) {
      expect_true(is.list(fit$posterior_mean_main),
        info = paste(ctx, "mixed posterior_mean_main should be a list")
      )
      expect_false(is.null(fit$posterior_mean_main$discrete),
        info = paste(ctx, "missing posterior_mean_main$discrete")
      )
      expect_false(is.null(fit$posterior_mean_main$continuous),
        info = paste(ctx, "missing posterior_mean_main$continuous")
      )
    } else {
      expect_false(is.null(fit$posterior_mean_main),
        info = paste(ctx, "missing posterior_mean_main")
      )
    }
  }
})

test_that("bgms fit objects have raw_samples for posterior-sample method", {
  for(spec in get_bgms_fixtures()) {
    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()

    expect_false(is.null(fit$raw_samples),
      info = paste(ctx, "missing raw_samples")
    )
    expect_false(is.null(fit$raw_samples$pairwise),
      info = paste(ctx, "missing raw_samples$pairwise")
    )
    expect_false(is.null(fit$raw_samples$main),
      info = paste(ctx, "missing raw_samples$main")
    )
    expect_true(is.list(fit$raw_samples$pairwise),
      info = paste(ctx, "raw_samples$pairwise not a list")
    )
    expect_true(is.list(fit$raw_samples$main),
      info = paste(ctx, "raw_samples$main not a list")
    )
  }
})



# ==============================================================================
# 4. Functional roundtrip tests
# ==============================================================================
# For every cached fixture type: simulate data <U+2192> predict on it <U+2192> verify
# structural soundness. This catches any mismatch between $arguments and
# the actual simulate/predict code paths after refactoring.
# ==============================================================================

test_that("simulate <U+2192> predict roundtrip works for all bgms fixtures", {
  for(spec in get_bgms_fixtures()) {
    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    n_sim = 20
    simulated = simulate(fit, nsim = n_sim, method = "posterior-mean", seed = 1)

    expect_true(is.matrix(simulated), info = paste(ctx, "simulate"))
    expect_equal(nrow(simulated), n_sim, info = paste(ctx, "nrow"))
    expect_equal(ncol(simulated), args$num_variables, info = paste(ctx, "ncol"))
    expect_equal(
      colnames(simulated), args$data_columnnames,
      info = paste(ctx, "colnames")
    )

    if(isTRUE(args$is_continuous)) {
      # GGM: predict returns list of mean/sd matrices
      colnames(simulated) = args$data_columnnames
      pred = predict(fit, newdata = simulated)
      expect_true(is.list(pred),
        info = paste(ctx, "predict type")
      )
      expect_equal(
        length(pred), args$num_variables,
        info = paste(ctx, "predict length")
      )
      for(j in seq_along(pred)) {
        expect_equal(nrow(pred[[j]]), n_sim,
          info = sprintf("%s predict var %d nrow", ctx, j)
        )
        expect_equal(ncol(pred[[j]]), 2,
          info = sprintf("%s predict var %d ncol", ctx, j)
        )
      }
    } else if(isTRUE(spec$is_mixed)) {
      # Mixed MRF: predict returns list with
      # discrete (probs) and continuous (mean/sd)
      probs = predict(
        fit,
        newdata = simulated,
        type = "probabilities"
      )
      expect_true(is.list(probs),
        info = paste(ctx, "predict type")
      )
      expect_equal(
        length(probs), args$num_variables,
        info = paste(ctx, "predict length")
      )

      for(j in seq_len(args$num_variables)) {
        vname = args$data_columnnames[j]
        expect_equal(nrow(probs[[j]]), n_sim,
          info = sprintf("%s predict %s nrow", ctx, vname)
        )
        expect_false(anyNA(probs[[j]]),
          info = sprintf("%s predict %s has NAs", ctx, vname)
        )

        if(args$variable_type[j] %in% c("ordinal", "blume-capel")) {
          # Discrete: probability rows sum to 1
          row_sums = rowSums(probs[[j]])
          expect_true(
            all(abs(row_sums - 1) < 1e-6),
            info = sprintf("%s predict %s probs don't sum to 1", ctx, vname)
          )
        } else {
          # Continuous: 2-column (mean, sd) matrix
          expect_equal(ncol(probs[[j]]), 2,
            info = sprintf("%s predict %s ncol", ctx, vname)
          )
        }
      }

      # type = "response" should return a matrix
      resp = predict(fit, newdata = simulated, type = "response")
      expect_true(is.matrix(resp), info = paste(ctx, "response matrix"))
      expect_equal(dim(resp), c(n_sim, args$num_variables),
        info = paste(ctx, "response dim")
      )
    } else {
      # OMRF: predict returns list of probability matrices
      probs = predict(
        fit,
        newdata = simulated,
        type = "probabilities"
      )
      expect_true(is.list(probs),
        info = paste(ctx, "predict type")
      )
      expect_equal(
        length(probs), args$num_variables,
        info = paste(ctx, "predict length")
      )
      for(j in seq_along(probs)) {
        expect_equal(nrow(probs[[j]]), n_sim,
          info = sprintf("%s predict var %d nrow", ctx, j)
        )
        # No NAs in probability output
        expect_false(anyNA(probs[[j]]),
          info = sprintf("%s predict var %d has NAs", ctx, j)
        )
        # Probability rows should sum to 1
        row_sums = rowSums(probs[[j]])
        expect_true(
          all(abs(row_sums - 1) < 1e-6),
          info = sprintf("%s predict var %d probs don't sum to 1", ctx, j)
        )
      }

      # type = "response" should work for all model types
      resp = predict(fit, newdata = simulated, type = "response")
      expect_true(is.matrix(resp), info = paste(ctx, "response matrix"))
      expect_equal(dim(resp), c(n_sim, args$num_variables),
        info = paste(ctx, "response dim")
      )
      expect_true(all(resp == round(resp)),
        info = paste(ctx, "response not integers")
      )
    }
  }
})

test_that("simulate <U+2192> predict roundtrip works for all bgmCompare fixtures", {
  for(spec in get_bgmcompare_fixtures()) {
    ctx = sprintf("[bgmCompare %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    n_sim = 20

    for(g in seq_len(args$num_groups)) {
      g_ctx = sprintf("%s group %d", ctx, g)

      simulated = simulate(fit,
        nsim = n_sim, group = g,
        method = "posterior-mean", seed = 1
      )

      expect_true(is.matrix(simulated), info = paste(g_ctx, "simulate"))
      expect_equal(nrow(simulated), n_sim, info = paste(g_ctx, "nrow"))
      expect_equal(
        ncol(simulated), args$num_variables,
        info = paste(g_ctx, "ncol")
      )
      expect_equal(colnames(simulated), args$data_columnnames,
        info = paste(g_ctx, "colnames")
      )

      # Values should be non-negative integers
      expect_true(all(simulated >= 0),
        info = paste(g_ctx, "negative values")
      )
      expect_true(all(simulated == round(simulated)),
        info = paste(g_ctx, "not integers")
      )

      # Predict
      probs = predict(
        fit,
        newdata = simulated,
        group = g, type = "probabilities"
      )
      expect_true(is.list(probs), info = paste(g_ctx, "predict type"))
      expect_equal(length(probs), args$num_variables,
        info = paste(g_ctx, "predict length")
      )

      for(j in seq_along(probs)) {
        expect_equal(nrow(probs[[j]]), n_sim,
          info = sprintf("%s predict var %d nrow", g_ctx, j)
        )
        row_sums = rowSums(probs[[j]], na.rm = TRUE)
        valid = !apply(probs[[j]], 1, function(x) any(is.na(x)))
        if(any(valid)) {
          expect_true(
            all(abs(row_sums[valid] - 1) < 1e-6),
            info = sprintf("%s predict var %d probs don't sum to 1", g_ctx, j)
          )
        }
      }

      resp = predict(fit, newdata = simulated, group = g, type = "response")
      expect_true(is.matrix(resp), info = paste(g_ctx, "response matrix"))
      expect_equal(dim(resp), c(n_sim, args$num_variables),
        info = paste(g_ctx, "response dim")
      )
    }
  }
})


# ==============================================================================
# 5. Posterior-sample method regression
# ==============================================================================
# The posterior-sample path reads raw_samples$pairwise / raw_samples$main.
# Verify it works for each model type and produces an sd attribute.
# ==============================================================================

test_that("simulate posterior-sample method works for all bgms fixtures", {
  for(spec in get_bgms_fixtures()) {
    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    n_draws = 2
    n_sim = 10
    sim = simulate(fit,
      nsim = n_sim, method = "posterior-sample",
      ndraws = n_draws, seed = 42,
      display_progress = "none"
    )

    # posterior-sample returns a list of matrices (one per draw)
    expect_true(is.list(sim), info = paste(ctx, "not a list"))
    expect_equal(length(sim), n_draws, info = paste(ctx, "wrong length"))

    for(d in seq_len(n_draws)) {
      expect_true(is.matrix(sim[[d]]),
        info = sprintf("%s draw %d not a matrix", ctx, d)
      )
      expect_equal(nrow(sim[[d]]), n_sim,
        info = sprintf("%s draw %d wrong nrow", ctx, d)
      )
      expect_equal(ncol(sim[[d]]), args$num_variables,
        info = sprintf("%s draw %d wrong ncol", ctx, d)
      )
    }
  }
})

test_that("predict posterior-sample method works for all bgms fixtures", {
  for(spec in get_bgms_fixtures()) {
    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)

    newdata = spec$get_prediction_data(n = 5)
    result = predict(fit,
      newdata = newdata, method = "posterior-sample",
      ndraws = 2, seed = 42
    )

    expect_true(is.list(result), info = paste(ctx, "not a list"))
    expect_equal(length(result), args$num_variables,
      info = paste(ctx, "wrong length")
    )

    sd_attr = attr(result, "sd")
    expect_false(is.null(sd_attr), info = paste(ctx, "missing sd attribute"))
    expect_equal(length(sd_attr), args$num_variables,
      info = paste(ctx, "sd wrong length")
    )
  }
})


# ==============================================================================
# 6. $arguments field type and value invariants
# ==============================================================================
# Verify that field types and ranges are what simulate/predict expect.
# ==============================================================================

test_that("bgms $arguments field types are correct for simulate/predict", {
  for(spec in get_bgms_fixtures()) {
    ctx = sprintf("[bgms %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)
    p = args$num_variables

    expect_true(
      is.numeric(args$num_variables) &&
        length(args$num_variables) == 1,
      info = paste(ctx, "num_variables")
    )
    expect_true(args$num_variables >= 1,
      info = paste(ctx, "num_variables >= 1")
    )

    expect_true(is.character(args$variable_type),
      info = paste(ctx, "variable_type character")
    )
    expect_true(
      all(args$variable_type %in%
        c("ordinal", "blume-capel", "continuous")),
      info = paste(ctx, "variable_type values")
    )

    expect_true(
      is.character(args$data_columnnames) &&
        length(args$data_columnnames) == p,
      info = paste(ctx, "data_columnnames length")
    )

    if(!isTRUE(spec$is_continuous) && !isTRUE(spec$is_mixed)) {
      # OMRF-only fields
      expect_true(
        is.numeric(args$num_categories) &&
          length(args$num_categories) == p,
        info = paste(ctx, "num_categories length")
      )
      expect_true(all(args$num_categories >= 1),
        info = paste(ctx, "num_categories >= 1")
      )
      expect_true(
        is.numeric(args$baseline_category) &&
          length(args$baseline_category) == p,
        info = paste(ctx, "baseline_category length")
      )
    }

    if(isTRUE(spec$is_mixed)) {
      pd = args$num_discrete
      qc = args$num_continuous
      expect_equal(
        pd + qc, p,
        info = paste(
          ctx, "num_discrete + num_continuous == p"
        )
      )
      expect_true(
        is.numeric(args$num_categories) &&
          length(args$num_categories) == pd,
        info = paste(
          ctx,
          "mixed num_categories length == num_discrete"
        )
      )
      expect_true(
        is.numeric(args$baseline_category) &&
          length(args$baseline_category) == pd,
        info = paste(
          ctx,
          "mixed baseline_category length == num_discrete"
        )
      )
      expect_true(
        is.numeric(args$discrete_indices) &&
          length(args$discrete_indices) == pd,
        info = paste(ctx, "discrete_indices length")
      )
      expect_true(
        is.numeric(args$continuous_indices) &&
          length(args$continuous_indices) == qc,
        info = paste(ctx, "continuous_indices length")
      )
    }
  }
})

test_that(paste(
  "bgmCompare $arguments field types",
  "are correct for simulate/predict"
), {
  for(spec in get_bgmcompare_fixtures()) {
    ctx = sprintf("[bgmCompare %s]", spec$label)
    fit = spec$get_fit()
    args = extract_arguments(fit)
    p = args$num_variables

    expect_true(is.numeric(args$num_groups) && args$num_groups >= 2,
      info = paste(ctx, "num_groups")
    )

    expect_true(is.numeric(args$num_variables) && args$num_variables >= 1,
      info = paste(ctx, "num_variables")
    )

    expect_true(
      is.numeric(args$num_categories) &&
        length(args$num_categories) == p,
      info = paste(ctx, "num_categories length")
    )

    expect_true(
      is.logical(args$is_ordinal_variable) &&
        length(args$is_ordinal_variable) == p,
      info = paste(ctx, "is_ordinal_variable")
    )

    expect_true(
      is.character(args$data_columnnames) &&
        length(args$data_columnnames) == p,
      info = paste(ctx, "data_columnnames")
    )

    expect_true(is.matrix(args$projection),
      info = paste(ctx, "projection is matrix")
    )
    expect_equal(nrow(args$projection), args$num_groups,
      info = paste(ctx, "projection nrow")
    )
  }
})


# ==============================================================================
# 7. simulate() returns data on the original category scale
# ==============================================================================
# Regression: simulate() emitted internal 0-based category codes while predict()
# recodes newdata from the original category values, so simulate() -> predict()
# mismatched (a "category values not observed" warning plus a silently miscoded
# prediction context) for any ordinal variable not already 0-based. simulate()
# now inverts the recode map, keeping the round trip on the original scale.
# ==============================================================================

test_that("simulate() returns ordinal data on the original category scale", {
  skip_on_cran()

  set.seed(1)
  n = 200
  x = cbind(
    a = sample(0:2, n, replace = TRUE), # already 0-based
    b = sample(1:3, n, replace = TRUE), # shifted (min = 1)
    c = sample(c(0L, 1L, 3L), n, replace = TRUE) # non-contiguous (gap at 2)
  )
  fit = bgm(
    x = x, variable_type = "ordinal",
    iter = 300, warmup = 300, chains = 1, cores = 1,
    update_method = "adaptive-metropolis",
    display_progress = "none", seed = 42
  )

  sim = simulate(fit, nsim = 50, method = "posterior-mean", seed = 7)

  # Every simulated value is one the model was trained on (original scale),
  # not an internal 0-based code.
  for(v in seq_len(ncol(x))) {
    expect_true(
      all(sim[, v] %in% sort(unique(x[, v]))),
      info = sprintf("variable %d off the original scale", v)
    )
  }

  # The round trip must not warn and must produce no NA prediction cells.
  expect_no_warning(
    probs <- predict(fit, newdata = sim, type = "probabilities")
  )
  expect_false(any(vapply(probs, anyNA, logical(1))))
})

test_that("simulate() leaves Blume-Capel scores unmapped for the round trip", {
  skip_on_cran()

  set.seed(2)
  n = 200
  # Blume-Capel carries no recode map, so its simulated scores stay on the
  # sampler scale that predict() reads directly (no back-mapping applied).
  x = cbind(
    bc = sample(0:3, n, replace = TRUE),
    o = sample(0:2, n, replace = TRUE)
  )
  fit = bgm(
    x = x, variable_type = c("blume-capel", "ordinal"),
    baseline_category = c(0L, 0L),
    iter = 300, warmup = 300, chains = 1, cores = 1,
    update_method = "adaptive-metropolis",
    display_progress = "none", seed = 43
  )

  sim = simulate(fit, nsim = 50, method = "posterior-mean", seed = 8)
  expect_no_warning(predict(fit, newdata = sim, type = "probabilities"))
})

test_that("mixed predict recodes discrete newdata to the original scale", {
  skip_on_cran()

  set.seed(3)
  n = 400
  d_raw = sample(1:3, n, replace = TRUE) # discrete on {1,2,3}
  cc = rnorm(n)
  x_shift = cbind(d = d_raw, c = cc)
  x_base = cbind(d = d_raw - 1L, c = cc) # same data recoded to {0,1,2}

  fit_it = function(x) {
    bgm(
      x = x, variable_type = c("ordinal", "continuous"),
      iter = 400, warmup = 400, chains = 1, cores = 1,
      update_method = "adaptive-metropolis",
      display_progress = "none", seed = 99
    )
  }
  fit_shift = fit_it(x_shift)
  fit_base = fit_it(x_base)

  # Relabeling {1,2,3} -> {0,1,2} leaves the internal model identical, so
  # predictions on corresponding original-scale newdata must match. Before the
  # fix the mixed path fed {1,2,3} to the sampler as codes, so they did not.
  probs_shift = predict(fit_shift, newdata = x_shift[1:20, ], type = "probabilities")
  probs_base = predict(fit_base, newdata = x_base[1:20, ], type = "probabilities")
  for(nm in names(probs_shift)) {
    expect_equal(probs_shift[[nm]], probs_base[[nm]],
      info = sprintf("relabel invariance %s", nm)
    )
  }

  # simulate() returns the original discrete scale and the round trip is clean.
  sim = simulate(fit_shift, nsim = 40, method = "posterior-mean", seed = 5)
  expect_true(all(sim[, "d"] %in% c(1, 2, 3)))
  expect_no_warning(predict(fit_shift, newdata = sim, type = "probabilities"))
})


# ==============================================================================
# 8. Sparse category codings (min > 0 AND gaps)
# ==============================================================================
# Regression for the recode path that CRAN 0.1.6.3 got wrong. There,
# recode_data_for_prediction() shifted ordinal newdata by its per-column
# minimum, so a sparse original coding such as {1,2,4,5} was mapped to
# {0,1,3,4} while the fit itself had collapsed those same values to {0,1,2,3}.
# Every value above a gap was silently attributed to the wrong category, and
# the top value landed outside the fitted category range altogether -- with no
# error and no warning. The stored recode map (arguments$category_levels, used
# by recode_data_for_prediction) makes the mapping absolute, so predict() on
# original-scale newdata and the simulate() -> predict() round trip both land
# on the categories the model was fitted on.
# ==============================================================================

test_that("sparse category codings recode to the fitted categories", {
  skip_on_cran()

  set.seed(20260802)
  n = 400
  p = 3
  latent = matrix(sample(0:3, n * p, replace = TRUE), n, p)
  sparse = c(1, 2, 4, 5) # gap at 3, and min > 0
  x_sparse = matrix(sparse[latent + 1L], n, p, dimnames = list(NULL, c("a", "b", "c")))
  x_dense = matrix(latent, n, p, dimnames = list(NULL, c("a", "b", "c")))

  fit_it = function(x) {
    bgm(
      x = x, variable_type = "ordinal",
      iter = 300, warmup = 300, chains = 1, cores = 1,
      update_method = "adaptive-metropolis",
      display_progress = "none", seed = 11
    )
  }
  fit_sparse = fit_it(x_sparse)
  fit_dense = fit_it(x_dense)

  args = extract_arguments(fit_sparse)

  # The fit stores the observed values, not a range: 4 categories, not 5.
  expect_equal(args$num_categories, rep(3L, p))
  for(v in seq_len(p)) expect_equal(as.numeric(args$category_levels[[v]]), sparse)

  # --- original-scale newdata is recoded through the map, not by min-shift ---
  newdata = matrix(rep(sparse, times = p), nrow = 4, ncol = p,
    dimnames = list(NULL, colnames(x_sparse)))
  recoded = bgms:::recode_data_for_prediction(
    newdata, args$num_categories, rep(TRUE, p),
    category_levels = args$category_levels,
    blume_capel_shift = args$blume_capel_shift
  )
  expect_equal(as.numeric(recoded), rep(0:3, times = p))
  # ... and specifically NOT the legacy answer, which mapped 4 -> 3 and 5 -> 4.
  expect_false(identical(as.numeric(recoded), as.numeric(newdata - 1)))

  # End-to-end: relabeling {1,2,4,5} -> {0,1,2,3} leaves the internal model
  # identical, so predictions on corresponding newdata must agree.
  probs_sparse = predict(fit_sparse, newdata = x_sparse[1:20, ], type = "probabilities")
  probs_dense = predict(fit_dense, newdata = x_dense[1:20, ], type = "probabilities")
  for(nm in names(probs_sparse)) {
    expect_equal(probs_sparse[[nm]], probs_dense[[nm]],
      info = sprintf("sparse/dense relabel invariance %s", nm)
    )
  }

  # --- simulate() -> predict() round trip stays on the sparse original scale ---
  sim = simulate(fit_sparse, nsim = 50, method = "posterior-mean", seed = 12)
  expect_true(all(sim %in% sparse))
  sim_recoded = bgms:::recode_data_for_prediction(
    sim, args$num_categories, rep(TRUE, p),
    category_levels = args$category_levels,
    blume_capel_shift = args$blume_capel_shift
  )
  expect_true(all(sim_recoded %in% 0:3))
  expect_no_warning(
    probs <- predict(fit_sparse, newdata = sim, type = "probabilities")
  )
  expect_false(any(vapply(probs, anyNA, logical(1))))
})


# ==============================================================================
# 9. simulate.bgmCompare() draws from the group it was asked for (F-073)
# ==============================================================================
# The only group-difference test simulate.bgmCompare() had described itself as
# soft. This one is numeric: data simulated from a group, scored back through
# predict() with the SAME group's parameters, must reproduce its own category
# margins -- E[1{X_v = c}] = E[P(X_v = c | X_-v)] holds for any fit that puts
# the same parameters into both paths, and fails if one of them scales,
# shifts, or selects the group differently from the other.
# ==============================================================================

test_that("simulated margins match predicted margins for the same group", {
  skip_on_cran()
  # Groups that differ a lot, so a group mix-up is not a rounding question.
  p = 4
  symmetric = function(values) {
    m = matrix(0, p, p)
    m[upper.tri(m)] = values
    m + t(m)
  }
  omega_1 = symmetric(c(0.6, 0.0, 0.1, 0.0, 0.5, 0.1))
  omega_2 = symmetric(c(0.0, 0.0, 0.6, 0.0, 0.0, 0.5))
  main = matrix(c(0, 0), nrow = p, ncol = 2, byrow = TRUE)

  draw = function(omega, seed) {
    simulate_mrf(
      500, p, num_categories = 2, pairwise = omega, main = main,
      variable_type = "ordinal", iter = 50, seed = seed
    )
  }
  x = rbind(draw(omega_1, 6), draw(omega_2, 106))
  # Full support WITHIN each group, not merely pooled: a category one group
  # never uses has its threshold set by the prior there, and the fit says so in
  # a warning. This block is about the numeric convention, so the data is chosen
  # to keep that question out of it rather than to silence the warning.
  for(g in 1:2) {
    rows = seq_len(500) + (g - 1L) * 500L
    expect_true(all(apply(x[rows, ], 2, function(z) length(unique(z))) == 3L),
      info = sprintf("group %d does not use every category", g))
  }

  fit = bgmCompare(
    x, group = rep(1:2, each = 500),
    iter = 300, warmup = 200, chains = 1, seed = 42,
    difference_selection = FALSE, display_progress = "none"
  )

  # simulate() returns the original scale and predict() returns one column per
  # INTERNAL category, and a sparsely observed category is collapsed into its
  # neighbour, so the two are not index-for-index the same set. The fit's own
  # level map says which original values a predicted column stands for; using
  # it keeps this a test of the numeric convention rather than of the coding.
  levels_of = extract_arguments(fit)$category_levels
  originals_for = function(v, k) {
    map = levels_of[[v]]
    as.numeric(names(map)[map == (k - 1L)])
  }

  nsim = 2000L
  for(g in 1:2) {
    ctx = sprintf("group %d", g)
    sim = simulate(fit, nsim = nsim, seed = 100L + g, group = g, iter = 500)
    probabilities = predict(fit, newdata = sim, group = g,
      type = "probabilities")

    residual = function(scored) {
      max(vapply(seq_len(p), function(v) {
        max(vapply(seq_len(ncol(scored[[v]])), function(k) {
          abs(mean(sim[, v] %in% originals_for(v, k)) - mean(scored[[v]][, k]))
        }, 0.0))
      }, 0.0))
    }

    for(v in seq_len(p)) {
      predicted = probabilities[[v]]
      for(k in seq_len(ncol(predicted))) {
        observed_rate = mean(sim[, v] %in% originals_for(v, k))
        predicted_rate = mean(predicted[, k])
        # The two estimate the same probability. The simulated margin is the
        # noisy one -- nsim near-independent draws -- so the band is its
        # binomial Monte-Carlo error, sqrt(pi (1 - pi) / nsim), at 4 standard
        # errors, with a floor of 0.01 so a category near 0 or 1 does not get
        # a vanishing band from its own vanishing variance. The predicted
        # margin averages conditional probabilities over the same rows and is
        # far tighter, so it contributes little to the spread.
        mc_error = sqrt(predicted_rate * (1 - predicted_rate) / nsim)
        expect_lt(
          abs(observed_rate - predicted_rate), max(4 * mc_error, 0.01),
          label = sprintf("%s, variable %d, category %d", ctx, v, k)
        )
      }
    }

    # And the agreement is a statement about THIS group: scoring the same
    # simulated data with the other group's parameters has to miss, or the
    # test above would pass for a fit that ignored `group` entirely.
    other = if(g == 1L) 2L else 1L
    cross = predict(fit, newdata = sim, group = other, type = "probabilities")
    expect_gt(residual(cross), 10 * residual(probabilities), label = ctx)
  }
})
