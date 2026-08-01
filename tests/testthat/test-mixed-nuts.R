# --------------------------------------------------------------------------- #
# Phase M.2 <U+2014> NUTS vs MH posterior agreement for Mixed MRF.
#
# Verifies that NUTS and adaptive-MH target the same pseudo-posterior by
# comparing posterior summaries across six conditions:
#
#   A  conditional PL, no edge selection, grouped variable order
#   B  marginal PL,    no edge selection, grouped variable order
#   C  conditional PL, edge selection,    grouped variable order
#   D  marginal PL,    edge selection,    grouped variable order
#   E  conditional PL, edge selection,    interleaved variable order
#   F  conditional PL, no edge selection, interleaved variable order
#
# Gated behind BGMS_RUN_SLOW_TESTS because each condition fits two
# chains x 5000 iterations for both samplers.
# --------------------------------------------------------------------------- #


# ---- Skip gate ---------------------------------------------------------------

skip_unless_slow_mixed = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run mixed NUTS vs MH tests"
  )
}


# ---- Data generator ----------------------------------------------------------

# Generate mixed MRF data with known truth (p=3 ordinal, q=2 continuous).
# Returns a data.frame with columns d1, d2, d3, c1, c2.
generate_mixed_test_data = function(seed = 2026, n = 200) {
  p = 3L
  q = 2L
  nc = c(2L, 2L, 2L)

  # True parameters in internal (block) order
  pairwise_disc = matrix(c(
    0, -0.2, 0.15,
    -0.2, 0, 0.0,
    0.15, 0.0, 0
  ), p, p, byrow = TRUE)

  pairwise_cross = matrix(c(
    0.25, 0.0,
    0.0, 0.2,
    -0.15, 0.1
  ), p, q, byrow = TRUE)

  pairwise_cont = matrix(c(
    -0.5, 0.0,
    0.0, -0.5
  ), q, q, byrow = TRUE)

  mux = matrix(0, p, max(nc) + 1)
  mux[1, 1:2] = c(0.0, 0.2)
  mux[2, 1:2] = c(-0.3, -0.1)
  mux[3, 1:2] = c(0.1, -0.15)

  muy = c(0.3, -0.2)

  sim = bgms:::sample_mixed_mrf_gibbs(
    num_states = as.integer(n),
    pairwise_disc_r = pairwise_disc,
    pairwise_cross_r = pairwise_cross,
    pairwise_cont_r = pairwise_cont,
    mux_r = mux, muy_r = muy,
    num_categories_r = nc,
    variable_type_r = rep("ordinal", p),
    baseline_category_r = rep(0L, p),
    iter = 2000L, seed = as.integer(seed)
  )

  df = data.frame(sim$x, sim$y)
  colnames(df) = c(paste0("d", seq_len(p)), paste0("c", seq_len(q)))
  df
}


# ---- Comparison helpers ------------------------------------------------------

ci_overlap = function(x, y) {
  ci_x = quantile(x, c(0.025, 0.975))
  ci_y = quantile(y, c(0.025, 0.975))
  overlap_lo = max(ci_x[1], ci_y[1])
  overlap_hi = min(ci_x[2], ci_y[2])
  if(overlap_lo >= overlap_hi) {
    return(0)
  }
  overlap_width = overlap_hi - overlap_lo
  avg_width = (diff(ci_x) + diff(ci_y)) / 2
  overlap_width / avg_width
}


compare_pairwise_samples = function(fit_nuts, fit_mh, label) {
  nuts_pw = do.call(rbind, fit_nuts$raw_samples$pairwise)
  mh_pw = do.call(rbind, fit_mh$raw_samples$pairwise)

  for(j in seq_len(ncol(nuts_pw))) {
    nuts_mean = mean(nuts_pw[, j])
    mh_mean = mean(mh_pw[, j])
    pooled_sd = sqrt((var(nuts_pw[, j]) + var(mh_pw[, j])) / 2)
    if(pooled_sd > 1e-10) {
      expect_lt(
        abs(nuts_mean - mh_mean) / pooled_sd, 3,
        label = paste0(label, " pairwise[", j, "] z < 3")
      )
    }

    ratio = sd(nuts_pw[, j]) / sd(mh_pw[, j])
    expect_gt(ratio, 0.6, label = paste0(label, " pairwise[", j, "] SD ratio > 0.6"))
    expect_lt(ratio, 1.7, label = paste0(label, " pairwise[", j, "] SD ratio < 1.7"))

    ov = ci_overlap(nuts_pw[, j], mh_pw[, j])
    expect_gt(ov, 0.7,
      label = paste0(label, " pairwise[", j, "] CI overlap > 0.7")
    )
  }
}


compare_main_samples = function(fit_nuts, fit_mh, label) {
  nuts_main = do.call(rbind, fit_nuts$raw_samples$main)
  mh_main = do.call(rbind, fit_mh$raw_samples$main)

  for(j in seq_len(ncol(nuts_main))) {
    nuts_mean = mean(nuts_main[, j])
    mh_mean = mean(mh_main[, j])
    pooled_sd = sqrt((var(nuts_main[, j]) + var(mh_main[, j])) / 2)
    if(pooled_sd > 1e-10) {
      expect_lt(
        abs(nuts_mean - mh_mean) / pooled_sd, 3,
        label = paste0(label, " main[", j, "] z < 3")
      )
    }
  }
}


# ---- Sampler configuration --------------------------------------------------

n_iter = 5000
n_warmup = 2000
n_chains = 2
pw_scale = 2.5
main_a = 0.5
main_b = 0.5


# ---- Condition A: conditional PL, no ES, grouped ----------------------------

test_that("M.2A: NUTS vs MH agree (conditional PL, no ES, grouped)", {
  skip_unless_slow_mixed()

  dat = generate_mixed_test_data(seed = 2026)
  vtype = c(rep("ordinal", 3), rep("continuous", 2))

  fit_nuts = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = FALSE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 101
  )

  fit_mh = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = FALSE, update_method = "adaptive-metropolis",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 201
  )

  compare_pairwise_samples(fit_nuts, fit_mh, "M.2A")
  compare_main_samples(fit_nuts, fit_mh, "M.2A")

  diff = abs(fit_nuts$posterior_mean_pairwise - fit_mh$posterior_mean_pairwise)
  expect_lt(max(diff), 0.1, label = "M.2A max abs diff in associations < 0.1")

  expect_equal(fit_nuts$nuts_diag$summary$total_divergences, 0)
})


# ---- Condition B: marginal PL, no ES, grouped -------------------------------

test_that("M.2B: NUTS vs MH agree (marginal PL, no ES, grouped)", {
  skip_unless_slow_mixed()

  dat = generate_mixed_test_data(seed = 2027)
  vtype = c(rep("ordinal", 3), rep("continuous", 2))

  fit_nuts = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = FALSE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 102
  )

  fit_mh = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = FALSE, update_method = "adaptive-metropolis",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 202
  )

  compare_pairwise_samples(fit_nuts, fit_mh, "M.2B")
  compare_main_samples(fit_nuts, fit_mh, "M.2B")

  expect_equal(fit_nuts$nuts_diag$summary$total_divergences, 0)
})


# ---- Condition C: conditional PL, ES, grouped -------------------------------

test_that("M.2C: NUTS vs MH agree (conditional PL, ES, grouped)", {
  skip_unless_slow_mixed()

  dat = generate_mixed_test_data(seed = 2028, n = 300)
  vtype = c(rep("ordinal", 3), rep("continuous", 2))

  fit_nuts = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = TRUE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 103
  )

  fit_mh = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = TRUE, update_method = "adaptive-metropolis",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 203
  )

  # PIP agreement
  pip_nuts = fit_nuts$posterior_mean_indicator
  pip_mh = fit_mh$posterior_mean_indicator
  pip_nuts_ut = pip_nuts[upper.tri(pip_nuts)]
  pip_mh_ut = pip_mh[upper.tri(pip_mh)]

  expect_lt(max(abs(pip_nuts_ut - pip_mh_ut)), 0.20,
    label = "M.2C max PIP diff < 0.20"
  )
  expect_gt(cor(pip_nuts_ut, pip_mh_ut), 0.95,
    label = "M.2C PIP correlation > 0.95"
  )

  expect_equal(fit_nuts$nuts_diag$summary$total_divergences, 0)
})


# ---- Condition D: marginal PL, ES, grouped ----------------------------------

test_that("M.2D: NUTS vs MH agree (marginal PL, ES, grouped)", {
  skip_unless_slow_mixed()

  dat = generate_mixed_test_data(seed = 2029, n = 300)
  vtype = c(rep("ordinal", 3), rep("continuous", 2))

  fit_nuts = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = TRUE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 104
  )

  fit_mh = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = TRUE, update_method = "adaptive-metropolis",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 204
  )

  pip_nuts = fit_nuts$posterior_mean_indicator
  pip_mh = fit_mh$posterior_mean_indicator
  pip_nuts_ut = pip_nuts[upper.tri(pip_nuts)]
  pip_mh_ut = pip_mh[upper.tri(pip_mh)]

  expect_lt(max(abs(pip_nuts_ut - pip_mh_ut)), 0.20,
    label = "M.2D max PIP diff < 0.20"
  )
  expect_gt(cor(pip_nuts_ut, pip_mh_ut), 0.95,
    label = "M.2D PIP correlation > 0.95"
  )

  expect_equal(fit_nuts$nuts_diag$summary$total_divergences, 0)
})


# ---- Condition E: conditional PL, ES, interleaved ----------------------------

test_that("M.2E: NUTS vs MH agree (conditional PL, ES, interleaved)", {
  skip_unless_slow_mixed()

  # Generate grouped data, then reorder to interleaved: d1, c1, d2, c2, d3
  dat_grouped = generate_mixed_test_data(seed = 2030, n = 300)
  dat = dat_grouped[, c("d1", "c1", "d2", "c2", "d3")]
  vtype = c("ordinal", "continuous", "ordinal", "continuous", "ordinal")

  fit_nuts = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = TRUE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 105
  )

  fit_mh = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = TRUE, update_method = "adaptive-metropolis",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 205
  )

  # PIP agreement
  pip_nuts = fit_nuts$posterior_mean_indicator
  pip_mh = fit_mh$posterior_mean_indicator
  pip_nuts_ut = pip_nuts[upper.tri(pip_nuts)]
  pip_mh_ut = pip_mh[upper.tri(pip_mh)]

  expect_lt(max(abs(pip_nuts_ut - pip_mh_ut)), 0.20,
    label = "M.2E max PIP diff < 0.20"
  )

  # Output ordering: names must be in user (interleaved) order
  expect_identical(rownames(pip_nuts), colnames(dat))
  expect_identical(colnames(pip_nuts), colnames(dat))
  expect_identical(rownames(pip_mh), colnames(dat))

  # Extractor consistency: indicator-weighted raw sample means must
  # match posterior_mean_pairwise at every position.
  pw_nuts = extract_pairwise_interactions(fit_nuts)
  ind_nuts = extract_indicators(fit_nuts)
  weighted_nuts = colMeans(pw_nuts * ind_nuts)
  expect_true(
    all(check_extractor_matrix_consistency(
      weighted_nuts, fit_nuts$posterior_mean_pairwise
    )),
    info = "M.2E NUTS extractor names match matrix positions"
  )

  pw_mh = extract_pairwise_interactions(fit_mh)
  ind_mh = extract_indicators(fit_mh)
  weighted_mh = colMeans(pw_mh * ind_mh)
  expect_true(
    all(check_extractor_matrix_consistency(
      weighted_mh, fit_mh$posterior_mean_pairwise
    )),
    info = "M.2E MH extractor names match matrix positions"
  )

  # Column name agreement between samplers
  expect_identical(
    colnames(extract_pairwise_interactions(fit_nuts)),
    colnames(extract_pairwise_interactions(fit_mh))
  )

  expect_equal(fit_nuts$nuts_diag$summary$total_divergences, 0)
})


# ---- Condition F: conditional PL, no ES, interleaved -------------------------

test_that("M.2F: NUTS vs MH agree (conditional PL, no ES, interleaved)", {
  skip_unless_slow_mixed()

  dat_grouped = generate_mixed_test_data(seed = 2031)
  dat = dat_grouped[, c("d1", "c1", "d2", "c2", "d3")]
  vtype = c("ordinal", "continuous", "ordinal", "continuous", "ordinal")

  fit_nuts = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = FALSE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 106
  )

  fit_mh = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = FALSE, update_method = "adaptive-metropolis",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 206
  )

  compare_pairwise_samples(fit_nuts, fit_mh, "M.2F")
  compare_main_samples(fit_nuts, fit_mh, "M.2F")

  # Output ordering checks
  assoc_nuts = fit_nuts$posterior_mean_pairwise
  expect_identical(rownames(assoc_nuts), colnames(dat))
  expect_identical(colnames(assoc_nuts), colnames(dat))

  pw_nuts_means = colMeans(extract_pairwise_interactions(fit_nuts))
  expect_true(
    all(check_extractor_matrix_consistency(
      pw_nuts_means, fit_nuts$posterior_mean_pairwise
    )),
    info = "M.2F NUTS extractor consistency"
  )

  expect_identical(
    colnames(extract_pairwise_interactions(fit_nuts)),
    colnames(extract_pairwise_interactions(fit_mh))
  )

  expect_equal(fit_nuts$nuts_diag$summary$total_divergences, 0)
})


# ---- Main effects comparison (grouped) --------------------------------------

test_that("M.2G: NUTS vs MH main effects agree (conditional PL, grouped)", {
  skip_unless_slow_mixed()

  dat = generate_mixed_test_data(seed = 2032)
  vtype = c(rep("ordinal", 3), rep("continuous", 2))

  fit_nuts = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = FALSE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 107
  )

  fit_mh = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = FALSE, update_method = "adaptive-metropolis",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 207
  )

  # Main effects via extractor
  me_nuts = extract_main_effects(fit_nuts)
  me_mh = extract_main_effects(fit_mh)

  # Discrete thresholds: rownames match
  expect_identical(rownames(me_nuts$discrete), rownames(me_mh$discrete))

  # Continuous means: rownames match
  expect_identical(rownames(me_nuts$continuous), rownames(me_mh$continuous))

  # Discrete threshold values agree
  for(i in seq_len(nrow(me_nuts$discrete))) {
    for(j in seq_len(ncol(me_nuts$discrete))) {
      if(!is.na(me_nuts$discrete[i, j]) && !is.na(me_mh$discrete[i, j])) {
        expect_lt(
          abs(me_nuts$discrete[i, j] - me_mh$discrete[i, j]), 0.3,
          label = paste0("M.2G discrete[", i, ",", j, "] diff < 0.3")
        )
      }
    }
  }

  # Continuous mean values agree
  for(i in seq_len(nrow(me_nuts$continuous))) {
    expect_lt(
      abs(me_nuts$continuous[i, 1] - me_mh$continuous[i, 1]), 0.3,
      label = paste0("M.2G continuous[", i, "] mean diff < 0.3")
    )
  }

  # Raw main samples z-test
  compare_main_samples(fit_nuts, fit_mh, "M.2G")
})


# ---- Downstream methods on NUTS fit -----------------------------------------

test_that("M.2H: coef/summary/simulate/predict work on mixed NUTS fit", {
  skip_unless_slow_mixed()

  dat = generate_mixed_test_data(seed = 2033, n = 150)
  vtype = c(rep("ordinal", 3), rep("continuous", 2))

  fit = bgm(
    dat,
    variable_type = vtype,
    iter = 2000, warmup = 1000, chains = 1,
    edge_selection = FALSE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 301
  )

  # coef
  coefs = coef(fit)
  expect_true(is.list(coefs))
  expect_true(!is.null(coefs$pairwise))
  expect_true(!is.null(coefs$main))
  expect_equal(nrow(coefs$pairwise), 5)
  expect_equal(ncol(coefs$pairwise), 5)

  # summary
  summ = summary(fit)
  expect_true(!is.null(summ))

  # simulate
  sim = simulate(fit, nsim = 10, seed = 1)
  expect_equal(ncol(sim), 5)
  expect_equal(nrow(sim), 10)
  expect_identical(colnames(sim), colnames(dat))

  # predict
  pred = predict(fit, newdata = dat[1:5, ])
  expect_true(is.list(pred))
  expect_equal(length(pred), 5)
  expect_identical(names(pred), colnames(dat))
})


# ---- Downstream methods on interleaved NUTS fit ------------------------------

test_that("M.2I: simulate/predict preserve interleaved column order", {
  skip_unless_slow_mixed()

  dat_grouped = generate_mixed_test_data(seed = 2034, n = 150)
  dat = dat_grouped[, c("d1", "c1", "d2", "c2", "d3")]
  vtype = c("ordinal", "continuous", "ordinal", "continuous", "ordinal")

  fit = bgm(
    dat,
    variable_type = vtype,
    iter = 2000, warmup = 1000, chains = 1,
    edge_selection = FALSE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    display_progress = "none", seed = 302
  )

  sim = simulate(fit, nsim = 10, seed = 2)
  expect_identical(colnames(sim), colnames(dat))

  pred = predict(fit, newdata = dat[1:5, ])
  expect_true(is.list(pred))
  expect_equal(length(pred), 5)
  expect_identical(names(pred), colnames(dat))
})


# ---- NUTS diagnostics -------------------------------------------------------

test_that("M.2J: NUTS diagnostics are clean for mixed MRF", {
  skip_unless_slow_mixed()

  dat = generate_mixed_test_data(seed = 2035)
  vtype = c(rep("ordinal", 3), rep("continuous", 2))

  fit = bgm(
    dat,
    variable_type = vtype,
    iter = 2000, warmup = 1000, chains = 2,
    edge_selection = FALSE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    delta = 0,
    display_progress = "none", seed = 303
  )

  expect_true(!is.null(fit$nuts_diag))
  expect_equal(fit$nuts_diag$summary$total_divergences, 0)
  expect_equal(fit$nuts_diag$summary$max_tree_depth_hits, 0)
  expect_gt(fit$nuts_diag$summary$min_ebfmi, 0.3)
})


# ---- Condition T: determinant tilt, conditional PL, no ES, grouped ----------

test_that("M.2T: NUTS vs MH agree under tilt (conditional PL, no ES, delta=1)", {
  skip_unless_slow_mixed()

  dat = generate_mixed_test_data(seed = 2036)
  vtype = c(rep("ordinal", 3), rep("continuous", 2))

  fit_nuts = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = FALSE, update_method = "nuts",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    delta = 1,
    display_progress = "none", seed = 104
  )

  fit_mh = bgm(
    dat,
    variable_type = vtype,
    iter = n_iter, warmup = n_warmup, chains = n_chains,
    edge_selection = FALSE, update_method = "adaptive-metropolis",
    pairwise_scale = pw_scale, main_alpha = main_a, main_beta = main_b,
    delta = 1,
    display_progress = "none", seed = 204
  )

  compare_pairwise_samples(fit_nuts, fit_mh, "M.2T")
  compare_main_samples(fit_nuts, fit_mh, "M.2T")

  diff = abs(fit_nuts$posterior_mean_pairwise - fit_mh$posterior_mean_pairwise)
  expect_lt(max(diff), 0.1,
    label = "M.2T max abs diff in associations < 0.1 (delta=1)"
  )

  expect_equal(fit_nuts$nuts_diag$summary$total_divergences, 0)
})
