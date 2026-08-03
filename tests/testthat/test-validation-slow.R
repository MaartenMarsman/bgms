# ==============================================================================
# Gated slow validation tests
# ==============================================================================
#
# Environment-gated wrappers around mixed MRF validation assertions. These run
# in nightly CI (where
# BGMS_RUN_SLOW_TESTS=true) and are skipped during local devtools::test()
# and CRAN checks.
#
# Each test fits a small mixed MRF network with enough iterations to
# produce meaningful posterior estimates, then checks a quantitative
# recovery or agreement criterion.
#
# Shared helper functions live in tests/testthat/helper-validation.R.
# ==============================================================================

# Load shared helpers once per file
helpers_path = file.path("tests", "testthat", "helper-validation.R")
if(!file.exists(helpers_path)) {
  # When running from testthat, the working directory is tests/testthat/
  helpers_path = file.path("helper-validation.R")
}
helpers_available = file.exists(helpers_path)


# ==============================================================================
# Gate: weekly certification tier (T2)
# ==============================================================================
#
# Parameter recovery, MH-vs-NUTS agreement and the estimate-simulate-re-estimate
# cycle are the heavy Monte-Carlo machinery, not a nightly heartbeat, so this
# file runs under BGMS_RUN_CERTIFICATION rather than BGMS_RUN_SLOW_TESTS.

skip_unless_certification_validation = function() {
  skip_unless_certification()
  skip_if(!helpers_available, "Validation helpers not found")
}

# ==============================================================================
# 1. Parameter recovery (adapted from group1)
# ==============================================================================

test_that("mixed MRF parameter recovery: cor > 0.8 (small network)", {
  skip_unless_certification_validation()
  source(helpers_path, local = TRUE)

  net = make_network(p = 2, q = 2, n_cat = c(1L, 2L), density = 1.0, seed = 101)
  dat = generate_data(net, n = 2000, source = "bgms", seed = 201)
  vt = c(rep("ordinal", 2), rep("continuous", 2))

  fit = bgm(dat,
    variable_type = vt,
    edge_selection = FALSE,
    iter = 10000, warmup = 5000, chains = 2, seed = 301,
    display_progress = "none"
  )

  true_blocks = list(
    mux = net$mux, muy = net$muy,
    pairwise_disc = net$pairwise_disc, pairwise_cross = net$pairwise_cross, pairwise_cont = net$pairwise_cont
  )
  est_blocks = extract_bgms_blocks(fit, net)

  true_flat = flatten_params(true_blocks)
  est_flat = flatten_params(est_blocks)

  r = cor(true_flat, est_flat)
  expect_gt(r, 0.8,
    label = sprintf("recovery correlation (%.4f)", r)
  )
})


# ==============================================================================
# 2. Metropolis vs NUTS agreement (adapted from group2)
# ==============================================================================

test_that("MH vs NUTS posterior agreement: cor > 0.95", {
  skip_unless_certification_validation()
  source(helpers_path, local = TRUE)

  net = make_network(p = 2, q = 2, n_cat = c(1L, 2L), density = 1.0, seed = 101)
  dat = generate_data(net, n = 2000, source = "bgms", seed = 201)
  vt = c(rep("ordinal", 2), rep("continuous", 2))

  fit_mh = bgm(dat,
    variable_type = vt,
    update_method = "adaptive-metropolis",
    edge_selection = FALSE,
    iter = 15000, warmup = 10000, chains = 2, seed = 401,
    display_progress = "none"
  )

  fit_nuts = bgm(dat,
    variable_type = vt,
    update_method = "nuts",
    edge_selection = FALSE,
    iter = 5000, warmup = 3000, chains = 2, seed = 402,
    display_progress = "none"
  )

  est_mh = flatten_params(extract_bgms_blocks(fit_mh, net))
  est_nuts = flatten_params(extract_bgms_blocks(fit_nuts, net))

  r = cor(est_mh, est_nuts)
  expect_gt(r, 0.95,
    label = sprintf("MH vs NUTS agreement (%.4f)", r)
  )
})


# ==============================================================================
# 3. Estimate-simulate-re-estimate cycle (mixed MRF)
# ==============================================================================

test_that("estimate-simulate-re-estimate cycle: cor > 0.7 (mixed MRF)", {
  skip_unless_certification_validation()
  source(helpers_path, local = TRUE)

  net = make_network(p = 2, q = 2, n_cat = c(1L, 2L), density = 1.0, seed = 101)
  dat = generate_data(net, n = 2000, source = "bgms", seed = 201)
  vt = c(rep("ordinal", 2), rep("continuous", 2))

  fit1 = bgm(dat,
    variable_type = vt,
    edge_selection = FALSE, iter = 5000, warmup = 2000,
    chains = 1, seed = 601,
    display_progress = "none"
  )

  simulated = simulate(fit1, nsim = 2000, seed = 701)

  fit2 = bgm(simulated,
    variable_type = vt,
    edge_selection = FALSE, iter = 5000, warmup = 2000,
    chains = 1, seed = 801,
    display_progress = "none"
  )

  pw1 = as.vector(fit1$posterior_mean_pairwise)
  pw2 = as.vector(fit2$posterior_mean_pairwise)

  r = cor(pw1, pw2)
  expect_gt(r, 0.7,
    label = sprintf("cycle pairwise correlation (%.4f)", r)
  )
})
