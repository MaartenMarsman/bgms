# ==============================================================================
# Conditional Distribution Correctness Tests
# ==============================================================================
#
# Verify that the C++ conditional prediction functions produce values
# matching hand-computed reference values for minimal networks.
#
# Covers compute_conditional_probs (OMRF), compute_conditional_ggm (GGM),
# and compute_conditional_mixed (mixed MRF).
#
# ==============================================================================


# ==============================================================================
# Test 1: OMRF binary conditional probabilities
# ==============================================================================
# Minimal example: p=2 binary ordinals, n=3.
#
# P(x_i = c | x_{-i}) proportional to:
#   c = 0:  1  (reference)
#   c = 1:  exp(main[i,1] + 1 * rest_i)
#
# rest_i = sum_{k != i} 2 * pairwise[k,i] * x_k
# (association scale: pairwise stores A = sigma/2)

test_that("OMRF binary conditional probabilities match hand computation", {
  n = 3L
  p = 2L

  observations = matrix(
    c(
      0L, 1L, 1L,
      1L, 0L, 1L
    ),
    nrow = n, ncol = p
  )

  pairwise = matrix(
    c(
      0.0, 0.15,
      0.15, 0.0
    ),
    nrow = p, byrow = TRUE
  )
  main = matrix(c(-0.5, 0.2), nrow = p, ncol = 1)
  num_categories = c(1L, 1L)
  variable_type = c("ordinal", "ordinal")
  baseline_category = c(0L, 0L)

  # --- Predict variable 0 (1st variable) ---
  # rest = 2 * pairwise[1,0] * x_2 = 2 * 0.15 * c(1, 0, 1) = c(0.3, 0, 0.3)
  rest_v0 = c(0.3, 0.0, 0.3)
  mu0 = main[1, 1] # -0.5

  # P(x_0 = 0) = 1 / (1 + exp(mu0 + rest))
  # P(x_0 = 1) = exp(mu0 + rest) / (1 + exp(mu0 + rest))
  logit_v0 = mu0 + rest_v0
  prob_v0_cat1 = exp(logit_v0) / (1 + exp(logit_v0))
  prob_v0_cat0 = 1 - prob_v0_cat1

  probs_cpp = compute_conditional_probs(
    observations = observations,
    predict_vars = 0L,
    pairwise = pairwise,
    main = main,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  expect_equal(probs_cpp[[1]][, 1], prob_v0_cat0, tolerance = 1e-10)
  expect_equal(probs_cpp[[1]][, 2], prob_v0_cat1, tolerance = 1e-10)

  # --- Predict variable 1 (2nd variable) ---
  # rest = 2 * pairwise[0,1] * x_1 = 2 * 0.15 * c(0, 1, 1) = c(0, 0.3, 0.3)
  rest_v1 = c(0.0, 0.3, 0.3)
  mu1 = main[2, 1] # 0.2

  logit_v1 = mu1 + rest_v1
  prob_v1_cat1 = exp(logit_v1) / (1 + exp(logit_v1))
  prob_v1_cat0 = 1 - prob_v1_cat1

  probs_cpp2 = compute_conditional_probs(
    observations = observations,
    predict_vars = 1L,
    pairwise = pairwise,
    main = main,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  expect_equal(probs_cpp2[[1]][, 1], prob_v1_cat0, tolerance = 1e-10)
  expect_equal(probs_cpp2[[1]][, 2], prob_v1_cat1, tolerance = 1e-10)

  # Probabilities sum to 1
  expect_equal(rowSums(probs_cpp[[1]]), rep(1, n), tolerance = 1e-10)
  expect_equal(rowSums(probs_cpp2[[1]]), rep(1, n), tolerance = 1e-10)
})


# ==============================================================================
# Test 2: Multi-category ordinal conditional probabilities
# ==============================================================================
# p=1 ordinal with 3 categories (0, 1, 2), n=4.
#
# P(x = c | ...) proportional to:
#   c = 0:  1
#   c = 1:  exp(main[1,1] + 1 * rest)
#   c = 2:  exp(main[1,2] + 2 * rest)

test_that("OMRF multi-category conditional probabilities match hand computation", {
  n = 4L
  p = 2L

  observations = matrix(
    c(
      0L, 1L, 2L, 1L,
      1L, 0L, 2L, 1L
    ),
    nrow = n, ncol = p
  )

  pairwise = matrix(
    c(
      0.0, 0.125,
      0.125, 0.0
    ),
    nrow = p, byrow = TRUE
  )
  # 2 categories (3 levels: 0, 1, 2) => 2 threshold columns
  main = matrix(
    c(
      -0.5, 0.1,
      0.2, -0.3
    ),
    nrow = p, ncol = 2, byrow = TRUE
  )
  num_categories = c(2L, 2L)
  variable_type = c("ordinal", "ordinal")
  baseline_category = c(0L, 0L)

  # Predict variable 0:
  # rest = 2 * pairwise[1,0] * x_2 = 2 * 0.125 * c(1, 0, 2, 1) = c(0.25, 0, 0.5, 0.25)
  rest = 2 * 0.125 * c(1, 0, 2, 1)

  # For each observation, compute unnormalized probabilities
  hand_probs = matrix(NA_real_, nrow = n, ncol = 3)
  for(i in seq_len(n)) {
    log_unnorm = c(
      0, # c = 0 (reference)
      main[1, 1] + 1 * rest[i], # c = 1
      main[1, 2] + 2 * rest[i] # c = 2
    )
    # Stable softmax
    max_val = max(log_unnorm)
    unnorm = exp(log_unnorm - max_val)
    hand_probs[i, ] = unnorm / sum(unnorm)
  }

  probs_cpp = compute_conditional_probs(
    observations = observations,
    predict_vars = 0L,
    pairwise = pairwise,
    main = main,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  expect_equal(as.matrix(probs_cpp[[1]]), hand_probs, tolerance = 1e-10)
  expect_equal(rowSums(probs_cpp[[1]]), rep(1, n), tolerance = 1e-10)
})


# ==============================================================================
# Test 3: GGM conditional mean and sd
# ==============================================================================
# For precision matrix Omega:
#   E[X_j | X_{-j}] = -(1/omega_jj) * sum_{k != j} omega_jk * x_k
#   SD[X_j | X_{-j}] = sqrt(1/omega_jj)
#
# Note: compute_conditional_ggm operates on centered data.

test_that("GGM conditional mean and sd match hand computation", {
  n = 4L
  p = 3L

  # Symmetric positive definite precision matrix
  omega = matrix(c(
    2.0, -0.5, 0.3,
    -0.5, 1.5, -0.2,
    0.3, -0.2, 1.0
  ), nrow = p, byrow = TRUE)

  # Centered observations
  x = matrix(c(
    1.0, -0.5, 0.3,
    0.2, 0.8, -0.4,
    -0.7, 0.1, 0.9,
    0.5, -0.3, 0.2
  ), nrow = n, byrow = TRUE)

  # Predict variable 0:
  # mean = -(1/2) * ((-0.5) * x[,2] + 0.3 * x[,3])
  hand_mean_v0 = -(1 / omega[1, 1]) *
    (omega[1, 2] * x[, 2] + omega[1, 3] * x[, 3])
  hand_sd_v0 = sqrt(1 / omega[1, 1])

  result = compute_conditional_ggm(
    observations = x,
    predict_vars = 0L,
    precision = omega
  )

  expect_equal(result[[1]][, 1], hand_mean_v0, tolerance = 1e-10)
  expect_equal(result[[1]][1, 2], hand_sd_v0, tolerance = 1e-10)

  # Predict variable 1:
  hand_mean_v1 = -(1 / omega[2, 2]) *
    (omega[2, 1] * x[, 1] + omega[2, 3] * x[, 3])
  hand_sd_v1 = sqrt(1 / omega[2, 2])

  result2 = compute_conditional_ggm(
    observations = x,
    predict_vars = 1L,
    precision = omega
  )

  expect_equal(result2[[1]][, 1], hand_mean_v1, tolerance = 1e-10)
  expect_equal(result2[[1]][1, 2], hand_sd_v1, tolerance = 1e-10)

  # Predict all variables at once
  result_all = compute_conditional_ggm(
    observations = x,
    predict_vars = c(0L, 1L, 2L),
    precision = omega
  )

  expect_equal(length(result_all), p)
  expect_equal(result_all[[1]][, 1], hand_mean_v0, tolerance = 1e-10)
  expect_equal(result_all[[2]][, 1], hand_mean_v1, tolerance = 1e-10)
})


# ==============================================================================
# Test 4: Mixed MRF discrete conditional probabilities
# ==============================================================================
# p=2 binary ordinals, q=1 continuous, n=3.
#
# For discrete variable s, the rest scores in the mixed MRF are:
#   rest = sum_{k != s} 2 * (x_k - ref_k) * pairwise_disc[k,s]
#        + sum_j 2 * pairwise_cross[s,j] * y_j
# (association scale: pairwise_disc stores A = sigma/2)
#
# Then P(x_s = c | rest) follows the same softmax as the pure OMRF.

test_that("mixed MRF discrete conditional probabilities match hand computation", {
  n = 3L
  p = 2L
  q = 1L

  x_obs = matrix(
    c(
      0L, 1L, 1L,
      1L, 0L, 1L
    ),
    nrow = n, ncol = p
  )
  y_obs = matrix(c(0.5, -0.3, 1.2), nrow = n, ncol = q)

  pairwise_disc = matrix(c(
    0.0, 0.15,
    0.15, 0.0
  ), nrow = p, byrow = TRUE)
  pairwise_cross = matrix(c(0.2, 0.4), nrow = p, ncol = q)
  pairwise_cont = matrix(-1.0, nrow = q, ncol = q)
  mux = matrix(c(-0.5, 0.2), nrow = p, ncol = 1)
  muy = c(0.1)

  num_categories = c(1L, 1L)
  variable_type = c("ordinal", "ordinal")
  baseline_category = c(0L, 0L)

  # --- Predict discrete variable 0 ---
  # rest_discrete = 2 * pairwise_disc[1,0] * (x_2 - 0) = 2 * 0.15 * c(1, 0, 1)
  # rest_continuous = 2 * pairwise_cross[0,0] * y = 2 * 0.2 * c(0.5, -0.3, 1.2)
  # rest = c(0.3, 0, 0.3) + c(0.2, -0.12, 0.48) = c(0.5, -0.12, 0.78)
  rest_v0 = c(0.5, -0.12, 0.78)
  mu0 = mux[1, 1] # -0.5

  logit_v0 = mu0 + rest_v0
  hand_prob_cat1 = exp(logit_v0) / (1 + exp(logit_v0))
  hand_prob_cat0 = 1 - hand_prob_cat1

  result = compute_conditional_mixed(
    x_observations = x_obs,
    y_observations = y_obs,
    predict_vars = 0L,
    pairwise_disc = pairwise_disc,
    pairwise_cross = pairwise_cross,
    pairwise_cont = pairwise_cont,
    mux = mux,
    muy = muy,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  expect_equal(result[[1]][, 1], hand_prob_cat0, tolerance = 1e-10)
  expect_equal(result[[1]][, 2], hand_prob_cat1, tolerance = 1e-10)
  expect_equal(rowSums(result[[1]]), rep(1, n), tolerance = 1e-10)

  # --- Predict discrete variable 1 ---
  # rest_discrete = 2 * pairwise_disc[0,1] * (x_1 - 0) = 2 * 0.15 * c(0, 1, 1)
  # rest_continuous = 2 * pairwise_cross[1,0] * y = 2 * 0.4 * c(0.5, -0.3, 1.2)
  # rest = c(0, 0.3, 0.3) + c(0.4, -0.24, 0.96) = c(0.4, 0.06, 1.26)
  rest_v1 = c(0.4, 0.06, 1.26)
  mu1 = mux[2, 1] # 0.2

  logit_v1 = mu1 + rest_v1
  hand_prob1_cat1 = exp(logit_v1) / (1 + exp(logit_v1))
  hand_prob1_cat0 = 1 - hand_prob1_cat1

  result2 = compute_conditional_mixed(
    x_observations = x_obs,
    y_observations = y_obs,
    predict_vars = 1L,
    pairwise_disc = pairwise_disc,
    pairwise_cross = pairwise_cross,
    pairwise_cont = pairwise_cont,
    mux = mux,
    muy = muy,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  expect_equal(result2[[1]][, 1], hand_prob1_cat0, tolerance = 1e-10)
  expect_equal(result2[[1]][, 2], hand_prob1_cat1, tolerance = 1e-10)
})


# ==============================================================================
# Test 5: Mixed MRF continuous conditional mean and sd
# ==============================================================================
# pairwise_cont is association-scale (negative-definite), precision = -2 * pairwise_cont.
# For continuous variable j in the mixed MRF:
#   theta_jj = -2 * pairwise_cont[j,j]
#   cond_var = 1 / theta_jj
#   cond_mean = muy_j + cond_var * (
#     sum_{k != j} 2 * pairwise_cont[j,k] * (y_k - muy_k)
#     + sum_s 2 * pairwise_cross[s,j] * (x_s - ref_s)
#   )
#   cond_sd = sqrt(cond_var)

test_that("mixed MRF continuous conditional matches hand computation", {
  n = 3L
  p = 2L
  q = 1L

  x_obs = matrix(
    c(
      0L, 1L, 1L,
      1L, 0L, 1L
    ),
    nrow = n, ncol = p
  )
  y_obs = matrix(c(0.5, -0.3, 1.2), nrow = n, ncol = q)

  pairwise_disc = matrix(c(
    0.0, 0.15,
    0.15, 0.0
  ), nrow = p, byrow = TRUE)
  pairwise_cross = matrix(c(0.2, 0.4), nrow = p, ncol = q)
  pairwise_cont = matrix(-1.0, nrow = q, ncol = q)
  mux = matrix(c(-0.5, 0.2), nrow = p, ncol = 1)
  muy = c(0.1)

  num_categories = c(1L, 1L)
  variable_type = c("ordinal", "ordinal")
  baseline_category = c(0L, 0L)

  # Predict continuous variable (index = p = 2, since 0-based: [0,1] are discrete)
  # theta_jj = -2 * pairwise_cont[1,1] = -2 * (-1.0) = 2.0
  # cond_var = 1/theta_jj = 0.5
  # lp_continuous = 0 (only 1 continuous variable, skip self)
  # lp_discrete = 2 * pairwise_cross[0,0] * (x_1 - 0) + 2 * pairwise_cross[1,0] * (x_2 - 0)
  #             = 2 * 0.2 * x_1 + 2 * 0.4 * x_2
  # For n=1: 0.4*0 + 0.8*1 = 0.8
  # For n=2: 0.4*1 + 0.8*0 = 0.4
  # For n=3: 0.4*1 + 0.8*1 = 1.2
  lp_discrete = c(0.8, 0.4, 1.2)

  theta_jj = -2 * pairwise_cont[1, 1]
  cond_var = 1 / theta_jj
  hand_mean = muy + cond_var * lp_discrete
  # = 0.1 + 0.5 * c(0.8, 0.4, 1.2) = c(0.5, 0.3, 0.7)
  hand_sd = sqrt(cond_var) # sqrt(0.5)

  result = compute_conditional_mixed(
    x_observations = x_obs,
    y_observations = y_obs,
    predict_vars = as.integer(p), # index 2 = first continuous variable
    pairwise_disc = pairwise_disc,
    pairwise_cross = pairwise_cross,
    pairwise_cont = pairwise_cont,
    mux = mux,
    muy = muy,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  expect_equal(result[[1]][, 1], hand_mean, tolerance = 1e-10)
  expect_equal(result[[1]][1, 2], hand_sd, tolerance = 1e-10)
})


# ==============================================================================
# Test 6: Mixed MRF with multiple continuous variables
# ==============================================================================
# p=1 binary ordinal, q=2 continuous, n=3.
# Tests cross-variable precision terms in the continuous conditional.

test_that("mixed MRF continuous conditional works for q > 1", {
  n = 3L
  p = 1L
  q = 2L

  x_obs = matrix(c(0L, 1L, 1L), nrow = n, ncol = p)
  y_obs = matrix(c(
    0.5, 0.2,
    -0.3, 0.8,
    1.2, -0.5
  ), nrow = n, ncol = q, byrow = TRUE)

  pairwise_disc = matrix(0.0, nrow = p, ncol = p)
  pairwise_cross = matrix(c(0.2, 0.1), nrow = p, ncol = q)
  pairwise_cont = matrix(c(-1.0, -0.15, -0.15, -0.75), nrow = q, byrow = TRUE)
  mux = matrix(-0.5, nrow = p, ncol = 1)
  muy = c(0.1, -0.2)

  num_categories = c(1L)
  variable_type = c("ordinal")
  baseline_category = c(0L)

  # --- Predict continuous variable 0 (internal index p = 1) ---
  # theta_jj = -2 * pairwise_cont[1,1] = -2 * (-1.0) = 2.0
  # cond_var_0 = 1/theta_jj = 0.5
  # lp_continuous = 2 * pairwise_cont[1,2] * (y[,2] - muy[2])
  #              = 2 * (-0.15) * (y[,2] + 0.2)
  #              = -0.3 * (y[,2] + 0.2)
  # For n=1: -0.3 * (0.2 + 0.2) = -0.12
  # For n=2: -0.3 * (0.8 + 0.2) = -0.3
  # For n=3: -0.3 * (-0.5 + 0.2) = 0.09
  lp_cont = c(-0.12, -0.3, 0.09)

  # lp_discrete = 2 * pairwise_cross[0,0] * (x - 0) = 0.4 * c(0, 1, 1)
  lp_disc = c(0.0, 0.4, 0.4)

  hand_mean_y0 = muy[1] + 0.5 * (lp_cont + lp_disc)
  hand_sd_y0 = sqrt(0.5)

  result = compute_conditional_mixed(
    x_observations = x_obs,
    y_observations = y_obs,
    predict_vars = as.integer(p), # index 1 = first continuous
    pairwise_disc = pairwise_disc,
    pairwise_cross = pairwise_cross,
    pairwise_cont = pairwise_cont,
    mux = mux,
    muy = muy,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  expect_equal(result[[1]][, 1], hand_mean_y0, tolerance = 1e-10)
  expect_equal(result[[1]][1, 2], hand_sd_y0, tolerance = 1e-10)

  # --- Predict continuous variable 1 (internal index p + 1 = 2) ---
  # theta_jj = -2 * pairwise_cont[2,2] = -2 * (-0.75) = 1.5
  # cond_var_1 = 1/1.5
  # lp_continuous = 2 * pairwise_cont[2,1] * (y[,1] - muy[1])
  #              = 2 * (-0.15) * (y[,1] - 0.1) = -0.3 * (y[,1] - 0.1)
  # For n=1: -0.3 * (0.5 - 0.1) = -0.12
  # For n=2: -0.3 * (-0.3 - 0.1) = 0.12
  # For n=3: -0.3 * (1.2 - 0.1) = -0.33
  lp_cont2 = c(-0.12, 0.12, -0.33)

  # lp_discrete = 2 * pairwise_cross[0,1] * (x - 0) = 0.2 * c(0, 1, 1)
  lp_disc2 = c(0.0, 0.2, 0.2)

  theta_jj_y1 = -2 * pairwise_cont[2, 2]
  hand_mean_y1 = muy[2] + (1 / theta_jj_y1) * (lp_cont2 + lp_disc2)
  hand_sd_y1 = sqrt(1 / theta_jj_y1)

  result2 = compute_conditional_mixed(
    x_observations = x_obs,
    y_observations = y_obs,
    predict_vars = as.integer(p + 1L), # index 2 = second continuous
    pairwise_disc = pairwise_disc,
    pairwise_cross = pairwise_cross,
    pairwise_cont = pairwise_cont,
    mux = mux,
    muy = muy,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  expect_equal(result2[[1]][, 1], hand_mean_y1, tolerance = 1e-10)
  expect_equal(result2[[1]][1, 2], hand_sd_y1, tolerance = 1e-10)
})


# ==============================================================================
# Test 7: Predicting all variables at once
# ==============================================================================
# Verify that predicting multiple variables yields the same results
# as predicting each one individually.

test_that("predicting all variables at once matches individual predictions", {
  n = 3L
  p = 2L
  q = 1L

  x_obs = matrix(
    c(
      0L, 1L, 1L,
      1L, 0L, 1L
    ),
    nrow = n, ncol = p
  )
  y_obs = matrix(c(0.5, -0.3, 1.2), nrow = n, ncol = q)

  pairwise_disc = matrix(c(0.0, 0.15, 0.15, 0.0), nrow = p, byrow = TRUE)
  pairwise_cross = matrix(c(0.2, 0.4), nrow = p, ncol = q)
  pairwise_cont = matrix(-1.0, nrow = q, ncol = q)
  mux = matrix(c(-0.5, 0.2), nrow = p, ncol = 1)
  muy = c(0.1)

  num_categories = c(1L, 1L)
  variable_type = c("ordinal", "ordinal")
  baseline_category = c(0L, 0L)

  # Predict all three variables at once
  result_all = compute_conditional_mixed(
    x_observations = x_obs,
    y_observations = y_obs,
    predict_vars = c(0L, 1L, as.integer(p)),
    pairwise_disc = pairwise_disc,
    pairwise_cross = pairwise_cross,
    pairwise_cont = pairwise_cont,
    mux = mux,
    muy = muy,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  # Predict each individually
  result_v0 = compute_conditional_mixed(
    x_observations = x_obs,
    y_observations = y_obs,
    predict_vars = 0L,
    pairwise_disc = pairwise_disc, pairwise_cross = pairwise_cross, pairwise_cont = pairwise_cont,
    mux = mux, muy = muy,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  result_v1 = compute_conditional_mixed(
    x_observations = x_obs,
    y_observations = y_obs,
    predict_vars = 1L,
    pairwise_disc = pairwise_disc, pairwise_cross = pairwise_cross, pairwise_cont = pairwise_cont,
    mux = mux, muy = muy,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  result_vy = compute_conditional_mixed(
    x_observations = x_obs,
    y_observations = y_obs,
    predict_vars = as.integer(p),
    pairwise_disc = pairwise_disc, pairwise_cross = pairwise_cross, pairwise_cont = pairwise_cont,
    mux = mux, muy = muy,
    num_categories = num_categories,
    variable_type = variable_type,
    baseline_category = baseline_category
  )

  expect_equal(length(result_all), 3L)
  expect_equal(as.matrix(result_all[[1]]), as.matrix(result_v0[[1]]),
    tolerance = 1e-10
  )
  expect_equal(as.matrix(result_all[[2]]), as.matrix(result_v1[[1]]),
    tolerance = 1e-10
  )
  expect_equal(as.matrix(result_all[[3]]), as.matrix(result_vy[[1]]),
    tolerance = 1e-10
  )
})
