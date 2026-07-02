# Tests for the per-sampler warmup staging in WarmupSchedule, via the
# test_warmup_schedule() C++ entry. Focus: the Gibbs settle/selection split
# and its scaling with the warmup budget; NUTS and adaptive-Metropolis
# boundaries unchanged.

ws = function(warmup, edge_selection, learn_sd, select_during_warmup,
              probes = integer(0)) {
  test_warmup_schedule(
    warmup, edge_selection, learn_sd, select_during_warmup,
    as.integer(probes)
  )
}

test_that("gibbs warmup splits into a settle window and selection-active warmup", {
  s = ws(2000,
    edge_selection = TRUE, learn_sd = FALSE,
    select_during_warmup = TRUE,
    probes = c(0, 299, 300, 1999, 2000)
  )
  # Settle window = first 15% of warmup; selection active for the remainder
  expect_equal(s$stage3c_start, 300)
  expect_equal(s$total_warmup, 2000)
  expect_equal(s$selection_enabled, c(FALSE, FALSE, TRUE, TRUE, TRUE))
})

test_that("gibbs settle window rescales with the warmup budget", {
  expect_equal(
    ws(400, TRUE, FALSE, TRUE)$stage3c_start,
    60
  )
  expect_equal(
    ws(10000, TRUE, FALSE, TRUE)$stage3c_start,
    1500
  )
  # warmup = 0: selection active from the first (sampling) iteration
  s0 = ws(0, TRUE, FALSE, TRUE, probes = c(0, 1))
  expect_equal(s0$stage3c_start, 0)
  expect_equal(s0$selection_enabled, c(TRUE, TRUE))
})

test_that("adaptive-Metropolis staging is unchanged: selection starts at sampling", {
  s = ws(2000,
    edge_selection = TRUE, learn_sd = FALSE,
    select_during_warmup = FALSE,
    probes = c(0, 1000, 1999, 2000)
  )
  expect_equal(s$stage3c_start, 2000)
  expect_equal(s$selection_enabled, c(FALSE, FALSE, FALSE, TRUE))
})

test_that("NUTS staging is unchanged: 85/10/5 split with selection in stage 3c", {
  s = ws(2000,
    edge_selection = TRUE, learn_sd = TRUE,
    select_during_warmup = FALSE,
    probes = c(1899, 1900, 2000)
  )
  expect_equal(s$stage3b_start, 1700)
  expect_equal(s$stage3c_start, 1900)
  expect_equal(s$selection_enabled, c(FALSE, TRUE, TRUE))
})

test_that("no edge selection means the gate never opens during warmup", {
  s = ws(2000,
    edge_selection = FALSE, learn_sd = FALSE,
    select_during_warmup = TRUE,
    probes = c(0, 300, 1999, 2000)
  )
  expect_equal(s$selection_enabled, c(FALSE, FALSE, FALSE, FALSE))
})

test_that("gibbs edge selection runs end to end with the staged warmup", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  set.seed(31)
  K_true = diag(5)
  K_true[1, 2] = K_true[2, 1] = 0.4
  Y = MASS::mvrnorm(300, mu = rep(0, 5), Sigma = solve(K_true))
  fit = bgm(
    as.data.frame(Y),
    variable_type = "continuous",
    update_method = "gibbs",
    edge_selection = TRUE,
    iter = 200, warmup = 100,
    chains = 1, cores = 1, seed = 5,
    display_progress = "none", verbose = FALSE
  )
  expect_s3_class(fit, "bgms")
  pips = colMeans(S7::prop(fit, "raw_samples")$indicator[[1L]])
  expect_true(all(pips >= 0 & pips <= 1))
})
