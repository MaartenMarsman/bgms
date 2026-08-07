# Regression coverage for check_warmup_complete().
#
# The function had no tests pinning its output: the compliance fixture skips
# nuts_diag$warmup_check entirely. These pin the return contract and the two
# reported autocorrelation fields, on synthetic traces with known behaviour.

test_that("check_warmup_complete returns one value per chain for every field", {
  set.seed(1)
  energy = t(replicate(3, as.numeric(arima.sim(list(ar = 0.5), n = 500))))
  out = check_warmup_complete(energy)

  expected = c(
    "warmup_incomplete", "energy_slope", "slope_significant",
    "ebfmi_first_half", "ebfmi_second_half", "var_ratio",
    "energy_tau", "slope_t_corrected"
  )
  expect_true(all(expected %in% names(out)))
  for(nm in expected) expect_length(out[[nm]], 3)

  expect_type(out$warmup_incomplete, "logical")
  expect_type(out$slope_significant, "logical")
  expect_type(out$energy_tau, "double")
  expect_type(out$slope_t_corrected, "double")
})


test_that("short or degenerate traces return NA diagnostics rather than erroring", {
  # Fewer than 20 draws overall.
  short = matrix(rnorm(2 * 10), nrow = 2)
  out = check_warmup_complete(short)
  expect_false(any(out$warmup_incomplete))
  expect_true(all(is.na(out$energy_tau)))

  # Enough columns, but one chain is almost entirely non-finite.
  set.seed(2)
  energy = t(replicate(2, as.numeric(arima.sim(list(ar = 0.5), n = 100))))
  energy[2, 5:100] = NA_real_
  out2 = check_warmup_complete(energy)
  expect_false(out2$warmup_incomplete[2])
  expect_true(is.na(out2$energy_tau[2]))
})


test_that("a stationary trace has tau near 3 and a corrected t below the threshold", {
  # rho = 0.5 gives an integrated autocorrelation time of about (1+rho)/(1-rho) = 3,
  # matching what real bgms energy traces show.
  set.seed(3)
  energy = matrix(as.numeric(arima.sim(list(ar = 0.5), n = 4000)), nrow = 1)
  out = check_warmup_complete(energy)

  expect_gt(out$energy_tau, 1.5)
  expect_lt(out$energy_tau, 6)
  expect_lt(abs(out$slope_t_corrected), 2.58)
})


test_that("a settling transient inflates the naive t and is caught by the flag", {
  # Exponential decay over the first ~600 draws, several energy SDs tall.
  set.seed(4)
  n = 2000
  noise = as.numeric(arima.sim(list(ar = 0.5), n = n))
  energy = matrix(noise + 6 * exp(-seq_len(n) / 300), nrow = 1)
  out = check_warmup_complete(energy)

  expect_true(out$slope_significant)
  expect_true(out$warmup_incomplete)
  expect_lt(out$energy_slope, 0) # settling drifts energy downward
  expect_gt(abs(out$slope_t_corrected), 2.58)
})


test_that("the reported fields do not change when the flag fires", {
  # energy_tau and slope_t_corrected are reported only. Removing the trend must
  # not alter warmup_incomplete, which is driven by the naive statistic.
  set.seed(5)
  energy = t(replicate(4, as.numeric(arima.sim(list(ar = 0.5), n = 1000))))
  out = check_warmup_complete(energy)

  naive_flag = abs(out$energy_slope) > 0 &
    out$slope_significant | out$ebfmi_first_half < 0.3 | out$var_ratio > 2.0
  expect_equal(out$warmup_incomplete, unname(naive_flag))
})


test_that("integrated_act recovers a known autocorrelation time", {
  set.seed(6)
  # tau = (1 + rho) / (1 - rho); rho = 0.8 gives tau = 9.
  x = as.numeric(arima.sim(list(ar = 0.8), n = 20000))
  expect_gt(integrated_act(x), 6)
  expect_lt(integrated_act(x), 13)

  # White noise has tau = 1.
  expect_lt(integrated_act(rnorm(20000)), 1.5)

  # Too short to estimate: documented to return 1.
  expect_equal(integrated_act(rnorm(5)), 1)
})


test_that("a constant half leaves the flag logical rather than NA", {
  # A half with zero variance makes ebfmi_first_half and var_ratio NaN and the
  # OLS standard error zero, so every criterion comparison returns NA. The old
  # `FALSE || NA` left warmup_incomplete at NA, which which() and any() then
  # skip silently: the chain vanished from the reported issues instead of
  # showing up in them.
  set.seed(21)
  tail_half = as.numeric(arima.sim(list(ar = 0.5), n = 100))
  energy = rbind(
    c(rep(3.5, 100), tail_half), # constant first half
    c(tail_half, rep(3.5, 100)) # constant second half
  )
  out = check_warmup_complete(energy)

  expect_type(out$warmup_incomplete, "logical")
  expect_false(anyNA(out$warmup_incomplete))
  expect_type(out$slope_significant, "logical")
  expect_false(anyNA(out$slope_significant))

  # The undefined criteria are still reported as the NaN they are; only the
  # flag derived from them is forced to TRUE/FALSE.
  expect_true(is.nan(out$ebfmi_first_half[1]))

  # And a fully constant trace, where all three criteria are undefined. The
  # trend criterion abstains on its own rather than dividing one rounding
  # error by another, so the trace passes quietly and the slope reads as the
  # zero it is.
  flat = expect_no_warning(
    check_warmup_complete(matrix(rep(2, 2 * 200), nrow = 2))
  )
  expect_equal(flat$warmup_incomplete, c(FALSE, FALSE))
  expect_equal(unname(flat$energy_slope), c(0, 0))
  expect_equal(flat$slope_significant, c(FALSE, FALSE))
})
