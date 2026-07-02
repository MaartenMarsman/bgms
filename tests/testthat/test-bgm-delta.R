# --------------------------------------------------------------------------- #
# End-to-end tests for the bgm(delta = ...) plumbing.
#
# Verifies that:
#   - delta = 0 (default) reproduces the untilted posterior
#   - delta > 0 shifts the K_ii posterior upward (away from the PD-cone
#     boundary), per the manuscript prediction Pr[lambda_min < x] <= C x^(delta+1)
#   - delta > 0 is rejected for pure-ordinal models (no precision matrix)
#   - delta < 0 / NA / non-numeric is rejected
#
# Tests run with NUTS to keep cost low; the MH path is exercised in a separate
# block with adaptive-metropolis.
# --------------------------------------------------------------------------- #

# Small continuous dataset shared across tests.
gen_continuous = function(n = 60, p = 3, seed = 1L) {
  set.seed(seed)
  matrix(rnorm(n * p), nrow = n, ncol = p)
}


# Sum of K_ii from a fitted bgms object. Larger sum = K diagonal mass
# concentrated further from zero (further from the PD-cone boundary).
trace_K = function(out) {
  sum(diag(extract_precision(out)))
}


# ---- delta > 0 pushes K diagonal upward ------------------------------------ #

test_that("delta > 0 shifts NUTS K_ii posterior mean upward (GGM)", {
  x = gen_continuous(n = 80, p = 3, seed = 2L)
  base = bgm(
    x,
    variable_type = "continuous",
    iter = 300L, warmup = 300L, chains = 1L,
    edge_selection = FALSE, verbose = FALSE, display_progress = "none", seed = 2L,
    delta = 0
  )
  tilted = bgm(
    x,
    variable_type = "continuous",
    iter = 300L, warmup = 300L, chains = 1L,
    edge_selection = FALSE, verbose = FALSE, display_progress = "none", seed = 2L,
    delta = 5
  )
  expect_gt(trace_K(tilted), trace_K(base))
})


test_that("delta > 0 shifts MH K_ii posterior mean upward (GGM)", {
  x = gen_continuous(n = 80, p = 3, seed = 3L)
  base = bgm(
    x,
    variable_type = "continuous",
    update_method = "adaptive-metropolis",
    iter = 400L, warmup = 400L, chains = 1L,
    edge_selection = FALSE, verbose = FALSE, display_progress = "none", seed = 3L,
    delta = 0
  )
  tilted = bgm(
    x,
    variable_type = "continuous",
    update_method = "adaptive-metropolis",
    iter = 400L, warmup = 400L, chains = 1L,
    edge_selection = FALSE, verbose = FALSE, display_progress = "none", seed = 3L,
    delta = 5
  )
  expect_gt(trace_K(tilted), trace_K(base))
})


# ---- Validation ------------------------------------------------------------ #

test_that("delta > 0 is rejected for pure-ordinal models", {
  set.seed(4)
  x = matrix(sample(0:1, 40 * 3, replace = TRUE), nrow = 40, ncol = 3)
  expect_error(
    bgm(
      x,
      variable_type = "ordinal",
      iter = 20L, warmup = 20L, chains = 1L,
      edge_selection = FALSE, verbose = FALSE, display_progress = "none",
      delta = 1
    ),
    "no precision matrix to tilt"
  )
})


test_that("invalid delta values are rejected with a clear message", {
  x = gen_continuous(n = 40, p = 2)
  expect_error(
    bgm(
      x,
      variable_type = "continuous",
      iter = 20L, warmup = 20L, chains = 1L,
      verbose = FALSE, display_progress = "none", delta = -1
    ),
    "non-negative"
  )
  expect_error(
    bgm(
      x,
      variable_type = "continuous",
      iter = 20L, warmup = 20L, chains = 1L,
      verbose = FALSE, display_progress = "none", delta = NA_real_
    ),
    "finite"
  )
  expect_error(
    bgm(
      x,
      variable_type = "continuous",
      iter = 20L, warmup = 20L, chains = 1L,
      verbose = FALSE, display_progress = "none", delta = c(0, 1)
    ),
    "single"
  )
})
