# --------------------------------------------------------------------------- #
# Tests for the C++ MCMC diagnostics: .compute_ess_cpp and .compute_rhat_cpp.
#
# ESS replaces coda::effectiveSize and is checked against it. Rhat is the
# classic split-Rhat (Gelman et al. 2013 / Stan; var_plus = (n-1)/n * W + B/n,
# Rhat = sqrt(var_plus / W)), NOT the df-adjusted coda::gelman.diag estimate, so
# it is checked against a pure-R implementation of that formula. The tests also
# verify safe behavior on pathological input (constant chains, NaN, Inf, very
# short chains, etc.).
# --------------------------------------------------------------------------- #

# Helper: build a 3D array [niter x nchains x nparam] from a matrix or vector.
make_array = function(x, niter, nchains, nparam = 1L) {
  array(x, dim = c(niter, nchains, nparam))
}

# Pure-R reference for classic split-Rhat, evaluated on the sub-chain array
# exactly as .compute_rhat_cpp sees it (splitting is done upstream by
# split_chains()). Mirrors the degenerate-case semantics: all sub-chains
# constant and equal -> NA; constant but unequal -> +Inf.
classic_rhat_ref = function(arr) {
  n = dim(arr)[1]
  m = dim(arr)[2]
  nparam = dim(arr)[3]
  vapply(seq_len(nparam), function(j) {
    X = matrix(arr[, , j], nrow = n, ncol = m)
    chain_means = colMeans(X)
    W = mean(apply(X, 2, stats::var)) # within: var() divides by n - 1
    B = n * stats::var(chain_means) # between: var() divides by m - 1
    if(W > 0) {
      var_plus = (n - 1) / n * W + B / n
      sqrt(var_plus / W)
    } else if(B > 0) {
      Inf
    } else {
      NA_real_
    }
  }, numeric(1))
}


# ---- Concordance with coda ------------------------------------------------- #

test_that("ESS matches coda::effectiveSize to machine precision", {
  skip_if_not_installed("coda")
  set.seed(42)
  niter = 500
  nchains = 2
  nparam = 5
  draws = array(rnorm(niter * nchains * nparam), dim = c(niter, nchains, nparam))

  ess_cpp = bgms:::.compute_ess_cpp(draws)

  ess_coda = numeric(nparam)
  for(j in seq_len(nparam)) {
    mcmc_list = coda::mcmc.list(
      lapply(seq_len(nchains), function(c) coda::mcmc(draws[, c, j]))
    )
    ess_coda[j] = coda::effectiveSize(mcmc_list)
  }
  expect_equal(ess_cpp, ess_coda, tolerance = 1e-10)
})

test_that("Rhat matches the pure-R classic split-Rhat on AR(1) chains", {
  set.seed(99)
  niter = 1000
  nchains = 4
  nparam = 6
  # Correlated draws with per-chain mean offsets so B and W are both non-trivial.
  draws = array(NA_real_, dim = c(niter, nchains, nparam))
  for(j in seq_len(nparam)) {
    for(c in seq_len(nchains)) {
      x = numeric(niter)
      x[1] = rnorm(1)
      for(i in 2:niter) x[i] = 0.8 * x[i - 1] + rnorm(1)
      draws[, c, j] = x + (c - 1) * 0.1
    }
  }
  split = bgms:::split_chains(draws)
  rhat_cpp = bgms:::.compute_rhat_cpp(split)
  rhat_ref = classic_rhat_ref(split)
  expect_equal(rhat_cpp, rhat_ref, tolerance = 1e-10)
})

test_that("ESS concordance holds for autocorrelated draws", {
  skip_if_not_installed("coda")
  set.seed(99)
  niter = 1000
  nchains = 3
  # AR(1) process with phi = 0.9
  draws = array(NA_real_, dim = c(niter, nchains, 1L))
  for(c in seq_len(nchains)) {
    x = numeric(niter)
    x[1] = rnorm(1)
    for(i in 2:niter) x[i] = 0.9 * x[i - 1] + rnorm(1)
    draws[, c, 1] = x
  }

  ess_cpp = bgms:::.compute_ess_cpp(draws)
  mcmc_list = coda::mcmc.list(
    lapply(seq_len(nchains), function(c) coda::mcmc(draws[, c, 1]))
  )
  ess_coda = unname(coda::effectiveSize(mcmc_list))
  expect_equal(ess_cpp, ess_coda, tolerance = 1e-10)
})


# ---- Single chain ---------------------------------------------------------- #

test_that("ESS works for single chain", {
  set.seed(1)
  draws = make_array(rnorm(200), niter = 200, nchains = 1)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_length(ess, 1)
  expect_true(is.finite(ess))
  expect_true(ess > 0)
})

test_that("Rhat returns NA for single chain", {
  set.seed(1)
  draws = make_array(rnorm(200), niter = 200, nchains = 1)
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_length(rhat, 1)
  expect_true(is.na(rhat))
})


# ---- Constant chains ------------------------------------------------------- #

test_that("ESS is NA for constant chain", {
  draws = make_array(rep(5.0, 200), niter = 100, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_true(is.na(ess))
})

test_that("Rhat is NA for constant chain", {
  draws = make_array(rep(5.0, 200), niter = 100, nchains = 2)
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_true(is.na(rhat))
})


# ---- Very short chains ----------------------------------------------------- #

test_that("ESS returns NA for single iteration", {
  draws = make_array(c(1.0, 2.0), niter = 1, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_true(is.na(ess))
})

test_that("Rhat returns NA for single iteration", {
  draws = make_array(c(1.0, 2.0), niter = 1, nchains = 2)
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_true(is.na(rhat))
})

test_that("ESS returns finite value for niter = 2", {
  set.seed(7)
  draws = make_array(rnorm(4), niter = 2, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_length(ess, 1)
  # May be NA if AR order saturates, but must not crash or be NaN
  expect_true(is.na(ess) || is.finite(ess))
})

test_that("ESS is finite for short chain (niter = 10)", {
  set.seed(3)
  draws = make_array(rnorm(20), niter = 10, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_true(is.finite(ess))
  expect_true(ess > 0)
})


# ---- NaN and Inf input ----------------------------------------------------- #

test_that("ESS returns NA when draws contain NaN", {
  draws = make_array(c(1, 2, NaN, 4, 5, 6, 7, 8, 9, 10), niter = 5, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_true(is.na(ess))
})

test_that("Rhat returns NA when draws contain NaN", {
  draws = make_array(c(1, 2, NaN, 4, 5, 6, 7, 8, 9, 10), niter = 5, nchains = 2)
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_true(is.na(rhat))
})

test_that("ESS returns NA when draws contain Inf", {
  draws = make_array(c(1, 2, Inf, 4, 5, 6, 7, 8, 9, 10), niter = 5, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_true(is.na(ess))
})

test_that("Rhat returns NA when draws contain Inf", {
  draws = make_array(c(1, 2, Inf, 4, 5, 6, 7, 8, 9, 10), niter = 5, nchains = 2)
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_true(is.na(rhat))
})

test_that("ESS returns NA when draws contain -Inf", {
  draws = make_array(c(1, 2, -Inf, 4, 5, 6, 7, 8, 9, 10), niter = 5, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_true(is.na(ess))
})

test_that("ESS returns NA when draws contain R's NA", {
  draws = make_array(c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10), niter = 5, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_true(is.na(ess))
})

test_that("Rhat returns NA when draws contain R's NA", {
  draws = make_array(c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10), niter = 5, nchains = 2)
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_true(is.na(rhat))
})

test_that("Rhat returns NA when draws contain -Inf", {
  draws = make_array(c(1, 2, -Inf, 4, 5, 6, 7, 8, 9, 10), niter = 5, nchains = 2)
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_true(is.na(rhat))
})

test_that("Rhat is finite for niter = 2", {
  set.seed(7)
  draws = make_array(rnorm(4), niter = 2, nchains = 2)
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_length(rhat, 1)
  expect_true(is.na(rhat) || is.finite(rhat))
})


# ---- Mixed pathological and good parameters -------------------------------- #

test_that("NaN in one parameter does not corrupt other parameters", {
  set.seed(5)
  good = rnorm(200)
  bad = c(rnorm(99), NaN, rnorm(100)) # NaN in chain 1
  draws = array(c(good, bad), dim = c(100, 2, 2))

  ess = bgms:::.compute_ess_cpp(draws)
  expect_length(ess, 2)
  expect_true(is.finite(ess[1])) # good parameter
  expect_true(ess[1] > 0)
  expect_true(is.na(ess[2])) # bad parameter

  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_length(rhat, 2)
  expect_true(is.finite(rhat[1]))
  expect_true(is.na(rhat[2]))
})


# ---- All-zero draws -------------------------------------------------------- #

test_that("ESS is NA for all-zero draws", {
  draws = make_array(rep(0.0, 200), niter = 100, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_true(is.na(ess))
})


# ---- Binary (0/1) draws ---------------------------------------------------- #

test_that("ESS is finite and positive for binary draws with variation", {
  set.seed(12)
  draws = make_array(sample(0:1, 200, replace = TRUE), niter = 100, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_true(is.finite(ess))
  expect_true(ess > 0)
})

test_that("Rhat is finite for binary draws with variation", {
  set.seed(12)
  draws = make_array(sample(0:1, 200, replace = TRUE), niter = 100, nchains = 2)
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_true(is.finite(rhat))
  expect_true(rhat > 0)
})


# ---- Multiple parameters --------------------------------------------------- #

test_that("batch ESS handles many parameters correctly", {
  set.seed(77)
  nparam = 50
  draws = array(rnorm(500 * 2 * nparam), dim = c(500, 2, nparam))
  ess = bgms:::.compute_ess_cpp(draws)
  expect_length(ess, nparam)
  expect_true(all(is.finite(ess)))
  expect_true(all(ess > 0))
})

test_that("batch Rhat handles many parameters correctly", {
  set.seed(77)
  nparam = 50
  draws = array(rnorm(500 * 2 * nparam), dim = c(500, 2, nparam))
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_length(rhat, nparam)
  expect_true(all(is.finite(rhat)))
  # Well-mixed iid draws should have Rhat close to 1
  expect_true(all(rhat > 0.95 & rhat < 1.05))
})


# ---- Near-constant draws (tiny variance) ----------------------------------- #

test_that("ESS handles near-constant draws without crashing", {
  set.seed(9)
  draws = make_array(1e10 + rnorm(200, sd = 1e-10), niter = 100, nchains = 2)
  ess = bgms:::.compute_ess_cpp(draws)
  expect_length(ess, 1)
  # Result may be NA (if variance is below threshold) or finite
  expect_true(is.na(ess) || is.finite(ess))
})


# ---- Output shape ---------------------------------------------------------- #

test_that("output length matches nparam dimension", {
  draws = array(rnorm(300), dim = c(50, 2, 3))
  expect_length(bgms:::.compute_ess_cpp(draws), 3)
  expect_length(bgms:::.compute_rhat_cpp(draws), 3)
})

test_that("ESS values are non-negative when finite", {
  set.seed(42)
  draws = array(rnorm(2000), dim = c(200, 2, 5))
  ess = bgms:::.compute_ess_cpp(draws)
  finite_ess = ess[is.finite(ess)]
  expect_true(all(finite_ess >= 0))
})


# ---- Rhat properties ------------------------------------------------------- #

test_that("Rhat is close to 1 for well-mixed iid chains", {
  set.seed(42)
  draws = array(rnorm(2000), dim = c(500, 2, 2))
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_true(all(abs(rhat - 1) < 0.05))
})

test_that("Rhat detects non-convergence (shifted chains)", {
  set.seed(42)
  niter = 500
  # Chain 1: mean 0, Chain 2: mean 10
  chain1 = rnorm(niter, mean = 0)
  chain2 = rnorm(niter, mean = 10)
  draws = array(c(chain1, chain2), dim = c(niter, 2, 1))
  rhat = bgms:::.compute_rhat_cpp(draws)
  expect_true(rhat > 1.5)
})


# ---- df-adjustment artifact fix (classic split-Rhat) ----------------------- #

test_that("a single brief excursion does not inflate Rhat", {
  # The typical near-saturated edge-indicator shape: 40,000 constant draws
  # across 4 chains with one 3-draw excursion in a single chain. The old
  # df-adjusted estimator returned ~1.2912 here regardless of the data; classic
  # split-Rhat sees this as essentially converged.
  niter = 10000
  nchains = 4
  x = array(1, dim = c(niter, nchains, 1))
  x[2000:2002, 2, 1] = 0
  rhat = bgms:::.compute_rhat_cpp(bgms:::split_chains(x))
  expect_lt(rhat, 1.01)
})

test_that("all-constant-and-equal sub-chains give NA", {
  x = array(1, dim = c(1000, 4, 1))
  rhat = bgms:::.compute_rhat_cpp(bgms:::split_chains(x))
  expect_true(is.na(rhat))
})

test_that("constant-but-unequal sub-chains give +Inf, not NA", {
  # Two chains stuck at 0, two stuck at 1: W == 0 but B > 0. Silence (NA) here
  # is the worst failure mode, so the estimator must raise a loud +Inf alarm.
  x = array(0, dim = c(1000, 4, 1))
  x[, 3:4, 1] = 1
  rhat = bgms:::.compute_rhat_cpp(bgms:::split_chains(x))
  expect_identical(rhat, Inf)
})

test_that("bgm RB inclusion Rhat is the classic split-Rhat on J draws, masked only for constant draws", {
  skip_on_cran()
  data = Wenchuan[, 1:6]
  fit = bgm(
    data,
    variable_type = "ordinal", chains = 2,
    iter = 400, warmup = 400, seed = 123,
    display_progress = "none", verbose = FALSE
  )
  summ = fit$posterior_summary_indicator
  reported = summ$Rhat

  # Independent classic split-Rhat from the Rao-Blackwellized (J) draws that the
  # summary reports, not the binary indicator draws.
  chains = fit$raw_samples$rb_inclusion
  nchains = length(chains)
  niter = nrow(chains[[1]])
  nparam = ncol(chains[[1]])
  arr = array(NA_real_, dim = c(niter, nchains, nparam))
  for(c in seq_len(nchains)) arr[, c, ] = chains[[c]]
  manual = classic_rhat_ref(bgms:::split_chains(arr))

  # Masking keys on the J draws being constant to double precision, not on the
  # indicator's flip count: where J varies there is an MCSE, an ESS and an Rhat
  # to report, and near-boundary edges in short runs often have zero flips.
  sds = apply(arr, 3, function(z) stats::sd(as.vector(z)))
  expect_true(all(is.na(reported[sds == 0])))
  expect_true(all(is.na(summ$n_eff[sds == 0])))
  expect_true(all(is.na(summ$mcse[sds == 0])))

  # Nothing that varies appreciably is masked; what masking remains beyond the
  # exactly-constant columns sits at the autocovariance kernel's numerical
  # floor, where the inclusion probability is at its bound.
  expect_true(all(sds[is.na(reported)] < 1e-6))

  # Edges whose indicator never flipped keep their numbers as long as J varies:
  # those are the near-boundary short-run edges where the MCSE is needed.
  zero_flip = summ[["n0->1"]] + summ[["n1->0"]] == 0
  varying = zero_flip & sds > 1e-6
  expect_true(any(varying))
  expect_true(all(is.finite(reported[varying])))
  expect_true(all(is.finite(summ$mcse[varying])))
  expect_true(all(is.finite(summ$n_eff[varying])))

  # On the unmasked edges the reported Rhat is exactly the classic split-Rhat on
  # the J draws, with no df adjustment.
  #
  # The tolerance sits above the cross-platform floating-point floor, not at the
  # local one: the two sides sum the same autocovariances in different orders,
  # so the last digits follow the machine's BLAS. Measured 3.1e-8 relative on
  # the Linux CI runner against a 1e-8 pin that holds on macOS (2026-08-02, run
  # 30719120127). What this guards against is a df adjustment or a wrong split,
  # which move Rhat by ~1e-2 -- four orders above this pin.
  keep = !is.na(reported)
  expect_equal(reported[keep], manual[keep], tolerance = 1e-6)
})


# ---- Indicator ESS -------------------------------------------------------- #

test_that("indicator ESS matches R reference implementation", {
  set.seed(42)
  niter = 500
  nchains = 2
  nparam = 3
  draws = array(
    sample(0:1, niter * nchains * nparam, replace = TRUE),
    dim = c(niter, nchains, nparam)
  )

  cpp_result = bgms:::.compute_indicator_ess_cpp(draws)

  # R reference. The mean pools every draw, but the transitions are counted
  # within each chain and summed: a chain boundary is not a transition (F-095).
  for(p in seq_len(nparam)) {
    vec = as.vector(draws[, , p])
    n_total = length(vec)
    p_hat = mean(vec)
    sd_r = sqrt(p_hat * (1 - p_hat))
    pairs = do.call(rbind, lapply(seq_len(nchains), function(cc) {
      ch = draws[, cc, p]
      cbind(curr = ch[-niter], nxt = ch[-1])
    }))
    g_curr = pairs[, "curr"]
    g_next = pairs[, "nxt"]
    n00 = sum(g_curr == 0 & g_next == 0)
    n01 = sum(g_curr == 0 & g_next == 1)
    n10 = sum(g_curr == 1 & g_next == 0)
    n11 = sum(g_curr == 1 & g_next == 1)
    a = n01 / (n00 + n01)
    b = n10 / (n10 + n11)
    tau_int = (2 - (a + b)) / (a + b)
    n_eff = n_total / tau_int
    mcse_r = sd_r / sqrt(n_eff)

    expect_equal(cpp_result[p, "mean"], p_hat, ignore_attr = TRUE)
    expect_equal(cpp_result[p, "sd"], sd_r, ignore_attr = TRUE)
    expect_equal(cpp_result[p, "mcse"], mcse_r, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(cpp_result[p, "n00"], n00, ignore_attr = TRUE)
    expect_equal(cpp_result[p, "n01"], n01, ignore_attr = TRUE)
    expect_equal(cpp_result[p, "n10"], n10, ignore_attr = TRUE)
    expect_equal(cpp_result[p, "n11"], n11, ignore_attr = TRUE)
    expect_equal(cpp_result[p, "n_eff_mixt"], n_eff, tolerance = 1e-12, ignore_attr = TRUE)
  }
})

test_that("indicator ESS returns correct column names", {
  draws = array(sample(0:1, 200, replace = TRUE), dim = c(50, 2, 2))
  result = bgms:::.compute_indicator_ess_cpp(draws)
  expect_equal(colnames(result), c("mean", "sd", "mcse", "n00", "n01", "n10", "n11", "n_eff_mixt"))
  expect_equal(nrow(result), 2)
})

test_that("indicator ESS handles all-zero draws (constant 0)", {
  draws = array(0, dim = c(100, 2, 1))
  result = bgms:::.compute_indicator_ess_cpp(draws)
  expect_equal(result[1, "mean"], 0, ignore_attr = TRUE)
  expect_equal(result[1, "sd"], 0, ignore_attr = TRUE)
  expect_true(is.na(result[1, "n_eff_mixt"]))
  expect_true(is.na(result[1, "mcse"]))
  expect_equal(result[1, "n01"], 0, ignore_attr = TRUE)
  expect_equal(result[1, "n10"], 0, ignore_attr = TRUE)
})

test_that("indicator ESS handles all-one draws (constant 1)", {
  draws = array(1, dim = c(100, 2, 1))
  result = bgms:::.compute_indicator_ess_cpp(draws)
  expect_equal(result[1, "mean"], 1, ignore_attr = TRUE)
  expect_equal(result[1, "sd"], 0, ignore_attr = TRUE)
  expect_true(is.na(result[1, "n_eff_mixt"]))
  expect_true(is.na(result[1, "mcse"]))
  expect_equal(result[1, "n01"], 0, ignore_attr = TRUE)
  expect_equal(result[1, "n10"], 0, ignore_attr = TRUE)
})

test_that("indicator ESS handles niter=1", {
  draws = array(c(1, 0), dim = c(1, 2, 1))
  result = bgms:::.compute_indicator_ess_cpp(draws)
  expect_true(is.na(result[1, "n_eff_mixt"]))
  expect_true(is.na(result[1, "mcse"]))
})

test_that("indicator ESS handles single chain", {
  set.seed(7)
  draws = array(sample(0:1, 200, replace = TRUE), dim = c(200, 1, 1))
  result = bgms:::.compute_indicator_ess_cpp(draws)
  expect_true(is.finite(result[1, "n_eff_mixt"]))
  expect_true(result[1, "n_eff_mixt"] > 0)
})

test_that("indicator ESS scales with multiple parameters", {
  set.seed(99)
  nparam = 10
  draws = array(
    sample(0:1, 500 * 2 * nparam, replace = TRUE),
    dim = c(500, 2, nparam)
  )
  result = bgms:::.compute_indicator_ess_cpp(draws)
  expect_equal(nrow(result), nparam)
  expect_true(all(result[, "n_eff_mixt"] > 0))
  # Transition counts sum to (niter - 1) per chain: one fewer pair than draws
  # in each chain, and no pair spanning a chain boundary (F-095).
  for(p in seq_len(nparam)) {
    expect_equal(
      unname(result[p, "n00"] + result[p, "n01"] + result[p, "n10"] + result[p, "n11"]),
      (500 - 1) * 2
    )
  }
})

test_that("indicator ESS returns all-NA row when draws contain NaN", {
  draws = array(c(1, 0, NaN, 1, 0, 1, 0, 1, 0, 1), dim = c(5, 2, 1))
  result = bgms:::.compute_indicator_ess_cpp(draws)
  expect_true(all(is.na(result[1, ])))
})

test_that("indicator ESS returns all-NA row when draws contain R's NA", {
  draws = array(c(1, 0, NA, 1, 0, 1, 0, 1, 0, 1), dim = c(5, 2, 1))
  result = bgms:::.compute_indicator_ess_cpp(draws)
  expect_true(all(is.na(result[1, ])))
})

test_that("indicator ESS returns all-NA row when draws contain Inf", {
  draws = array(c(1, 0, Inf, 1, 0, 1, 0, 1, 0, 1), dim = c(5, 2, 1))
  result = bgms:::.compute_indicator_ess_cpp(draws)
  expect_true(all(is.na(result[1, ])))
})

test_that("NaN in one indicator parameter does not corrupt others", {
  set.seed(5)
  good = sample(0:1, 200, replace = TRUE)
  bad = c(sample(0:1, 99, replace = TRUE), NaN, sample(0:1, 100, replace = TRUE))
  draws = array(c(good, bad), dim = c(100, 2, 2))

  result = bgms:::.compute_indicator_ess_cpp(draws)
  expect_true(is.finite(result[1, "mean"])) # good parameter
  expect_true(is.finite(result[1, "n_eff_mixt"]))
  expect_true(all(is.na(result[2, ]))) # bad parameter
})


# ---- NUTS diagnostics return-type contract --------------------------------- #

test_that("summarize_nuts_diagnostics honors the integer matrix contract", {
  mk_chain = function() {
    list(
      treedepth__ = c(2, 3, 2, 4),
      divergent__ = c(0, 0, 1, 0),
      energy__ = c(1.2, 0.8, 1.5, 0.9),
      accept_prob__ = c(0.9, 0.8, 0.95, 0.7)
    )
  }
  res = bgms:::summarize_nuts_diagnostics(
    list(mk_chain(), mk_chain()),
    nuts_max_depth = 4, verbose = FALSE
  )
  # Count fields are integer matrices; real-valued fields are double.
  expect_true(is.integer(res$treedepth))
  expect_true(is.integer(res$divergent))
  expect_true(is.double(res$energy))
  expect_true(is.double(res$accept_prob))
  # Exact integer tree-depth comparison: depth 4 hit once per chain.
  expect_equal(res$summary$max_tree_depth_hits, 2L)
})


# ---- Compare fallback parameter labels ------------------------------------- #

test_that("summarize_manual_compare brackets fallback parameter labels", {
  arr = array(rnorm(5 * 2 * 3), dim = c(5, 2, 3))
  res = bgms:::summarize_manual_compare(arr, "main_samples", param_names = NULL)
  expect_equal(res$parameter, c("param [1]", "param [2]", "param [3]"))
})


# ---- Transitions do not cross chain boundaries ------------------------------ #

test_that("the transition scan restarts at every chain boundary", {
  # Chain 1 is all zeros, chain 2 all ones. Within either chain nothing ever
  # flips, so both directional counts are zero. Scanning straight through the
  # pooled buffer would read the 0 -> 1 step at the boundary as a transition
  # and report n0->1 = 1 (F-095).
  draws = make_array(c(rep(0, 50), rep(1, 50)), niter = 50, nchains = 2)
  res = bgms:::.compute_indicator_ess_cpp(draws)

  expect_equal(res[1, "n01"], 0, ignore_attr = TRUE)
  expect_equal(res[1, "n10"], 0, ignore_attr = TRUE)
  expect_equal(res[1, "n00"], 49, ignore_attr = TRUE) # 49 within-chain steps in chain 1
  expect_equal(res[1, "n11"], 49, ignore_attr = TRUE) # 49 within-chain steps in chain 2
  # No flip in either direction, so there is no transition ESS to report.
  expect_true(is.na(res[1, "n_eff_mixt"]))
  expect_true(is.na(res[1, "mcse"]))
  # Every draw is still counted once for the pooled mean.
  expect_equal(res[1, "mean"], 0.5, ignore_attr = TRUE)

  # The counts are the within-chain totals summed, for any number of chains:
  # three chains that each alternate 0,1,0,1 give 3 * 1 of each per pair.
  alt = make_array(rep(c(0, 1, 0, 1), 3), niter = 4, nchains = 3)
  alt_res = bgms:::.compute_indicator_ess_cpp(alt)
  expect_equal(alt_res[1, "n01"], 6, ignore_attr = TRUE) # 2 per chain, 3 chains
  expect_equal(alt_res[1, "n10"], 3, ignore_attr = TRUE) # 1 per chain, 3 chains
  expect_equal(unname(alt_res[1, "n00"] + alt_res[1, "n11"]), 0)
})


# ---- E-BFMI on partly non-finite energy traces ------------------------------ #

test_that("E-BFMI is computed from the finite energy draws and NA is reported", {
  mk_chain = function(energy) {
    n = length(energy)
    list(
      treedepth__ = rep(2, n),
      divergent__ = rep(0, n),
      energy__ = energy,
      accept_prob__ = rep(0.9, n)
    )
  }
  set.seed(31)
  clean = as.numeric(arima.sim(list(ar = 0.5), n = 200))

  # One draw of chain 2 is non-finite. check_warmup_complete() already assessed
  # the finite draws only; compute_ebfmi() saw the raw trace, so the whole chain
  # came out NA -- which then made min_ebfmi NA for the run and dropped the
  # chain out of `which(ebfmi < 0.2)`, because NA < 0.2 is NA.
  holed = clean
  holed[100] = NA_real_
  res = bgms:::summarize_nuts_diagnostics(
    list(mk_chain(clean), mk_chain(holed)),
    nuts_max_depth = 10, verbose = FALSE
  )
  expect_false(anyNA(res$ebfmi))
  expect_true(is.finite(res$summary$min_ebfmi))
  # The single hole barely moves the statistic.
  expect_equal(res$ebfmi[[2]], res$ebfmi[[1]], tolerance = 0.05)

  # A chain with no usable energy has no E-BFMI. It must not erase the
  # minimum for the chains that do, and it must be named in the issues.
  res2 = bgms:::summarize_nuts_diagnostics(
    list(mk_chain(clean), mk_chain(rep(NA_real_, 200))),
    nuts_max_depth = 10, verbose = FALSE
  )
  expect_true(is.na(res2$ebfmi[[2]]))
  expect_equal(res2$summary$min_ebfmi, res2$ebfmi[[1]], ignore_attr = TRUE)
  expect_true(res2$has_issues)

  withr::local_options(bgms.verbose = TRUE)
  msg = capture.output(bgms:::summarize_nuts_diagnostics(
    list(mk_chain(clean), mk_chain(rep(NA_real_, 200))),
    nuts_max_depth = 10, verbose = TRUE
  ))
  expect_true(any(grepl("not computable in chain 2", msg, fixed = TRUE)))
})


# ---- Slab summary carries no structurally-NA Rhat --------------------------- #

test_that("summarize_slab reports no Rhat column", {
  # The included-only draws of one parameter are pooled into a single vector,
  # so .compute_rhat_cpp() returns NA for every row by contract. The column was
  # NA in every fit and read by nobody.
  set.seed(32)
  niter = 200
  nchains = 2
  nparam = 3
  pw = array(rnorm(niter * nchains * nparam), dim = c(niter, nchains, nparam))
  ind = array(rbinom(niter * nchains * nparam, 1, 0.7), dim = c(niter, nchains, nparam))
  pw[ind == 0] = 0

  slab = bgms:::summarize_slab(NULL, array3d = pw, array3d_ind = ind)
  expect_equal(colnames(slab), c("parameter", "mean", "mcse", "sd", "n_eff"))
  expect_false("Rhat" %in% colnames(slab))
  expect_true(all(is.finite(slab$n_eff)))

  # summarize_pair(), the only consumer, reads mean and mcse and still produces
  # its own Rhat from the full effect chain.
  pair = bgms:::summarize_pair(
    list(list(pairwise_samples = matrix(0, niter, nparam))),
    summ_slab = slab, array3d_id = ind, array3d_pw = pw
  )
  expect_true("Rhat" %in% colnames(pair))
  expect_true(all(is.finite(pair$Rhat)))
})


# ---- Adaptive-Metropolis acceptance vs move rate ---------------------------- #

test_that("am_diag separates the tuner acceptance probability from the move rate", {
  mk_chain = function(seed) {
    set.seed(seed)
    niter = 50
    list(
      main_samples = matrix(rnorm(niter * 2), niter, 2),
      # An excluded pairwise effect is pinned at 0 and never moves, so its move
      # rate is 0 no matter how the sampler is tuned.
      pairwise_samples = cbind(rnorm(niter), rep(0, niter)),
      am_accept_prob__ = rep(0.44, niter)
    )
  }
  res = bgms:::summarize_am_diagnostics(
    list(mk_chain(1), mk_chain(2)),
    names_main = c("m1", "m2"),
    names_pairwise = c("p1", "p2")
  )

  # The acceptance rate is the sampler's own per-sweep trace, comparable with
  # target_accept; the move rate keeps its own name and its per-parameter shape.
  expect_equal(dim(res$accept_prob), c(2L, 50L))
  expect_equal(rownames(res$accept_prob), c("chain 1", "chain 2"))
  expect_equal(res$summary$mean_accept_prob, 0.44)

  expect_equal(dim(res$move_rate), c(4L, 2L))
  expect_equal(rownames(res$move_rate), c("m1", "m2", "p1", "p2"))
  # The pinned parameter is the point: its move rate is 0 while the sampler was
  # accepting at the target rate the whole time, so the old reading of the move
  # rate as "the acceptance rate" understated it without bound.
  expect_equal(unname(res$move_rate["p2", ]), c(0, 0))
  expect_equal(res$summary$mean_move_rate, mean(res$move_rate))
  expect_false(isTRUE(all.equal(
    res$summary$mean_move_rate, res$summary$mean_accept_prob
  )))

  # A fit without the trace reports no acceptance rate rather than the move
  # rate under its name.
  bare = lapply(list(mk_chain(1), mk_chain(2)), function(chain) {
    chain$am_accept_prob__ = NULL
    chain
  })
  res_bare = bgms:::summarize_am_diagnostics(
    bare,
    names_main = c("m1", "m2"), names_pairwise = c("p1", "p2")
  )
  expect_null(res_bare$accept_prob)
  expect_true(is.na(res_bare$summary$mean_accept_prob))
  expect_equal(res_bare$move_rate, res$move_rate)
})


test_that("a real adaptive-metropolis fit carries both am_diag quantities", {
  skip_on_cran()
  fit = get_bgms_fit_adaptive_metropolis()

  expect_true(is.matrix(fit$am_diag$accept_prob))
  expect_equal(nrow(fit$am_diag$accept_prob), 2L)
  expect_true(all(fit$am_diag$accept_prob >= 0 & fit$am_diag$accept_prob <= 1))
  expect_true(is.matrix(fit$am_diag$move_rate))
  expect_equal(ncol(fit$am_diag$move_rate), 2L)
  expect_true(is.finite(fit$am_diag$summary$mean_accept_prob))
  expect_true(is.finite(fit$am_diag$summary$mean_move_rate))
})
