# Regression tests pinning fixed defects in Blume-Capel handling, GGM
# prediction and simulation scales, split-Rhat, and group-difference priors.

# ---------------------------------------------------------------------------
# mixed model with a single continuous variable
# ---------------------------------------------------------------------------

test_that("mixed model with one continuous variable does not crash", {
  skip_on_cran()

  set.seed(2)
  n = 60
  x = cbind(
    d1 = sample(0:2, n, replace = TRUE),
    d2 = sample(0:2, n, replace = TRUE),
    c1 = rnorm(n)
  )

  expect_no_error(
    bgm(
      x,
      variable_type = c("ordinal", "ordinal", "continuous"),
      edge_selection = FALSE, update_method = "adaptive-metropolis",
      iter = 30, warmup = 30, chains = 1, seed = 1, display_progress = "none"
    )
  )
})

# ---------------------------------------------------------------------------
# OMRF SBM allocation samples carry node names, not edge names
# ---------------------------------------------------------------------------

test_that("OMRF SBM allocation samples carry node names", {
  skip_on_cran()

  set.seed(3)
  p = 4
  x = matrix(sample(0:2, 60 * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("V", seq_len(p))

  fit = bgm(
    x,
    variable_type = "ordinal", edge_selection = TRUE,
    edge_prior = sbm_prior(),
    iter = 30, warmup = 30, chains = 1, seed = 1, display_progress = "none"
  )

  alloc_names = fit$raw_samples$parameter_names$allocations
  expect_equal(alloc_names, colnames(x))
})

# ---------------------------------------------------------------------------
# explicit integer baseline_category = 0L is honoured
# ---------------------------------------------------------------------------

test_that("explicit baseline_category = 0L works for Blume-Capel", {
  skip_on_cran()

  set.seed(4)
  p = 3
  x = matrix(sample(0:3, 60 * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("V", seq_len(p))

  # An explicit integer 0L must be accepted as a provided baseline, not
  # read as missing.
  expect_no_error(
    bgm(
      x,
      variable_type = "blume-capel", baseline_category = 0L,
      edge_selection = FALSE, update_method = "adaptive-metropolis",
      iter = 20, warmup = 20, chains = 1, seed = 1, display_progress = "none"
    )
  )
})

# ---------------------------------------------------------------------------
# Blume-Capel prediction/simulation on the original category scale
# ---------------------------------------------------------------------------

test_that("Blume-Capel recode round-trips through the stored shift", {
  # Data coded 2..5 shifts to internal 0..3; the reverse map restores 2..5.
  x = matrix(c(2, 3, 4, 5), ncol = 1)
  shift = c(2)

  internal = recode_data_for_prediction(
    x,
    num_categories = 3L, is_ordinal = FALSE,
    category_levels = NULL, blume_capel_shift = shift
  )
  expect_equal(as.numeric(internal), c(0, 1, 2, 3))

  original = recode_simulated_to_original(
    internal,
    category_levels = list(NULL), blume_capel_shift = shift
  )
  expect_equal(as.numeric(original), c(2, 3, 4, 5))
})

test_that("simulate() returns Blume-Capel draws on the original scale", {
  skip_on_cran()

  set.seed(7)
  p = 3
  x = matrix(sample(1:5, 150 * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("B", seq_len(p))

  fit = bgm(
    x,
    variable_type = "blume-capel", baseline_category = 3,
    edge_selection = FALSE, update_method = "adaptive-metropolis",
    iter = 200, warmup = 200, chains = 1, seed = 1, display_progress = "none"
  )

  # The shift back to the original 1-based coding is stored per variable.
  expect_equal(extract_arguments(fit)$blume_capel_shift, rep(1, p))

  sim = simulate(fit, nsim = 200, seed = 2)
  # Original scale is 1..5; the internal 0-based scale would include 0.
  expect_gte(min(sim), 1)
  expect_lte(max(sim), 5)
})

# ---------------------------------------------------------------------------
# GGM prediction centers newdata on the training means
# ---------------------------------------------------------------------------

test_that("GGM predict uses stored training means, not newdata means", {
  skip_on_cran()

  set.seed(8)
  p = 4
  n = 120
  x = sweep(matrix(rnorm(n * p), n, p), 2, c(5, -3, 10, 0), "+")
  colnames(x) = paste0("V", seq_len(p))

  fit = bgm(
    x,
    variable_type = "continuous", update_method = "adaptive-metropolis",
    edge_selection = FALSE, iter = 300, warmup = 300, chains = 1, seed = 1,
    display_progress = "none"
  )

  # Training means are stored so prediction centers on the fitted scale.
  expect_equal(extract_arguments(fit)$column_means, colMeans(x),
    ignore_attr = TRUE
  )

  # With one row, centering by newdata's own mean gave back the input value
  # (a variable predicting itself). Centering on the training means does not.
  row1 = x[1, , drop = FALSE]
  pred = predict(fit, newdata = row1)
  pred_means = vapply(pred, function(m) m[, "mean"], numeric(1))
  expect_false(isTRUE(all.equal(unname(pred_means), unname(row1[1, ]))))
})

# ---------------------------------------------------------------------------
# Rhat splits each chain so it detects within-chain drift (split-Rhat)
# ---------------------------------------------------------------------------

test_that("split_chains halves each chain", {
  a = array(seq_len(10 * 2 * 3), dim = c(10, 2, 3))
  s = split_chains(a)
  expect_equal(dim(s), c(5L, 4L, 3L))
  # sub-chains 1-2 are the two halves of the original chain 1
  expect_equal(s[, 1, ], a[1:5, 1, ])
  expect_equal(s[, 2, ], a[6:10, 1, ])
  # an odd iteration count drops the middle draw
  odd = array(1:7, dim = c(7, 1, 1))
  expect_equal(as.vector(split_chains(odd)), c(1, 2, 3, 5, 6, 7))
})

test_that("split-Rhat flags a drifting chain that classic Rhat misses", {
  set.seed(17)
  # A single chain that drifts: classic Gelman-Rubin needs >1 chain and returns
  # NA, but split-Rhat compares the two halves and sees the trend.
  drift = array(cumsum(rnorm(2000, mean = 0.03)) + rnorm(2000), dim = c(2000, 1, 1))
  expect_true(is.na(.compute_rhat_cpp(drift)))
  expect_gt(.compute_rhat_cpp(split_chains(drift)), 1.1)

  # A stationary chain still passes.
  stat = array(rnorm(2000), dim = c(2000, 1, 1))
  expect_lt(.compute_rhat_cpp(split_chains(stat)), 1.05)
})

# ---------------------------------------------------------------------------
# bgmCompare difference prior family is independent of the interaction
# prior, defaulting to Normal (F-119)
# ---------------------------------------------------------------------------

test_that("difference_family selects the difference prior independently", {
  skip_on_cran()

  set.seed(8)
  ng = 2
  npg = 40
  p = 3
  x = matrix(sample(0:2, ng * npg * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("V", seq_len(p))
  g = rep(seq_len(ng), each = npg)

  run = function(...) {
    bgmCompare(
      x,
      group_indicator = g, ...,
      difference_selection = FALSE, iter = 80, warmup = 80, chains = 1,
      seed = 3, display_progress = "none"
    )
  }

  cauchy = run(difference_family = "Cauchy")
  normal = run(difference_family = "Normal")

  # Same seed, different prior family -> different difference draws.
  expect_false(isTRUE(all.equal(
    cauchy$raw_samples$pairwise[[1]], normal$raw_samples$pairwise[[1]]
  )))

  # The default is Normal (F-119). Omitting difference_family must reproduce
  # the explicit "Normal" run draw for draw at the same seed, and must not
  # reproduce the "Cauchy" one -- the previous version of this block passed
  # difference_family = "Cauchy" explicitly and so asserted nothing about the
  # default. Exact equality: same seed, same sampler path, no tolerance.
  default = run()
  expect_equal(
    default$raw_samples$pairwise[[1]], normal$raw_samples$pairwise[[1]]
  )
  expect_false(isTRUE(all.equal(
    default$raw_samples$pairwise[[1]], cauchy$raw_samples$pairwise[[1]]
  )))

  # The baseline interaction prior defaults to the Normal too, mirroring bgm()
  # (F-119). Same idiom: the default fit must reproduce the explicit Normal run
  # and diverge from the explicit Cauchy one. extract_arguments() does not
  # surface the slab family for compare fits, so the draws are the evidence.
  expect_equal(
    default$raw_samples$pairwise[[1]],
    run(interaction_prior = normal_prior(scale = 1))$raw_samples$pairwise[[1]]
  )
  expect_false(isTRUE(all.equal(
    default$raw_samples$pairwise[[1]],
    run(interaction_prior = cauchy_prior(scale = 1))$raw_samples$pairwise[[1]]
  )))

  # An unknown family is rejected.
  expect_error(
    bgmCompare(
      x,
      group_indicator = g, difference_family = "laplace",
      difference_selection = FALSE, iter = 10, warmup = 10, chains = 1,
      seed = 1, display_progress = "none"
    )
  )
})

# ---------------------------------------------------------------------------
# GGM summary with edge selection keeps the mixture columns,
# on the association scale
# ---------------------------------------------------------------------------

test_that("GGM selection summary is a mixture summary on the association scale", {
  skip_on_cran()

  set.seed(9)
  p = 4
  n = 200
  x = matrix(rnorm(n * p), n, p)
  colnames(x) = paste0("V", seq_len(p))

  fit = bgm(
    x,
    variable_type = "continuous", update_method = "adaptive-metropolis",
    edge_selection = TRUE, iter = 300, warmup = 300, chains = 1, seed = 2,
    display_progress = "none"
  )

  pw = fit$posterior_summary_pairwise
  # Same summarizer as the discrete models: the derived composite ESS and its
  # inclusion-share (bottleneck) columns are present for the weight row.
  expect_true(all(c("n_eff", "share_incl") %in% colnames(pw)))
  # And the mean is on the association scale coef() reports.
  cf = coef(fit)$pairwise
  expect_equal(sort(pw$mean), sort(cf[upper.tri(cf)]), tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# GGM simulate() returns data on the original scale
# ---------------------------------------------------------------------------

test_that("GGM simulate() draws are on the training data scale", {
  skip_on_cran()

  set.seed(10)
  p = 3
  n = 150
  shift = c(5, -3, 10)
  x = sweep(matrix(rnorm(n * p), n, p), 2, shift, "+")
  colnames(x) = paste0("V", seq_len(p))

  fit = bgm(
    x,
    variable_type = "continuous", update_method = "adaptive-metropolis",
    edge_selection = FALSE, iter = 200, warmup = 200, chains = 1, seed = 1,
    display_progress = "none"
  )

  sim = simulate(fit, nsim = 2000, seed = 4)
  # Column means of the simulated data sit near the training means, not zero.
  expect_equal(unname(colMeans(sim)), unname(colMeans(x)), tolerance = 0.5)
})

# ---------------------------------------------------------------------------
# bgmCompare Blume-Capel prediction/simulation on the original
# category scale
# ---------------------------------------------------------------------------

test_that("bgmCompare handles Blume-Capel variables coded 1-5", {
  skip_on_cran()

  set.seed(11)
  ng = 2
  npg = 60
  p = 3
  x = matrix(sample(1:5, ng * npg * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("B", seq_len(p))
  g = rep(seq_len(ng), each = npg)

  fit = bgmCompare(
    x,
    group_indicator = g, variable_type = "blume-capel",
    baseline_category = 3, difference_selection = FALSE,
    iter = 100, warmup = 100, chains = 1, seed = 2,
    display_progress = "none"
  )

  # The shift back to the original 1-based coding is stored per variable.
  expect_equal(extract_arguments(fit)$blume_capel_shift, rep(1, p))

  # simulate() returns the original 1..5 coding, not the internal 0..4 one.
  sim = simulate(fit, nsim = 200, group = 1, seed = 3)
  expect_gte(min(sim), 1)
  expect_lte(max(sim), 5)

  # predict() on original-scale newdata runs and returns finite probabilities.
  pred = predict(fit, newdata = x[1:5, ], group = 1)
  expect_true(all(vapply(pred, function(m) all(is.finite(m)), logical(1))))
})

# ---------------------------------------------------------------------------
# parallel dispatch keeps the R API on the main thread
#
# Chains run under a helper-thread parallelFor while the main thread polls
# for interrupts, progress, and the R callback. Chain seeding is independent
# of the dispatch, so serial and parallel runs must produce identical draws,
# and the callback must fire during a parallel run.
# ---------------------------------------------------------------------------

test_that("serial and parallel dispatch produce identical draws", {
  skip_on_cran()

  # RcppParallel >= 6.0.0 bundles oneTBB 2022, whose scheduler does not preserve
  # the bitwise identity of the second chain between serial and parallel dispatch
  # on Windows. Each chain still samples the same posterior (posterior means
  # agree to within Monte Carlo error), so this is a reproducibility limit, not a
  # correctness bug; the bitwise guarantee still holds on the other platforms and
  # on earlier RcppParallel.
  if(Sys.info()[["sysname"]] == "Windows" &&
    utils::packageVersion("RcppParallel") >= "6.0.0") {
    skip("bitwise serial/parallel identity not preserved on Windows under RcppParallel >= 6.0.0")
  }

  set.seed(11)
  p = 4
  x = matrix(sample(0:2, 120 * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("V", seq_len(p))

  serial = bgm(
    x,
    iter = 100, warmup = 100, chains = 2, cores = 1, seed = 7,
    display_progress = "none"
  )
  parallel = bgm(
    x,
    iter = 100, warmup = 100, chains = 2, cores = 2, seed = 7,
    display_progress = "none"
  )

  expect_identical(serial$raw_samples$pairwise, parallel$raw_samples$pairwise)
  expect_identical(serial$raw_samples$main, parallel$raw_samples$main)
})

test_that("the progress callback fires during a parallel run", {
  skip_on_cran()

  set.seed(12)
  p = 3
  x = matrix(sample(0:2, 100 * p, replace = TRUE), ncol = p)
  colnames(x) = paste0("V", seq_len(p))

  calls = 0L
  fit = bgm(
    x,
    iter = 600, warmup = 600, chains = 2, cores = 2, seed = 1,
    progress_callback = function(done, total) calls <<- calls + 1L
  )

  # At least the final finish() report; runs longer than the poll throttle
  # also report from the main-thread loop mid-run.
  expect_gte(calls, 1L)
  expect_s3_class(fit, "bgms")
})
