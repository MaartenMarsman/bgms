# ==============================================================================
# Unit tests for validate_sampler()
# Phase A.7 of the R scaffolding refactor.
# ==============================================================================

# Helper: minimal valid call with sensible defaults
vs = function(...) {
  defaults = list(
    update_method     = c("nuts", "adaptive-metropolis", "gibbs"),
    target_accept     = NULL,
    iter              = 1000L,
    warmup            = 250L,
    nuts_max_depth    = 10L,
    learn_mass_matrix = TRUE,
    chains            = 2L,
    cores             = 2L,
    seed              = 42L,
    display_progress  = "none",
    is_continuous     = FALSE,
    edge_selection    = FALSE,
    verbose           = FALSE
  )
  args = modifyList(defaults, list(...))
  do.call(validate_sampler, args)
}


# ==============================================================================
# 1. update_method  — match.arg
# ==============================================================================

test_that("default pair resolves to 'nuts'", {
  res = vs()
  expect_equal(res$update_method, "nuts")
})

test_that("explicit 'adaptive-metropolis' passes through", {
  res = vs(update_method = "adaptive-metropolis")
  expect_equal(res$update_method, "adaptive-metropolis")
})

test_that("invalid update_method errors", {
  expect_error(vs(update_method = "bogus"), "arg")
})

test_that("removed 'hamiltonian-mc' errors", {
  expect_error(vs(update_method = "hamiltonian-mc"), "arg")
})


# ==============================================================================
# 2. GGM guard  (is_continuous = TRUE)
# ==============================================================================

test_that("GGM defaults → nuts silently", {
  res = vs(is_continuous = TRUE)
  expect_equal(res$update_method, "nuts")
})

test_that("GGM + explicit 'adaptive-metropolis' OK", {
  res = vs(is_continuous = TRUE, update_method = "adaptive-metropolis")
  expect_equal(res$update_method, "adaptive-metropolis")
})

test_that("GGM + explicit 'nuts' OK", {
  res = vs(is_continuous = TRUE, update_method = "nuts")
  expect_equal(res$update_method, "nuts")
})

test_that("GGM + explicit 'gibbs' OK", {
  res = vs(is_continuous = TRUE, update_method = "gibbs")
  expect_equal(res$update_method, "gibbs")
})

test_that("gibbs is rejected for non-continuous data", {
  expect_error(
    vs(is_continuous = FALSE, update_method = "gibbs"),
    "available only for the Gaussian"
  )
})


# ==============================================================================
# 3. target_accept  — defaults and clamping
# ==============================================================================

test_that("NULL target_accept → 0.44 for adaptive-metropolis", {
  res = vs(update_method = "adaptive-metropolis", target_accept = NULL)
  expect_equal(res$target_accept, 0.44)
})

test_that("NULL target_accept → 0.80 for nuts", {
  res = vs(update_method = "nuts", target_accept = NULL)
  expect_equal(res$target_accept, 0.80)
})

test_that("user-supplied target_accept passes through", {
  res = vs(target_accept = 0.8)
  expect_equal(res$target_accept, 0.8)
})

test_that("target_accept clamped below epsilon", {
  res = vs(target_accept = 0)
  expect_gt(res$target_accept, 0)
})

test_that("target_accept clamped above 1 - epsilon", {
  res = vs(target_accept = 1)
  expect_lt(res$target_accept, 1)
})


# ==============================================================================
# 4. iter / warmup
# ==============================================================================

test_that("valid iter accepted", {
  res = vs(iter = 500L)
  expect_equal(res$iter, 500L)
})

test_that("zero iter errors", {
  expect_error(vs(iter = 0L), "iter")
})

test_that("negative iter errors", {
  expect_error(vs(iter = -1L), "iter")
})

test_that("valid warmup accepted (including 0)", {
  res = vs(warmup = 0L)
  expect_equal(res$warmup, 0L)
})

test_that("negative warmup errors", {
  expect_error(vs(warmup = -5L), "warmup")
})


# ==============================================================================
# 5. warmup warnings  (verbose = TRUE)
# ==============================================================================

test_that("no-edge-selection: warmup < 20 warns (NUTS)", {
  expect_warning(
    vs(
      update_method = "nuts", warmup = 10L,
      edge_selection = FALSE, verbose = TRUE
    ),
    "no mass matrix"
  )
})

test_that("no-edge-selection: 20 <= warmup < 150 warns (NUTS)", {
  expect_warning(
    vs(
      update_method = "nuts", warmup = 50L,
      edge_selection = FALSE, verbose = TRUE
    ),
    "proportional allocation"
  )
})

test_that("no-edge-selection: warmup >= 150 no warning (NUTS)", {
  expect_silent(
    vs(
      update_method = "nuts", warmup = 200L,
      edge_selection = FALSE, verbose = TRUE
    )
  )
})

test_that("edge-selection: warmup < 50 warns (NUTS)", {
  expect_warning(
    vs(
      update_method = "nuts", warmup = 30L,
      edge_selection = TRUE, verbose = TRUE
    ),
    "very short"
  )
})

test_that("edge-selection: 50 <= warmup < 200 warns (NUTS)", {
  expect_warning(
    vs(
      update_method = "nuts", warmup = 100L,
      edge_selection = TRUE, verbose = TRUE
    ),
    "proposal SD tuning skipped"
  )
})

test_that("edge-selection: 200 <= warmup < 300 warns (NUTS)", {
  expect_warning(
    vs(
      update_method = "nuts", warmup = 250L,
      edge_selection = TRUE, verbose = TRUE
    ),
    "limited proposal SD tuning"
  )
})

test_that("edge-selection: warmup >= 300 no warning (NUTS)", {
  expect_silent(
    vs(
      update_method = "nuts", warmup = 300L,
      edge_selection = TRUE, verbose = TRUE
    )
  )
})

test_that("adaptive-metropolis never fires warmup warnings", {
  expect_silent(
    vs(
      update_method = "adaptive-metropolis", warmup = 5L,
      edge_selection = TRUE, verbose = TRUE
    )
  )
})

test_that("verbose = FALSE suppresses all warmup warnings", {
  expect_silent(
    vs(
      update_method = "nuts", warmup = 5L,
      edge_selection = TRUE, verbose = FALSE
    )
  )
})


# ==============================================================================
# 6. nuts_max_depth
# ==============================================================================

test_that("nuts_max_depth passes through", {
  res = vs(nuts_max_depth = 8L)
  expect_equal(res$nuts_max_depth, 8L)
})

test_that("nuts_max_depth clamped to >= 1", {
  res = vs(nuts_max_depth = 1L)
  expect_equal(res$nuts_max_depth, 1L)
})


# ==============================================================================
# 7. learn_mass_matrix
# ==============================================================================

test_that("learn_mass_matrix TRUE passes through", {
  res = vs(learn_mass_matrix = TRUE)
  expect_true(res$learn_mass_matrix)
})

test_that("learn_mass_matrix FALSE passes through", {
  res = vs(learn_mass_matrix = FALSE)
  expect_false(res$learn_mass_matrix)
})


# ==============================================================================
# 8. chains / cores
# ==============================================================================

test_that("valid chains accepted", {
  res = vs(chains = 4L)
  expect_equal(res$chains, 4L)
})

test_that("zero chains errors", {
  expect_error(vs(chains = 0L), "chains")
})

test_that("valid cores accepted", {
  res = vs(cores = 1L)
  expect_equal(res$cores, 1L)
})

test_that("zero cores errors", {
  expect_error(vs(cores = 0L), "cores")
})


# ==============================================================================
# 9. seed
# ==============================================================================

test_that("integer seed passes through", {
  res = vs(seed = 123L)
  expect_equal(res$seed, 123L)
})

test_that("NULL seed generates a random integer", {
  res = vs(seed = NULL)
  expect_true(is.integer(res$seed))
  expect_length(res$seed, 1L)
})

test_that("negative seed errors", {
  expect_error(vs(seed = -1L), "seed")
})

test_that("NA seed errors", {
  expect_error(vs(seed = NA), "seed")
})


# ==============================================================================
# 10. display_progress  →  progress_type
# ==============================================================================

test_that("'per-chain' → 2L", {
  res = vs(display_progress = "per-chain")
  expect_equal(res$progress_type, 2L)
})

test_that("'total' → 1L", {
  res = vs(display_progress = "total")
  expect_equal(res$progress_type, 1L)
})

test_that("'none' → 0L", {
  res = vs(display_progress = "none")
  expect_equal(res$progress_type, 0L)
})

test_that("TRUE → 2L (per-chain)", {
  res = vs(display_progress = TRUE)
  expect_equal(res$progress_type, 2L)
})

test_that("FALSE → 0L (none)", {
  res = vs(display_progress = FALSE)
  expect_equal(res$progress_type, 0L)
})


# ==============================================================================
# 11. Full return structure
# ==============================================================================

test_that("return list has all expected elements", {
  res = vs()
  expected_names = c(
    "update_method", "target_accept", "iter", "warmup",
    "nuts_max_depth", "learn_mass_matrix",
    "chains", "cores", "seed", "progress_type", "progress_callback",
    "verbose"
  )
  expect_named(res, expected_names)
})


# ==============================================================================
# 12. R CMD check's two-core limit
# ==============================================================================

test_that("_R_CHECK_LIMIT_CORES_ caps the resolved sampler cores at two", {
  detected = as.integer(parallel::detectCores())
  skip_if(detected < 3L, "machine has fewer than three cores")

  withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = "TRUE")
  expect_equal(vs(cores = detected)$cores, 2L)
  expect_equal(vs(cores = 1L)$cores, 1L)
  expect_equal(normalize_parallel_cores(detected), 2L)
  expect_true(check_limit_cores())
})

test_that("the core cap is off without the check environment variable", {
  detected = as.integer(parallel::detectCores())
  skip_if(detected < 3L, "machine has fewer than three cores")

  withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = NA)
  expect_false(check_limit_cores())
  expect_equal(vs(cores = detected)$cores, detected)
  # The literal string "false" is what R CMD check writes when it does not
  # want the limit, and it is not a request for two cores.
  withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = "false")
  expect_false(check_limit_cores())
  expect_equal(normalize_parallel_cores(detected), detected)
})

test_that("resolved cores never exceed the machine and stay integer", {
  detected = as.integer(parallel::detectCores())
  withr::local_envvar(`_R_CHECK_LIMIT_CORES_` = NA)
  expect_equal(normalize_parallel_cores(detected + 100L), detected)
  expect_identical(normalize_parallel_cores(2), 2L)
  expect_identical(vs(cores = 2)$cores, 2L)
})
