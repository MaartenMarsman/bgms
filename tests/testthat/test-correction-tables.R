# --------------------------------------------------------------------------- #
# Tests for the hierarchical-graph-prior correction tables.
#
# The table math has closed forms when the per-edge tilt slope is constant c:
#   C(theta) = (1 - theta + theta * exp(c))^E,
# so the per-pair curve is f(theta) = log(1 - theta + theta * exp(c)) up to
# an additive constant and the implied local slope is c everywhere. c = 0 is
# the untilted case where every curve is identically zero.
#
# The assertions are on the DEPLOYED quantities -- (fprime_density, fprime) as
# the C++ SBMCorrection reader consumes them -- rather than on anything the
# table merely happens to carry.
# --------------------------------------------------------------------------- #

# R mirror of SBMCorrection::fprime_at (src/priors/edge_prior_correction.h):
# linear interpolation on the tabulated slope curve, constant extension past
# either end. Written as its own function so the tests below read the table the
# way the deployed chain does, and so the mirror can be checked against the
# compiled reader itself (see the test interface call further down).
fprime_at_mirror = function(tab, d) {
  x = tab$fprime_density
  y = tab$fprime
  n = length(x)
  vapply(d, function(dd) {
    if(dd <= x[1]) {
      return(y[1])
    }
    if(dd >= x[n]) {
      return(y[n])
    }
    lo = max(which(x <= dd))
    t = (dd - x[lo]) / (x[lo + 1] - x[lo])
    y[lo] + t * (y[lo + 1] - y[lo])
  }, numeric(1))
}

test_that("theta grid is dense at both ends and stays inside (0, 1)", {
  g = ggm_correction_theta_grid(120L)

  expect_false(is.unsorted(g))
  expect_true(all(g > 0 & g < 1))
  expect_gte(sum(g < 0.1), 25)
  expect_gte(sum(g > 0.9), 25)
})

test_that("untilted sweep gives identically zero curves", {
  theta = ggm_correction_theta_grid(80L)
  tab = correction_table_from_edens(theta,
    edens_raw = theta,
    num_pairs = 1e6
  )

  expect_equal(tab$f, rep(0, length(theta)))
  expect_equal(tab$logC, rep(0, length(theta)))
  expect_equal(tab$fprime, rep(0, length(tab$fprime_density)),
    tolerance = 1e-12
  )
  expect_equal(tab$num_repaired, 0)
})

test_that("constant-slope tilt recovers the closed form", {
  c0 = -0.3
  theta = ggm_correction_theta_grid(120L)
  edens = theta * exp(c0) / (1 - theta + theta * exp(c0))
  tab = correction_table_from_edens(theta, edens, num_pairs = 1e6)

  expect_equal(tab$fprime, rep(c0, length(tab$fprime_density)),
    tolerance = 1e-10
  )
  expect_lte(length(tab$fprime_density), sum(theta <= 0.90))

  f_exact = log(1 - theta + theta * exp(c0))
  expect_equal(tab$f, f_exact - f_exact[1], tolerance = 2e-4)

  # Read through the deployed lookup: a constant-slope cell must return c0 at
  # any density, inside the tabulated range and past both ends.
  probe = c(
    min(tab$fprime_density) - 0.1, min(tab$fprime_density),
    mean(range(tab$fprime_density)), max(tab$fprime_density),
    max(tab$fprime_density) + 0.1
  )
  expect_equal(fprime_at_mirror(tab, probe), rep(c0, length(probe)),
    tolerance = 1e-10
  )
})

test_that("isotonic repair lifts a crash dip and reports it", {
  theta = ggm_correction_theta_grid(60L)
  edens = theta * exp(-0.3) / (1 - theta + theta * exp(-0.3))
  edens_dipped = edens
  edens_dipped[50] = 0.01

  tab = correction_table_from_edens(theta, edens_dipped, num_pairs = 1e6)

  expect_gte(tab$num_repaired, 1)
  expect_false(is.unsorted(tab$edens))
})

test_that("a repaired region deploys the ties-mean slope, not the run's last", {
  # The repair pools grid points 43:45 onto one density while theta keeps
  # rising, so the raw slope curve carries three distinct values at a single
  # abscissa. The C++ reader brackets by binary search and would land on the
  # last -- the lowest -- so the table must be aggregated before it is handed
  # over.
  theta = ggm_correction_theta_grid(60L)
  edens = theta * exp(-0.3) / (1 - theta + theta * exp(-0.3))
  edens[50] = 0.01
  tab = correction_table_from_edens(theta, edens, num_pairs = 1e6)

  # The contract the reader needs: a strictly increasing abscissa, no ties.
  expect_gt(length(tab$fprime_density), 1)
  expect_true(all(diff(tab$fprime_density) > 0))
  expect_length(tab$fprime, length(tab$fprime_density))

  # Rebuild the un-aggregated curve the way the table used to return it, find
  # the tied run, and check the table now carries that run's mean.
  keep = theta <= tab$fprime_theta_cap &
    tab$edens > 0.5 / tab$num_pairs & tab$edens < 1 - 0.5 / tab$num_pairs
  d_raw = tab$edens[keep]
  fp_raw = log((d_raw / (1 - d_raw)) * (1 - theta[keep]) / theta[keep])
  runs = rle(d_raw)
  tied = which(runs$lengths > 1)
  expect_gte(length(tied), 1)

  ends = cumsum(runs$lengths)
  for(k in tied) {
    hi = ends[k]
    lo = hi - runs$lengths[k] + 1L
    d_tied = runs$values[k]
    run_mean = mean(fp_raw[lo:hi])
    run_last = fp_raw[hi]
    # The bias is real, not a rounding artefact: the run must actually fall.
    expect_lt(run_last, run_mean)
    expect_equal(fprime_at_mirror(tab, d_tied), run_mean, tolerance = 1e-12)
  }
})

test_that("the compiled reader and the R mirror agree on the deployed table", {
  # SBMCorrection::fprime_at has no direct R entry, but compute_ce_sbm's fixed
  # point exposes it: with one cluster every expected degree density equals the
  # common edge probability e, so the returned slope is exactly f'(e) and e is
  # recoverable from the slope. Any disagreement between the compiled bracket
  # and the mirror above shows up as a mismatch here.
  theta = ggm_correction_theta_grid(60L)
  edens = theta * exp(-0.3) / (1 - theta + theta * exp(-0.3))
  edens[50] = 0.01
  tab = correction_table_from_edens(theta, edens, num_pairs = 1e6)

  quad_theta = seq(0.0025, 0.9975, length.out = 200L)
  quad_f = stats::approx(tab$theta, tab$f, quad_theta, rule = 2)$y
  q = 5L

  for(th in c(0.05, 0.3, 0.6, 0.8, 0.95)) {
    ce = test_sbm_compute_ce(
      cluster_assign = rep(1L, q),
      block_probs = matrix(th, 1L, 1L),
      fprime_density = tab$fprime_density,
      fprime = tab$fprime,
      quad_theta = quad_theta,
      quad_f = quad_f
    )
    c_hat = ce[1, 2]
    # The fixed point the C++ converged to, read back out of its own answer.
    e_hat = th * exp(c_hat) / (1 - th + th * exp(c_hat))
    expect_equal(c_hat, fprime_at_mirror(tab, e_hat),
      tolerance = 1e-8,
      info = sprintf("theta = %g", th)
    )
  }
})

test_that("a single pair yields the logC curve without slope pieces", {
  c0 = -0.3
  theta = ggm_correction_theta_grid(80L)
  edens = theta * exp(c0) / (1 - theta + theta * exp(c0))

  # num_pairs = 1: the resolvable density window (0.5, 0.5) is empty, so
  # the slope curve cannot be tabulated; the integrated curve still can.
  tab = correction_table_from_edens(theta, edens, num_pairs = 1)

  expect_null(tab$fprime)
  expect_null(tab$fprime_density)
  expect_true(all(is.finite(tab$logC)))
  f_exact = log(1 - theta + theta * exp(c0))
  expect_equal(tab$f, f_exact - f_exact[1], tolerance = 2e-4)
})

test_that("table build with cache round-trips and reuses the file", {
  cache_dir = file.path(tempdir(), "bgms-ctable-test")
  unlink(cache_dir, recursive = TRUE)
  old = options(bgms.correction_cache_dir = cache_dir)
  on.exit(options(old), add = TRUE)

  tab1 = ggm_correction_table(
    p = 4, n_grid = 12L, n_samples = 100L, n_warmup = 100L, n_seeds = 1L,
    update_method = "gibbs"
  )
  files = list.files(cache_dir)
  expect_length(files, 1)

  tab2 = ggm_correction_table(
    p = 4, n_grid = 12L, n_samples = 100L, n_warmup = 100L, n_seeds = 1L,
    update_method = "gibbs"
  )
  expect_identical(tab1, tab2)

  expect_identical(tab1$cell$q, 4L)
  expect_identical(tab1$cell$slab_family, "cauchy")
  expect_equal(tab1$cell$eta, 1)
  expect_equal(tab1$cell$delta, 0.5 * log(4))
})

test_that("the cache key is scoped to the package version", {
  cache_dir = file.path(tempdir(), "bgms-ctable-version-test")
  unlink(cache_dir, recursive = TRUE)
  old = options(bgms.correction_cache_dir = cache_dir)
  on.exit(options(old), add = TRUE)

  version = as.character(utils::packageVersion("bgms"))
  cell = ggm_correction_cell(
    4, 0.5 * log(4), cauchy_prior(scale = 2.5), exponential_prior(eta = 1)
  )
  key = ggm_correction_table_key(cell, 12L, 100L, 100L, 1L, "gibbs", 1L)
  expect_match(key, version, fixed = TRUE)

  # The sweep that builds the table is code, so a table another release wrote
  # into the shared cache directory must not be served to this one. Two
  # versions' tables for the same cell were observed to differ by up to 0.042
  # in edge density under the unversioned key.
  stale = file.path(cache_dir, sub(version, "0.0.0.0", key, fixed = TRUE))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(sentinel = TRUE), stale)

  table = ggm_correction_table(
    p = 4, n_grid = 12L, n_samples = 100L, n_warmup = 100L, n_seeds = 1L,
    update_method = "gibbs"
  )
  expect_null(table$sentinel)
  expect_identical(table$cell$q, 4L)
  expect_true(file.exists(file.path(cache_dir, key)))
})

test_that("the cache key is scoped to the sweep seed and the table version", {
  cell = ggm_correction_cell(
    4, 0.5 * log(4), cauchy_prior(scale = 2.5), exponential_prior(eta = 1)
  )
  key_of = function(base_seed) {
    ggm_correction_table_key(cell, 12L, 100L, 100L, 1L, "gibbs", base_seed)
  }

  # base_seed offsets every chain in the sweep, so two base_seeds are two
  # different Monte-Carlo tables for one cell and must not share a file.
  expect_false(identical(key_of(1L), key_of(2L)))

  # The reader's contract changed twice over -- aggregated fprime abscissae,
  # no fed grid -- so v1 files must never be served to this reader again.
  expect_match(key_of(1L), "^ggm_ctable_v2_")
})

test_that("a table that cannot be cached is still returned, with a warning", {
  # The read side already tolerates an unusable cache; the write side must too,
  # or a read-only cache directory throws away a sweep that has already run.
  cache_dir = file.path(tempdir(), "bgms-ctable-unwritable")
  unlink(cache_dir, recursive = TRUE)
  # A regular FILE where the cache directory should be: dir.create cannot
  # replace it and saveRDS cannot open a path underneath it.
  writeLines("not a directory", cache_dir)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  withr::local_options(bgms.correction_cache_dir = cache_dir)

  table = NULL
  warned = character(0)
  capture.output(
    withCallingHandlers(
      table <- ggm_correction_table(
        p = 4, n_grid = 12L, n_samples = 100L, n_warmup = 100L, n_seeds = 1L,
        update_method = "gibbs"
      ),
      warning = function(w) {
        warned <<- c(warned, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  )
  # Exactly one warning, and it says what the failure means for the caller --
  # saveRDS's own "cannot open compressed file" must not escape on its own.
  expect_length(warned, 1L)
  expect_match(warned, "could not be cached")
  expect_identical(table$cell$q, 4L)
  expect_true(all(is.finite(table$logC)))
})

test_that("a two-variable Stochastic-Block cell skips the sweep it cannot use", {
  # One tilted pair leaves the resolvable density window empty, so the slope
  # curve is NULL for certain. That is knowable up front and the sweep costs
  # minutes, so the warning must arrive without one having run.
  cache_dir = file.path(tempdir(), "bgms-ctable-sbm-precheck")
  unlink(cache_dir, recursive = TRUE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  withr::local_options(bgms.correction_cache_dir = cache_dir)

  prior = list(
    edge_selection = TRUE, edge_prior = "Stochastic-Block",
    interaction_prior_type = "cauchy", pairwise_scale = 2.5,
    scale_prior_type = "exponential", scale_rate = 1, scale_shape = 1,
    delta = 0.5 * log(2)
  )
  sampler = list(cores = 1L, verbose = FALSE, progress_type = 0L)

  expect_warning(
    out <- ggm_edge_prior_correction(prior, sampler, num_variables = 2L),
    "slope curve is not resolvable"
  )
  expect_null(out)
  # No sweep ran: nothing was built, so nothing was cached.
  expect_false(dir.exists(cache_dir))
})

test_that("the progress bar renders the label, counts, and percentage", {
  pb = new_correction_progress(120L, prefix = "Correction table")
  mid = paste(capture.output(pb$update(70L)), collapse = "")
  expect_match(mid, "Correction table:", fixed = TRUE)
  expect_match(mid, "70/120", fixed = TRUE)
  expect_match(mid, "58.3%", fixed = TRUE)
  full = paste(capture.output(pb$update(120L)), collapse = "")
  expect_match(full, "120/120 (100.0%)", fixed = TRUE)
})

test_that("the parallel sweep matches the serial sweep cell for cell", {
  skip_on_cran()
  skip_on_os("windows")

  args = list(
    p = 4, theta = c(0.2, 0.5, 0.8), delta = 0.5 * log(4),
    interaction_prior = cauchy_prior(scale = 2.5),
    precision_scale_prior = gamma_prior(shape = 1, eta = 1),
    n_samples = 200L, n_warmup = 100L, n_seeds = 2L, update_method = "gibbs"
  )
  serial = do.call(sweep_prior_edge_density, c(args, cores = 1L))
  parallel = do.call(sweep_prior_edge_density, c(args, cores = 2L))

  expect_equal(parallel$edens_raw, serial$edens_raw)
})

test_that("the build announces itself once; a cache hit is silent", {
  cache_dir = file.path(tempdir(), "bgms-ctable-msg-test")
  unlink(cache_dir, recursive = TRUE)
  old = options(bgms.correction_cache_dir = cache_dir)
  on.exit(options(old), add = TRUE)

  build_args = list(
    p = 4, n_grid = 12L, n_samples = 100L, n_warmup = 100L, n_seeds = 1L,
    update_method = "gibbs", verbose = TRUE
  )
  # capture.output silences the progress bar (stdout); messages pass through.
  capture.output(
    expect_message(
      do.call(ggm_correction_table, build_args),
      "Building the edge-selection prior correction table"
    )
  )
  capture.output(
    expect_no_message(do.call(ggm_correction_table, build_args))
  )
})

test_that("gibbs and adaptive-metropolis sweeps agree on edge density", {
  skip_on_cran()

  theta = c(0.2, 0.5, 0.8)
  args = list(
    p = 6, theta = theta, delta = 0.5 * log(6),
    interaction_prior = cauchy_prior(scale = 2.5),
    precision_scale_prior = gamma_prior(shape = 1, eta = 1),
    n_samples = 3000L, n_warmup = 500L, n_seeds = 2L
  )
  sweep_gibbs = do.call(
    sweep_prior_edge_density,
    c(args, update_method = "gibbs")
  )
  sweep_am = do.call(
    sweep_prior_edge_density,
    c(args, update_method = "adaptive-metropolis")
  )

  expect_lt(max(abs(sweep_gibbs$edens_raw - sweep_am$edens_raw)), 0.02)
})
