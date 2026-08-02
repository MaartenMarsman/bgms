# The extrapolation notice reports the retained sweeps, not the whole run. The
# sampler initializes from a complete graph, so during warmup every mediating
# block is ~q and lands past the hull; quoting that share would report an
# initial transient as if it described the posterior.
#
# The size-cap block fits past the surface's trained hull and costs ~7.5 min on
# the 2-core CI runner -- most of this file. That is a heavy surface case, not
# a nightly heartbeat, so it runs in the weekly certification tier (T2,
# BGMS_RUN_CERTIFICATION); the counter arithmetic around it stays local (T0).

chain_with = function(...) {
  counters = c(
    n_hit = 0, n_miss = 0, n_pred = 0, n_add = 0, cache_size = 0,
    n_extrap = 0, max_extrap_size = 0, n_slope_floor = 0,
    n_pred_retained = 0, n_extrap_retained = 0, max_extrap_size_retained = 0
  )
  over = c(...)
  counters[names(over)] = over
  list(list(zratio = list(counters = as.list(counters))))
}

test_that("no extrapolation stays silent", {
  expect_no_message(
    out <- bgms:::zratio_extrapolation_notice(chain_with(n_pred = 5000))
  )
  expect_false(out)
})

test_that("a warmup-only transient is reported as one", {
  chains = chain_with(
    n_pred = 5000, n_extrap = 600, max_extrap_size = 198,
    n_pred_retained = 2000, n_extrap_retained = 0,
    max_extrap_size_retained = 0
  )
  expect_message(
    bgms:::zratio_extrapolation_notice(chains),
    "during warmup only"
  )
  # 600 of the 3000 warmup evaluations, and the stored draws are called out as
  # unaffected rather than the run being reported at 12%.
  expect_message(bgms:::zratio_extrapolation_notice(chains), "20.0% of warmup")
  expect_message(bgms:::zratio_extrapolation_notice(chains), "stored draws are unaffected")
})

test_that("a retained share is quoted with warmup in brackets", {
  chains = chain_with(
    n_pred = 5000, n_extrap = 900, max_extrap_size = 198,
    n_pred_retained = 2000, n_extrap_retained = 100,
    max_extrap_size_retained = 96
  )
  # Retained: 100 / 2000 = 5.0%. Warmup: (900 - 100) / (5000 - 2000) = 26.7%.
  expect_message(bgms:::zratio_extrapolation_notice(chains), "5.0% of the hierarchical")
  expect_message(bgms:::zratio_extrapolation_notice(chains), "warmup 26.7%")
  # The retained maximum, not the warmup one, sizes the reported block.
  expect_message(bgms:::zratio_extrapolation_notice(chains), "largest 96 variables")
  # The measured band statistics replace the old vague "can slightly reduce
  # accuracy", and the notice quotes a median and a maximum rather than a bare
  # figure a reader would take for a bound.
  expect_message(bgms:::zratio_extrapolation_notice(chains), "median of 0.0006 nats")
  expect_message(bgms:::zratio_extrapolation_notice(chains), "at most 0.0060")
})

test_that("counters missing the retained split claim no phase", {
  # Output from a chain that predates the split carries no retained counters, so
  # the notice reports the whole-run share. Claiming "warmup only" there would
  # reassure the reader on evidence the output does not contain.
  legacy = list(list(zratio = list(counters = list(
    n_pred = 5000, n_extrap = 600, max_extrap_size = 198
  ))))
  # One notice, read twice: asserting the absence of "warmup" against a second
  # call would let the notice itself through to the console, since a filtered
  # expectation only intercepts the messages it matches.
  msg = capture_messages(bgms:::zratio_extrapolation_notice(legacy))
  expect_match(msg, "12.0% of the hierarchical", all = FALSE)
  expect_no_match(msg, "warmup")
})

test_that("the retained split reaches R from a fit past the size cap", {
  skip_on_cran()
  skip_unless_certification()
  # Only a fit on more variables than the cap can extrapolate at all: at or
  # below it the cap tier anchors the hull past every reachable block.
  set.seed(6)
  q = 90
  y = matrix(rnorm(120 * q), 120, q)
  msgs = character(0)
  withCallingHandlers(
    bgm(
      x = y, variable_type = "continuous", iter = 50, warmup = 50,
      interaction_prior = normal_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 1, rate = 2),
      precision_graph_prior = "hierarchical",
      update_method = "gibbs", chains = 1, cores = 1, seed = 3,
      display_progress = "none", verbose = TRUE
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  # The complete-graph init guarantees warmup extrapolation past the cap, so the
  # notice fires and its phase split proves the counters reached R.
  hit = grep("anchored size range", msgs, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "warmup")
})
