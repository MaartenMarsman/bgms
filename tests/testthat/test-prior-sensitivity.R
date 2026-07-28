# Tests for the random interaction-slab-scale hyperprior and the refit-based
# prior_sensitivity_check().
#
# The hyperprior plumbing (scale storage/summary/extractor, mean-1 and scope
# validation, prior-only recovery) is unchanged. prior_sensitivity_check() now
# refits the model at a grid of fixed scales and classifies every edge; the old
# single-fit conditional-density reweighting curve has been removed.

test_that("the sampled scale is stored, summarized, and extractable", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:6],
    interaction_scale_prior = gamma_prior(2, 2),
    iter = 1000, warmup = 500, chains = 2, seed = 1,
    update_method = "adaptive-metropolis", display_progress = "none"
  )

  raw = fit$raw_samples
  expect_false(is.null(raw$interaction_scale))
  expect_equal(length(raw$interaction_scale), 2L)
  expect_equal(length(raw$interaction_scale[[1]]), 1000L)
  expect_true(all(unlist(raw$interaction_scale) > 0))

  draws = extract_scale_draws(fit)
  expect_equal(ncol(draws), 1L)
  expect_equal(colnames(draws), "interaction_scale")
  expect_true(all(draws > 0))

  s = summary(fit)
  expect_false(is.null(s$interaction_scale))
  expect_true(all(c("2.5%", "97.5%", "Rhat", "n_eff") %in%
    colnames(s$interaction_scale)))
})


test_that("fixed-scale fits expose no scale draws and error helpfully", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:5],
    iter = 500, warmup = 500, chains = 1, seed = 2,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  expect_null(fit$raw_samples$interaction_scale)
  expect_null(summary(fit)$interaction_scale)
  expect_error(extract_scale_draws(fit), "fixed interaction slab scale")
})


test_that("interaction_scale_prior validates the mean-1 and normal/omrf scope", {
  data("Wenchuan", package = "bgms")
  x = Wenchuan[, 1:4]

  expect_error(
    bgm(x,
      interaction_scale_prior = gamma_prior(2, 3),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"
    ),
    "mean 1"
  )
  expect_error(
    bgm(x,
      interaction_scale_prior = gamma_prior(shape = 2),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"
    ),
    "raw frame"
  )
  expect_error(
    bgm(x,
      interaction_prior = cauchy_prior(2.5),
      interaction_scale_prior = gamma_prior(2, 2),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"
    ),
    "normal slab"
  )
  expect_no_error(
    bgm(x,
      interaction_scale_prior = exponential_prior(rate = 1),
      iter = 50, warmup = 50, chains = 1, seed = 4,
      update_method = "adaptive-metropolis", display_progress = "none"
    )
  )
})


test_that("with the edge held out the scale draws recover the mean-1 hyperprior", {
  set.seed(42)
  n = 400
  x = cbind(sample(0:2, n, TRUE), sample(0:2, n, TRUE))
  fit = bgm(
    x,
    interaction_scale_prior = gamma_prior(2, 2),
    edge_prior = bernoulli_prior(1e-3),
    iter = 6000, warmup = 1500, chains = 2, seed = 3, edge_selection = TRUE,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  raw = fit$raw_samples
  s = do.call(c, raw$interaction_scale)
  g = do.call(rbind, raw$indicator)[, 1]
  expect_lt(mean(g), 0.05)
  u = s[g == 0]
  expect_equal(mean(u), 1, tolerance = 0.1)
  expect_equal(stats::var(u), 0.5, tolerance = 0.2)
})


# ---- refit-based prior_sensitivity_check() ----------------------------------

test_that("the wobble q95 pools threshold-relevant edges only", {
  # Two near-saturated edges (|log10 BF| > 3 at s0) with huge replicate
  # spread must not inflate the yardstick; the median and per-edge spread
  # keep covering every edge.
  lbf_s0 = c(0.2, -1.5, 2.9, 8.0, -12.0)
  lbf_rep = c(0.3, -1.3, 2.7, 9.3, -10.7)
  d = abs(lbf_s0 - lbf_rep)
  w = wobble_yardstick(lbf_s0, lbf_rep)
  expect_equal(w$per_edge, d)
  expect_equal(w$median, stats::median(d))
  expect_equal(w$q95, stats::quantile(d[1:3], 0.95, names = FALSE))
  expect_lt(w$q95, min(d[4:5]))
  # with every edge saturated the yardstick is undefined
  expect_true(is.na(wobble_yardstick(c(5, -7), c(6, -8))$q95))
})


test_that("prior_sensitivity_check needs edge selection", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:5],
    edge_selection = FALSE,
    iter = 300, warmup = 300, chains = 1, seed = 7,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  expect_error(prior_sensitivity_check(fit), "edge selection")
})


test_that("prior_sensitivity_check refits a fixed-scale fit and returns the grid object", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:6],
    chains = 2, iter = 1500, warmup = 1000, seed = 21,
    display_progress = "none"
  )
  ps = suppressWarnings(suppressMessages(prior_sensitivity_check(
    fit,
    scale_multipliers = c(0.5, 1, 2.5), seed = 21
  )))

  expect_s3_class(ps, "bgms_prior_sensitivity")
  # one grid row per multiplier plus the s0 replicate
  expect_equal(sum(!ps$grid$replicate), 3L)
  expect_equal(sum(ps$grid$replicate), 1L)
  expect_equal(nrow(ps$edges), 15L) # choose(6, 2)
  # verdicts and movers take only the documented levels
  expect_true(all(unlist(ps$verdict) %in%
    c("presence", "undecided", "absence", NA)))
  expect_true(all(ps$edges$mover %in%
    c("stable", "indistinguishable-from-wobble", "moved-beyond-wobble")))
  # chosen-scale verdict comes from the s0 refit column
  expect_equal(ps$edges$chosen_scale_verdict, ps$verdict[ps$chosen_index, ])
  # the data-preferred scale is reported without a refit
  expect_true(is.finite(ps$preferred_scale$s_hat))
  # print and plot run
  expect_output(print(ps), "data prefer")
  expect_silent({
    pdf(tempfile())
    plot(ps)
    dev.off()
  })
})


test_that("the chosen-scale verdict uses the per-edge prior odds, not 1/2", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:5],
    edge_prior = bernoulli_prior(0.2),
    chains = 2, iter = 1500, warmup = 1000, seed = 5, display_progress = "none"
  )
  ps = suppressWarnings(suppressMessages(prior_sensitivity_check(fit, seed = 5)))
  expect_true(all(abs(ps$edges$prior_inclusion_probability - 0.2) < 1e-8))
  # log10 BF divides posterior odds by the 0.2/0.8 prior odds
  p = ps$edges$chosen_scale_pip
  expected = log10((p / (1 - p)) / (0.2 / 0.8))
  expect_equal(ps$edges$chosen_scale_log10_bf, expected, tolerance = 1e-8)
})


test_that("warm-started short refits agree with cold full refits within wobble", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:8],
    chains = 4, iter = 2000, warmup = 1000, seed = 11,
    display_progress = "none"
  )
  ws = extract_warm_state(fit)
  s0 = get_fit_spec(fit)$prior$pairwise_scale
  pipv = function(f) {
    m = extract_posterior_inclusion_probabilities(f)
    m[upper.tri(m)]
  }
  warm = suppressMessages(refit_at_scale(fit, s0, ws,
    warmup = 500, iter = 1000,
    seed = 101, cores = 4L, sampler = "nuts"
  ))
  cold = suppressMessages(refit_at_scale(fit, s0,
    warm_state = NULL,
    warmup = 1000, iter = 2000, seed = 202, cores = 4L, sampler = "nuts"
  ))
  # agreement well within the s0-replicate wobble yardstick (~0.05 PIP q95)
  expect_lt(stats::median(abs(pipv(warm) - pipv(cold))), 0.03)
})


test_that("learned-scale extras appear only for hyperprior fits", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fixed = bgm(Wenchuan[, 1:6],
    chains = 2, iter = 1200, warmup = 800, seed = 6,
    display_progress = "none"
  )
  ps_fixed = suppressWarnings(suppressMessages(prior_sensitivity_check(fixed, seed = 6)))
  expect_null(ps_fixed$learned)

  learned = bgm(Wenchuan[, 1:6],
    interaction_scale_prior = gamma_prior(2, 2),
    chains = 2, iter = 1200, warmup = 800, seed = 6,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  ps_learned = suppressWarnings(suppressMessages(prior_sensitivity_check(learned, seed = 6)))
  expect_false(is.null(ps_learned$learned))
  # the scale-averaged verdict counts colMeans(indicator) through the prior odds
  gpool = do.call(rbind, learned$raw_samples$indicator)
  marg_pip = unname(colMeans(gpool))
  po = ps_learned$edges$prior_inclusion_probability
  po = po / (1 - po)
  expect_equal(
    unname(as.integer(ps_learned$learned$counts["presence"])),
    sum(log10((marg_pip / (1 - marg_pip)) / po) >=
      log10(ps_learned$evidence_threshold))
  )
})


test_that("prior_sensitivity_check runs for GGM and mixed fits (cold refits)", {
  skip_on_cran()
  set.seed(32)
  xg = matrix(rnorm(220 * 6), 220, 6)
  fg = bgm(xg,
    variable_type = "continuous",
    iter = 1000, warmup = 800, chains = 2, seed = 8,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  psg = suppressWarnings(suppressMessages(prior_sensitivity_check(fg,
    iter = 800, warmup = 500,
    seed = 8
  )))
  expect_s3_class(psg, "bgms_prior_sensitivity")
  expect_equal(psg$model_type, "ggm")
  expect_false(psg$warm) # continuous fits refit cold

  n = 250
  xm = data.frame(
    a = sample(0:2, n, TRUE), b = sample(0:2, n, TRUE),
    c = sample(0:1, n, TRUE), y1 = rnorm(n), y2 = rnorm(n)
  )
  vt = c("ordinal", "ordinal", "ordinal", "continuous", "continuous")
  fm = bgm(xm,
    variable_type = vt, iter = 1000, warmup = 700, chains = 2, seed = 9,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  psm = suppressWarnings(suppressMessages(prior_sensitivity_check(fm,
    iter = 700, warmup = 500,
    seed = 9
  )))
  expect_equal(psm$model_type, "mixed_mrf")
})


test_that("a warm-start list with the wrong length errors", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:5],
    chains = 2, iter = 300, warmup = 300, seed = 3,
    display_progress = "none"
  )
  ws = extract_warm_state(fit)
  spec = get_fit_spec(fit)
  spec$sampler$warmup = 100L
  spec$sampler$iter = 100L
  spec$sampler$cores = 2L
  spec$sampler$display_progress = "none"
  spec$sampler$progress_type = 0L
  spec$sampler$progress_callback = NULL
  # three parameter vectors but only two chains
  spec$initial_state = list(parameters = c(ws$parameters, ws$parameters[1]))
  expect_error(run_sampler(spec), "one entry per chain")
})


test_that("the random scale is fenced for the hierarchical graph prior", {
  set.seed(33)
  x = matrix(rnorm(200 * 4), 200, 4)
  expect_error(
    bgm(x,
      variable_type = "continuous",
      precision_graph_prior = "hierarchical",
      interaction_scale_prior = gamma_prior(2, 2),
      iter = 10, warmup = 10, chains = 1, display_progress = "none"
    ),
    "joint"
  )
})
