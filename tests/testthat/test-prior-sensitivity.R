# Tests for the refit-based prior_sensitivity_check().
#
# prior_sensitivity_check() refits the model at a grid of fixed scales and
# classifies every edge; it works on any bgm() fit with edge selection.
# Structural tests use small fits and reduced anchor grids: every refit is a
# full MCMC run, so the anchor set is the dominant cost. The warm-vs-cold
# agreement certification runs in the BGMS_RUN_SLOW_TESTS tier.

skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run prior-sensitivity certifications"
  )
}

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


test_that("prior_sensitivity_check builds the anchored curve object", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:6],
    chains = 2, iter = 400, warmup = 300, seed = 21,
    display_progress = "none"
  )
  ps = suppressWarnings(suppressMessages(prior_sensitivity_check(fit,
    anchors = c(0.5, 1, 2), iter = 400, warmup = 300, ess_floor = 100,
    seed = 21
  )))

  expect_s3_class(ps, "bgms_prior_sensitivity")
  # one grid row per anchor (the 1x row is the original fit) plus the replicate
  expect_equal(sum(!ps$grid$replicate), 3L)
  expect_equal(sum(ps$grid$replicate), 1L)
  expect_equal(sum(ps$grid$original_fit), 1L)
  expect_equal(ps$grid$seconds[ps$grid$original_fit], 0)
  expect_equal(nrow(ps$edges), 15L) # choose(6, 2)
  # the curve is dense and the anchors lie exactly on it
  expect_gte(length(ps$multipliers), 41L)
  expect_equal(ps$multipliers[ps$anchor_index], ps$anchors)
  expect_equal(ps$multipliers[ps$chosen_index], 1)
  expect_equal(dim(ps$log10_bf), c(length(ps$multipliers), 15L))
  # Exactness lives on the anchor fits, not the pooled curve: the 1x verdict
  # column and every chosen-scale quantity are the original fit's own RB
  # analysis. The reported PIP is the fit's RB inclusion, and the log10 BF is
  # its exact prior-odds transform (checked on the pip scale, since the log10
  # BF derivative amplifies a 1e-10 pip agreement near saturation).
  rb = extract_posterior_inclusion_probabilities(fit)
  rbv = unname(rb[upper.tri(rb)][order_upper_tri_rowmajor(6)])
  po = ps$edges$prior_inclusion_probability /
    (1 - ps$edges$prior_inclusion_probability)
  expect_equal(ps$edges$chosen_scale_pip, rbv, tolerance = 1e-10)
  pc = ps$edges$chosen_scale_pip
  expect_equal(
    ps$edges$chosen_scale_log10_bf,
    log10((pc / (1 - pc)) / po),
    tolerance = 1e-12
  )
  expect_equal(
    ps$edges$verdict_x1,
    verdict_from_lbf(ps$edges$chosen_scale_log10_bf, log10(ps$evidence_threshold)),
    tolerance = 1e-10
  )
  # the pooled curve at 1x tracks the fit's own raw 1x inclusion proportions:
  # dominated by the 1x anchor (an identity reweight) but pooled with the
  # neighbours, so on the log10 BF scale it agrees on median within a modest
  # margin. The exact 1x analysis lives in the chosen-scale columns above; this
  # is a tracking check on absolute log10 BF, robust to the tiny per-point MCSE
  # of near-saturated edges that makes a normalised ratio platform-unstable.
  g = do.call(rbind, fit$raw_samples$indicator)
  raw_lbf = unname(log10((colMeans(g) / (1 - colMeans(g))) / po))
  fin = is.finite(ps$log10_bf[ps$chosen_index, ]) & is.finite(raw_lbf)
  expect_lt(
    stats::median(abs(ps$log10_bf[ps$chosen_index, ] - raw_lbf)[fin]),
    0.5
  )
  # verdicts and movers take only the documented levels
  expect_true(all(unlist(ps$verdict) %in%
    c("presence", "undecided", "absence", NA)))
  expect_true(all(ps$edges$mover %in%
    c("stable", "indistinguishable-from-wobble", "moved-beyond-wobble")))
  # pooling keeps the curve finite even where an edge saturates at some scale
  expect_false(any(is.infinite(ps$log10_bf)))
  # the data-preferred scale is reported without a refit
  expect_true(is.finite(ps$preferred_scale$s_hat))
  # print and plot run
  expect_output(print(ps), "are the edge verdicts robust")
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
    chains = 2, iter = 500, warmup = 400, seed = 5, display_progress = "none"
  )
  # every asserted quantity is the original fit's own RB analysis, so the
  # anchor grid is minimal
  ps = suppressWarnings(suppressMessages(prior_sensitivity_check(fit,
    anchors = c(1, 2), iter = 300, warmup = 300, seed = 5
  )))
  expect_true(all(abs(ps$edges$prior_inclusion_probability - 0.2) < 1e-8))
  # log10 BF divides posterior odds by the 0.2/0.8 prior odds
  p = ps$edges$chosen_scale_pip
  expected = log10((p / (1 - p)) / (0.2 / 0.8))
  expect_equal(ps$edges$chosen_scale_log10_bf, expected, tolerance = 1e-8)
})


test_that("warm-started short refits agree with cold full refits within wobble", {
  skip_on_cran()
  skip_unless_slow()
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


test_that("anchors must include a multiplier other than 1", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:4],
    chains = 1, iter = 200, warmup = 200, seed = 3,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  expect_error(prior_sensitivity_check(fit, anchors = 1), "other than 1")
  expect_error(prior_sensitivity_check(fit, anchors = c(-1, 2)), "positive")
  expect_error(prior_sensitivity_check(fit, ess_floor = 0), "positive")
})


test_that("prior_sensitivity_check runs for GGM and mixed fits (cold refits)", {
  skip_on_cran()
  set.seed(32)
  xg = matrix(rnorm(180 * 5), 180, 5)
  fg = bgm(xg,
    variable_type = "continuous",
    iter = 300, warmup = 250, chains = 2, seed = 8,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  psg = suppressWarnings(suppressMessages(prior_sensitivity_check(fg,
    anchors = c(1, 2), iter = 250, warmup = 200, ess_floor = 100,
    seed = 8
  )))
  expect_s3_class(psg, "bgms_prior_sensitivity")
  expect_equal(psg$model_type, "ggm")
  expect_false(psg$warm) # continuous fits refit cold
  # curve points that clear the ESS floor carry a finite BF; masked ones are NA
  finite_pts = !is.na(psg$curve$anchor_used)
  expect_true(any(finite_pts))
  expect_true(all(is.finite(psg$log10_bf[finite_pts, 1])))
  expect_true(all(is.na(psg$log10_bf[!finite_pts, ])))

  n = 180
  xm = data.frame(
    a = sample(0:2, n, TRUE), b = sample(0:2, n, TRUE),
    c = sample(0:1, n, TRUE), y1 = rnorm(n), y2 = rnorm(n)
  )
  vt = c("ordinal", "ordinal", "ordinal", "continuous", "continuous")
  fm = bgm(xm,
    variable_type = vt, iter = 300, warmup = 250, chains = 2, seed = 9,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  psm = suppressWarnings(suppressMessages(prior_sensitivity_check(fm,
    anchors = c(1, 2), iter = 250, warmup = 200, ess_floor = 100,
    seed = 9
  )))
  expect_equal(psm$model_type, "mixed_mrf")
})


test_that("a warm-start list with the wrong length errors", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:5],
    chains = 2, iter = 150, warmup = 200, seed = 3,
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
