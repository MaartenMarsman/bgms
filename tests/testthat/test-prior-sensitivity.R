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
  expect_equal(w$censored, 0L)
})


test_that("an edge saturating in one refit does not carry the yardstick to Inf", {
  # A threshold-relevant edge whose replicate saturates has a censored spread.
  # Pooling the Inf would make the yardstick infinite, and every verdict move
  # would then read as run-to-run noise.
  lbf_s0 = c(0.2, -1.5, 2.9, -1.0)
  lbf_rep = c(0.3, -1.3, 2.7, -Inf)
  w = wobble_yardstick(lbf_s0, lbf_rep)

  expect_true(is.finite(w$q95))
  expect_equal(w$q95, stats::quantile(abs(lbf_s0 - lbf_rep)[1:3], 0.95, names = FALSE))
  expect_equal(w$censored, 1L)
  # The per-edge column still reports the censored edge as such.
  expect_true(is.infinite(w$per_edge[4]))
  # The median covers every edge with a measurable spread.
  expect_equal(w$median, stats::median(abs(lbf_s0 - lbf_rep)[1:3]))

  # Both refits saturating leaves nothing measurable, and the yardstick drops
  # out rather than becoming Inf.
  none = wobble_yardstick(c(-1.0, 0.5), c(-Inf, Inf))
  expect_true(is.na(none$q95))
  expect_equal(none$censored, 2L)
})


test_that("vary resolves against the frame the fit used", {
  standardized = list(prior = list(
    pairwise_scale = 2, scale_rate = 0.5, scale_eta = 1
  ))
  raw_frame = list(prior = list(
    pairwise_scale = 2, scale_rate = 1, scale_eta = NA_real_
  ))
  discrete = list(prior = list(pairwise_scale = 1, scale_rate = NA_real_))

  # "auto" follows the frame the diagonal prior was written in.
  expect_equal(resolve_vary(standardized, "auto")$mode, "slab-and-diagonal")
  expect_equal(resolve_vary(raw_frame, "auto")$mode, "slab")
  # A model with no precision diagonal has nothing to tie the slab to.
  expect_equal(resolve_vary(discrete, "auto")$mode, "none")
  expect_equal(resolve_vary(discrete, "slab-and-diagonal")$mode, "none")

  # An explicit mode overrides the frame; a raw-frame fit asked for the joint
  # sweep gets the eta its own scale implies (rate * s = 1 * 2).
  expect_equal(resolve_vary(raw_frame, "slab-and-diagonal")$eta, 2)
  expect_equal(resolve_vary(standardized, "slab-and-diagonal")$eta, 1)
  expect_true(is.na(resolve_vary(standardized, "slab")$eta))

  expect_error(resolve_vary(standardized, "diagonal"), "should be one of")
})


test_that("only slab-and-diagonal moves the raw diagonal rate", {
  standardized = list(prior = list(
    pairwise_scale = 1, scale_rate = 1, scale_eta = 1
  ))
  joint = resolve_vary(standardized, "slab-and-diagonal")
  slab = resolve_vary(standardized, "slab")

  # eta / s: a wider slab lowers the raw rate, so the prior's shape is fixed.
  expect_equal(vary_diagonal_rate(joint, 2.5), 0.4)
  expect_equal(vary_diagonal_rate(joint, 0.4), 2.5)
  # NULL leaves the fit's own rate in place, which is what refit_at_scale
  # already does with a spec that stores the rate already resolved.
  expect_null(vary_diagonal_rate(slab, 2.5))
  expect_null(vary_diagonal_rate(resolve_vary(
    list(prior = list(pairwise_scale = 1, scale_rate = NA_real_)), "auto"
  ), 2.5))
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
  skip_unless_slow()
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
  skip_unless_slow()
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
  # A default GGM specifies its diagonal in the standardized frame, so "auto"
  # holds eta fixed and the report says which prior moved.
  expect_equal(psg$vary$mode, "slab-and-diagonal")
  expect_match(
    paste(utils::capture.output(print(psg)), collapse = "\n"),
    "holding the standardized rate eta"
  )
  psg_slab = suppressWarnings(suppressMessages(prior_sensitivity_check(fg,
    vary = "slab", anchors = c(1, 2), iter = 250, warmup = 200,
    ess_floor = 100, seed = 8
  )))
  expect_equal(psg_slab$vary$mode, "slab")
  expect_match(
    paste(utils::capture.output(print(psg_slab)), collapse = "\n"),
    "interaction slab scale alone"
  )
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
  expect_equal(psm$vary$mode, "slab-and-diagonal")
})


test_that("compare_anchor_draws aligns every gated difference with its indicator", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgmCompare(
    x = Wenchuan[1:100, 1:4], group_indicator = rep(1:2, each = 50),
    iter = 200, warmup = 200, chains = 2, seed = 13,
    difference_selection = TRUE, display_progress = "none"
  )
  d = anchor_draws(fit)
  names_all = get_raw_samples(fit)$parameter_names

  # 6 pairwise differences + 4 variables x 4 thresholds of main difference,
  # against 4 + 6 = 10 indicators.
  expect_equal(ncol(d$theta[[1]]), 6L + 16L)
  expect_equal(ncol(d$indicator[[1]]), 10L)
  expect_equal(dim(d$gamma[[1]]), dim(d$theta[[1]]))
  expect_equal(d$family, "cauchy")

  # Each parameter carries its own indicator's draws, not a neighbour's. The
  # pairwise differences come first, in the order the indicators name them.
  pair_indicator = which(!grepl("(main)", names_all$indicator, fixed = TRUE))
  expect_equal(
    d$gamma[[1]][, seq_along(pair_indicator)],
    d$indicator[[1]][, pair_indicator],
    ignore_attr = TRUE
  )
  # A main-effect difference block repeats its variable's indicator, one column
  # per threshold.
  main_indicator = which(grepl("(main)", names_all$indicator, fixed = TRUE))
  expect_equal(
    d$gamma[[1]][, 6L + seq_len(4L)],
    d$indicator[[1]][, rep(main_indicator[1], 4L)],
    ignore_attr = TRUE
  )

  # The theta columns are the difference draws themselves, not the baselines.
  raw = get_raw_samples(fit)
  expect_equal(
    d$theta[[1]][, seq_len(6L)], raw$pairwise[[1]][, 7:12],
    ignore_attr = TRUE
  )
})


test_that("prior_sensitivity_check traces bgmCompare difference verdicts", {
  skip_on_cran()
  skip_unless_slow()
  data("Wenchuan", package = "bgms")
  fit = bgmCompare(
    x = Wenchuan[, 1:5], group_indicator = rep(1:2, length.out = nrow(Wenchuan)),
    iter = 800, warmup = 800, chains = 2, seed = 21,
    difference_selection = TRUE, display_progress = "none"
  )
  ps = suppressWarnings(suppressMessages(prior_sensitivity_check(fit,
    anchors = c(0.5, 1, 2), iter = 500, warmup = 500, ess_floor = 100, seed = 3
  )))

  expect_s3_class(ps, "bgms_prior_sensitivity")
  expect_equal(ps$edges$edge, get_raw_samples(fit)$parameter_names$indicator)
  # It sweeps the difference scale, not the interaction one.
  expect_equal(ps$unit$scale_field, "difference_scale")
  expect_equal(ps$chosen_scale, get_fit_spec(fit)$prior$difference_scale)
  # There is no data-preferred difference scale to compare against.
  expect_true(is.na(ps$preferred_scale$s_hat))

  out = paste(utils::capture.output(print(ps)), collapse = "\n")
  expect_match(out, "are the difference verdicts robust to the difference scale")
  # main_difference_selection is FALSE by default, so those indicators were
  # never updated; the report says so instead of counting them as undecided.
  expect_match(out, "never updated by the sampler and carry no verdict")

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(ps))
})


test_that("prior_sensitivity_check needs difference selection", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgmCompare(
    x = Wenchuan[1:80, 1:4], group_indicator = rep(1:2, each = 40),
    iter = 150, warmup = 150, chains = 2, seed = 4,
    difference_selection = FALSE, display_progress = "none"
  )
  expect_error(prior_sensitivity_check(fit), "difference selection")
})


test_that("the difference-scale reweighting reproduces a refit at that scale", {
  skip_on_cran()
  skip_unless_slow()
  data("Wenchuan", package = "bgms")
  x = Wenchuan[, 1:6]
  g = rep(1:2, length.out = nrow(x))
  compare_at = function(scale, seed) {
    bgmCompare(
      x = x, group_indicator = g, difference_scale = scale,
      iter = 3000, warmup = 2000, chains = 4, seed = seed,
      difference_selection = TRUE, display_progress = "none"
    )
  }
  pip_of = function(f) rowMeans(sapply(get_raw_samples(f)$rb_inclusion, colMeans))

  f1 = compare_at(1, 11)
  f2 = compare_at(2, 12)
  f2b = compare_at(2, 13)

  rw = anchor_reweight(anchor_draws(f1), s_a = 1, s_grid = c(1, 2))
  nm = get_raw_samples(f1)$parameter_names$indicator
  pairwise = !grepl("(main)", nm, fixed = TRUE)

  # Reweighting to the anchor's own scale is the identity up to Monte Carlo.
  expect_lt(max(abs(rw$pip[1, pairwise] - pip_of(f1)[pairwise])), 0.02)

  # A doubling is a real extrapolation; it stays usable and lands on the refit
  # to within a small multiple of the refit's own run-to-run spread. Without
  # this the whole curve would be reweighting an untested density.
  expect_gt(rw$ess[2], 400)
  noise = max(abs(pip_of(f2)[pairwise] - pip_of(f2b)[pairwise]))
  expect_lt(max(abs(rw$pip[2, pairwise] - pip_of(f2)[pairwise])), 4 * noise)
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
