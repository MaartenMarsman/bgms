# Tests for the refit-based prior_sensitivity_check().
#
# prior_sensitivity_check() refits the model at a grid of fixed scales and
# classifies every edge; it works on any bgm() fit with edge selection.
# Structural tests use small fits and reduced anchor grids: every refit is a
# full MCMC run, so the anchor set is the dominant cost. The warm-vs-cold
# agreement certification and the other refit cross-validations run in the
# weekly certification tier (T2, BGMS_RUN_CERTIFICATION); the verdict and
# compare-trace checks run nightly (T1, BGMS_RUN_SLOW_TESTS).

skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run prior-sensitivity certifications"
  )
}

# ---- refit-based prior_sensitivity_check() ----------------------------------

test_that("the wobble q95 pools threshold-relevant edges only", {
  # Two near-saturated edges (|log BF| > 3 log(10) at s0) with huge replicate
  # spread must not inflate the yardstick; the median and per-edge spread
  # keep covering every edge.
  lbf_s0 = log(10) * c(0.2, -1.5, 2.9, 8.0, -12.0)
  lbf_rep = log(10) * c(0.3, -1.3, 2.7, 9.3, -10.7)
  d = abs(lbf_s0 - lbf_rep)
  w = wobble_yardstick(lbf_s0, lbf_rep)
  expect_equal(w$per_edge, d)
  expect_equal(w$median, stats::median(d))
  expect_equal(w$q95, stats::quantile(d[1:3], 0.95, names = FALSE))
  expect_lt(w$q95, min(d[4:5]))
  # with every edge saturated the yardstick is undefined
  expect_true(is.na(wobble_yardstick(log(10) * c(5, -7), log(10) * c(6, -8))$q95))
  expect_equal(w$censored, 0L)
})


test_that("an edge saturating in one refit does not carry the yardstick to Inf", {
  # A threshold-relevant edge whose replicate saturates has a censored spread.
  # Pooling the Inf would make the yardstick infinite, and every verdict move
  # would then read as run-to-run noise.
  lbf_s0 = log(10) * c(0.2, -1.5, 2.9, -1.0)
  lbf_rep = log(10) * c(0.3, -1.3, 2.7, -Inf)
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
  none = wobble_yardstick(log(10) * c(-1.0, 0.5), c(-Inf, Inf))
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
  expect_equal(dim(ps$log_bf), c(length(ps$multipliers), 15L))
  # Exactness lives on the anchor fits, not the pooled curve: the 1x verdict
  # column and every chosen-scale quantity are the original fit's own RB
  # analysis. The reported PIP is the fit's RB inclusion, and the log BF is
  # its exact prior-odds transform (checked on the pip scale, since the log
  # BF derivative amplifies a 1e-10 pip agreement near saturation).
  rb = extract_posterior_inclusion_probabilities(fit)
  rbv = unname(rb[upper.tri(rb)][order_upper_tri_rowmajor(6)])
  po = ps$edges$prior_inclusion_probability /
    (1 - ps$edges$prior_inclusion_probability)
  expect_equal(ps$edges$chosen_scale_pip, rbv, tolerance = 1e-10)
  pc = ps$edges$chosen_scale_pip
  expect_equal(
    ps$edges$chosen_scale_log_bf,
    log((pc / (1 - pc)) / po),
    tolerance = 1e-12
  )
  expect_equal(
    ps$edges$verdict_x1,
    verdict_from_lbf(ps$edges$chosen_scale_log_bf, log(ps$evidence_threshold)),
    tolerance = 1e-10
  )
  # the pooled curve at 1x tracks the fit's own raw 1x inclusion proportions:
  # dominated by the 1x anchor (an identity reweight) but pooled with the
  # neighbours, so on the log BF scale it agrees on median within a modest
  # margin. The exact 1x analysis lives in the chosen-scale columns above; this
  # is a tracking check on absolute log BF, robust to the tiny per-point MCSE
  # of near-saturated edges that makes a normalised ratio platform-unstable.
  g = do.call(rbind, fit$raw_samples$indicator)
  raw_lbf = unname(log((colMeans(g) / (1 - colMeans(g))) / po))
  fin = is.finite(ps$log_bf[ps$chosen_index, ]) & is.finite(raw_lbf)
  expect_lt(
    stats::median(abs(ps$log_bf[ps$chosen_index, ] - raw_lbf)[fin]),
    0.5 * log(10)
  )
  # verdicts and movers take only the documented levels
  expect_true(all(unlist(ps$verdict) %in%
    c("presence", "undecided", "absence", NA)))
  expect_true(all(ps$edges$mover %in%
    c("stable", "indistinguishable-from-wobble", "moved-beyond-wobble")))
  # pooling keeps the curve finite even where an edge saturates at some scale
  expect_false(any(is.infinite(ps$log_bf)))
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
  # the log BF divides posterior odds by the 0.2/0.8 prior odds
  p = ps$edges$chosen_scale_pip
  expected = log((p / (1 - p)) / (0.2 / 0.8))
  expect_equal(ps$edges$chosen_scale_log_bf, expected, tolerance = 1e-8)
})


test_that("warm-started short refits agree with cold full refits within wobble", {
  skip_on_cran()
  skip_unless_certification()
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
  skip_unless_certification()
  set.seed(32)
  xg = matrix(rnorm(180 * 5), 180, 5)
  fg = bgm(xg,
    variable_type = "continuous",
    iter = 300, warmup = 250, chains = 2, seed = 8,
    # joint pinned: this smoke's anchors were recorded under the joint spec,
    # before the hierarchical default (F-010)
    precision_graph_prior = "joint",
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
  expect_true(all(is.finite(psg$log_bf[finite_pts, 1])))
  expect_true(all(is.na(psg$log_bf[!finite_pts, ])))

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
  # anchor_draws() reports tolower(spec$prior$difference_prior_type), which is
  # tolower(difference_family). The fit is at defaults, and that default is
  # "Normal" as of F-119, so the family is "normal". Exact, not a tolerance.
  expect_equal(d$family, "normal")

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


test_that("the compare noise yardstick runs per indicator, not per gated parameter", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgmCompare(
    x = Wenchuan[1:100, 1:4], group_indicator = rep(1:2, each = 50),
    iter = 200, warmup = 200, chains = 2, seed = 13,
    difference_selection = TRUE, display_progress = "none"
  )
  # 22 gated difference parameters against 10 indicators: averaging the gamma
  # draws instead of the indicator draws recycles the per-indicator prior odds
  # across the parameter columns (with a length warning) and hands the wobble
  # rule a 22-long yardstick for a 10-row edge table.
  warns = character(0)
  ps = withCallingHandlers(
    suppressMessages(prior_sensitivity_check(
      fit,
      anchors = c(0.5, 1, 2), iter = 200, warmup = 200,
      ess_floor = 50, seed = 13
    )),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("longer object length", warns)))
  expect_equal(length(ps$wobble$per_edge), nrow(ps$edges))
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


# Weekly certification (T2): eleven independent 4-chain fits.
#
# The gate asks whether the reweighted curve at a doubled difference scale
# reproduces a refit at that scale, in units of refit-to-refit spread. Every
# term in that ratio used to rest on a SINGLE fit, so the test flipped between
# seeds and went red on the 2026-08-01 nightly. Report 06's 20-seed study
# separated the two effects: a real systematic reweighting bias of about 0.01
# in inclusion probability at the extrapolation end, under a pass/fail decided
# by seed luck (4/20 tripped the x4 gate).
#
# All three terms are now pooled:
#
#   prediction   the mean over THREE anchor fits at the chosen scale. This is
#                the term that mattered: the reweighting prediction inherits
#                its anchor's Monte Carlo error whole, and pooling the other
#                two sides without this one leaves the gate seed-fragile.
#   reference    the mean of EIGHT refits at the target scale, so the
#                numerator is the reweighting bias rather than the bias plus
#                one refit's noise.
#   yardstick    four independent refit pairs, each contributing its
#                max-over-edges deviation, averaged.
#
# The gate itself is the maintainer's, unchanged: gap < 4 x pooled_noise.
#
# Measured over 8 seed bases (s0 in 11, 31, 51, 71, 91, 111, 131, 151), ratio
# = gap / pooled_noise:
#
#   1 anchor,  1 refit   (the original)   0.61 1.45 1.48 2.70 3.31 3.33 4.58 4.94
#                                         median 3.00, 2/8 over the gate
#   1 anchor,  8 refits                   0.78 1.00 1.53 2.01 2.89 3.52 3.87 5.44
#                                         median 2.45, 1/8 over the gate
#   3 anchors, 8 refits  (this block)     0.41 0.47 0.92 1.41 1.69 1.81 2.12 2.83
#                                         median 1.55, 0/8 over the gate
#
# The pooled gap over those bases is 0.0018-0.0088 in inclusion probability,
# inside the ~0.01 bound ?prior_sensitivity_check documents. Report 12 carries
# the table.
test_that("the difference-scale reweighting reproduces a refit at that scale", {
  skip_on_cran()
  skip_unless_certification()
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

  # Three anchor fits at the chosen scale and four independent refit pairs at
  # the doubled scale, all on deterministic seeds.
  anchors = lapply(c(11, 1011, 2011), function(s) compare_at(1, s))
  pair_seeds = list(c(12, 13), c(14, 15), c(16, 17), c(18, 19))
  pips = lapply(pair_seeds, function(s) {
    list(
      pip_of(compare_at(2, s[1])),
      pip_of(compare_at(2, s[2]))
    )
  })

  nm = get_raw_samples(anchors[[1]])$parameter_names$indicator
  pairwise = !grepl("(main)", nm, fixed = TRUE)
  npw = sum(pairwise)
  rws = lapply(anchors, function(f) {
    anchor_reweight(anchor_draws(f), s_a = 1, s_grid = c(1, 2))
  })

  # Reweighting to an anchor's own scale is the identity up to Monte Carlo, and
  # that has to hold for every anchor the prediction pools.
  for(i in seq_along(anchors)) {
    expect_lt(
      max(abs(rws[[i]]$pip[1, pairwise] - pip_of(anchors[[i]])[pairwise])), 0.02
    )
    expect_gt(rws[[i]]$ess[2], 400)
  }

  # Yardstick: four independent two-refit deviations, pooled.
  pooled_noise = mean(vapply(
    pips, function(p) max(abs(p[[1]][pairwise] - p[[2]][pairwise])), numeric(1)
  ))
  # Prediction and reference, each pooled over its own fits.
  prediction = rowMeans(vapply(rws, function(r) r$pip[2, pairwise], numeric(npw)))
  reference = rowMeans(vapply(
    unlist(pips, recursive = FALSE), function(p) p[pairwise], numeric(npw)
  ))

  # A doubling is a real extrapolation; it stays usable and lands on the refit
  # to within a small multiple of the refit's own run-to-run spread. Without
  # this the whole curve would be reweighting an untested density.
  gap = max(abs(prediction - reference))
  expect_lt(gap, 4 * pooled_noise)
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

test_that("the refit gate abstains on an all-NA warmup check without warning", {
  # A degenerate source fit leaves every energy diagnostic NA. min()/max() over
  # nothing warns and returns an infinity, which reads as a criterion rather
  # than as an absent one.
  expect_equal(finite_reduce(c(NA_real_, NA_real_), min, Inf), Inf)
  expect_equal(finite_reduce(numeric(0), max, 0), 0)
  expect_equal(finite_reduce(c(NA_real_, 2, 5), max, 0), 5)
  expect_identical(finite_reduce(c(NA_real_, NaN), stats::median, NA_real_), NA_real_)
  expect_silent(finite_reduce(rep(NA_real_, 3L), min, Inf))
  expect_silent(finite_reduce(rep(NA_real_, 3L), max, 0))

  # Infinities are not finite values to reduce either: an unusable diagnostic
  # abstains whichever way it is unusable.
  expect_equal(finite_reduce(c(-Inf, Inf), min, Inf), Inf)
})


# ------------------------------------------------------------------------------
# Coherence between the reported verdict and the stability interval
# ------------------------------------------------------------------------------

# A minimal bgms_prior_sensitivity object, enough for print(). Only the fields
# print() reads are filled; the curve machinery is not re-run.
fake_sensitivity = function(chosen_verdict, curve_verdict, multipliers,
                            chosen_idx, anchors) {
  n_anchor = length(anchors)
  stab = stability_interval(
    curve_verdict, multipliers, chosen_idx,
    target = chosen_verdict
  )
  edges = data.frame(
    edge = "A-B",
    prior_inclusion_probability = 0.5,
    chosen_scale_pip = 0.6,
    chosen_scale_log_bf = 0.4,
    chosen_scale_mcse = 0.01,
    chosen_scale_verdict = chosen_verdict,
    stability_lower = stab[1],
    stability_upper = stab[2],
    mover = "stable",
    insufficient = FALSE,
    insufficient_noisy = FALSE,
    insufficient_disagree = FALSE,
    saturated = FALSE,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      edges = edges,
      grid = data.frame(
        multiplier = anchors, replicate = rep(FALSE, n_anchor),
        original_fit = multipliers[chosen_idx] == anchors,
        usable = rep(TRUE, n_anchor), forced = rep(FALSE, n_anchor)
      ),
      anchors = anchors,
      multipliers = multipliers,
      chosen_scale = 1,
      chosen_index = chosen_idx,
      anchor_verdict = matrix(chosen_verdict, n_anchor, 1L),
      curve = list(ess_floor = 400),
      wobble = list(q95 = 0.3, anchor = anchors[n_anchor], censored = 0L),
      preferred_scale = list(s_hat = NA_real_),
      unit = list(
        noun = "edge", nouns = "edges",
        headline = "are the edge verdicts robust to the slab scale?"
      ),
      vary = list(mode = "none"),
      evidence_threshold = 10,
      tolerance = 0.5 * log(10),
      refit_sampler = "nuts",
      warm = TRUE,
      runtime_seconds = 3
    ),
    class = "bgms_prior_sensitivity"
  )
}


test_that("the stability interval is anchored on the verdict the table reports", {
  # COHERENCE GATE. The table's chosen_scale_verdict is the original fit's own
  # Rao-Blackwellized verdict; the pooled curve is a different estimator of the
  # same quantity and can disagree at 1x on a borderline edge. Taking the
  # interval's target from the curve made the bounds describe a verdict the
  # table never showed.
  multipliers = c(0.5, 0.7, 1, 1.4, 2)
  chosen_idx = 3L

  # Agreement: unchanged behavior, the whole grid.
  agree = rep("presence", 5L)
  expect_equal(
    stability_interval(agree, multipliers, chosen_idx, target = "presence"),
    c(0.5, 2)
  )
  expect_equal(
    stability_interval(agree, multipliers, chosen_idx),
    stability_interval(agree, multipliers, chosen_idx, target = "presence")
  )

  # Disagreement at 1x: the reported verdict holds nowhere around the chosen
  # scale, so the bounds are missing rather than describing the curve's own
  # verdict across the full range.
  disagree = rep("undecided", 5L)
  expect_equal(
    stability_interval(disagree, multipliers, chosen_idx, target = "presence"),
    c(NA_real_, NA_real_)
  )
  # What the old, curve-anchored call would have reported.
  expect_equal(
    stability_interval(disagree, multipliers, chosen_idx),
    c(0.5, 2)
  )

  # Partial agreement still walks out from the chosen scale only.
  mixed = c("undecided", "presence", "presence", "presence", "absence")
  expect_equal(
    stability_interval(mixed, multipliers, chosen_idx, target = "presence"),
    c(0.7, 1.4)
  )

  # print() describes the same verdict the table does: with the two estimators
  # disagreeing, the edge is an exception, not a verdict that "holds".
  ps = fake_sensitivity("presence", disagree, multipliers, chosen_idx,
    anchors = c(0.5, 1, 2)
  )
  expect_true(is.na(ps$edges$stability_lower))
  expect_true(is.na(ps$edges$stability_upper))
  expect_equal(ps$edges$chosen_scale_verdict, "presence")
  out = paste(utils::capture.output(print(ps)), collapse = "\n")
  expect_match(out, "0 of 1 verdicts hold across the whole")
  expect_false(grepl("All 1 verdicts hold", out, fixed = TRUE))

  # And with the two agreeing it is still reported as holding throughout.
  ok = fake_sensitivity("presence", agree, multipliers, chosen_idx,
    anchors = c(0.5, 1, 2)
  )
  expect_equal(unname(unlist(ok$edges[c("stability_lower", "stability_upper")])), c(0.5, 2))
  expect_match(
    paste(utils::capture.output(print(ok)), collapse = "\n"),
    "All 1 verdicts hold"
  )
})


test_that("a forced 1x anchor keeps the gate's own record in the grid", {
  # The curve still uses the original fit when it fails its gate (it is the
  # analysis under check), but $grid$usable must stay the gate's verdict and
  # $grid$forced must say the curve overrode it. print() then names it as kept,
  # not as excluded.
  ps = fake_sensitivity("presence", rep("presence", 5L), c(0.5, 0.7, 1, 1.4, 2),
    3L,
    anchors = c(0.5, 1, 2)
  )
  ps$grid$usable[2] = FALSE
  ps$grid$forced[2] = TRUE
  out = paste(utils::capture.output(print(ps)), collapse = "\n")
  expect_match(out, "did not pass its convergence check")
  expect_match(out, "still reported")
  expect_false(grepl("excluded from the verdicts", out, fixed = TRUE))

  # A genuinely excluded refit still reads as excluded.
  ps2 = ps
  ps2$grid$forced[2] = FALSE
  out2 = paste(utils::capture.output(print(ps2)), collapse = "\n")
  expect_match(out2, "did not converge and is excluded from the verdicts")
})


test_that("a fit whose swept prior carries no scale errors before any refit", {
  # A beta-prime slab is scale-free, so bgm() records pairwise_scale = NA and
  # the whole multiplier grid is NA. The old code ran every refit first and
  # failed afterwards.
  data("Wenchuan", package = "bgms")
  fit = suppressWarnings(bgm(Wenchuan[, 1:4],
    interaction_prior = beta_prime_prior(),
    iter = 200, warmup = 200, chains = 1, seed = 11,
    update_method = "adaptive-metropolis", display_progress = "none"
  ))
  expect_true(is.na(get_fit_spec(fit)$prior$pairwise_scale))
  t0 = Sys.time()
  expect_error(
    prior_sensitivity_check(fit),
    "no usable pairwise_scale"
  )
  # Fast means fast: no refit had time to run.
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 5)
})


test_that("prior_sensitivity_check runs end to end on a single-indicator fit", {
  skip_on_cran()
  # A two-variable network has one edge, and apply()/vapply() over a one-column
  # matrix collapse to a vector; the curve and the verdict columns lost their
  # grid-by-edge shape.
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:2],
    chains = 2, iter = 300, warmup = 300, seed = 31,
    display_progress = "none"
  )
  ps = suppressWarnings(suppressMessages(prior_sensitivity_check(fit,
    anchors = c(0.5, 1, 2), iter = 300, warmup = 300, ess_floor = 50,
    seed = 31
  )))
  expect_s3_class(ps, "bgms_prior_sensitivity")
  expect_equal(nrow(ps$edges), 1L)
  expect_equal(dim(ps$log_bf), c(length(ps$multipliers), 1L))
  expect_equal(dim(ps$log_bf_mcse), c(length(ps$multipliers), 1L))
  expect_equal(dim(ps$verdict), c(length(ps$multipliers), 1L))
  expect_equal(dim(ps$anchor_verdict), c(length(ps$anchors), 1L))
  expect_true(all(paste0("verdict_x", ps$anchors) %in% names(ps$edges)))
  # The gate record and the forced column travel with the grid. Whether this
  # fit's own draws clear the gate is a property of the sampler on the machine
  # running the test, not of the reshape under test, so what is asserted here
  # is the structure of the record rather than its verdict: only the original
  # fit can ever be forced, and a forced row carries the gate's own rejection
  # rather than the override. The override's behavior is pinned exactly, on a
  # synthetic grid, in "a forced 1x anchor keeps the gate's own record in the
  # grid".
  expect_true("forced" %in% names(ps$grid))
  expect_false(any(ps$grid$forced & !ps$grid$original_fit))
  expect_true(all(!ps$grid$usable[ps$grid$forced]))
  expect_output(print(ps), "are the edge verdicts robust")
})
