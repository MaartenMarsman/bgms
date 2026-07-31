# Tests for Rao-Blackwellized inclusion probabilities.
#
# The sampler records, per indicator update, the one-step RB draw
#   J_t = gamma_t + (1 - 2 gamma_t) alpha_t
# with gamma_t the pre-move state and alpha_t the birth/death acceptance
# probability. Averaging J_t is a lower-variance, boundary-stable estimator of
# the inclusion probability. These tests check the plumbing, the alignment of
# the RB draws with the raw indicator draws, agreement with the raw PIP for
# well-mixed edges, and the boundary behaviour for saturated edges.
#
# The OMRF and saturated-edge fits are session-cached: several tests assert
# different properties of the same posterior, so each config is fit once. The
# Monte-Carlo-saturation boundary behaviour needs a large, well-mixed fit to
# realise machine 0/1 edges, so those two tests run in the BGMS_RUN_SLOW_TESTS
# tier; the plumbing, alignment, and well-mixed agreement stay local.

skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run the RB saturation-boundary tests"
  )
}

rb_cache = new.env(parent = emptyenv())

rb_omrf_fit = function() {
  if(is.null(rb_cache$omrf)) {
    data("Wenchuan", package = "bgms", envir = environment())
    rb_cache$omrf = bgm(
      Wenchuan[, 1:8],
      iter = 1000, warmup = 300, chains = 2, seed = 4242,
      edge_selection = TRUE, display_progress = "none"
    )
  }
  rb_cache$omrf
}

rb_saturated_fit = function() {
  if(is.null(rb_cache$saturated)) {
    data("Wenchuan", package = "bgms", envir = environment())
    rb_cache$saturated = bgm(
      Wenchuan[, 1:10],
      iter = 800, warmup = 400, chains = 2, seed = 909,
      edge_selection = TRUE, display_progress = "none"
    )
  }
  rb_cache$saturated
}

test_that("rb_inclusion draws are exposed and aligned for bgm (OMRF)", {
  fit = rb_omrf_fit()

  raw = fit$raw_samples
  expect_false(is.null(raw$rb_inclusion))
  expect_equal(length(raw$rb_inclusion), length(raw$indicator))
  # Same shape as the raw indicator draws (iterations x edges).
  expect_equal(dim(raw$rb_inclusion[[1]]), dim(raw$indicator[[1]]))

  rb_all = do.call(rbind, raw$rb_inclusion)
  # Every recorded draw is a probability in [0, 1].
  expect_true(all(is.finite(rb_all)))
  expect_true(all(rb_all >= 0 & rb_all <= 1))
})

test_that("RB draw aligns with the empirical flip behaviour per edge", {
  # For each edge, the mean RB draw restricted to iterations whose pre-move
  # state was 1 should approximate 1 minus the empirical 1 -> 0 flip rate
  # (J = 1 - alpha when gamma = 1). A scrambled edge index would break this.
  fit = rb_omrf_fit()
  raw = fit$raw_samples
  n_chains = length(raw$indicator)
  n_edges = ncol(raw$indicator[[1]])

  stay_rate = numeric(n_edges)
  mean_rb_at_1 = numeric(n_edges)
  n_at_1 = integer(n_edges)

  for(e in seq_len(n_edges)) {
    n1 = 0L
    stays = 0L
    rb_at_1 = numeric(0)
    for(c in seq_len(n_chains)) {
      ind_e = raw$indicator[[c]][, e]
      rb_e = raw$rb_inclusion[[c]][, e]
      if(length(ind_e) < 2) next
      prior_state = ind_e[-length(ind_e)]
      next_state = ind_e[-1]
      rb_next = rb_e[-1]
      at_1 = prior_state == 1
      n1 = n1 + sum(at_1)
      stays = stays + sum(at_1 & next_state == 1)
      rb_at_1 = c(rb_at_1, rb_next[at_1])
    }
    n_at_1[e] = n1
    if(n1 > 0) {
      stay_rate[e] = stays / n1
      mean_rb_at_1[e] = mean(rb_at_1)
    }
  }

  active = n_at_1 >= 50
  expect_true(sum(active) >= 2)
  # Rao-Blackwellization reduces variance, so mean(J | state 1) tracks the
  # empirical stay rate but is smoother; require a loose per-edge match.
  expect_true(max(abs(stay_rate[active] - mean_rb_at_1[active])) < 0.1)
})

test_that("RB inclusion matrix matches raw PIP for well-mixed edges", {
  fit = rb_omrf_fit()
  rb = extract_posterior_inclusion_probabilities(fit, estimator = "rb")
  pip = extract_posterior_inclusion_probabilities(fit)

  expect_equal(dim(rb), dim(pip))
  expect_equal(rownames(rb), rownames(pip))

  lt = lower.tri(rb)
  rbv = rb[lt]
  pv = pip[lt]
  # Well-mixed edges: raw PIP away from the boundary. RB and raw both estimate
  # the same quantity, so they agree within Monte Carlo error.
  mixed = pv > 0.1 & pv < 0.9
  if(any(mixed)) {
    expect_true(max(abs(rbv[mixed] - pv[mixed])) < 0.05)
  }
})

test_that("RB is interior for Monte-Carlo-saturated edges and never worse than raw", {
  skip_unless_slow()
  fit = rb_saturated_fit()
  rb = extract_posterior_inclusion_probabilities(fit, estimator = "rb")
  pip = extract_posterior_inclusion_probabilities(fit)
  lt = lower.tri(rb)
  rbv = rb[lt]
  pv = pip[lt]

  # RB never saturates more edges than the raw estimator: it converts
  # Monte-Carlo-induced 0/1 estimates into interior, finite-Bayes-factor
  # values. (Edges with overwhelming per-iteration evidence can still reach
  # machine 0/1 because alpha underflows; RB is never worse than raw.)
  n_sat_raw = sum(pv %in% c(0, 1))
  n_sat_rb = sum(rbv %in% c(0, 1))
  expect_true(n_sat_rb <= n_sat_raw)

  # Any edge whose raw PIP saturated but that RB rescued must be interior.
  rescued = (pv %in% c(0, 1)) & !(rbv %in% c(0, 1))
  if(any(rescued)) {
    expect_true(all(rbv[rescued] > 0 & rbv[rescued] < 1))
  }
})

test_that("rb_inclusion is exposed, interior, and edge-aligned for GGM", {
  set.seed(1)
  n = 200
  p = 6
  x = matrix(stats::rnorm(n * p), n, p)
  x[, 2] = x[, 1] + stats::rnorm(n, sd = 0.5) # induce one strong edge

  fit = bgm(
    x,
    variable_type = "continuous",
    iter = 800, warmup = 400, chains = 2, seed = 21,
    edge_selection = TRUE, display_progress = "none"
  )
  raw = fit$raw_samples
  expect_false(is.null(raw$rb_inclusion))
  expect_equal(dim(raw$rb_inclusion[[1]]), dim(raw$indicator[[1]]))

  rb = extract_posterior_inclusion_probabilities(fit, estimator = "rb")
  vals = rb[lower.tri(rb)]
  expect_true(all(is.finite(vals)))
  expect_true(all(vals >= 0 & vals <= 1))

  # Flip-alignment on the GGM path: its rb storage uses a row-major-with-diagonal
  # index, so a mismatch would scramble edge labels. mean(J | pre-state 1) must
  # track the empirical stay-from-1 rate per edge.
  n_chains = length(raw$indicator)
  n_edges = ncol(raw$indicator[[1]])
  stay_rate = numeric(n_edges)
  mean_rb_at_1 = numeric(n_edges)
  n_at_1 = integer(n_edges)
  for(e in seq_len(n_edges)) {
    n1 = 0L
    stays = 0L
    rb_at_1 = numeric(0)
    for(c in seq_len(n_chains)) {
      ind_e = raw$indicator[[c]][, e]
      rb_e = raw$rb_inclusion[[c]][, e]
      if(length(ind_e) < 2) next
      prior_state = ind_e[-length(ind_e)]
      at_1 = prior_state == 1
      n1 = n1 + sum(at_1)
      stays = stays + sum(at_1 & ind_e[-1] == 1)
      rb_at_1 = c(rb_at_1, rb_e[-1][at_1])
    }
    n_at_1[e] = n1
    if(n1 > 0) {
      stay_rate[e] = stays / n1
      mean_rb_at_1[e] = mean(rb_at_1)
    }
  }
  active = n_at_1 >= 50
  expect_true(sum(active) >= 1)
  expect_true(max(abs(stay_rate[active] - mean_rb_at_1[active])) < 0.1)
})

test_that("extract_inclusion_bf is finite for saturated edges and matches the RB odds", {
  skip_unless_slow()
  fit = rb_saturated_fit()
  logbf = extract_inclusion_bf(fit, log = TRUE)
  rb = extract_posterior_inclusion_probabilities(fit, estimator = "rb")
  pip = extract_posterior_inclusion_probabilities(fit)

  lt = lower.tri(logbf)
  bfv = logbf[lt]
  rbv = rb[lt]
  pv = pip[lt]

  # Every edge is proposed every sweep, so no NA. The accumulator odds are
  # finite even where the naive RB average saturates, except for edges whose
  # death acceptance underflows all the way below exp(-745): those are +Inf
  # honestly (an exactly-zero denominator), never NaN or -Inf-by-cancellation.
  expect_true(all(!is.na(bfv)))
  expect_true(all(bfv > -Inf)) # no spurious -Inf from cancellation
  expect_true(any(pv == 1)) # there ARE saturated edges in this fit

  # The accumulator rescues boundary edges the naive RB average cannot: more
  # edges get a finite Bayes factor than are strictly interior in the RB
  # probability matrix (i.e. at least one saturated edge becomes finite).
  n_interior_prob = sum(rbv > 0 & rbv < 1)
  n_finite_bf = sum(is.finite(bfv))
  expect_true(n_finite_bf > n_interior_prob)

  # Where the RB probability is interior, the BF must equal its log odds
  # (default prior is 0.5, so posterior odds == Bayes factor).
  interior = rbv > 0.02 & rbv < 0.98
  if(any(interior)) {
    expect_equal(bfv[interior], log(rbv[interior] / (1 - rbv[interior])),
      tolerance = 1e-6
    )
  }

  # The default (log = FALSE) is the Bayes factor itself. Saturated edges whose
  # log-scale value is +Inf stay +Inf on the Bayes factor scale, as do edges
  # whose finite log exceeds the double-precision ceiling of about 709.78 nats.
  bf = extract_inclusion_bf(fit)
  expect_equal(bf, exp(logbf))
  overflow = !is.na(logbf) & (logbf == Inf | logbf > 709.79)
  expect_true(all(bf[overflow] == Inf))
})

test_that("extract_inclusion_bf returns Bayes factors by default and logs on request", {
  fit = rb_omrf_fit()

  bf = extract_inclusion_bf(fit)
  logbf = extract_inclusion_bf(fit, log = TRUE)

  expect_equal(dim(bf), dim(logbf))
  expect_equal(dimnames(bf), dimnames(logbf))
  expect_equal(bf, exp(logbf))
  expect_true(all(bf[!is.na(bf)] >= 0))

  # log = TRUE is the pre-argument behaviour: the log odds of the RB inclusion
  # probability, at the default prior inclusion probability of 1/2.
  rb = extract_posterior_inclusion_probabilities(fit, estimator = "rb")
  lt = lower.tri(rb)
  interior = rb[lt] > 0.02 & rb[lt] < 0.98
  expect_true(any(interior))
  expect_equal(logbf[lt][interior],
    log(rb[lt][interior] / (1 - rb[lt][interior])),
    tolerance = 1e-6
  )

  expect_error(extract_inclusion_bf(fit, log = NA), "single logical value")
  expect_error(extract_inclusion_bf(fit, log = "yes"), "single logical value")
})

test_that("rb_bf_scale maps the boundary values as documented", {
  rb_bf_scale = bgms:::rb_bf_scale
  x = matrix(c(-Inf, Inf, NA_real_, 0, 800, log(3)), nrow = 2)

  expect_identical(rb_bf_scale(x, TRUE), x)

  bf = rb_bf_scale(x, FALSE)
  expect_identical(bf[1, 1], 0) # -Inf: no inclusion evidence remains
  expect_identical(bf[2, 1], Inf) # +Inf: no exclusion evidence remains
  expect_true(is.na(bf[1, 2]))
  expect_identical(bf[2, 2], 1)
  expect_identical(bf[1, 3], Inf) # 800 nats overflows double precision
  expect_equal(bf[2, 3], 3)
})

test_that("rb_inclusion is exposed and interior for a mixed MRF", {
  set.seed(2)
  n = 200
  cont = matrix(stats::rnorm(n * 2), n, 2)
  disc = matrix(sample(0:2, n * 2, replace = TRUE), n, 2)
  dat = data.frame(cont, disc)

  fit = bgm(
    dat,
    variable_type = c("continuous", "continuous", "ordinal", "ordinal"),
    iter = 600, warmup = 400, chains = 2, seed = 33,
    edge_selection = TRUE, display_progress = "none"
  )
  raw = fit$raw_samples
  expect_false(is.null(raw$rb_inclusion))
  expect_equal(dim(raw$rb_inclusion[[1]]), dim(raw$indicator[[1]]))

  rb = extract_posterior_inclusion_probabilities(fit, estimator = "rb")
  vals = rb[lower.tri(rb)]
  expect_true(all(is.finite(vals)))
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("extract_inclusion_bf pins the bgmCompare interleaved flattening", {
  data("Wenchuan", package = "bgms")
  x = Wenchuan[1:80, 1:5]
  group_ind = rep(1:2, each = 40)

  fit = bgmCompare(
    x = x, group_indicator = group_ind,
    iter = 600, warmup = 400, chains = 2, seed = 77,
    difference_selection = TRUE, main_difference_selection = TRUE,
    display_progress = "none"
  )

  logbf = extract_inclusion_bf(fit, log = TRUE)
  rb = extract_posterior_inclusion_probabilities(fit, estimator = "rb")

  # The bgmCompare method takes the same log argument; the default is the Bayes
  # factor scale.
  expect_equal(extract_inclusion_bf(fit), exp(logbf))

  # Same VxV shape and names as the RB probability matrix, and symmetric.
  expect_equal(dim(logbf), dim(rb))
  expect_equal(dimnames(logbf), dimnames(rb))
  expect_equal(logbf[lower.tri(logbf)], t(logbf)[lower.tri(logbf)])

  # The diagonal carries the main-effect difference Bayes factors, so the
  # interleaved (i, j >= i) flattening must place a finite value there when the
  # main difference was updated. A scrambled flatten would drop these.
  expect_true(any(is.finite(diag(logbf))))
  # No NaN anywhere (only NA / finite / +-Inf are legitimate).
  expect_false(any(is.nan(logbf)))

  # Default difference prior is Bernoulli(0.5), so the Bayes factor equals the
  # posterior odds: interior RB probabilities must match log(rb / (1 - rb)) on
  # and off the diagonal.
  interior = is.finite(rb) & rb > 0.02 & rb < 0.98
  if(any(interior)) {
    expect_equal(logbf[interior], log(rb[interior] / (1 - rb[interior])),
      tolerance = 1e-6
    )
  }
})

test_that("estimator = 'rb' errors without edge selection", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:5],
    iter = 200, warmup = 150, chains = 1, seed = 7,
    edge_selection = FALSE, display_progress = "none"
  )
  expect_error(
    extract_posterior_inclusion_probabilities(fit, estimator = "rb"),
    "edge_selection = TRUE"
  )
})

test_that("rb_inclusion is exposed for bgmCompare difference selection", {
  data("Wenchuan", package = "bgms")
  x = Wenchuan[1:80, 1:5]
  group_ind = rep(1:2, each = 40)

  fit = bgmCompare(
    x = x, group_indicator = group_ind,
    iter = 600, warmup = 400, chains = 2, seed = 13,
    difference_selection = TRUE, display_progress = "none"
  )

  raw = fit$raw_samples
  expect_false(is.null(raw$rb_inclusion))
  expect_equal(length(raw$rb_inclusion), length(raw$indicator))
  expect_equal(dim(raw$rb_inclusion[[1]]), dim(raw$indicator[[1]]))

  rb = extract_posterior_inclusion_probabilities(fit, estimator = "rb")
  pip = extract_posterior_inclusion_probabilities(fit)
  expect_equal(dim(rb), dim(pip))

  # Selected difference indicators receive finite RB estimates in [0, 1].
  finite_rb = rb[is.finite(rb)]
  expect_true(length(finite_rb) > 0)
  expect_true(all(finite_rb >= 0 & finite_rb <= 1))
})

test_that("the printed bgm summary reports the RB inclusion estimate", {
  fit = rb_omrf_fit()
  ind = summary(fit)$indicator

  # The RB precision columns sit beside the raw directional flip counts, which
  # stay retrievable so their asymmetry survives. The transition ESS is gone.
  expect_true(all(c("mean", "mcse", "sd", "n_eff", "Rhat", "n0->1", "n1->0") %in%
    colnames(ind)))
  expect_false("n_eff_mixt" %in% colnames(ind))
  # ESS/Rhat are finite or NA (constant columns), never NaN.
  expect_false(any(is.nan(ind$n_eff)))
  expect_false(any(is.nan(ind$Rhat)))
  # Directional flip counts are non-negative integers, retrievable per edge.
  expect_true(all(ind[["n0->1"]] >= 0 & ind[["n1->0"]] >= 0))

  # The summary mean equals the RB inclusion probability, edge for edge.
  rb = extract_posterior_inclusion_probabilities(fit, estimator = "rb")
  expect_equal(ind$mean, rb[lower.tri(rb)], tolerance = 1e-8)
})

test_that("the bgmCompare summary RB inclusion has NA rows for unselected mains", {
  data("Wenchuan", package = "bgms")
  x = Wenchuan[1:80, 1:5]
  group_ind = rep(1:2, each = 40)

  fit = bgmCompare(
    x = x, group_indicator = group_ind,
    iter = 600, warmup = 400, chains = 2, seed = 13,
    difference_selection = TRUE, main_difference_selection = FALSE,
    display_progress = "none"
  )
  ind = summary(fit)$indicator

  expect_true(all(c("mean", "mcse", "sd", "n_eff", "Rhat", "n0->1", "n1->0") %in%
    colnames(ind)))
  expect_false("n_eff_mixt" %in% colnames(ind))
  expect_false(any(is.nan(ind$n_eff)))
  # Unselected main-effect differences were never updated, so some rows are NA;
  # the selected pairwise differences are finite.
  expect_true(any(is.na(ind$mean)))
  expect_true(any(is.finite(ind$mean)))
})
