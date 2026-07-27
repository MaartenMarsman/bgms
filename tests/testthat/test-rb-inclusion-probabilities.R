# Tests for Rao-Blackwellized inclusion probabilities.
#
# The sampler records, per indicator update, the one-step RB draw
#   J_t = gamma_t + (1 - 2 gamma_t) alpha_t
# with gamma_t the pre-move state and alpha_t the birth/death acceptance
# probability. Averaging J_t is a lower-variance, boundary-stable estimator of
# the inclusion probability. These tests check the plumbing, the alignment of
# the RB draws with the raw indicator draws, agreement with the raw PIP for
# well-mixed edges, and the boundary behaviour for saturated edges.

test_that("rb_inclusion draws are exposed and aligned for bgm (OMRF)", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:8],
    iter = 2000, warmup = 500, chains = 2, seed = 20260727,
    edge_selection = TRUE, display_progress = "none"
  )

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
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:8],
    iter = 2000, warmup = 500, chains = 2, seed = 11,
    edge_selection = TRUE, display_progress = "none"
  )
  raw = fit$raw_samples
  n_chains = length(raw$indicator)
  n_edges = ncol(raw$indicator[[1]])

  stay_rate = numeric(n_edges)
  mean_rb_at_1 = numeric(n_edges)
  n_at_1 = integer(n_edges)

  for(e in seq_len(n_edges)) {
    n1 = 0L; stays = 0L; rb_at_1 = numeric(0)
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
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:8],
    iter = 2000, warmup = 500, chains = 2, seed = 4242,
    edge_selection = TRUE, display_progress = "none"
  )
  rb = extract_rb_inclusion_probabilities(fit)
  pip = extract_posterior_inclusion_probabilities(fit)

  expect_equal(dim(rb), dim(pip))
  expect_equal(rownames(rb), rownames(pip))

  lt = lower.tri(rb)
  rbv = rb[lt]; pv = pip[lt]
  # Well-mixed edges: raw PIP away from the boundary. RB and raw both estimate
  # the same quantity, so they agree within Monte Carlo error.
  mixed = pv > 0.1 & pv < 0.9
  if(any(mixed)) {
    expect_true(max(abs(rbv[mixed] - pv[mixed])) < 0.05)
  }
})

test_that("RB is interior for Monte-Carlo-saturated edges and never worse than raw", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:12],
    iter = 2000, warmup = 500, chains = 2, seed = 909,
    edge_selection = TRUE, display_progress = "none"
  )
  rb = extract_rb_inclusion_probabilities(fit)
  pip = extract_posterior_inclusion_probabilities(fit)
  lt = lower.tri(rb)
  rbv = rb[lt]; pv = pip[lt]

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

test_that("rb_inclusion is exposed and interior for GGM (continuous)", {
  set.seed(1)
  n = 200; p = 6
  x = matrix(stats::rnorm(n * p), n, p)
  x[, 2] = x[, 1] + stats::rnorm(n, sd = 0.5)   # induce one strong edge

  fit = bgm(
    x, variable_type = "continuous",
    iter = 1500, warmup = 500, chains = 2, seed = 21,
    edge_selection = TRUE, display_progress = "none"
  )
  raw = fit$raw_samples
  expect_false(is.null(raw$rb_inclusion))
  expect_equal(dim(raw$rb_inclusion[[1]]), dim(raw$indicator[[1]]))

  rb = extract_rb_inclusion_probabilities(fit)
  vals = rb[lower.tri(rb)]
  expect_true(all(is.finite(vals)))
  expect_true(all(vals >= 0 & vals <= 1))
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
    iter = 1500, warmup = 500, chains = 2, seed = 33,
    edge_selection = TRUE, display_progress = "none"
  )
  raw = fit$raw_samples
  expect_false(is.null(raw$rb_inclusion))
  expect_equal(dim(raw$rb_inclusion[[1]]), dim(raw$indicator[[1]]))

  rb = extract_rb_inclusion_probabilities(fit)
  vals = rb[lower.tri(rb)]
  expect_true(all(is.finite(vals)))
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("extract_rb_inclusion_probabilities errors without edge selection", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[, 1:5],
    iter = 500, warmup = 200, chains = 1, seed = 7,
    edge_selection = FALSE, display_progress = "none"
  )
  expect_error(
    extract_rb_inclusion_probabilities(fit),
    "edge_selection = TRUE"
  )
})

test_that("rb_inclusion is exposed for bgmCompare difference selection", {
  data("Wenchuan", package = "bgms")
  x = Wenchuan[1:80, 1:5]
  group_ind = rep(1:2, each = 40)

  fit = bgmCompare(
    x = x, group_indicator = group_ind,
    iter = 1500, warmup = 500, chains = 2, seed = 13,
    difference_selection = TRUE, display_progress = "none"
  )

  raw = fit$raw_samples
  expect_false(is.null(raw$rb_inclusion))
  expect_equal(length(raw$rb_inclusion), length(raw$indicator))
  expect_equal(dim(raw$rb_inclusion[[1]]), dim(raw$indicator[[1]]))

  rb = extract_rb_inclusion_probabilities(fit)
  pip = extract_posterior_inclusion_probabilities(fit)
  expect_equal(dim(rb), dim(pip))

  # Selected difference indicators receive finite RB estimates in [0, 1].
  finite_rb = rb[is.finite(rb)]
  expect_true(length(finite_rb) > 0)
  expect_true(all(finite_rb >= 0 & finite_rb <= 1))
})
