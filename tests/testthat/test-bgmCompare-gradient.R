# --------------------------------------------------------------------------- #
# Tests for the bgmCompare log-pseudoposterior and gradient (F-112).
#
# The compare path had no [[Rcpp::export]] at all, so its log-posterior could
# not be called from R and could not be checked against anything.
# bgmCompare_test_logp_and_gradient() exposes it, mirroring
# ggm_test_logp_and_gradient() and mixed_test_logp_and_gradient().
#
# Two independent checks:
#   (a) the VALUE against a from-scratch R implementation of the multi-group
#       pseudolikelihood plus priors, written from the model definition rather
#       than from the C++ (compare_reference_logp() below);
#   (b) the GRADIENT against central finite differences of the C++ value.
#
# Both are run on a group-differing-support dataset under the union semantics
# of F-075 --- a retained category with zero observations in one group --- to
# check that a structural zero does not break either.
# --------------------------------------------------------------------------- #

# ---- Building the C++ call's arguments from a data set ---------------------- #

compare_setup = function(x, group, num_categories, inclusion_indicator = NULL) {
  num_variables = ncol(x)
  num_groups = length(unique(group))
  is_ordinal = rep(1L, num_variables)
  baseline_category = rep(0L, num_variables)

  # Group-sorted observations and the group row ranges the C++ expects.
  ord = order(group)
  obs = x[ord, , drop = FALSE]
  g_sorted = group[ord]
  group_indices = t(vapply(seq_len(num_groups), function(g) {
    rows = which(g_sorted == g)
    c(min(rows) - 1L, max(rows) - 1L)
  }, numeric(2)))
  storage.mode(group_indices) = "integer"

  counts = lapply(seq_len(num_groups), function(g) {
    m = matrix(0L, nrow = max(num_categories), ncol = num_variables)
    for(v in seq_len(num_variables)) {
      for(cat in seq_len(num_categories[v])) {
        m[cat, v] = sum(obs[g_sorted == g, v] == cat)
      }
    }
    m
  })
  bc_stats = lapply(seq_len(num_groups), function(g) {
    matrix(0L, nrow = 2, ncol = num_variables)
  })
  pair_stats = lapply(seq_len(num_groups), function(g) {
    o = obs[g_sorted == g, , drop = FALSE]
    t(o) %*% o
  })

  main_effect_indices = matrix(NA_integer_, num_variables, 2)
  for(v in seq_len(num_variables)) {
    main_effect_indices[v, 1] = if(v > 1) {
      1L + main_effect_indices[v - 1, 2]
    } else {
      0L
    }
    main_effect_indices[v, 2] = main_effect_indices[v, 1] +
      num_categories[v] - 1L
  }

  pairwise_effect_indices = matrix(NA_integer_, num_variables, num_variables)
  tel = 0L
  for(v1 in seq_len(num_variables - 1)) {
    for(v2 in seq(v1 + 1, num_variables)) {
      pairwise_effect_indices[v1, v2] = tel
      pairwise_effect_indices[v2, v1] = tel
      tel = tel + 1L
    }
  }

  if(is.null(inclusion_indicator)) {
    inclusion_indicator = matrix(1L, num_variables, num_variables)
  }

  one = matrix(1, num_groups, num_groups)
  V = diag(num_groups) - one / num_groups
  projection = eigen(V)$vectors[, -num_groups, drop = FALSE]
  if(num_groups == 2) projection = projection / sqrt(2)

  list(
    observations = matrix(as.integer(obs), nrow(obs), ncol(obs)),
    group_indices = group_indices,
    num_groups = num_groups,
    counts_per_category = counts,
    blume_capel_stats = bc_stats,
    pairwise_stats = pair_stats,
    num_categories = as.integer(num_categories),
    is_ordinal_variable = as.integer(is_ordinal),
    baseline_category = baseline_category,
    main_effect_indices = main_effect_indices,
    pairwise_effect_indices = pairwise_effect_indices,
    inclusion_indicator = inclusion_indicator,
    projection = projection
  )
}

compare_call = function(params, s, ...) {
  bgmCompare_test_logp_and_gradient(
    params, s$observations, s$group_indices, s$num_groups,
    s$counts_per_category, s$blume_capel_stats, s$pairwise_stats,
    s$num_categories, s$is_ordinal_variable, s$baseline_category,
    s$main_effect_indices, s$pairwise_effect_indices,
    s$inclusion_indicator, s$projection, ...
  )
}

# Number of active parameters in the flat vector, i.e. its required length.
compare_param_length = function(s) {
  nm = sum(s$num_categories)
  np = ncol(s$observations) * (ncol(s$observations) - 1) / 2
  nmd = sum(s$num_categories[diag(s$inclusion_indicator) == 1L]) *
    (s$num_groups - 1)
  npd = sum(s$inclusion_indicator[upper.tri(s$inclusion_indicator)]) *
    (s$num_groups - 1)
  nm + np + nmd + npd
}

# ---- An independent R reference for the log-pseudoposterior ----------------- #
#
# Written from the model: for group g with contrast row proj_g, the
# group-specific threshold of variable v at category c is
#   mu[v, c, g] = M[v, c, 1] + sum_k proj_g[k] * M[v, c, k + 1]
# and likewise for the pairwise weights. The pseudolikelihood is the sum over
# persons and variables of the ordinal conditional
#   mu[v, x_iv, g] + 2 * x_iv * sum_{u != v} w[v, u, g] * x_iu - log Z_iv
# with mu[v, 0, g] = 0. Priors: threshold_prior on the overall thresholds,
# interaction_prior on the overall pairwise effects, difference_prior on every
# active difference.
compare_reference_logp = function(params, s,
                                  pairwise_scale = 1, difference_scale = 1,
                                  main_alpha = 1, main_beta = 1) {
  V = ncol(s$observations)
  G = s$num_groups
  nc = s$num_categories
  nm = sum(nc)
  np = V * (V - 1) / 2

  # --- unpack the flat vector into (row, group) matrices ---
  M = matrix(0, nm, G)
  P = matrix(0, np, G)
  M[, 1] = params[seq_len(nm)]
  P[, 1] = params[nm + seq_len(np)]
  off = nm + np
  for(v in seq_len(V)) {
    if(s$inclusion_indicator[v, v] != 1L) next
    rows = (s$main_effect_indices[v, 1] + 1):(s$main_effect_indices[v, 2] + 1)
    for(r in rows) {
      for(g in 2:G) {
        off = off + 1
        M[r, g] = params[off]
      }
    }
  }
  for(v1 in seq_len(V - 1)) {
    for(v2 in seq(v1 + 1, V)) {
      if(s$inclusion_indicator[v1, v2] != 1L) next
      r = s$pairwise_effect_indices[v1, v2] + 1
      for(g in 2:G) {
        off = off + 1
        P[r, g] = params[off]
      }
    }
  }

  lp = 0
  for(g in seq_len(G)) {
    proj = s$projection[g, ]
    rows = (s$group_indices[g, 1] + 1):(s$group_indices[g, 2] + 1)
    xg = s$observations[rows, , drop = FALSE]

    # group-specific thresholds and weights
    mu = matrix(0, V, max(nc))
    for(v in seq_len(V)) {
      r0 = s$main_effect_indices[v, 1] + 1
      for(cat in seq_len(nc[v])) {
        mu[v, cat] = M[r0 + cat - 1, 1] +
          if(s$inclusion_indicator[v, v] == 1L) {
            sum(proj * M[r0 + cat - 1, -1])
          } else {
            0
          }
      }
    }
    W = matrix(0, V, V)
    for(v1 in seq_len(V - 1)) {
      for(v2 in seq(v1 + 1, V)) {
        r = s$pairwise_effect_indices[v1, v2] + 1
        w = P[r, 1] +
          if(s$inclusion_indicator[v1, v2] == 1L) sum(proj * P[r, -1]) else 0
        W[v1, v2] = w
        W[v2, v1] = w
      }
    }

    for(i in seq_len(nrow(xg))) {
      for(v in seq_len(V)) {
        rest = 2 * sum(W[v, ] * xg[i, ])
        scores = 0:nc[v]
        lin = c(0, mu[v, seq_len(nc[v])]) + scores * rest
        lp = lp + lin[xg[i, v] + 1] - (max(lin) + log(sum(exp(lin - max(lin)))))
      }
    }
  }

  # --- priors ---
  for(v in seq_len(V)) {
    r0 = s$main_effect_indices[v, 1] + 1
    for(cat in seq_len(nc[v])) {
      xv = M[r0 + cat - 1, 1]
      lp = lp + xv * main_alpha - log1p(exp(xv)) * (main_alpha + main_beta)
      if(s$inclusion_indicator[v, v] == 1L) {
        lp = lp + sum(dcauchy(M[r0 + cat - 1, -1], 0, difference_scale, TRUE))
      }
    }
  }
  for(v1 in seq_len(V - 1)) {
    for(v2 in seq(v1 + 1, V)) {
      r = s$pairwise_effect_indices[v1, v2] + 1
      lp = lp + dcauchy(P[r, 1], 0, pairwise_scale, TRUE)
      if(s$inclusion_indicator[v1, v2] == 1L) {
        lp = lp + sum(dcauchy(P[r, -1], 0, difference_scale, TRUE))
      }
    }
  }

  lp
}

compare_fd_gradient = function(params, s, eps = 1e-5) {
  vapply(seq_along(params), function(k) {
    pp = params
    pm = params
    pp[k] = pp[k] + eps
    pm[k] = pm[k] - eps
    (compare_call(pp, s)$value - compare_call(pm, s)$value) / (2 * eps)
  }, numeric(1))
}

compare_max_rel_error = function(analytic, fd) {
  max(abs(analytic - fd) / pmax(abs(analytic), abs(fd), 1))
}


# ---- Tests ------------------------------------------------------------------ #

test_that("compare log-posterior matches an independent R reference", {
  set.seed(11)
  n = 40
  V = 3
  x = matrix(sample(0:2, 2 * n * V, replace = TRUE), 2 * n, V)
  group = rep(1:2, each = n)
  s = compare_setup(x, group, num_categories = rep(2L, V))

  set.seed(12)
  params = rnorm(compare_param_length(s), sd = 0.4)

  got = compare_call(params, s)$value
  want = compare_reference_logp(params, s)
  expect_equal(got, want, tolerance = 1e-8)
})

test_that("compare gradient matches central finite differences", {
  set.seed(11)
  n = 40
  V = 3
  x = matrix(sample(0:2, 2 * n * V, replace = TRUE), 2 * n, V)
  group = rep(1:2, each = n)
  s = compare_setup(x, group, num_categories = rep(2L, V))

  set.seed(12)
  params = rnorm(compare_param_length(s), sd = 0.4)

  analytic = compare_call(params, s)$gradient
  expect_length(analytic, length(params))
  expect_lt(compare_max_rel_error(analytic, compare_fd_gradient(params, s)), 1e-5)
})

test_that("inactive differences drop out of the parameter vector", {
  set.seed(13)
  n = 30
  V = 3
  x = matrix(sample(0:2, 2 * n * V, replace = TRUE), 2 * n, V)
  group = rep(1:2, each = n)

  ind = matrix(1L, V, V)
  ind[1, 2] = ind[2, 1] = 0L # this edge's difference is excluded
  ind[3, 3] = 0L # this variable's threshold differences are excluded
  s = compare_setup(x, group, num_categories = rep(2L, V), inclusion_indicator = ind)

  set.seed(14)
  params = rnorm(compare_param_length(s), sd = 0.4)

  got = compare_call(params, s)
  expect_equal(got$value, compare_reference_logp(params, s), tolerance = 1e-8)
  expect_lt(compare_max_rel_error(got$gradient, compare_fd_gradient(params, s)), 1e-5)

  # A vector of the wrong length is rejected rather than read past its end.
  expect_error(compare_call(params[-1], s), "active parameter vector")
})

test_that("a category with no observations in one group is handled", {
  # Union semantics (F-075): group 1 spans {0,1,2}, group 2 spans {1,2,3}, so
  # all four categories are retained and each group has one empty cell. The
  # log-posterior and its gradient must still agree with the reference.
  g1 = cbind(rep(0:2, each = 12), rep(0:3, times = 9), rep(0:2, times = 12))
  g2 = cbind(rep(1:3, each = 12), rep(0:3, times = 9), rep(0:2, times = 12))
  x = rbind(g1, g2)
  group = rep(1:2, each = 36)
  s = compare_setup(x, group, num_categories = c(3L, 3L, 2L))

  # The structural zeros are real: category 3 is empty in group 1, category 0
  # in group 2.
  expect_equal(s$counts_per_category[[1]][3, 1], 0L)
  expect_equal(sum(s$observations[seq_len(36), 1] == 0), 12L)

  set.seed(15)
  params = rnorm(compare_param_length(s), sd = 0.4)

  got = compare_call(params, s)
  expect_true(is.finite(got$value))
  expect_true(all(is.finite(got$gradient)))
  expect_equal(got$value, compare_reference_logp(params, s), tolerance = 1e-8)
  expect_lt(compare_max_rel_error(got$gradient, compare_fd_gradient(params, s)), 1e-5)
})

test_that("the hook works for three groups", {
  set.seed(16)
  n = 25
  V = 3
  x = matrix(sample(0:2, 3 * n * V, replace = TRUE), 3 * n, V)
  group = rep(1:3, each = n)
  s = compare_setup(x, group, num_categories = rep(2L, V))

  set.seed(17)
  params = rnorm(compare_param_length(s), sd = 0.4)

  got = compare_call(params, s)
  expect_equal(got$value, compare_reference_logp(params, s), tolerance = 1e-8)
  expect_lt(compare_max_rel_error(got$gradient, compare_fd_gradient(params, s)), 1e-5)
})
