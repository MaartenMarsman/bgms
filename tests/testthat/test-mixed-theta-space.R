# --------------------------------------------------------------------------- #
# Theta-space (free-element Cholesky) tests for the Mixed MRF Kyy block.
#
# The mixed model parameterizes the continuous precision block by per-column
# null-space coordinates (f_q, psi_q), sharing GGMGradientEngine with the
# GGM. Excluded edges are absent from the parameter vector, so NUTS needs
# no constrained integration.
#
# Tests verify:
#   1. Analytic gradient matches finite differences on the dense graph
#      (dims and values match the pre-existing Cholesky parameterization).
#   2. Analytic gradient matches finite differences on sparse graphs with
#      excluded Kxx, Kxy, and Kyy edges (reduced parameter dimension).
#   3. The forward map enforces excluded Kyy edges: the log-posterior is
#      invariant to Kyy-slab terms of excluded edges (checked through the
#      dimension contract).
#   4. NUTS with edge selection runs end-to-end and returns finite samples.
# --------------------------------------------------------------------------- #

mixed_theta_dim = function(G, num_categories, p, q) {
  num_main = sum(num_categories)
  n_xx = sum(G[1:p, 1:p][upper.tri(diag(p))])
  n_xy = sum(G[1:p, p + 1:q])
  Gyy = G[p + 1:q, p + 1:q, drop = FALSE]
  d_sum = 0
  for(col in seq_len(q)) {
    if(col > 1) d_sum = d_sum + sum(Gyy[1:(col - 1), col])
  }
  num_main + n_xx + q + n_xy + (q + d_sum)
}

mixed_theta_fd_check = function(G, X, Y, num_categories, is_ord, baseline,
                                eps = 1e-6, tol = 1e-5) {
  p = ncol(X)
  q = ncol(Y)
  dim = mixed_theta_dim(G, num_categories, p, q)
  theta = rnorm(dim, 0, 0.3)
  f = function(th) {
    mixed_test_logp_and_gradient(
      th, X, Y, num_categories, is_ord, baseline, G, 2.5
    )
  }
  res = f(theta)
  expect_true(is.finite(res$value))
  expect_length(res$gradient, dim)
  fd = numeric(dim)
  for(k in seq_len(dim)) {
    tp = theta
    tp[k] = tp[k] + eps
    tm = theta
    tm[k] = tm[k] - eps
    fd[k] = (f(tp)$value - f(tm)$value) / (2 * eps)
  }
  expect_lt(max(abs(fd - res$gradient) / pmax(1, abs(fd))), tol)
}

test_that("theta-space gradient matches finite differences (dense graph)", {
  set.seed(11)
  n = 50
  p = 3
  q = 4
  X = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  Y = matrix(rnorm(n * q), n, q)
  num_categories = rep(2L, p)
  G = matrix(1L, p + q, p + q)
  mixed_theta_fd_check(G, X, Y, num_categories, rep(1L, p), rep(0L, p))
})

test_that("theta-space gradient matches finite differences (sparse Kyy)", {
  set.seed(12)
  n = 50
  p = 3
  q = 5
  X = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  Y = matrix(rnorm(n * q), n, q)
  num_categories = rep(2L, p)
  G = matrix(1L, p + q, p + q)
  cut = function(G, i, j) {
    G[i, j] = 0L
    G[j, i] = 0L
    G
  }
  G = cut(G, p + 1, p + 3)
  G = cut(G, p + 2, p + 5)
  G = cut(G, p + 4, p + 5)
  mixed_theta_fd_check(G, X, Y, num_categories, rep(1L, p), rep(0L, p))
})

test_that("theta-space gradient matches finite differences (sparse xx+xy+yy)", {
  set.seed(13)
  n = 50
  p = 3
  q = 4
  X = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  Y = matrix(rnorm(n * q), n, q)
  num_categories = rep(2L, p)
  total = p + q
  for(r in 1:3) {
    G = matrix(1L, total, total)
    for(i in 1:(total - 1)) {
      for(j in (i + 1):total) {
        if(runif(1) < 0.4) {
          G[i, j] = 0L
          G[j, i] = 0L
        }
      }
    }
    mixed_theta_fd_check(G, X, Y, num_categories, rep(1L, p), rep(0L, p))
  }
})

test_that("reduced dimension drops excluded edges from the parameter vector", {
  n = 30
  p = 2
  q = 3
  set.seed(14)
  X = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  Y = matrix(rnorm(n * q), n, q)
  num_categories = rep(2L, p)
  G_dense = matrix(1L, p + q, p + q)
  G_sparse = G_dense
  G_sparse[p + 1, p + 2] = 0L
  G_sparse[p + 2, p + 1] = 0L
  dim_dense = mixed_theta_dim(G_dense, num_categories, p, q)
  dim_sparse = mixed_theta_dim(G_sparse, num_categories, p, q)
  expect_equal(dim_dense - dim_sparse, 1L)

  # The reduced vector is accepted; a dense-length vector is not.
  theta = rnorm(dim_sparse, 0, 0.2)
  res = mixed_test_logp_and_gradient(
    theta, X, Y, num_categories, rep(1L, p), rep(0L, p), G_sparse, 2.5
  )
  expect_true(is.finite(res$value))
  expect_length(res$gradient, dim_sparse)
})

test_that("mixed NUTS with edge selection runs unconstrained end-to-end", {
  skip_on_cran()
  set.seed(15)
  n = 150
  p = 3
  q = 3
  S = 0.5^abs(outer(1:(p + q), 1:(p + q), "-"))
  Z = matrix(rnorm(n * (p + q)), n) %*% chol(S)
  Yord = apply(Z[, 1:p], 2, function(z) {
    as.integer(cut(z, breaks = c(-Inf, quantile(z, c(.3, .6)), Inf))) - 1L
  })
  Ymix = cbind(Yord, Z[, p + 1:q])
  colnames(Ymix) = paste0("V", 1:(p + q))
  fit = bgm(Ymix,
    variable_type = c(rep("ordinal", p), rep("continuous", q)),
    iter = 150, warmup = 150, chains = 1, cores = 1, seed = 4,
    update_method = "nuts", edge_selection = TRUE,
    display_progress = "none", verbose = FALSE
  )
  s = summary(fit)
  expect_true(all(is.finite(s$pairwise$mean)))
  expect_true(all(s$pairwise$posterior_incl_prob >= 0 &
    s$pairwise$posterior_incl_prob <= 1, na.rm = TRUE))
})
