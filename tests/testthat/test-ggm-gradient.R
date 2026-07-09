# --------------------------------------------------------------------------- #
# Tests for the GGM free-element Cholesky gradient engine.
#
# Validates the C++ logp_and_gradient against central finite differences
# across a grid of test configurations.
# --------------------------------------------------------------------------- #

# ---- Helpers ----------------------------------------------------------------

make_edge_matrix = function(p, included_edges) {
  E = matrix(0L, nrow = p, ncol = p)
  if(length(included_edges) > 0) {
    for(k in seq_along(included_edges)) {
      ij = included_edges[[k]]
      i = ij[1]
      j = ij[2]
      E[i, j] = 1L
      E[j, i] = 1L
    }
  }
  E
}

theta_dim = function(edge_mat) {
  p = nrow(edge_mat)
  # p diagonals + number of included edges
  p + sum(edge_mat[upper.tri(edge_mat)] == 1L)
}

fd_gradient = function(theta, suf_stat, n, edge_mat, pairwise_scale, eps = 1e-6) {
  g = numeric(length(theta))
  for(k in seq_along(theta)) {
    t_plus = theta
    t_minus = theta
    t_plus[k] = t_plus[k] + eps
    t_minus[k] = t_minus[k] - eps
    f_plus = ggm_test_logp_and_gradient(t_plus, suf_stat, n, edge_mat, pairwise_scale)$value
    f_minus = ggm_test_logp_and_gradient(t_minus, suf_stat, n, edge_mat, pairwise_scale)$value
    g[k] = (f_plus - f_minus) / (2 * eps)
  }
  g
}

check_gradient = function(p, edge_mat, n = 200, seed = 1, pairwise_scale = 1,
                          eps = 1e-6, tol = 1e-4) {
  set.seed(seed)
  X = matrix(rnorm(n * p), nrow = n, ncol = p)
  S = t(X) %*% X

  d = theta_dim(edge_mat)
  theta = rnorm(d, sd = 0.2)

  ag = ggm_test_logp_and_gradient(theta, S, n, edge_mat, pairwise_scale)
  fd = fd_gradient(theta, S, n, edge_mat, pairwise_scale, eps = eps)

  denom = pmax(abs(ag$gradient), abs(fd), 1)
  rel_err = abs(ag$gradient - fd) / denom
  max(rel_err)
}


# ---- Forward map tests ------------------------------------------------------ #

test_that("forward map produces valid upper-triangular Phi and symmetric K", {
  p = 4
  all_edges = list()
  idx = 1L
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      all_edges[[idx]] = c(i, j)
      idx = idx + 1L
    }
  }
  edge_mat = make_edge_matrix(p, all_edges)
  d = theta_dim(edge_mat)
  set.seed(42)
  theta = rnorm(d, sd = 0.3)

  result = ggm_test_forward_map(theta, edge_mat)
  Phi = result$Phi
  K = result$K

  # Phi is upper triangular
  expect_true(all(abs(Phi[lower.tri(Phi)]) < 1e-14))
  # Phi has positive diagonal
  expect_true(all(diag(Phi) > 0))
  # K = Phi^T Phi
  expect_equal(K, t(Phi) %*% Phi, tolerance = 1e-12)
  # K is symmetric
  expect_equal(K, t(K), tolerance = 1e-14)
  # K is positive definite
  expect_true(all(eigen(K, only.values = TRUE)$values > 0))
})

test_that("forward map enforces zero entries for excluded edges", {
  p = 4
  edge_mat = make_edge_matrix(p, list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(3, 4)))
  # Edge (2,4) is excluded
  d = theta_dim(edge_mat)
  set.seed(42)
  theta = rnorm(d, sd = 0.3)

  result = ggm_test_forward_map(theta, edge_mat)
  K = result$K

  # K[2,4] and K[4,2] should be zero (1-based indexing)
  expect_equal(K[2, 4], 0, tolerance = 1e-12)
  expect_equal(K[4, 2], 0, tolerance = 1e-12)
})

test_that("Jacobian matches analytical formula for complete graph", {
  p = 4
  all_edges = list()
  idx = 1L
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      all_edges[[idx]] = c(i, j)
      idx = idx + 1L
    }
  }
  edge_mat = make_edge_matrix(p, all_edges)
  d = theta_dim(edge_mat)
  set.seed(42)
  theta = rnorm(d, sd = 0.3)

  result = ggm_test_forward_map(theta, edge_mat)
  psi = result$psi

  # For complete graph: log|det J| = p*log(2) + sum_q (2 + (p-q)) * psi_q
  # (0-based: 2 + (p-1-q) for q < p-1, and 2 for q = p-1)
  expected = p * log(2) + sum((p - seq_len(p) + 2) * psi)

  expect_equal(result$log_det_jacobian, expected, tolerance = 1e-10)
})


# ---- Gradient validation tests ----------------------------------------------- #

test_that("gradient matches FD for complete graph (p=4)", {
  p = 4
  all_edges = list()
  idx = 1L
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      all_edges[[idx]] = c(i, j)
      idx = idx + 1L
    }
  }
  edge_mat = make_edge_matrix(p, all_edges)
  max_err = check_gradient(p = 4, edge_mat = edge_mat)
  expect_lt(max_err, 1e-4)
})

test_that("gradient matches FD for p=4 with missing edge (2,4)", {
  edge_mat = make_edge_matrix(
    4,
    list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(3, 4))
  )
  max_err = check_gradient(p = 4, edge_mat = edge_mat)
  expect_lt(max_err, 1e-4)
})

test_that("gradient matches FD for p=4 with empty graph", {
  edge_mat = make_edge_matrix(4, list())
  max_err = check_gradient(p = 4, edge_mat = edge_mat)
  expect_lt(max_err, 1e-4)
})

test_that("gradient matches FD for p=6 with sparse graph", {
  set.seed(99)
  edges_6 = list()
  idx = 1L
  for(i in 1:5) {
    for(j in (i + 1):6) {
      if(runif(1) < 0.6) {
        edges_6[[idx]] = c(i, j)
        idx = idx + 1L
      }
    }
  }
  if(length(edges_6) == 0) edges_6 = list(c(1, 2))
  edge_mat = make_edge_matrix(6, edges_6)
  max_err = check_gradient(p = 6, edge_mat = edge_mat)
  expect_lt(max_err, 1e-4)
})

test_that("gradient matches FD for p=6 with complete graph", {
  p = 6
  all_edges = list()
  idx = 1L
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      all_edges[[idx]] = c(i, j)
      idx = idx + 1L
    }
  }
  edge_mat = make_edge_matrix(p, all_edges)
  max_err = check_gradient(p = 6, edge_mat = edge_mat)
  expect_lt(max_err, 1e-4)
})

test_that("gradient matches FD for p=8 stress test", {
  p = 8
  set.seed(7)
  edges = list()
  idx = 1L
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      if(runif(1) < 0.5) {
        edges[[idx]] = c(i, j)
        idx = idx + 1L
      }
    }
  }
  if(length(edges) == 0) edges = list(c(1, 2))
  edge_mat = make_edge_matrix(p, edges)
  max_err = check_gradient(p = p, edge_mat = edge_mat)
  expect_lt(max_err, 1e-4)
})
