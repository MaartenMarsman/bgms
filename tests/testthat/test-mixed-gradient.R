# --------------------------------------------------------------------------- #
# Tests for the mixed MRF Cholesky gradient engine.
#
# Validates logp_and_gradient against central finite differences
# for the mixed MRF model with precision in the NUTS vector.
# --------------------------------------------------------------------------- #

# ---- Helpers ----------------------------------------------------------------

mixed_fd_gradient = function(params, x, y, num_cats, is_ord, base_cat,
                             edge_ind, scale, eps = 1e-5) {
  n_total = length(params)
  fd = numeric(n_total)
  for(k in seq_len(n_total)) {
    p_plus = params
    p_minus = params
    p_plus[k] = p_plus[k] + eps
    p_minus[k] = p_minus[k] - eps
    fp = mixed_test_logp_and_gradient(
      p_plus, x, y, num_cats, as.integer(is_ord),
      base_cat, edge_ind, scale
    )$value
    fm = mixed_test_logp_and_gradient(
      p_minus, x, y, num_cats, as.integer(is_ord),
      base_cat, edge_ind, scale
    )$value
    fd[k] = (fp - fm) / (2 * eps)
  }
  fd
}

mixed_check_gradient = function(params, x, y, num_cats, is_ord, base_cat,
                                edge_ind, scale, eps = 1e-5) {
  res = mixed_test_logp_and_gradient(
    params, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, scale
  )
  fd = mixed_fd_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, scale, eps
  )
  ag = res$gradient
  denom = pmax(abs(ag), abs(fd), 1)
  rel_err = abs(ag - fd) / denom
  max(rel_err)
}


# ---- Ordinal-only tests ----------------------------------------------------- #

test_that("gradient matches FD for ordinal-only (p=3, q=2, conditional)", {
  set.seed(42)
  n = 80
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  n_main = sum(num_cats)
  n_pw = p * (p - 1) / 2
  n_chol = q * (q + 1) / 2
  set.seed(123)
  params = rnorm(n_main + n_pw + q + p * q + n_chol, sd = 0.3)
  err = mixed_check_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, 2.5
  )
  expect_lt(err, 1e-5)
})

test_that("gradient matches FD for ordinal-only (p=3, q=2, marginal)", {
  set.seed(42)
  n = 80
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  n_main = sum(num_cats)
  n_pw = p * (p - 1) / 2
  n_chol = q * (q + 1) / 2
  set.seed(123)
  params = rnorm(n_main + n_pw + q + p * q + n_chol, sd = 0.3)
  err = mixed_check_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, 2.5
  )
  expect_lt(err, 1e-5)
})


# ---- Blume-Capel tests ------------------------------------------------------- #

test_that("gradient matches FD for mixed ord+BC (p=3, q=2, conditional)", {
  set.seed(99)
  n = 70
  p = 3
  q = 2
  x = matrix(sample(0:3, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = c(3L, 3L, 3L)
  is_ord = c(1L, 0L, 1L)
  base_cat = c(0L, 1L, 0L)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  n_main = 3 + 2 + 3
  n_pw = 3
  n_chol = q * (q + 1) / 2
  set.seed(2)
  params = rnorm(n_main + n_pw + q + p * q + n_chol, sd = 0.3)
  err = mixed_check_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, 2.5
  )
  expect_lt(err, 1e-5)
})

test_that("gradient matches FD for mixed ord+BC (p=3, q=2, marginal)", {
  set.seed(99)
  n = 70
  p = 3
  q = 2
  x = matrix(sample(0:3, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = c(3L, 3L, 3L)
  is_ord = c(1L, 0L, 1L)
  base_cat = c(0L, 1L, 0L)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  n_main = 3 + 2 + 3
  n_pw = 3
  n_chol = q * (q + 1) / 2
  set.seed(2)
  params = rnorm(n_main + n_pw + q + p * q + n_chol, sd = 0.3)
  err = mixed_check_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, 2.5
  )
  expect_lt(err, 1e-5)
})


# ---- Larger q tests --------------------------------------------------------- #

test_that("gradient matches FD for ordinal (p=2, q=3, conditional)", {
  set.seed(42)
  n = 60
  p = 2
  q = 3
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  n_main = sum(num_cats)
  n_pw = p * (p - 1) / 2
  n_chol = q * (q + 1) / 2
  set.seed(1)
  params = rnorm(n_main + n_pw + q + p * q + n_chol, sd = 0.3)
  err = mixed_check_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, 2.5
  )
  expect_lt(err, 1e-5)
})

test_that("gradient matches FD for ordinal (p=4, q=4, conditional)", {
  set.seed(77)
  n = 100
  p = 4
  q = 4
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  n_main = sum(num_cats)
  n_pw = p * (p - 1) / 2
  n_chol = q * (q + 1) / 2
  set.seed(3)
  params = rnorm(n_main + n_pw + q + p * q + n_chol, sd = 0.2)
  err = mixed_check_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, 2.5
  )
  expect_lt(err, 1e-5)
})


# ---- Sparse edges test ------------------------------------------------------- #

test_that("gradient matches FD with sparse edges (p=3, q=2, conditional)", {
  set.seed(42)
  n = 60
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(0L, total, total)
  edge_ind[1, 2] = edge_ind[2, 1] = 1L
  edge_ind[1, 4] = edge_ind[4, 1] = 1L
  edge_ind[2, 5] = edge_ind[5, 2] = 1L
  n_main = sum(num_cats)
  n_pw = 1
  n_cross = 2
  # Kyy theta block: q log-diagonals + one f_q coordinate per included yy
  # edge; the (4,5) edge is excluded here, so the block is just the psis.
  n_chol = q
  set.seed(4)
  params = rnorm(n_main + n_pw + q + n_cross + n_chol, sd = 0.3)
  err = mixed_check_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, 2.5
  )
  expect_lt(err, 1e-5)
})

test_that("gradient matches FD with sparse edges (p=3, q=2, marginal)", {
  set.seed(42)
  n = 60
  p = 3
  q = 2
  x = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(0L, total, total)
  edge_ind[1, 2] = edge_ind[2, 1] = 1L
  edge_ind[1, 4] = edge_ind[4, 1] = 1L
  edge_ind[2, 5] = edge_ind[5, 2] = 1L
  n_main = sum(num_cats)
  n_pw = 1
  n_cross = 2
  # Kyy theta block: q log-diagonals + one f_q coordinate per included yy
  # edge; the (4,5) edge is excluded here, so the block is just the psis.
  n_chol = q
  set.seed(4)
  params = rnorm(n_main + n_pw + q + n_cross + n_chol, sd = 0.3)
  err = mixed_check_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, 2.5
  )
  expect_lt(err, 1e-5)
})
