# --------------------------------------------------------------------------- #
# Finite-difference checks for MixedMRFModel::logp_and_gradient_full.
#
# Targets the Cholesky-to-K Jacobian + per-column Pfaffian split that ports
# GGM Fix 2 to the mixed-MRF Kyy block. The key code paths exercised are:
#   - the unified ldj = q*log(2) + sum_j (q+1-j) psi_j Jacobian
#   - the per-column Pfaffian 0.5 * sum_qq log det(A_qq diag(M_qq^-1) A_qq^T)
#   - the matching Pfaffian adjoint that flows into R_bar at excluded-edge
#     positions, then into the diagonal psi gradient through the exp() chain
#
# Full Kyy graph: A_qq are empty so Pfaffian = 0 — this verifies the
# refactor doesn't change anything in the no-constraint case.
# Sparse Kyy:    Pfaffian is non-trivial; tests both identity and random
#                positive diagonal mass (the latter is what NUTS uses after
#                warmup).
# --------------------------------------------------------------------------- #

# ---- Helpers ---------------------------------------------------------------- #

mixed_full_fd_gradient = function(params, x, y, num_cats, is_ord, base_cat,
                                  edge_ind, scale, inv_mass = NULL,
                                  delta = 0, eps = 1e-5) {
  n_total = length(params)
  fd = numeric(n_total)
  for(k in seq_len(n_total)) {
    p_plus = params
    p_minus = params
    p_plus[k] = p_plus[k] + eps
    p_minus[k] = p_minus[k] - eps
    fp = mixed_test_logp_and_gradient_full(
      p_plus, x, y, num_cats, as.integer(is_ord),
      base_cat, edge_ind, scale,
      inv_mass_diag = inv_mass, delta = delta
    )$value
    fm = mixed_test_logp_and_gradient_full(
      p_minus, x, y, num_cats, as.integer(is_ord),
      base_cat, edge_ind, scale,
      inv_mass_diag = inv_mass, delta = delta
    )$value
    fd[k] = (fp - fm) / (2 * eps)
  }
  fd
}

mixed_full_check_gradient = function(params, x, y, num_cats, is_ord, base_cat,
                                     edge_ind, scale, inv_mass = NULL,
                                     delta = 0, eps = 1e-5) {
  res = mixed_test_logp_and_gradient_full(
    params, x, y, num_cats, as.integer(is_ord),
    base_cat, edge_ind, scale,
    inv_mass_diag = inv_mass, delta = delta
  )
  fd = mixed_full_fd_gradient(
    params, x, y, num_cats, is_ord, base_cat,
    edge_ind, scale, inv_mass, delta, eps
  )
  ag = res$gradient
  denom = pmax(abs(ag), abs(fd), 1)
  max(abs(ag - fd) / denom)
}

# Build a small mixed model with q continuous and p discrete variables,
# returning everything needed by the gradient harness in one shot.
mixed_full_make_setup = function(p = 2, q = 3, n = 60, seed = 42) {
  set.seed(seed)
  x = matrix(sample(0:1, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)
  num_cats = rep(2L, p)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L

  n_main = sum(num_cats)
  n_pw_xx = p * (p - 1L) / 2L
  n_means = q
  n_xy = p * q
  n_chol = q * (q + 1L) / 2L
  total_dim = n_main + n_pw_xx + n_means + n_xy + n_chol

  set.seed(seed + 1)
  params = rnorm(total_dim, sd = 0.3)

  list(
    x = x, y = y, num_cats = num_cats, is_ord = is_ord, base_cat = base_cat,
    edge_ind = edge_ind, params = params,
    p = p, q = q, total = total, total_dim = total_dim
  )
}


# ---- Full Kyy graph: Pfaffian = 0 ------------------------------------------- #

test_that("gradient matches FD on full Kyy graph (identity mass)", {
  s = mixed_full_make_setup()
  err = mixed_full_check_gradient(
    s$params, s$x, s$y, s$num_cats, s$is_ord, s$base_cat,
    s$edge_ind, 2.5
  )
  expect_lt(err, 1e-5)
})


# ---- Sparse Kyy graph: Pfaffian is non-trivial ------------------------------ #

test_that("gradient matches FD with one excluded Kyy edge (identity mass)", {
  s = mixed_full_make_setup()
  # Exclude one Kyy edge: between continuous-1 and continuous-3.
  s$edge_ind[s$p + 1L, s$p + 3L] = 0L
  s$edge_ind[s$p + 3L, s$p + 1L] = 0L
  err = mixed_full_check_gradient(
    s$params, s$x, s$y, s$num_cats, s$is_ord, s$base_cat,
    s$edge_ind, 2.5
  )
  expect_lt(err, 1e-5)
})


test_that("gradient matches FD with one excluded Kyy edge (random mass)", {
  s = mixed_full_make_setup()
  s$edge_ind[s$p + 1L, s$p + 3L] = 0L
  s$edge_ind[s$p + 3L, s$p + 1L] = 0L
  set.seed(7)
  inv_mass = runif(s$total_dim, min = 0.5, max = 2.0)
  err = mixed_full_check_gradient(
    s$params, s$x, s$y, s$num_cats, s$is_ord, s$base_cat,
    s$edge_ind, 2.5,
    inv_mass = inv_mass
  )
  expect_lt(err, 1e-5)
})


test_that("gradient matches FD with two excluded Kyy edges (random mass)", {
  s = mixed_full_make_setup(q = 4)
  # Exclude (p+1, p+3) and (p+2, p+4): two columns carry constraints.
  s$edge_ind[s$p + 1L, s$p + 3L] = 0L
  s$edge_ind[s$p + 3L, s$p + 1L] = 0L
  s$edge_ind[s$p + 2L, s$p + 4L] = 0L
  s$edge_ind[s$p + 4L, s$p + 2L] = 0L
  set.seed(13)
  inv_mass = runif(s$total_dim, min = 0.5, max = 2.0)
  err = mixed_full_check_gradient(
    s$params, s$x, s$y, s$num_cats, s$is_ord, s$base_cat,
    s$edge_ind, 2.5,
    inv_mass = inv_mass
  )
  expect_lt(err, 1e-5)
})


# ---- Determinant tilt on the Kyy block -------------------------------------- #

test_that("delta = 0 is a no-op (bit-identical to no-tilt call)", {
  s = mixed_full_make_setup()
  ref = mixed_test_logp_and_gradient_full(
    s$params, s$x, s$y, s$num_cats, as.integer(s$is_ord),
    s$base_cat, s$edge_ind, 2.5
  )
  with_zero = mixed_test_logp_and_gradient_full(
    s$params, s$x, s$y, s$num_cats, as.integer(s$is_ord),
    s$base_cat, s$edge_ind, 2.5,
    delta = 0
  )
  expect_equal(with_zero$value, ref$value)
  expect_equal(with_zero$gradient, ref$gradient)
})


test_that("gradient matches FD with delta > 0 on full Kyy", {
  s = mixed_full_make_setup()
  err = mixed_full_check_gradient(
    s$params, s$x, s$y, s$num_cats, s$is_ord, s$base_cat,
    s$edge_ind, 2.5,
    delta = 1.5
  )
  expect_lt(err, 1e-5)
})


test_that("gradient matches FD with delta > 0 on sparse Kyy (random mass)", {
  s = mixed_full_make_setup()
  s$edge_ind[s$p + 1L, s$p + 3L] = 0L
  s$edge_ind[s$p + 3L, s$p + 1L] = 0L
  set.seed(21)
  inv_mass = runif(s$total_dim, min = 0.5, max = 2.0)
  err = mixed_full_check_gradient(
    s$params, s$x, s$y, s$num_cats, s$is_ord, s$base_cat,
    s$edge_ind, 2.5,
    inv_mass = inv_mass, delta = 2.0
  )
  expect_lt(err, 1e-5)
})


test_that("delta shifts logp by exactly delta * log|Kyy|", {
  # Analytic check: tilting by |Kyy|^delta adds delta * log|Kyy| to logp.
  s = mixed_full_make_setup()
  base = mixed_test_logp_and_gradient_full(
    s$params, s$x, s$y, s$num_cats, as.integer(s$is_ord),
    s$base_cat, s$edge_ind, 2.5,
    delta = 0
  )$value
  tilted = mixed_test_logp_and_gradient_full(
    s$params, s$x, s$y, s$num_cats, as.integer(s$is_ord),
    s$base_cat, s$edge_ind, 2.5,
    delta = 3.0
  )$value

  # log|Kyy| from psi entries (Cholesky-of-precision diagonals are exp(psi)).
  n_main = sum(s$num_cats)
  n_pw_xx = s$p * (s$p - 1L) / 2L
  n_means = s$q
  n_xy = s$p * s$q
  chol_offset = n_main + n_pw_xx + n_means + n_xy
  psi_positions = chol_offset + cumsum(seq_len(s$q))
  log_det_Kyy = 2 * sum(s$params[psi_positions])

  expect_equal(tilted - base, 3.0 * log_det_Kyy, tolerance = 1e-12)
})
