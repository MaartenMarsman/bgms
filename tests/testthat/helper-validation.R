# ===========================================================================
# Shared helpers for the mixed MRF validation test suite
# ===========================================================================
# Provides (all four are used by test-validation-slow.R):
#   - make_network()       : generate reproducible true parameter sets
#   - generate_data()      : simulate data from true parameters via bgms
#   - extract_bgms_blocks(): pull (mux, muy, pairwise_disc, pairwise_cross, pairwise_cont) from bgms fit
#   - flatten_params()     : flatten all blocks to a single named vector
# ===========================================================================

# ------------------------------------------------------------------
# make_network
# ------------------------------------------------------------------
# Build a reproducible mixed MRF parameter set.
#
# @param p              Number of discrete variables.
# @param q              Number of continuous variables.
# @param n_cat          Integer vector of length p: number of categories per
#                       discrete variable (bgms convention = max index,
#                       e.g. binary = 1).
# @param variable_type  Character vector of length p: "ordinal" or
#                       "blume-capel" per discrete variable. Default: all
#                       ordinal.
# @param baseline_category  Integer vector of length p: baseline category
#                       for Blume-Capel variables. Ignored for ordinal.
#                       Default: all zeros.
# @param density        Approximate fraction of non-zero edges.
# @param seed           Random seed for reproducibility.
#
# Returns: named list with mux, muy, pairwise_disc, pairwise_cross, pairwise_cont, n_cat, p, q,
#   variable_type, baseline_category.
# ------------------------------------------------------------------
make_network = function(p, q, n_cat, variable_type = rep("ordinal", p),
                        baseline_category = rep(0L, p),
                        density = 0.5, seed = 42) {
  set.seed(seed)

  max_cat = max(n_cat)
  # BC variables use 2 columns (alpha, beta); ensure mux is wide enough
  mux_cols = max(max_cat, 2L)

  # --- Main effects (mux) ---
  # Ordinal: threshold parameters, NA-padded.
  # Blume-Capel: alpha (linear) and beta (quadratic), rest NA.
  mux = matrix(NA_real_, nrow = p, ncol = mux_cols)
  for(i in seq_len(p)) {
    if(variable_type[i] == "blume-capel") {
      mux[i, 1] = round(runif(1, -0.3, 0.3), 2) # alpha
      mux[i, 2] = round(runif(1, -0.4, -0.05), 2) # beta (negative = peaked)
    } else {
      vals = sort(round(seq(-0.5, 0.5, length.out = n_cat[i]), 2))
      if(n_cat[i] == 1) vals = round(runif(1, -0.3, 0.3), 2)
      mux[i, seq_len(n_cat[i])] = vals
    }
  }

  # --- Continuous means ---
  muy = round(runif(q, -0.5, 0.5), 2)

  # --- pairwise_disc: association-scale discrete interactions (A = <U+03C3>/2, symmetric, zero diagonal) ---
  n_edges_xx = p * (p - 1) / 2
  mask_xx = rbinom(n_edges_xx, 1, density)
  vals_xx = mask_xx * round(runif(n_edges_xx, 0.15, 0.4) * sample(c(-1, 1), n_edges_xx, replace = TRUE), 2)
  pairwise_disc = matrix(0, p, p)
  pairwise_disc[upper.tri(pairwise_disc)] = vals_xx
  pairwise_disc = pairwise_disc + t(pairwise_disc)
  pairwise_disc = 0.5 * pairwise_disc

  # --- pairwise_cont: association-scale continuous block (negative-definite; precision = -2 * pairwise_cont) ---
  # Build a positive-definite precision matrix first, then convert to association scale.
  n_edges_yy = q * (q - 1) / 2
  mask_yy = rbinom(n_edges_yy, 1, density)
  vals_yy = mask_yy * round(runif(n_edges_yy, 0.05, 0.2) * sample(c(-1, 1), n_edges_yy, replace = TRUE), 2)
  pairwise_cont = matrix(0, q, q)
  pairwise_cont[upper.tri(pairwise_cont)] = vals_yy
  pairwise_cont = pairwise_cont + t(pairwise_cont)
  diag(pairwise_cont) = abs(rowSums(pairwise_cont)) + runif(q, 1.2, 1.8)
  pairwise_cont = -0.5 * pairwise_cont

  # --- pairwise_cross: ordinal-continuous cross interactions ---
  n_cross = p * q
  mask_xy = rbinom(n_cross, 1, density)
  vals_xy = mask_xy * round(runif(n_cross, 0.1, 0.3) * sample(c(-1, 1), n_cross, replace = TRUE), 2)
  pairwise_cross = matrix(vals_xy, p, q)

  list(
    mux = mux, muy = muy,
    pairwise_disc = pairwise_disc, pairwise_cross = pairwise_cross, pairwise_cont = pairwise_cont,
    n_cat = n_cat, p = p, q = q,
    variable_type = variable_type,
    baseline_category = as.integer(baseline_category)
  )
}

# ------------------------------------------------------------------
# generate_data
# ------------------------------------------------------------------
# Simulate data from true network parameters.
#
# @param net        Output from make_network().
# @param n          Number of observations.
# @param source     "bgms" (the only simulator wired up here).
# @param iter       Gibbs burn-in iterations.
# @param seed       Random seed.
#
# Returns: data.frame (n x (p+q)) with ordinal columns first, then continuous.
# ------------------------------------------------------------------
generate_data = function(net, n, source = "bgms", iter = 1000L, seed = 1) {
  if(source != "bgms") {
    stop('generate_data(): source must be "bgms", not "', source, '".')
  }
  # simulate_mrf() does not support mixed types; use the C++ Gibbs
  # sampler directly via sample_mixed_mrf_gibbs().
  sim = bgms:::sample_mixed_mrf_gibbs(
    num_states = as.integer(n),
    pairwise_disc_r = net$pairwise_disc,
    pairwise_cross_r = net$pairwise_cross,
    pairwise_cont_r = net$pairwise_cont,
    mux_r = net$mux,
    muy_r = net$muy,
    num_categories_r = as.integer(net$n_cat),
    variable_type_r = net$variable_type,
    baseline_category_r = net$baseline_category,
    iter = as.integer(iter),
    seed = as.integer(seed)
  )
  df = as.data.frame(cbind(sim$x, sim$y))
  names(df) = c(paste0("X", seq_len(net$p)), paste0("Y", seq_len(net$q)))
  df
}

# ------------------------------------------------------------------
# extract_bgms_blocks
# ------------------------------------------------------------------
# Extract parameter blocks from a bgms fit object.
#
# @param fit  A bgms object.
# @param net  The true network (for dimension reference).
#
# Returns: list(mux, muy, pairwise_disc, pairwise_cross, pairwise_cont).
# ------------------------------------------------------------------
extract_bgms_blocks = function(fit, net) {
  pm = coef(fit)
  mux = pm$main$discrete # p x max_cat
  muy_vec = pm$main$continuous[, "mean"] # length q

  pw = pm$pairwise # (p+q) x (p+q)
  p = net$p
  q = net$q
  pairwise_disc = pw[seq_len(p), seq_len(p)]
  pairwise_cross = pw[seq_len(p), p + seq_len(q)]
  # pairwise_cont off-diagonal from pairwise, diagonal converted from residual variance
  pairwise_cont = pw[p + seq_len(q), p + seq_len(q)]
  diag(pairwise_cont) = -1 / (2 * fit$posterior_mean_residual_variance)

  list(mux = mux, muy = muy_vec, pairwise_disc = pairwise_disc, pairwise_cross = pairwise_cross, pairwise_cont = pairwise_cont)
}

# ------------------------------------------------------------------
# flatten_params
# ------------------------------------------------------------------
# Flatten parameter blocks to a single named vector for comparison.
# Excludes NA entries in mux.
#
# @param blocks  list(mux, muy, pairwise_disc, pairwise_cross, pairwise_cont).
# @param prefix  Optional prefix for names.
#
# Returns: named numeric vector.
# ------------------------------------------------------------------
flatten_params = function(blocks, prefix = "") {
  mux_vals = as.vector(t(blocks$mux))
  mux_keep = !is.na(mux_vals)
  mux_named = mux_vals[mux_keep]
  names(mux_named) = paste0(prefix, "mux_", which(mux_keep))

  muy_named = blocks$muy
  names(muy_named) = paste0(prefix, "muy_", seq_along(muy_named))

  # pairwise_disc upper triangle
  disc_ut = blocks$pairwise_disc[upper.tri(blocks$pairwise_disc)]
  names(disc_ut) = paste0(prefix, "disc_", seq_along(disc_ut))

  # pairwise_cross full
  cross_vals = as.vector(blocks$pairwise_cross)
  names(cross_vals) = paste0(prefix, "cross_", seq_along(cross_vals))

  # pairwise_cont upper triangle (includes diagonal)
  cont_ut = blocks$pairwise_cont[upper.tri(blocks$pairwise_cont, diag = TRUE)]
  names(cont_ut) = paste0(prefix, "cont_", seq_along(cont_ut))

  c(mux_named, muy_named, disc_ut, cross_vals, cont_ut)
}
