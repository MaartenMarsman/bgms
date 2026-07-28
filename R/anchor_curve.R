# ==============================================================================
# Anchored sensitivity curve: importance reweighting between fixed-scale fits.
# ==============================================================================
#
# A fit at fixed slab scale s_a is reweighted to a nearby scale s with per-draw
# slab ratios over the currently included edges; the likelihood cancels and no
# normalizer is involved (spike terms cancel, and a fixed-scale anchor has no
# scale posterior to integrate). Self-normalized importance sampling then gives
# the inclusion probability, Bayes factor, and importance ESS at any target
# scale within the anchor's usable radius; log-spaced anchors make the radii
# overlap and the stitched curve continuous. Design and validation gates:
# dev/audit/2026-07-28-anchored-curve-spec.md.
# ==============================================================================


# ------------------------------------------------------------------
# anchor_draws
# ------------------------------------------------------------------
# Per-chain draws of a fit in the slab frame, for reweighting.
#
# @param fit  A fitted bgms object.
#
# Returns: list(theta = per-chain iter x edges matrices in the frame the
#   slab prior applies to (GGM pairwise values carry the -0.5 factor),
#   gamma = per-chain iter x edges indicator matrices,
#   family = "normal" or "cauchy").
# ------------------------------------------------------------------
anchor_draws = function(fit) {
  spec = get_fit_spec(fit)
  raw = get_raw_samples(fit)
  slab_factor = if(identical(spec$model_type, "ggm")) -0.5 else 1
  theta = lapply(raw$pairwise, function(m) slab_factor * m)
  gamma = raw$indicator
  family = spec$prior$interaction_prior_type %||% "normal"
  list(theta = theta, gamma = gamma, family = tolower(family))
}


# ------------------------------------------------------------------
# anchor_log_weights
# ------------------------------------------------------------------
# Per-draw log importance weights from anchor scale s_a to each target
# scale, for one chain. Only included edges contribute; the weight is
# the slab-density ratio summed over them.
#
# @param theta    iter x edges matrix, slab frame.
# @param gamma    iter x edges indicator matrix.
# @param family   "normal" or "cauchy".
# @param s_a      Anchor slab scale.
# @param s_grid   Target slab scales.
#
# Returns: iter x length(s_grid) matrix of log weights.
# ------------------------------------------------------------------
anchor_log_weights = function(theta, gamma, family, s_a, s_grid) {
  incl = gamma == 1
  m_t = rowSums(incl)
  if(identical(family, "normal")) {
    a_t = rowSums((theta * incl)^2)
    lw = vapply(s_grid, function(s) {
      -m_t * log(s / s_a) - 0.5 * a_t * (1 / s^2 - 1 / s_a^2)
    }, numeric(length(m_t)))
  } else {
    base = log1p((theta / s_a)^2)
    lw = vapply(s_grid, function(s) {
      -m_t * log(s / s_a) -
        rowSums(incl * (log1p((theta / s)^2) - base))
    }, numeric(length(m_t)))
  }
  if(is.null(dim(lw))) lw = matrix(lw, nrow = 1L)
  lw
}


# ------------------------------------------------------------------
# anchor_reweight
# ------------------------------------------------------------------
# Reweighted inclusion probabilities and importance ESS at each target
# scale, from one anchor fit's draws: pooled over chains and per chain.
#
# @param draws    Output of anchor_draws() for the anchor fit.
# @param s_a      Anchor slab scale.
# @param s_grid   Target slab scales (length P).
#
# Returns: list(pip = P x E pooled, ess = length-P pooled importance ESS,
#   chain_pip = list of P x E per chain, chain_ess = C x P).
# ------------------------------------------------------------------
anchor_reweight = function(draws, s_a, s_grid) {
  n_chain = length(draws$theta)
  lw = lapply(seq_len(n_chain), function(c) {
    anchor_log_weights(draws$theta[[c]], draws$gamma[[c]], draws$family, s_a, s_grid)
  })
  n_grid = length(s_grid)
  n_edge = ncol(draws$gamma[[1]])
  chain_pip = vector("list", n_chain)
  chain_ess = matrix(NA_real_, n_chain, n_grid)
  pip = matrix(NA_real_, n_grid, n_edge)
  ess = numeric(n_grid)
  gamma_all = do.call(rbind, draws$gamma)
  for(p in seq_len(n_grid)) {
    lw_all = unlist(lapply(lw, function(m) m[, p]))
    w_all = exp(lw_all - max(lw_all))
    ess[p] = sum(w_all)^2 / sum(w_all^2)
    pip[p, ] = colSums(w_all * gamma_all) / sum(w_all)
  }
  for(c in seq_len(n_chain)) {
    cp = matrix(NA_real_, n_grid, n_edge)
    for(p in seq_len(n_grid)) {
      w = exp(lw[[c]][, p] - max(lw[[c]][, p]))
      chain_ess[c, p] = sum(w)^2 / sum(w^2)
      cp[p, ] = colSums(w * draws$gamma[[c]]) / sum(w)
    }
    chain_pip[[c]] = cp
  }
  list(pip = pip, ess = ess, chain_pip = chain_pip, chain_ess = chain_ess)
}


# ------------------------------------------------------------------
# assemble_curve
# ------------------------------------------------------------------
# Stitch the display curve by precision-weighted pooling of the anchors.
# At each display point every anchor whose importance ESS clears the floor
# contributes its reweighted PIP, weighted by inverse variance on the PIP
# scale (w = ESS / (p(1-p)), the reciprocal of p(1-p)/ESS). Pooling on the
# PIP scale, then transforming to the log10 BF, removes the staircase seams
# and the infinities that winner-take-all selection produced at anchor
# switch points (a capped-edge anchor no longer hands off discontinuously to
# a finite one). Points where no anchor clears the floor are masked. The
# per-chain curves pool with the same weights, so the downstream per-point
# MCSE and unanimity checks carry over unchanged. This is the MBAR-lite /
# multistate-bridge direction; full self-consistent MBAR weights are future
# work. Exactness at anchors lives on the anchor fits' own RB statistics
# (verdict columns, chosen-scale quantities), not on these pooled rows.
#
# @param reweights     List over anchors of anchor_reweight() output.
# @param usable        Logical per anchor (convergence-gated).
# @param ess_floor     Minimum importance ESS for an anchor to contribute.
# @param anchor_index  Grid position of each anchor (for anchor_used tagging).
#
# Returns: list(pip = P x E, chain_pip = list C of P x E, anchor_used =
#   dominant anchor per point (NA if masked), ess = length-P best ESS).
# ------------------------------------------------------------------
assemble_curve = function(reweights, usable, ess_floor, anchor_index = NULL) {
  n_grid = length(reweights[[1]]$ess)
  n_edge = ncol(reweights[[1]]$pip)
  n_chain = length(reweights[[1]]$chain_pip)
  ess_mat = vapply(reweights, `[[`, numeric(n_grid), "ess")
  if(is.null(dim(ess_mat))) ess_mat = matrix(ess_mat, nrow = 1L)
  ess_mat[, !usable] = 0 # an unusable anchor never contributes

  pip = matrix(NA_real_, n_grid, n_edge)
  chain_pip = replicate(n_chain, matrix(NA_real_, n_grid, n_edge), simplify = FALSE)
  anchor_used = rep(NA_integer_, n_grid)
  best_ess = rep(NA_real_, n_grid)

  for(p in seq_len(n_grid)) {
    contrib = which(ess_mat[p, ] >= ess_floor)
    if(length(contrib) == 0L) next # masked: no anchor reaches this scale
    best_ess[p] = max(ess_mat[p, contrib])
    anchor_used[p] = contrib[which.max(ess_mat[p, contrib])]
    wsum = numeric(n_edge)
    pnum = numeric(n_edge)
    cnum = replicate(n_chain, numeric(n_edge), simplify = FALSE)
    for(a in contrib) {
      pa = reweights[[a]]$pip[p, ]
      # Inverse-variance weight on the PIP scale (var ~ p(1-p)/ESS), clamped
      # so a saturated edge does not divide by zero.
      wa = ess_mat[p, a] / pmax(pa * (1 - pa), 1e-6)
      wsum = wsum + wa
      pnum = pnum + wa * pa
      for(cc in seq_len(n_chain)) {
        cnum[[cc]] = cnum[[cc]] + wa * reweights[[a]]$chain_pip[[cc]][p, ]
      }
    }
    pip[p, ] = pnum / wsum
    for(cc in seq_len(n_chain)) chain_pip[[cc]][p, ] = cnum[[cc]] / wsum
  }

  # Tag each pooled point with its dominant (highest-ESS) anchor for
  # reporting; at an anchor's own grid point that is the anchor itself.
  for(a in seq_along(anchor_index)) {
    if(usable[a] && !is.na(anchor_used[anchor_index[a]])) {
      anchor_used[anchor_index[a]] = a
    }
  }
  list(pip = pip, chain_pip = chain_pip, anchor_used = anchor_used, ess = best_ess)
}
