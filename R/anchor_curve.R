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
# Returns: list(theta = per-chain iter x parameters matrices in the frame the
#   slab prior applies to (GGM pairwise values carry the -0.5 factor),
#   gamma = per-chain iter x parameters indicator matrices aligned with theta,
#   indicator = per-chain iter x indicators matrices, the unit the inclusion
#   probability is reported on (identical to gamma for bgm()),
#   family = "normal" or "cauchy").
# ------------------------------------------------------------------
anchor_draws = function(fit) {
  spec = get_fit_spec(fit)
  raw = get_raw_samples(fit)
  if(identical(spec$model_type, "compare")) {
    return(compare_anchor_draws(fit, spec, raw))
  }
  slab_factor = if(identical(spec$model_type, "ggm")) -0.5 else 1
  theta = lapply(raw$pairwise, function(m) slab_factor * m)
  gamma = raw$indicator
  family = spec$prior$interaction_prior_type %||% "normal"
  list(
    theta = theta, gamma = gamma, indicator = gamma,
    family = tolower(family),
    diagonal = anchor_diagonal(spec, raw)
  )
}


# ------------------------------------------------------------------
# anchor_diagonal
# ------------------------------------------------------------------
# The Gamma prior on the precision diagonal, in the frame the weight ratio
# needs it, or NULL for a model that has no such diagonal.
#
# refit_at_scale() rewrites the raw diagonal rate whenever the resolved vary
# mode moves it (rate = eta / s; see vary_diagonal_rate()), so a fit reweighted
# under that mode differs from its target in the diagonal prior as well as in
# the slab, and the weight has to carry both ratios. The sufficient statistic
# the rate multiplies is the sum of the values the prior is evaluated at, and
# the sampler evaluates it at K_jj / 2 (ggm_model.cpp:320, mixed_mrf_
# gradient.cpp:559), not at K_jj; the raw main draws store K_jj itself, the
# same columns the per-draw predict path reconstructs the precision matrix
# from (build_precision_from_draw()).
#
# @param spec  The fit's spec.
# @param raw   Its raw samples.
#
# Returns: list(sum = per-chain vectors of sum_j K_jj / 2, n = number of
#   diagonal elements, shape = the Gamma shape alpha), or NULL.
# ------------------------------------------------------------------
anchor_diagonal = function(spec, raw) {
  cols = if(identical(spec$model_type, "ggm")) {
    seq_len(ncol(raw$main[[1]]))
  } else if(identical(spec$model_type, "mixed_mrf")) {
    which(endsWith(raw$parameter_names[["main"]], "(precision diag)"))
  } else {
    integer(0)
  }
  if(length(cols) == 0L) {
    return(NULL)
  }
  list(
    sum = lapply(raw$main, function(m) 0.5 * rowSums(m[, cols, drop = FALSE])),
    n = length(cols),
    shape = spec$prior$scale_shape %||% 1
  )
}


# ------------------------------------------------------------------
# compare_anchor_draws
# ------------------------------------------------------------------
# Anchor draws for a bgmCompare() fit, where the swept prior is the difference
# slab and one indicator can gate several parameters.
#
# A pairwise difference indicator gates that pair's difference in every group
# contrast; a main-effect difference indicator gates the whole block of that
# variable's threshold differences, again in every contrast. The prior is
# evaluated once per gated parameter, so the importance weight sums over
# parameters while the inclusion probability is reported per indicator.
#
# @param fit   A fitted bgmCompare object.
# @param spec  Its spec.
# @param raw   Its raw samples.
#
# Returns: as anchor_draws().
# ------------------------------------------------------------------
compare_anchor_draws = function(fit, spec, raw) {
  arguments = extract_arguments(fit)
  names_all = raw$parameter_names
  num_variables = as.integer(
    arguments[["num_variables"]] %||% arguments[["no_variables"]]
  )

  num_pairs = length(names_all$pairwise_baseline)
  num_contrasts = length(names_all$pairwise_diff) / num_pairs
  num_main_baseline = length(names_all$main_baseline)

  # Main parameters per variable: one per threshold for an ordinal variable,
  # a linear and a quadratic term for a Blume-Capel one. The field is named
  # per model type -- is_ordinal_variable on a compare fit, is_ordinal on the
  # single-network ones -- and is read exactly, not by partial match.
  is_ordinal = arguments[["is_ordinal_variable"]] %||% arguments[["is_ordinal"]]
  block = ifelse(is_ordinal, as.integer(arguments[["num_categories"]]), 2L)

  # Indicators run over the upper triangle with the diagonal: (v, v) is that
  # variable's main-effect difference, (i, j) the pair's.
  idx = compare_indicator_index(num_variables)
  main_indicator = which(idx[, 1] == idx[, 2])
  pair_indicator = which(idx[, 1] != idx[, 2])

  # Which indicator gates each difference parameter, in the column order the
  # draws carry: pairwise differences contrast by contrast, then main ones.
  pairwise_owner = rep(pair_indicator, times = num_contrasts)
  main_owner = rep(rep(main_indicator, times = block), times = num_contrasts)
  owner = c(pairwise_owner, main_owner)

  pairwise_diff_cols = num_pairs + seq_len(num_pairs * num_contrasts)
  main_diff_cols = num_main_baseline + seq_len(num_main_baseline * num_contrasts)

  theta = vector("list", length(raw$pairwise))
  gamma = vector("list", length(raw$pairwise))
  for(c in seq_along(raw$pairwise)) {
    theta[[c]] = cbind(
      raw$pairwise[[c]][, pairwise_diff_cols, drop = FALSE],
      raw$main[[c]][, main_diff_cols, drop = FALSE]
    )
    gamma[[c]] = raw$indicator[[c]][, owner, drop = FALSE]
  }

  family = spec$prior$difference_prior_type %||% "cauchy"
  list(
    theta = theta, gamma = gamma, indicator = raw$indicator,
    family = tolower(family)
  )
}


# ------------------------------------------------------------------
# anchor_log_weights
# ------------------------------------------------------------------
# Per-draw log importance weights from anchor scale s_a to each target
# scale, for one chain. Only included edges contribute; the weight is
# the slab-density ratio summed over them.
#
# Under vary = "slab-and-diagonal" the refits also move the Gamma prior on the
# precision diagonal, rate = eta / s, so the target differs from the anchor
# there too and the weight carries that ratio as well. For Gamma(alpha,
# eta / s) on each of the n diagonal values x_j = K_jj / 2, the per-draw term is
#
#   n * alpha * log((eta / s) / (eta / s_a)) - (eta / s - eta / s_a) * S_t
#     = -n * alpha * log(s / s_a) - eta * (1 / s - 1 / s_a) * S_t,
#
# with S_t = sum_j K_jj / 2 at draw t. The leading term is the same for every
# draw and so cancels in the self-normalized average; the S_t term is what
# actually reweights. Omitting it left the curve reweighting only the slab
# while the refits it is compared against had moved both priors.
#
# @param theta    iter x edges matrix, slab frame.
# @param gamma    iter x edges indicator matrix.
# @param family   "normal" or "cauchy".
# @param s_a      Anchor slab scale.
# @param s_grid   Target slab scales.
# @param diagonal The anchor's anchor_diagonal() record, or NULL; `sum` is
#                 this chain's per-draw S_t.
# @param eta      Standardized diagonal rate held fixed by the sweep, or NULL
#                 when the sweep leaves the diagonal prior alone.
#
# Returns: iter x length(s_grid) matrix of log weights.
# ------------------------------------------------------------------
anchor_log_weights = function(theta, gamma, family, s_a, s_grid,
                              diagonal = NULL, eta = NULL) {
  # The reweighting identity is written per slab family; a scale-free family
  # (beta-prime) has no s to sweep and any other one has a different density
  # ratio. Silently applying the Cauchy formula would return a curve, so the
  # fence stops instead of guessing.
  if(!family %in% c("normal", "cauchy")) {
    stop(
      "The anchored sensitivity curve supports a normal or Cauchy slab on ",
      "the swept parameters; this fit uses a '", family, "' slab, for which ",
      "there is no slab scale to reweight across."
    )
  }
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
  # Only when the sweep actually moved the diagonal prior.
  if(!is.null(eta) && !is.null(diagonal)) {
    s_t = diagonal$sum
    dlw = vapply(s_grid, function(s) {
      -diagonal$n * diagonal$shape * log(s / s_a) -
        eta * (1 / s - 1 / s_a) * s_t
    }, numeric(length(s_t)))
    if(is.null(dim(dlw))) dlw = matrix(dlw, nrow = 1L)
    lw = lw + dlw
  }
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
# @param eta      Standardized diagonal rate the sweep holds fixed, or NULL
#                 when the sweep leaves the precision diagonal alone. The
#                 reported importance ESS follows from the same weights.
#
# Returns: list(pip = P x E pooled, ess = length-P pooled importance ESS,
#   chain_pip = list of P x E per chain, chain_ess = C x P).
# ------------------------------------------------------------------
anchor_reweight = function(draws, s_a, s_grid, eta = NULL) {
  n_chain = length(draws$theta)
  diag = draws$diagonal
  lw = lapply(seq_len(n_chain), function(c) {
    anchor_log_weights(
      draws$theta[[c]], draws$gamma[[c]], draws$family, s_a, s_grid,
      diagonal = if(is.null(diag)) NULL else list(sum = diag$sum[[c]], n = diag$n, shape = diag$shape),
      eta = eta
    )
  })
  # The weight sums over gated parameters; the inclusion probability is
  # reported on the indicators, which are the same thing for bgm() and a
  # coarser unit for bgmCompare().
  indicator = draws$indicator %||% draws$gamma
  n_grid = length(s_grid)
  n_edge = ncol(indicator[[1]])
  chain_pip = vector("list", n_chain)
  chain_ess = matrix(NA_real_, n_chain, n_grid)
  pip = matrix(NA_real_, n_grid, n_edge)
  ess = numeric(n_grid)
  gamma_all = do.call(rbind, indicator)
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
      cp[p, ] = colSums(w * indicator[[c]]) / sum(w)
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
# PIP scale, then transforming to the log BF, removes the staircase seams
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
      ea = ess_mat[p, a]
      # Inverse-variance weight on the PIP scale (var ~ p(1-p)/ESS). The
      # plug-in variance is computed on a smoothed proportion -- the Bayes
      # estimator under a Jeffreys-like half-count, (ess * pa + 0.5)/(ess + 1)
      # -- so an anchor whose reweighted PIP saturates to 0 or 1 gets the
      # variance its own ESS supports rather than a near-zero one; on the raw
      # plug-in a distant, low-ESS, saturated anchor outweighed every other
      # anchor at the same grid point. The estimate being pooled is still the
      # unsmoothed pa; only its weight is smoothed, and the 1e-6 floor stays
      # as a backstop.
      pa_s = (ea * pa + 0.5) / (ea + 1)
      wa = ea / pmax(pa_s * (1 - pa_s), 1e-6)
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
