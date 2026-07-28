# ==============================================================================
# prior_sensitivity_check: scale-robustness of edge inclusion verdicts
# ==============================================================================
#
# Refit the model at a small grid of interaction slab scales (multipliers of the
# fit's chosen scale) and classify every edge at every scale, with a Monte Carlo
# wobble guard calibrated from a replicate refit at the chosen scale. Warm starts
# and the refit engine live in refit_engine.R. See the package vignette
# "Checking prior sensitivity" and Bartos et al. (arXiv 2604.21596).
#
# The previous single-fit conditional-density reweighting curve was removed here:
# on real network data the scale posterior collapses and its importance ESS is ~1
# across the advertised band (phase-B review, 2026-07-28). Git history and
# dev/audit/2026-07-28-phase-b-review.md are the archive if it is ever revived.
# ==============================================================================


# ------------------------------------------------------------------
# verdict_from_bf
# ------------------------------------------------------------------
# Map an inclusion Bayes factor to an evidence verdict.
#
# @param bf         Numeric inclusion Bayes factor(s).
# @param threshold  Upper evidence threshold; lower is 1 / threshold.
#
# Returns: character vector of "presence" / "absence" / "undecided".
# ------------------------------------------------------------------
verdict_from_bf = function(bf, threshold) {
  out = rep("undecided", length(bf))
  out[bf >= threshold] = "presence"
  out[bf <= 1 / threshold] = "absence"
  out[is.na(bf)] = NA_character_
  out
}


#' @title Prior Sensitivity of Edge Inclusion Verdicts
#'
#' @description
#' Checks whether each edge's inclusion verdict is robust to the interaction
#' (pairwise) slab scale, by refitting the model at a small grid of scales and
#' classifying every edge at every scale. Works on any \code{\link{bgm}()} fit
#' with edge selection --- no random-scale hyperprior is required. A built-in
#' Monte Carlo wobble guard, calibrated from a replicate refit at the chosen
#' scale, keeps the check from crying wolf on boundary edges.
#'
#' @details
#' The check refits the model at scales \code{scale_multipliers * s0}, where
#' \eqn{s_0} is the chosen scale from \code{interaction_prior}. Multipliers
#' (not absolute scales) are used so one default serves \code{normal_prior(1)}
#' and \code{cauchy_prior(2.5)} alike. The multiplier \code{1} is not redundant:
#' its refit is a replicate of the original analysis, and the spread of
#' \eqn{|\Delta \log_{10} \mathrm{BF}|} between the chosen-scale refit and this
#' replicate is the empirical wobble yardstick. All cross-scale comparisons run
#' through refits only, so a single identical pipeline produces every verdict.
#'
#' \strong{Warm starts.} For ordinal (omrf) fits each refit starts from the
#' original fit's per-chain final state, and a NUTS refit additionally carries
#' the adapted step size and diagonal mass matrix, so a short warmup suffices
#' and the whole check costs about one original fit. Because the warm starts sit
#' near the chosen-scale posterior, cross-chain dispersion is reduced by
#' construction, which weakens split-\eqn{\hat R} as a between-chain diagnostic;
#' the refit gate therefore leans on per-chain verdict agreement and the
#' indicator transition ESS, not on \eqn{\hat R} alone. Continuous (GGM) and
#' mixed-MRF fits refit cold (full warmup), costing about one fit per scale.
#'
#' \strong{The mover rule.} An edge is flagged scale-sensitive only if its
#' verdict differs somewhere on the grid \emph{and} its \eqn{\log_{10}}
#' Bayes-factor change across scales exceeds
#' \code{max(tolerance, 2 * MCSE, wobble)}, with the MCSE from the
#' Rao-Blackwellized machinery and the wobble the 95th percentile of the
#' s0-replicate spread. Each edge is reported as \code{stable},
#' \code{moved-beyond-wobble}, or \code{indistinguishable-from-wobble}; a bare
#' verdict flip inside the replicate noise is never reported as a move.
#'
#' \strong{Edge-level sufficiency.} An edge whose per-chain verdicts disagree, or
#' whose between-chain-inflated \eqn{\log_{10}} BF band straddles a verdict
#' threshold, is marked \code{insufficient} at that scale: the refit cannot
#' certify its verdict. This errs toward caution --- disagreement widens the
#' band rather than vanishing into a pooled estimate.
#'
#' \strong{Two-level gate.} Each refit passes a refit-level gate (continuous
#' split-\eqn{\hat R} below 1.01, bulk RB-inclusion median \eqn{\hat R} below
#' 1.01, E-BFMI and energy variance-ratio for NUTS) before its verdicts are
#' used; a refit that fails is reported as unusable rather than silently pooled.
#'
#' A fit run with an \code{interaction_scale_prior} is a different analysis
#' (scale learning); for such a fit the report adds the realized-versus-nominal
#' scale interval and the scale-averaged verdict, clearly labeled as the
#' learned-scale analysis, not a robustness average.
#'
#' @param bgms_object A fitted \code{bgms} object from \code{\link{bgm}()} run
#'   with \code{edge_selection = TRUE}.
#' @param scale_multipliers Numeric vector of positive multipliers of the chosen
#'   scale at which to refit. Default \code{c(0.5, 1, 2.5)} (narrow / chosen /
#'   wide). The multiplier \code{1} is always included.
#' @param evidence_threshold Positive numeric. Inclusion Bayes factor threshold
#'   for a presence verdict; \code{1 / evidence_threshold} is the absence
#'   threshold. Default: \code{10}.
#' @param refit_sampler One of \code{"same-as-fit"} (default; inherit the
#'   original fit's update method) or an explicit \code{"nuts"},
#'   \code{"adaptive-metropolis"}, or \code{"gibbs"}. NUTS refits of an ordinal
#'   fit carry the adapted metric and run fastest.
#' @param iter,warmup Integer sampling and warmup iterations per refit, or
#'   \code{NULL} (default) to use the validated short schedule for warm NUTS
#'   refits and inherit the original fit's schedule otherwise.
#' @param tolerance Numeric. Minimum \eqn{\log_{10}} BF change across scales for
#'   a verdict flip to count as a move, before the MCSE and wobble floors.
#'   Default: \code{0.5}.
#' @param include_preferred_scale Logical. Add an extra grid point at the
#'   data-preferred scale \eqn{\hat s}. Default: \code{FALSE}.
#' @param cores Integer thread count for each refit's chains. Default: the
#'   original fit's core count.
#' @param seed Integer base seed for the refits. Default: \code{1}.
#' @param keep_fits Logical. Retain the full refit objects in the result (for
#'   power users); the default keeps only per-scale summaries. Default:
#'   \code{FALSE}.
#'
#' @return An object of class \code{"bgms_prior_sensitivity"}: a list with the
#'   per-edge \code{edges} table (per-scale verdicts and \eqn{\log_{10}} BFs,
#'   chosen-scale verdict, stability range, mover category, insufficiency flag),
#'   a \code{grid} data frame (one row per refit with its convergence gate),
#'   the \code{log10_bf} scale-by-edge matrix, the \code{wobble} yardstick, the
#'   data-\code{preferred_scale}, and the settings used.
#'
#' @seealso \code{\link{bgm}()}, \code{\link{extract_posterior_inclusion_probabilities}()}
#' @family diagnostics
#' @references \insertRef{bartos2026}{bgms}
#' @export
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:6], chains = 2)
#' ps = prior_sensitivity_check(fit)
#' ps
#' plot(ps)
#' }
prior_sensitivity_check = function(bgms_object,
                                   scale_multipliers = c(0.5, 1, 2.5),
                                   evidence_threshold = 10,
                                   refit_sampler = "same-as-fit",
                                   iter = NULL,
                                   warmup = NULL,
                                   tolerance = 0.5,
                                   include_preferred_scale = FALSE,
                                   cores = NULL,
                                   seed = 1L,
                                   keep_fits = FALSE) {
  if(!inherits(bgms_object, "bgms")) {
    stop("prior_sensitivity_check() requires a fit from bgm().")
  }
  if(!is.numeric(evidence_threshold) || length(evidence_threshold) != 1L ||
     evidence_threshold <= 1) {
    stop("'evidence_threshold' must be a single number greater than 1.")
  }
  if(!is.numeric(scale_multipliers) || any(scale_multipliers <= 0)) {
    stop("'scale_multipliers' must be positive numbers.")
  }

  spec = get_fit_spec(bgms_object)
  if(is.null(spec) || !isTRUE(spec$prior$edge_selection)) {
    stop("prior_sensitivity_check() needs edge selection. Refit with ",
         "edge_selection = TRUE.")
  }

  chosen_scale = spec$prior$pairwise_scale
  lthr = log10(evidence_threshold)
  if(is.null(cores)) cores = spec$sampler$chains

  # --- Sampler resolution + speed recommendation ------------------------------
  rs = resolve_refit_sampler(bgms_object, refit_sampler)
  if(rs$recommend_nuts) {
    message("prior_sensitivity_check(): refits inherit the fit's '", rs$method,
            "' sampler and cost about one original fit each. For a faster ",
            "check pass refit_sampler = \"nuts\" (licensed because the check ",
            "compares refits only to refits).")
  }

  # --- Warm state (omrf) or cold refits (other models) ------------------------
  is_omrf = identical(spec$model_type, "omrf")
  warm_state = if(is_omrf) extract_warm_state(bgms_object) else NULL
  warm_metric = is_omrf && identical(rs$method, "nuts") &&
    !is.null(warm_state$inv_mass)
  if(!is_omrf) {
    message("prior_sensitivity_check(): warm starts are implemented for ",
            "ordinal models; this ", spec$model_type, " fit refits cold ",
            "(full warmup, about one fit per scale).")
  }

  # --- Data-preferred scale (no refit) ----------------------------------------
  preferred = data_preferred_scale(bgms_object)

  # --- Scale grid: multipliers (with 1) plus an s0 replicate ------------------
  multipliers = sort(unique(c(scale_multipliers, 1)))
  if(include_preferred_scale && is.finite(preferred$s_hat)) {
    multipliers = sort(unique(c(multipliers, preferred$s_hat / chosen_scale)))
  }
  scales = multipliers * chosen_scale
  rl = refit_run_length(bgms_object, rs$method, warm_metric, warmup, iter)

  # One refit per multiplier, plus a replicate at multiplier 1 (the wobble).
  jobs = data.frame(multiplier = c(multipliers, 1),
                    replicate = c(rep(FALSE, length(multipliers)), TRUE))
  refit_cores = min(as.integer(cores), spec$sampler$chains)

  refits = vector("list", nrow(jobs))
  walls = numeric(nrow(jobs))
  for(i in seq_len(nrow(jobs))) {
    t0 = Sys.time()
    refits[[i]] = refit_at_scale(
      bgms_object, scale = jobs$multiplier[i] * chosen_scale,
      warm_state = warm_state, warmup = rl$warmup, iter = rl$iter,
      seed = seed + i, cores = refit_cores, sampler = rs$method)
    walls[i] = as.numeric(Sys.time() - t0, units = "secs")
  }

  gates = lapply(refits, refit_convergence_gate)
  stats = lapply(refits, refit_edge_stats, evidence_threshold = evidence_threshold)

  grid_idx = seq_len(length(multipliers))          # non-replicate grid points
  rep_idx = length(multipliers) + 1L               # the s0 replicate
  s0_idx = which(multipliers == 1)                 # chosen-scale grid point
  usable = vapply(gates, `[[`, logical(1), "usable")
  edge_names = stats[[1]]$edge
  n_edges = length(edge_names)

  # --- Scale-by-edge log10 BF and verdict matrices (grid points only) ---------
  lbf_mat = t(vapply(stats[grid_idx], `[[`, numeric(n_edges), "lbf"))
  mcse_mat = t(vapply(stats[grid_idx], `[[`, numeric(n_edges), "mcse_lbf"))
  verdict_mat = t(vapply(stats[grid_idx], `[[`, character(n_edges), "verdict"))
  rownames(lbf_mat) = rownames(verdict_mat) = paste0("m", multipliers)

  # A refit that failed its gate must not contribute a verdict anywhere. Mask its
  # grid row to NA so it stays out of the count table, the verdict_x* columns, the
  # stability range, and the mover rule (which restricts to usable rows anyway).
  gate_fail = which(!usable[grid_idx])
  if(length(gate_fail) > 0) {
    lbf_mat[gate_fail, ] = NA_real_
    mcse_mat[gate_fail, ] = NA_real_
    verdict_mat[gate_fail, ] = NA_character_
  }
  s0_usable = usable[s0_idx]
  rep_usable = usable[rep_idx]

  # --- Wobble yardstick from the s0 refit vs its replicate --------------------
  # Only defined when both the chosen-scale refit and its replicate are usable.
  if(s0_usable && rep_usable) {
    d_wobble = abs(stats[[s0_idx]]$lbf - stats[[rep_idx]]$lbf)
    wobble_q95 = stats::quantile(d_wobble, 0.95, names = FALSE, na.rm = TRUE)
    wobble_med = stats::median(d_wobble, na.rm = TRUE)
  } else {
    warning("The chosen-scale (s0) refit or its replicate failed the ",
            "convergence gate; the wobble yardstick and chosen-scale verdicts ",
            "are reported as NA. Increase iter/warmup and rerun.", call. = FALSE)
    d_wobble = rep(NA_real_, n_edges)
    wobble_q95 = NA_real_
    wobble_med = NA_real_
  }

  # --- Per-edge chosen-scale sufficiency (from the s0 refit) ------------------
  # A near-constant edge has NA MCSE; treat missing uncertainty as zero (it is a
  # saturated/zero-flip edge, masked below anyway). If the s0 refit is unusable,
  # its verdicts are already NA-masked and sufficiency is undefined (NA).
  s0 = stats[[s0_idx]]
  if(s0_usable) {
    band = s0$band_half; band[is.na(band)] = 0
    mcse = s0$mcse_lbf; mcse[is.na(mcse)] = 0
    hw = pmax(band, 2 * mcse)
    straddle = (!s0$zeroflip) &
      (((s0$lbf - hw < lthr) & (s0$lbf + hw > lthr)) |
       ((s0$lbf - hw < -lthr) & (s0$lbf + hw > -lthr)))
    insufficient = (!s0$zeroflip) & (!s0$unanimous | straddle)
    insufficient[is.na(insufficient)] = FALSE
  } else {
    insufficient = rep(NA, n_edges)
  }

  # --- Mover category + stability range over the usable grid ------------------
  usable_grid = grid_idx[usable[grid_idx]]
  um = match(usable_grid, grid_idx)
  mover = character(n_edges)
  stability_lower = rep(NA_real_, n_edges)
  stability_upper = rep(NA_real_, n_edges)
  s0_col = match(s0_idx, grid_idx)
  for(e in seq_len(n_edges)) {
    v = verdict_mat[um, e]
    moved = length(unique(v[!is.na(v)])) > 1L
    lbf_e = lbf_mat[um, e]; lbf_e = lbf_e[is.finite(lbf_e)]
    dlbf = if(length(lbf_e)) diff(range(lbf_e)) else 0
    mcse_e = mcse_mat[um, e]; mcse_e = mcse_e[is.finite(mcse_e)]
    move_thr = max(c(tolerance, if(length(mcse_e)) 2 * max(mcse_e) else 0,
                     wobble_q95), na.rm = TRUE)
    mover[e] = if(!moved) "stable"
      else if(dlbf > move_thr) "moved-beyond-wobble"
      else "indistinguishable-from-wobble"
    si = stability_interval(verdict_mat[, e], multipliers, s0_col)
    stability_lower[e] = si[1]; stability_upper[e] = si[2]
  }

  # --- Edges table ------------------------------------------------------------
  # Chosen-scale columns read from the (masked) s0 grid row, so an unusable s0
  # refit leaves them NA rather than reporting uncertified verdicts.
  edges = data.frame(edge = edge_names,
                     prior_inclusion_probability = s0$prior_odds / (1 + s0$prior_odds),
                     chosen_scale_pip = if(s0_usable) s0$pip else NA_real_,
                     chosen_scale_log10_bf = lbf_mat[s0_col, ],
                     chosen_scale_mcse = mcse_mat[s0_col, ],
                     chosen_scale_verdict = verdict_mat[s0_col, ],
                     stability_lower = stability_lower,
                     stability_upper = stability_upper,
                     mover = mover,
                     insufficient = insufficient,
                     saturated = if(s0_usable) s0$zeroflip else NA,
                     stringsAsFactors = FALSE, row.names = NULL)
  vcols = as.data.frame(t(verdict_mat), stringsAsFactors = FALSE)
  names(vcols) = paste0("verdict_x", multipliers)
  edges = cbind(edges, vcols)

  # --- Grid table (per refit convergence gate) --------------------------------
  grid = data.frame(
    multiplier = jobs$multiplier,
    scale = jobs$multiplier * chosen_scale,
    replicate = jobs$replicate,
    usable = usable,
    rhat_continuous = vapply(gates, `[[`, numeric(1), "rhat_cont"),
    rhat_continuous_max = vapply(gates, `[[`, numeric(1), "rhat_cont_max"),
    ess_continuous = vapply(gates, `[[`, numeric(1), "ess_cont"),
    indicator_pair_ess = vapply(gates, `[[`, numeric(1), "pair_ess"),
    rb_median_rhat = vapply(gates, `[[`, numeric(1), "rb_med_rhat"),
    seconds = walls,
    row.names = NULL)

  # --- Learned-scale extras (hyperprior fits only) ----------------------------
  learned = NULL
  if(!is.null(spec$prior$interaction_scale_prior_type)) {
    learned = learned_scale_analysis(bgms_object, s0$prior_odds, evidence_threshold)
  }

  structure(
    list(
      edges = edges,
      grid = grid,
      multipliers = multipliers,
      scales = scales,
      chosen_scale = chosen_scale,
      chosen_index = s0_col,
      log10_bf = lbf_mat,
      log10_bf_mcse = mcse_mat,
      verdict = verdict_mat,
      wobble = list(q95 = wobble_q95, median = wobble_med, per_edge = d_wobble),
      preferred_scale = preferred,
      evidence_threshold = evidence_threshold,
      tolerance = tolerance,
      refit_sampler = rs$method,
      warm = warm_metric,
      model_type = spec$model_type,
      runtime_seconds = sum(walls),
      learned = learned,
      edge_names = edge_names,
      fits = if(keep_fits) refits else NULL
    ),
    class = "bgms_prior_sensitivity"
  )
}


# ------------------------------------------------------------------
# learned_scale_analysis
# ------------------------------------------------------------------
# For a fit run with an interaction_scale_prior, the scale-learning extras: the
# realized-vs-nominal multiplier interval and the scale-averaged verdict counts.
# Labeled in the report as a different analysis from the robustness grid.
# ------------------------------------------------------------------
learned_scale_analysis = function(fit, prior_odds, evidence_threshold) {
  spec = get_fit_spec(fit)
  raw = get_raw_samples(fit)
  shape = spec$prior$interaction_scale_shape
  rate = spec$prior$interaction_scale_rate
  chosen_scale = spec$prior$pairwise_scale
  scale_draws = do.call(c, raw$interaction_scale)
  u = scale_draws / chosen_scale
  probs = c(0.025, 0.5, 0.975)
  scale_prior_check = data.frame(
    quantity = c("mean", "2.5%", "50%", "97.5%"),
    nominal = c(shape / rate, stats::qgamma(probs, shape = shape, rate = rate)),
    realized = c(mean(u), stats::quantile(u, probs, names = FALSE)),
    row.names = NULL)
  gamma_pooled = do.call(rbind, raw$indicator)
  marg_pip = colMeans(gamma_pooled)
  marg_bf = (marg_pip / (1 - marg_pip)) / prior_odds
  marg_verdict = verdict_from_bf(marg_bf, evidence_threshold)
  list(scale_prior_check = scale_prior_check,
       marginalized_verdict = marg_verdict,
       counts = table(factor(marg_verdict,
         levels = c("presence", "undecided", "absence"))))
}


# ------------------------------------------------------------------
# order_upper_tri_rowmajor
# ------------------------------------------------------------------
# Permutation that reorders the column-major upper-triangle extraction
# matrix[upper.tri(matrix)] into the row-major edge order used by the
# raw pairwise / indicator columns and edge names.
#
# @param p  Number of variables.
#
# Returns: integer permutation of length p(p-1)/2.
# ------------------------------------------------------------------
order_upper_tri_rowmajor = function(p) {
  # Column-major upper triangle visits (i, j) with j outer, i < j inner.
  cm = which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  # Row-major target order sorts by row then column.
  order(cm[, "row"], cm[, "col"])
}


# ------------------------------------------------------------------
# stability_interval
# ------------------------------------------------------------------
# Maximal contiguous relative-scale range around the chosen scale over
# which the verdict is unchanged, using only reliable (non-NA) grid
# points. NA verdicts break the run.
#
# @param verdict     Per-grid-point verdict (may contain NA).
# @param relative    Relative scale at each grid point.
# @param chosen_idx  Grid index of the chosen scale.
#
# Returns: c(lower, upper) relative scale, or c(NA, NA).
# ------------------------------------------------------------------
stability_interval = function(verdict, relative, chosen_idx) {
  target = verdict[chosen_idx]
  if(is.na(target)) {
    return(c(NA_real_, NA_real_))
  }
  n = length(verdict)
  lo = chosen_idx
  while(lo - 1L >= 1L && !is.na(verdict[lo - 1L]) && verdict[lo - 1L] == target) {
    lo = lo - 1L
  }
  hi = chosen_idx
  while(hi + 1L <= n && !is.na(verdict[hi + 1L]) && verdict[hi + 1L] == target) {
    hi = hi + 1L
  }
  c(relative[lo], relative[hi])
}


# ==============================================================================
# Methods for bgms_prior_sensitivity
# ==============================================================================

#' @title Print a Prior Sensitivity Check
#'
#' @description
#' Prints the data-preferred scale, the per-scale verdict counts, the mover
#' summary with its wobble yardstick, and any unusable refits from a
#' \code{\link{prior_sensitivity_check}()} result.
#'
#' @param x A \code{bgms_prior_sensitivity} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{prior_sensitivity_check}()}
#' @family diagnostics
#' @export
print.bgms_prior_sensitivity = function(x, ...) {
  edges = x$edges
  n_edges = nrow(edges)

  cat("Prior sensitivity of edge inclusion verdicts\n\n")

  # The single most informative line: what scale the data prefer.
  ps = x$preferred_scale
  if(is.finite(ps$s_hat)) {
    cat(sprintf("Chosen scale %.3g; the data prefer approximately %.3g [%.3g, %.3g].\n",
                x$chosen_scale, ps$s_hat, ps$lo, ps$hi))
  } else {
    cat(sprintf("Chosen scale %.3g.\n", x$chosen_scale))
  }
  cat(sprintf("Evidence thresholds: presence >= %g, absence <= %.3g.\n\n",
              x$evidence_threshold, 1 / x$evidence_threshold))

  # Per-scale verdict counts, one column per multiplier (grid points only).
  cat("Verdict counts by slab scale (multiplier x chosen scale):\n")
  vt = apply(x$verdict, 1, function(v)
    c(presence = sum(v == "presence", na.rm = TRUE),
      undecided = sum(v == "undecided", na.rm = TRUE),
      absence = sum(v == "absence", na.rm = TRUE)))
  colnames(vt) = sprintf("%.2gx", x$multipliers)
  print(vt)
  cat("\n")

  # Mover summary with the wobble yardstick.
  mv = table(factor(edges$mover,
    levels = c("stable", "indistinguishable-from-wobble", "moved-beyond-wobble")))
  cat("Scale sensitivity (wobble q95 = ", format(x$wobble$q95, digits = 2),
      " log10 BF from the s0 replicate):\n", sep = "")
  cat("  stable:                        ", mv[["stable"]], "\n")
  cat("  indistinguishable-from-wobble: ", mv[["indistinguishable-from-wobble"]], "\n")
  cat("  moved-beyond-wobble:           ", mv[["moved-beyond-wobble"]], "\n")
  n_insuf = sum(edges$insufficient, na.rm = TRUE)
  if(n_insuf > 0) {
    cat(sprintf("  %d edge(s) insufficient at the chosen scale (chains disagree or\n", n_insuf))
    cat("    the Bayes-factor band straddles a threshold): cannot certify from this refit.\n")
  }
  cat("\n")

  # Unusable refits, if any.
  bad = x$grid[!x$grid$usable, ]
  if(nrow(bad) > 0) {
    cat("Unusable refit(s) (failed the convergence gate; excluded from verdicts):\n")
    for(i in seq_len(nrow(bad))) {
      cat(sprintf("  multiplier %.2g: median continuous Rhat %.3f, RB median Rhat %.3f\n",
                  bad$multiplier[i], bad$rhat_continuous[i], bad$rb_median_rhat[i]))
    }
    cat("\n")
  }

  cat(sprintf("Refit sampler: %s%s. Whole check: %.0f s (%d refits).\n",
              x$refit_sampler, if(x$warm) " (warm-started)" else "",
              x$runtime_seconds, nrow(x$grid)))

  # Learned-scale block (hyperprior fits only), clearly labeled.
  if(!is.null(x$learned)) {
    spc = x$learned$scale_prior_check
    lc = x$learned$counts
    cat("\nLearned-scale analysis (this fit used interaction_scale_prior; a\n")
    cat("different question from robustness):\n")
    cat(sprintf("  realized multiplier u = s/s0: mean %.2f (nominal %.2f), 95%% [%.2f, %.2f]\n",
                spc$realized[spc$quantity == "mean"], spc$nominal[spc$quantity == "mean"],
                spc$realized[spc$quantity == "2.5%"], spc$realized[spc$quantity == "97.5%"]))
    cat(sprintf("  scale-averaged verdicts: presence %d, undecided %d, absence %d\n",
                lc[["presence"]], lc[["undecided"]], lc[["absence"]]))
  }

  cat("\nUse plot() for the Bayes-factor-vs-scale lines, and $edges for the ",
      "per-edge table.\n", sep = "")
  invisible(x)
}


# ------------------------------------------------------------------
# verdict_palette
# ------------------------------------------------------------------
# Fixed colors for the three evidence verdicts.
# ------------------------------------------------------------------
verdict_palette = function() {
  c(presence = "#1b7837", undecided = "#7f7f7f", absence = "#b2182b")
}


# ------------------------------------------------------------------
# verdict_color
# ------------------------------------------------------------------
# Color for a verdict, defaulting to grey for NA / unknown (e.g. an
# edge whose chosen-scale grid point was grayed out at low ESS).
# ------------------------------------------------------------------
verdict_color = function(v, pal = verdict_palette()) {
  if(length(v) != 1L || is.na(v) || !v %in% names(pal)) {
    return("#7f7f7f")
  }
  pal[[v]]
}

#' @title Plot a Prior Sensitivity Check
#'
#' @description
#' Draws the \eqn{\log_{10}} inclusion-Bayes-factor lines against the relative
#' slab scale (multiplier of the chosen scale) on a log axis. The top panel
#' draws one line per edge through the grid points, with the undecided band
#' shaded, the chosen scale marked, the s0-replicate wobble band shown, and
#' moved-beyond-wobble edges at full opacity over faint stable edges. The bottom
#' panel is a forest of scale-stability ranges for the mover edges.
#'
#' @param x A \code{bgms_prior_sensitivity} object.
#' @param max_labels Integer. Maximum mover edges to label in the forest panel.
#'   Default: \code{25}.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly. Called for the side effect of drawing.
#'
#' @seealso \code{\link{prior_sensitivity_check}()}
#' @family diagnostics
#' @export
plot.bgms_prior_sensitivity = function(x, max_labels = 25L, ...) {
  pal = verdict_palette()
  rel = x$multipliers
  thr = log10(x$evidence_threshold)
  edges = x$edges
  n_edges = nrow(edges)
  movers = which(edges$mover == "moved-beyond-wobble" & !(edges$saturated %in% TRUE))

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(2, 1), mar = c(4.2, 4.2, 2.2, 1))

  # --- Top panel: log10 BF lines through the grid points ----------------------
  finite_bf = x$log10_bf[is.finite(x$log10_bf)]
  ylim = if(length(finite_bf)) range(c(finite_bf, -thr, thr)) else c(-thr, thr) * 2
  graphics::plot(NA, NA,
    xlim = range(rel), ylim = ylim, log = "x",
    xlab = "relative slab scale (multiplier of the chosen scale)",
    ylab = expression(log[10] * " inclusion Bayes factor"),
    main = "Inclusion Bayes factor vs. slab scale")
  graphics::rect(min(rel), -thr, max(rel), thr,
    col = grDevices::adjustcolor(pal[["undecided"]], 0.12), border = NA)
  graphics::abline(h = c(-thr, 0, thr), col = "grey70", lty = c(2, 1, 2))
  graphics::abline(v = 1, col = "grey40", lty = 3)

  # Stable edges faint, movers solid; x$log10_bf is scale-by-edge.
  for(e in seq_len(n_edges)) {
    if(isTRUE(edges$saturated[e])) next
    is_mover = edges$mover[e] == "moved-beyond-wobble"
    col_e = grDevices::adjustcolor(
      verdict_color(edges$chosen_scale_verdict[e], pal), if(is_mover) 0.9 else 0.2)
    graphics::lines(rel, x$log10_bf[, e], col = col_e, lwd = if(is_mover) 2 else 1)
  }
  n_saturated = sum(edges$saturated, na.rm = TRUE)
  if(n_saturated > 0) {
    graphics::mtext(paste0(n_saturated,
      " edge(s): presence or absence at every scale considered"),
      side = 3, line = 0.2, cex = 0.75, adj = 1, col = "grey40")
  }

  # --- Bottom panel: forest of stability ranges for movers --------------------
  if(length(movers) == 0L) {
    graphics::plot.new()
    graphics::text(0.5, 0.5,
      "No moved-beyond-wobble edges: every verdict is scale-robust.", cex = 1)
    return(invisible(x))
  }
  ord = movers[order(edges$stability_upper[movers] - edges$stability_lower[movers])]
  show = utils::head(ord, max_labels)
  yy = seq_along(show)
  xr = range(c(edges$stability_lower[show], edges$stability_upper[show], 1), na.rm = TRUE)
  graphics::plot(NA, NA,
    xlim = xr, ylim = c(0.5, length(show) + 0.5),
    log = "x", yaxt = "n", xlab = "relative slab scale",
    ylab = "", main = "Scale-stability ranges (mover edges)")
  graphics::abline(v = 1, col = "grey40", lty = 3)
  for(k in seq_along(show)) {
    e = show[k]
    graphics::segments(edges$stability_lower[e], yy[k],
      edges$stability_upper[e], yy[k],
      col = verdict_color(edges$chosen_scale_verdict[e], pal), lwd = 3)
    graphics::points(1, yy[k], pch = 18,
      col = verdict_color(edges$chosen_scale_verdict[e], pal))
  }
  graphics::axis(2, at = yy, labels = edges$edge[show], las = 1, cex.axis = 0.7)
  invisible(x)
}
