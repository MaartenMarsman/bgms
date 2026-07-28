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
#' s0-replicate spread over threshold-relevant edges
#' (\eqn{|\log_{10} \mathrm{BF}| \le 3} at s0; near-saturated edges would
#' inflate it). The \code{$edges$mover} column stores \code{stable},
#' \code{indistinguishable-from-wobble}, or \code{moved-beyond-wobble}; the
#' printed report shows the same categories in plain language ("robust",
#' "changed, within run-to-run noise", "changed, beyond run-to-run noise",
#' with edges failing the sufficiency check below printed as "not
#' certifiable"). A bare verdict flip inside the replicate noise is never
#' reported as a move.
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
#'   fit carry the adapted metric and run fastest; when an inherited slower
#'   sampler makes the check cost more than about a minute, a message suggests
#'   the switch. Refitting with a different sampler than the original fit is
#'   sound because every comparison the check makes runs refit-against-refit
#'   under one identical pipeline; the original fit never enters a comparison.
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
    stop(
      "prior_sensitivity_check() needs edge selection. Refit with ",
      "edge_selection = TRUE."
    )
  }

  chosen_scale = spec$prior$pairwise_scale
  lthr = log10(evidence_threshold)
  if(is.null(cores)) cores = spec$sampler$chains

  # --- Sampler resolution -----------------------------------------------------
  rs = resolve_refit_sampler(bgms_object, refit_sampler)

  # --- Warm state (omrf) or cold refits (other models) ------------------------
  is_omrf = identical(spec$model_type, "omrf")
  warm_state = if(is_omrf) extract_warm_state(bgms_object) else NULL
  warm_metric = is_omrf && identical(rs$method, "nuts") &&
    !is.null(warm_state$inv_mass)

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
  jobs = data.frame(
    multiplier = c(multipliers, 1),
    replicate = c(rep(FALSE, length(multipliers)), TRUE)
  )
  refit_cores = min(as.integer(cores), spec$sampler$chains)

  refits = vector("list", nrow(jobs))
  walls = numeric(nrow(jobs))
  for(i in seq_len(nrow(jobs))) {
    t0 = Sys.time()
    refits[[i]] = refit_at_scale(
      bgms_object,
      scale = jobs$multiplier[i] * chosen_scale,
      warm_state = warm_state, warmup = rl$warmup, iter = rl$iter,
      seed = seed + i, cores = refit_cores, sampler = rs$method
    )
    walls[i] = as.numeric(Sys.time() - t0, units = "secs")
    # Suggest the faster sampler only when the cost is material (projected
    # total above ~60 s); see the man page for why the switch is sound.
    if(i == 1L && rs$recommend_nuts && walls[1] * nrow(jobs) > 60) {
      message(
        "Refits use the fit's ", rs$method, " sampler (~",
        round(walls[1]), " s each); refit_sampler = \"nuts\" is usually faster."
      )
    }
  }

  gates = lapply(refits, refit_convergence_gate)
  stats = lapply(refits, refit_edge_stats, evidence_threshold = evidence_threshold)

  grid_idx = seq_along(multipliers) # non-replicate grid points
  rep_idx = length(multipliers) + 1L # the s0 replicate
  s0_idx = which(multipliers == 1) # chosen-scale grid point
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
    wob = wobble_yardstick(stats[[s0_idx]]$lbf, stats[[rep_idx]]$lbf)
    d_wobble = wob$per_edge
    wobble_q95 = wob$q95
    wobble_med = wob$median
  } else {
    warning("The chosen-scale (s0) refit or its replicate failed the ",
      "convergence gate; the wobble yardstick and chosen-scale verdicts ",
      "are reported as NA. Increase iter/warmup and rerun.",
      call. = FALSE
    )
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
    band = s0$band_half
    band[is.na(band)] = 0
    mcse = s0$mcse_lbf
    mcse[is.na(mcse)] = 0
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
    lbf_e = lbf_mat[um, e]
    lbf_e = lbf_e[is.finite(lbf_e)]
    dlbf = if(length(lbf_e)) diff(range(lbf_e)) else 0
    mcse_e = mcse_mat[um, e]
    mcse_e = mcse_e[is.finite(mcse_e)]
    move_thr = max(c(
      tolerance, if(length(mcse_e)) 2 * max(mcse_e) else 0,
      wobble_q95
    ), na.rm = TRUE)
    mover[e] = if(!moved) {
      "stable"
    } else if(dlbf > move_thr) {
      "moved-beyond-wobble"
    } else {
      "indistinguishable-from-wobble"
    }
    si = stability_interval(verdict_mat[, e], multipliers, s0_col)
    stability_lower[e] = si[1]
    stability_upper[e] = si[2]
  }

  # --- Edges table ------------------------------------------------------------
  # Chosen-scale columns read from the (masked) s0 grid row, so an unusable s0
  # refit leaves them NA rather than reporting uncertified verdicts.
  edges = data.frame(
    edge = edge_names,
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
    stringsAsFactors = FALSE, row.names = NULL
  )
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
    row.names = NULL
  )

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
      edge_names = edge_names,
      fits = if(keep_fits) refits else NULL
    ),
    class = "bgms_prior_sensitivity"
  )
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
# wobble_yardstick
# ------------------------------------------------------------------
# Monte Carlo wobble from the chosen-scale refit and its replicate. The
# q95 pools threshold-relevant edges only (|log10 BF| <= 3 at s0):
# near-saturated edges have an exploding log-odds derivative, and their
# replicate spread would inflate the yardstick past any genuine
# scale-driven move. The per-edge spread and its median cover all edges.
#
# @param lbf_s0   Per-edge log10 BF from the chosen-scale refit.
# @param lbf_rep  Per-edge log10 BF from the s0 replicate.
#
# Returns: list(q95, median, per_edge); q95 is NA when no edge is
# threshold-relevant.
# ------------------------------------------------------------------
wobble_yardstick = function(lbf_s0, lbf_rep) {
  per_edge = abs(lbf_s0 - lbf_rep)
  relevant = which(abs(lbf_s0) <= 3)
  list(
    q95 = stats::quantile(per_edge[relevant], 0.95, names = FALSE, na.rm = TRUE),
    median = stats::median(per_edge, na.rm = TRUE),
    per_edge = per_edge
  )
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
#' Prints the answer first: how many edge verdicts are robust to the slab
#' scale, which edges genuinely depend on it (by name), which cannot be
#' certified, the verdict counts per scale, and how the chosen scale compares
#' with the size of the estimated interactions. Machinery (sampler, refit
#' cost, the run-to-run noise band) is confined to a closing Details line.
#'
#' @param x A \code{bgms_prior_sensitivity} object.
#' @param max_rows Integer. Maximum edges to name in the scale-dependent
#'   table; the rest are counted and left to \code{$edges}. Default: \code{10}.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{prior_sensitivity_check}()}
#' @family diagnostics
#' @export
print.bgms_prior_sensitivity = function(x, max_rows = 10L, ...) {
  edges = x$edges
  n_edges = nrow(edges)
  mlab = sprintf("%.2gx", x$multipliers)
  mlist = paste0(
    paste(mlab[-length(mlab)], collapse = ", "), ", and ", mlab[length(mlab)]
  )

  cat("Prior sensitivity check: are the edge verdicts robust to the slab scale?\n")
  cat(sprintf("Refits at %s the chosen scale; %d edges.\n\n", mlist, n_edges))

  # Without a converged chosen-scale refit and replicate there is no yardstick
  # and no certified verdict; say so instead of tabulating NAs.
  grid_rows = x$grid[!x$grid$replicate, ]
  s0_ok = grid_rows$usable[x$chosen_index] && any(x$grid$usable[x$grid$replicate])
  if(!s0_ok) {
    cat("The chosen-scale refit or its repeat run did not converge, so no\n")
    cat("verdict could be certified from this run. Increase iter/warmup\n")
    cat("(or pass refit_sampler = \"nuts\") and rerun.\n\n")
    cat(sprintf(
      "Details: %s refits (%s), %d refits in %.0f s.\n",
      x$refit_sampler,
      if(x$warm) "warm-started from the fit" else "cold, full warmup",
      nrow(x$grid), x$runtime_seconds
    ))
    return(invisible(x))
  }

  # One category per edge: not-certifiable wins, then the mover category.
  uncert = edges$insufficient %in% TRUE
  robust = !uncert & edges$mover == "stable"
  within = !uncert & edges$mover == "indistinguishable-from-wobble"
  beyond = !uncert & edges$mover == "moved-beyond-wobble"
  counts = c(sum(robust), sum(within), sum(beyond), sum(uncert))
  labels = c(
    "robust (same verdict at every scale)",
    "changed, within run-to-run noise",
    "changed, beyond run-to-run noise",
    "not certifiable (chains disagree)"
  )
  keep = counts > 0 | seq_along(counts) <= 3L
  cat(sprintf("  %-38s %4d\n", labels[keep], counts[keep]), sep = "")
  cat("\n")

  # Name the edges whose verdict genuinely depends on the scale.
  if(sum(beyond) > 0) {
    cat(if(sum(beyond) == 1L) {
      "1 edge's verdict genuinely depends on the scale:\n"
    } else {
      sprintf("%d edges' verdicts genuinely depend on the scale:\n", sum(beyond))
    })
    idx = which(beyond)
    show = utils::head(idx, max_rows)
    vmat = t(x$verdict[, show, drop = FALSE])
    tab = cbind(edge = edges$edge[show], vmat)
    colnames(tab) = c("edge", mlab)
    rownames(tab) = rep("", nrow(tab))
    print(tab, quote = FALSE, print.gap = 2)
    if(length(idx) > length(show)) {
      cat(sprintf("  ...and %d more; see $edges.\n", length(idx) - length(show)))
    }
    cat("\n")
  }

  # Name the edges the refits cannot certify.
  if(sum(uncert) > 0) {
    nm = edges$edge[uncert]
    shown = utils::head(nm, max_rows)
    cat(sprintf(
      "%d edge%s cannot be certified from this run: %s%s\n\n",
      length(nm), if(length(nm) == 1L) "" else "s", paste(shown, collapse = ", "),
      if(length(nm) > length(shown)) {
        sprintf(", and %d more (see $edges).", length(nm) - length(shown))
      } else {
        "."
      }
    ))
  }

  # Verdict counts by scale.
  cat("Verdict counts by scale:\n")
  vt = apply(x$verdict, 1, function(v) {
    c(
      presence = sum(v == "presence", na.rm = TRUE),
      undecided = sum(v == "undecided", na.rm = TRUE),
      absence = sum(v == "absence", na.rm = TRUE)
    )
  })
  colnames(vt) = mlab
  out = utils::capture.output(print(vt))
  cat(paste0("  ", out, collapse = "\n"), "\n", sep = "")
  ab = vt["absence", ]
  if(ab[length(ab)] > ab[1]) {
    cat("More absence at wider scales is expected: a wider slab strengthens\n")
    cat("evidence against borderline edges.\n")
  }
  cat("\n")

  # Refits that failed their convergence check.
  bad = x$grid[!x$grid$usable & !x$grid$replicate, ]
  if(nrow(bad) > 0) {
    bl = sprintf("%.2gx", bad$multiplier)
    bl = if(length(bl) > 1L) {
      paste0(paste(bl[-length(bl)], collapse = ", "), " and ", bl[length(bl)])
    } else {
      bl
    }
    cat(sprintf(
      "The %s refit%s did not converge and %s excluded from the verdicts.\n\n",
      bl, if(nrow(bad) == 1L) "" else "s", if(nrow(bad) == 1L) "is" else "are"
    ))
  }

  # How the chosen scale compares with the estimated interactions.
  ps = x$preferred_scale
  if(is.finite(ps$s_hat)) {
    ratio = x$chosen_scale / ps$s_hat
    if(ratio > 2) {
      cat(sprintf(
        "Note: the chosen scale (%.3g) is much wider than the estimated\ninteractions (about %.3g [%.3g, %.3g]); absence verdicts in\nparticular depend on this choice.\n\n",
        x$chosen_scale, ps$s_hat, ps$lo, ps$hi
      ))
    } else if(ratio < 0.5) {
      cat(sprintf(
        "Note: the chosen scale (%.3g) is much narrower than the estimated\ninteractions (about %.3g [%.3g, %.3g]); presence verdicts in\nparticular depend on this choice.\n\n",
        x$chosen_scale, ps$s_hat, ps$lo, ps$hi
      ))
    } else {
      cat(sprintf(
        "The chosen scale (%.3g) matches the size of the estimated interactions (about %.3g).\n\n",
        x$chosen_scale, ps$s_hat
      ))
    }
  }

  secs = if(x$runtime_seconds < 10) {
    sprintf("%.1f s", x$runtime_seconds)
  } else {
    sprintf("%.0f s", x$runtime_seconds)
  }
  cat(sprintf(
    "Details: %s refits (%s), %d refits in %s;\nrun-to-run noise band %.2g log10 BF (95th pct of the repeated\nchosen-scale refit). ?prior_sensitivity_check for how to read this.\n",
    x$refit_sampler,
    if(x$warm) "warm-started from the fit" else "cold, full warmup",
    nrow(x$grid), secs, x$wobble$q95
  ))
  invisible(x)
}


# ------------------------------------------------------------------
# mover_palette
# ------------------------------------------------------------------
# Fixed-order categorical colors for the named mover lines (Okabe-Ito,
# colorblind-safe; warm/cool alternating for adjacent separation). Every
# mover line is also name-labeled, so identity never rides on color alone.
# ------------------------------------------------------------------
mover_palette = function() {
  c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#882255", "#44AA99", "#661100", "#999933"
  )
}


# ------------------------------------------------------------------
# spread_labels
# ------------------------------------------------------------------
# Nudge label y-positions apart until adjacent labels are at least `gap`
# apart, preserving order. Returns the adjusted positions.
# ------------------------------------------------------------------
spread_labels = function(y, gap) {
  ord = order(y)
  ys = y[ord]
  for(k in seq_along(ys)[-1]) {
    if(ys[k] - ys[k - 1] < gap) ys[k] = ys[k - 1] + gap
  }
  out = y
  out[ord] = ys
  out
}


#' @title Plot a Prior Sensitivity Check
#'
#' @description
#' One panel, answer first: the title states how many edge verdicts depend on
#' the slab scale. Each edge's \eqn{\log_{10}} inclusion Bayes factor (the
#' evidence for the edge) is drawn across the refit scales. Edges whose
#' verdict genuinely depends on the scale are colored and labeled by name;
#' all other edges are the muted background. The shaded band is the undecided
#' zone between the evidence thresholds; the zones are labeled in the right
#' margin.
#'
#' @param x A \code{bgms_prior_sensitivity} object.
#' @param max_labels Integer. Maximum scale-dependent edges to color and
#'   label by name; the rest are counted in a corner note. Default: \code{10}.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly. Called for the side effect of drawing.
#'
#' @seealso \code{\link{prior_sensitivity_check}()}
#' @family diagnostics
#' @export
plot.bgms_prior_sensitivity = function(x, max_labels = 10L, ...) {
  rel = x$multipliers
  thr = log10(x$evidence_threshold)
  edges = x$edges
  n_edges = nrow(edges)
  mlab = sprintf("%.2gx", rel)
  ink = "grey25"
  muted = "grey55"
  faint = grDevices::adjustcolor("grey55", 0.35)

  # Same partition as print(): an edge the refits cannot certify is not
  # reported as scale-dependent.
  movers = which(
    edges$mover == "moved-beyond-wobble" &
      !(edges$saturated %in% TRUE) & !(edges$insufficient %in% TRUE)
  )
  # Label the largest evidence swings first.
  swing = apply(x$log10_bf, 2, function(v) {
    v = v[is.finite(v)]
    if(length(v)) diff(range(v)) else 0
  })
  movers = movers[order(swing[movers], decreasing = TRUE)]
  named = utils::head(movers, max_labels)

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mar = c(5.1, 4.4, 3.6, 7.5), mgp = c(2.6, 0.7, 0),
    col.axis = ink, col.lab = ink, col.main = ink
  )

  # The y-range covers the decision band and every named line; edges that
  # live entirely outside it never change verdict and clip silently.
  named_bf = x$log10_bf[, named, drop = FALSE]
  named_bf = named_bf[is.finite(named_bf)]
  ylim = range(c(-1.5 * thr, 1.5 * thr, named_bf))
  ylim = ylim + c(-0.4, 0.4)

  n_dep = length(movers)
  title = if(n_dep == 0L) {
    "No edge verdict depends on the slab scale"
  } else if(n_dep == 1L) {
    "1 edge verdict depends on the slab scale"
  } else {
    sprintf("%d edge verdicts depend on the slab scale", n_dep)
  }

  graphics::plot(NA, NA,
    xlim = range(rel), ylim = ylim, log = "x", axes = FALSE,
    xlab = "slab scale (relative to the chosen scale)",
    ylab = expression("evidence for the edge (" * log[10] * " Bayes factor)"),
    main = title
  )
  graphics::axis(1, at = rel, labels = mlab, col = muted, col.ticks = muted)
  graphics::axis(2, col = muted, col.ticks = muted, las = 1)

  # Verdict zones: shaded undecided band, dashed thresholds, margin labels.
  usr = graphics::par("usr")
  graphics::rect(10^usr[1], -thr, 10^usr[2], thr,
    col = grDevices::adjustcolor("grey60", 0.12), border = NA
  )
  graphics::abline(h = c(-thr, thr), col = muted, lty = 2)
  graphics::abline(v = 1, col = muted, lty = 3)
  graphics::mtext("chosen", side = 3, at = 1, line = 0.1, cex = 0.75, col = muted)
  zone_x = 10^(usr[1] + 0.015 * diff(usr[1:2]))
  pad = 0.05 * diff(ylim)
  graphics::text(zone_x, thr + pad, "presence", adj = c(0, 0), cex = 0.8, col = muted)
  graphics::text(zone_x, thr - pad, "undecided", adj = c(0, 1), cex = 0.8, col = muted)
  graphics::text(zone_x, -thr - pad, "absence", adj = c(0, 1), cex = 0.8, col = muted)

  # Background: every other edge in one muted color.
  for(e in setdiff(seq_len(n_edges), named)) {
    if(isTRUE(edges$saturated[e])) next
    graphics::lines(rel, x$log10_bf[, e], col = faint, lwd = 1)
  }

  # Foreground: the named movers, colored in fixed order and name-labeled.
  if(length(named)) {
    pal = mover_palette()
    # Anchor each name at the line's last finite point (0 if none is finite).
    anchor = vapply(named, function(e) {
      v = x$log10_bf[, e]
      v = v[is.finite(v)]
      if(length(v)) v[length(v)] else 0
    }, numeric(1))
    end_y = spread_labels(anchor, gap = 0.05 * diff(ylim))
    for(k in seq_along(named)) {
      e = named[k]
      graphics::lines(rel, x$log10_bf[, e], col = pal[k], lwd = 2.5)
      graphics::points(rel, x$log10_bf[, e], col = pal[k], pch = 16, cex = 0.9)
      graphics::segments(
        max(rel), anchor[k],
        10^(usr[2] + 0.005 * diff(usr[1:2])), end_y[k],
        col = grDevices::adjustcolor(pal[k], 0.5), lwd = 0.8, xpd = NA
      )
      graphics::text(
        10^(usr[2] + 0.015 * diff(usr[1:2])), end_y[k], edges$edge[e],
        xpd = NA, adj = 0, cex = 0.75, col = ink
      )
    }
  }

  # Corner notes: unnamed movers and off-scale edges, counted not drawn.
  notes = character(0)
  if(n_dep > length(named)) {
    notes = c(notes, sprintf("and %d more; see $edges", n_dep - length(named)))
  }
  n_off = sum(vapply(seq_len(n_edges), function(e) {
    v = x$log10_bf[, e]
    v = v[is.finite(v)]
    length(v) == 0L || all(v > ylim[2]) || all(v < ylim[1])
  }, logical(1)))
  if(n_off > 0) {
    notes = c(notes, sprintf(
      "%d edges beyond the plot range keep their verdict at every scale", n_off
    ))
  }
  if(length(notes)) {
    graphics::mtext(paste(notes, collapse = "; "),
      side = 1, line = 3.9, adj = 0, cex = 0.75, col = muted
    )
  }
  invisible(x)
}
