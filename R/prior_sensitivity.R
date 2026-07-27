# ==============================================================================
# prior_sensitivity_check: scale-robustness of edge inclusion verdicts
# ==============================================================================
#
# From a single fit made with a random interaction slab scale
# (bgm(..., interaction_scale_prior = )), report for every edge (a) the
# inclusion Bayes factor already marginalized over the scale posterior and (b)
# the whole inclusion-Bayes-factor-versus-scale curve, recovered post hoc by
# conditional-density reweighting. See the package vignette
# "Checking prior sensitivity" and Bartos et al. (arXiv 2604.21596).
# ==============================================================================


# ------------------------------------------------------------------
# slab_log_density
# ------------------------------------------------------------------
# Log slab density at a scale, switched on the slab family so Cauchy
# support is a one-line addition.
#
# @param theta   Numeric (vector or matrix) of interaction values.
# @param scale   Positive slab scale.
# @param family  "normal" or "cauchy".
#
# Returns: log density, same shape as theta.
# ------------------------------------------------------------------
slab_log_density = function(theta, scale, family) {
  switch(family,
    normal = stats::dnorm(theta, 0, scale, log = TRUE),
    cauchy = stats::dcauchy(theta, 0, scale, log = TRUE),
    stop(
      "prior_sensitivity_check() supports a normal or cauchy slab; got '",
      family, "'."
    )
  )
}


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
#' Assesses whether each edge's inclusion verdict is robust to the interaction
#' (pairwise) slab scale, from a single model fitted with a random-scale
#' hyperprior (\code{bgm(..., interaction_scale_prior = )}). Reports, per edge,
#' the inclusion Bayes factor marginalized over the scale posterior (the
#' headline scale-robust verdict) and the whole inclusion-Bayes-factor-versus-
#' scale curve, recovered by conditional-density reweighting of the single fit.
#'
#' @details
#' The marginalized inclusion probability is the mean of the indicator draws,
#' which already integrates over the scale posterior because the edge prior
#' does not involve the scale. The inclusion Bayes factor divides the posterior
#' inclusion odds by the prior inclusion odds from
#' \code{\link{extract_prior_inclusion_probabilities}()} (never a hard-coded
#' \eqn{1/2}), so a non-\eqn{1/2} edge prior is handled correctly.
#'
#' The fixed-scale curve at scale \eqn{s} is estimated by reweighting each draw
#' by the conditional density of \eqn{s} given that draw's included
#' interactions, \eqn{p(s \mid \theta, \gamma) \propto \pi(s) \prod_{j\
#' included} \mathrm{slab}(\theta_j; s)}. Grid points whose importance ESS
#' \eqn{(\sum w)^2 / \sum w^2} falls below \code{ess_floor} are unreliable and
#' returned as \code{NA} (the reliability window). Edges that are always or
#' never selected are reported as "presence at every scale considered" (or
#' absence), never as numbers.
#'
#' @param bgms_object A fitted \code{bgms} object from \code{\link{bgm}()} run
#'   with an \code{interaction_scale_prior}.
#' @param evidence_threshold Positive numeric. Inclusion Bayes factor threshold
#'   for a presence verdict; \code{1 / evidence_threshold} is the absence
#'   threshold. Default: \code{10}.
#' @param grid Optional numeric vector of absolute slab scales at which to
#'   evaluate the curve. Default \code{NULL} uses \code{n_grid} log-spaced
#'   points spanning the central 95\% of the scale posterior, with the chosen
#'   scale inserted.
#' @param n_grid Integer. Number of grid points when \code{grid} is \code{NULL}.
#'   Default: \code{25}.
#' @param ess_floor Numeric. Minimum importance ESS for a grid point to be
#'   reported; points below it are returned as \code{NA}. Default: \code{400}.
#'
#' @return An object of class \code{"bgms_prior_sensitivity"}: a list with
#'   \describe{
#'     \item{edges}{A data frame, one row per edge, with the marginalized
#'       inclusion probability and Bayes factor, the marginalized and
#'       chosen-scale verdicts, the scale-stability interval (relative to the
#'       chosen scale), and \code{mover} / \code{saturated} flags.}
#'     \item{grid}{A data frame with the absolute \code{scale}, the
#'       \code{relative} scale (scale / chosen scale), and the importance
#'       \code{ess} at each grid point.}
#'     \item{log10_bf}{A grid-by-edge matrix of \eqn{\log_{10}} inclusion Bayes
#'       factors, \code{NA} where the ESS is below the floor or the edge is
#'       saturated.}
#'     \item{log10_bf_by_chain}{A list of the same matrix per chain, for the
#'       Monte Carlo fragility check.}
#'     \item{chosen_scale}{The scale from \code{interaction_prior}.}
#'     \item{evidence_threshold, ess_floor}{The settings used.}
#'   }
#'
#' @seealso \code{\link{bgm}()}, \code{\link{extract_scale_draws}()},
#'   \code{\link{extract_prior_inclusion_probabilities}()}
#' @family diagnostics
#' @references \insertRef{bartos2026}{bgms}
#' @export
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:6],
#'   interaction_scale_prior = gamma_prior(2, 2),
#'   chains = 2
#' )
#' ps = prior_sensitivity_check(fit)
#' ps
#' plot(ps)
#' }
prior_sensitivity_check = function(bgms_object,
                                   evidence_threshold = 10,
                                   grid = NULL,
                                   n_grid = 25L,
                                   ess_floor = 400) {
  if(!inherits(bgms_object, "bgms")) {
    stop("prior_sensitivity_check() requires a fit from bgm().")
  }
  if(!is.numeric(evidence_threshold) || length(evidence_threshold) != 1L ||
    evidence_threshold <= 1) {
    stop("'evidence_threshold' must be a single number greater than 1.")
  }

  spec = get_fit_spec(bgms_object)
  if(is.null(spec) || is.null(spec$prior$interaction_scale_prior_type)) {
    stop(
      "This fit used a fixed interaction slab scale, so there is nothing to ",
      "check. Refit with a random scale, for example bgm(..., ",
      "interaction_scale_prior = gamma_prior(shape = 2, rate = 2))."
    )
  }
  if(!isTRUE(spec$prior$edge_selection)) {
    stop(
      "prior_sensitivity_check() needs the indicator draws from edge ",
      "selection. Refit with edge_selection = TRUE."
    )
  }

  raw = get_raw_samples(bgms_object)
  if(is.null(raw$interaction_scale) || is.null(raw$indicator)) {
    stop(
      "This fit does not store the scale and indicator draws needed for the ",
      "check. Refit with the current bgms version, edge_selection = TRUE, and ",
      "an interaction_scale_prior."
    )
  }

  # --- Fit-level quantities ---------------------------------------------------
  family = spec$prior$interaction_prior_type
  chosen_scale = spec$prior$pairwise_scale
  shape = spec$prior$interaction_scale_shape
  rate = spec$prior$interaction_scale_rate
  edge_names = raw$parameter_names$indicator %||% raw$parameter_names$pairwise

  # Per-edge prior inclusion odds (never a hard-coded 1/2). The edge prior does
  # not involve the scale, so one value per edge serves the whole curve.
  prior_pip_mat = extract_prior_inclusion_probabilities(bgms_object)
  num_vars = nrow(prior_pip_mat)
  prior_q = prior_pip_mat[upper.tri(prior_pip_mat)][
    order_upper_tri_rowmajor(num_vars)
  ]
  prior_odds = prior_q / (1 - prior_q)

  # --- Evaluation grid --------------------------------------------------------
  scale_draws = do.call(c, raw$interaction_scale)
  if(is.null(grid)) {
    span = stats::quantile(scale_draws, c(0.025, 0.975), names = FALSE)
    grid = exp(seq(log(span[1]), log(span[2]), length.out = n_grid))
  }
  grid = sort(unique(c(grid, chosen_scale)))
  n_grid = length(grid)

  # --- Reweighting per chain and pooled ---------------------------------------
  per_chain = lapply(seq_along(raw$indicator), function(ci) {
    reweight_curve(
      theta = raw$pairwise[[ci]],
      gamma = raw$indicator[[ci]],
      grid = grid, chosen_scale = chosen_scale,
      shape = shape, rate = rate, family = family,
      prior_odds = prior_odds
    )
  })

  theta_pooled = do.call(rbind, raw$pairwise)
  gamma_pooled = do.call(rbind, raw$indicator)
  pooled = reweight_curve(
    theta = theta_pooled, gamma = gamma_pooled,
    grid = grid, chosen_scale = chosen_scale,
    shape = shape, rate = rate, family = family,
    prior_odds = prior_odds
  )

  # --- Marginalized (scale-averaged) verdict ----------------------------------
  marg_pip = colMeans(gamma_pooled)
  marg_bf = (marg_pip / (1 - marg_pip)) / prior_odds
  marg_verdict = verdict_from_bf(marg_bf, evidence_threshold)

  # Saturated edges: always or never selected. The curve is uninformative and
  # the verdict is scale-robust by construction.
  saturated = (marg_pip <= 0) | (marg_pip >= 1)

  # --- Gray out low-ESS grid points and saturated edges -----------------------
  reliable = pooled$ess >= ess_floor
  log10_bf = pooled$log10_bf
  log10_bf[!reliable, ] = NA_real_
  log10_bf[, saturated] = NA_real_

  # --- Per-edge verdicts, stability interval, mover flag ----------------------
  chosen_idx = which.min(abs(grid - chosen_scale))
  relative = grid / chosen_scale

  n_edges = length(edge_names)
  chosen_verdict = character(n_edges)
  stability_lower = rep(NA_real_, n_edges)
  stability_upper = rep(NA_real_, n_edges)
  mover = logical(n_edges)

  for(e in seq_len(n_edges)) {
    if(saturated[e]) {
      chosen_verdict[e] = if(marg_pip[e] >= 1) "presence" else "absence"
      next
    }
    v = verdict_from_bf(log10_to_bf(log10_bf[, e]), evidence_threshold)
    chosen_verdict[e] = v[chosen_idx]
    si = stability_interval(v, relative, chosen_idx)
    stability_lower[e] = si[1]
    stability_upper[e] = si[2]
    reliable_v = v[!is.na(v)]
    mover[e] = (length(unique(reliable_v)) > 1L) ||
      (!is.na(chosen_verdict[e]) && chosen_verdict[e] != marg_verdict[e])
  }

  edges = data.frame(
    edge = edge_names,
    prior_inclusion_probability = prior_q,
    marginalized_inclusion_probability = marg_pip,
    marginalized_bf = marg_bf,
    marginalized_verdict = marg_verdict,
    chosen_scale_verdict = chosen_verdict,
    stability_lower = stability_lower,
    stability_upper = stability_upper,
    mover = mover,
    saturated = saturated,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  structure(
    list(
      edges = edges,
      grid = data.frame(scale = grid, relative = relative, ess = pooled$ess),
      log10_bf = log10_bf,
      log10_bf_by_chain = lapply(per_chain, function(pc) pc$log10_bf),
      chosen_scale = chosen_scale,
      chosen_index = chosen_idx,
      evidence_threshold = evidence_threshold,
      ess_floor = ess_floor,
      slab_family = family,
      edge_names = edge_names
    ),
    class = "bgms_prior_sensitivity"
  )
}


# ------------------------------------------------------------------
# reweight_curve
# ------------------------------------------------------------------
# Conditional-density reweighting of a set of draws onto a scale grid.
# For each draw t the conditional density of the scale given that
# draw's included interactions is p(s | theta_t, gamma_t) proportional
# to pi(s) prod_{included} slab(theta_j; s); normalizing it over the
# grid (log-space trapezoid) gives a per-draw distribution p_t over the
# grid. The conditional inclusion probability at grid point x is
# sum_t gamma_{e,t} p_t(x) / sum_t p_t(x).
#
# @param theta        niter x n_edges matrix of interaction draws.
# @param gamma        niter x n_edges 0/1 indicator draws.
# @param grid         Absolute scale grid.
# @param chosen_scale Centering scale s0.
# @param shape,rate   Mean-1 gamma hyperprior on u = s / s0.
# @param family       Slab family ("normal"/"cauchy").
# @param prior_odds   Per-edge prior inclusion odds.
#
# Returns: list(log10_bf = n_grid x n_edges, ess = n_grid).
# ------------------------------------------------------------------
reweight_curve = function(theta, gamma, grid, chosen_scale,
                          shape, rate, family, prior_odds) {
  niter = nrow(gamma)
  n_grid = length(grid)

  # Log unnormalized conditional density of the scale per draw and grid point,
  # including the mean-1 hyperprior on u = s / s0 and the log-space quadrature
  # weight log(s). Constants that do not depend on the draw cancel in the ratio.
  lw = matrix(NA_real_, nrow = niter, ncol = n_grid)
  log_hyper = stats::dgamma(grid / chosen_scale,
    shape = shape, rate = rate,
    log = TRUE
  )
  for(x in seq_len(n_grid)) {
    log_slab = slab_log_density(theta, grid[x], family)
    slab_sum = rowSums(gamma * log_slab)
    lw[, x] = log_hyper[x] + slab_sum + log(grid[x])
  }

  # Per-draw normalization over the grid (softmax across grid points).
  row_max = apply(lw, 1, max)
  p = exp(lw - row_max)
  p = p / rowSums(p)

  # Marginal weight mass and importance ESS per grid point.
  weight_mass = colSums(p)
  ess = weight_mass^2 / colSums(p^2)

  # Conditional inclusion probability: sum_t gamma p / sum_t p.
  num = crossprod(gamma, p) # n_edges x n_grid
  pip = sweep(num, 2, weight_mass, "/")
  pip = pmin(pmax(pip, 0), 1)

  post_odds = pip / (1 - pip)
  bf = sweep(post_odds, 1, prior_odds, "/")
  log10_bf = t(log10(bf))

  list(log10_bf = log10_bf, ess = ess)
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
# log10_to_bf
# ------------------------------------------------------------------
# Undo the log10 transform, mapping NA through.
# ------------------------------------------------------------------
log10_to_bf = function(x) 10^x


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
#' Prints the headline scale-robust verdict counts, the number of mover edges,
#' and the reliability window from a \code{\link{prior_sensitivity_check}()}
#' result.
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
  counts = table(factor(edges$marginalized_verdict,
    levels = c("presence", "undecided", "absence")
  ))
  n_movers = sum(edges$mover, na.rm = TRUE)
  n_saturated = sum(edges$saturated)
  reliable = sum(x$grid$ess >= x$ess_floor)

  cat("Prior sensitivity of edge inclusion verdicts\n\n")
  cat("Chosen interaction slab scale:", format(x$chosen_scale, digits = 3), "\n")
  cat(
    "Evidence thresholds: presence >=", x$evidence_threshold,
    " absence <=", format(1 / x$evidence_threshold, digits = 3), "\n\n"
  )

  cat("Marginalized verdicts (averaged over the scale posterior):\n")
  cat("  evidence of presence: ", counts[["presence"]], "\n")
  cat("  undecided:            ", counts[["undecided"]], "\n")
  cat("  evidence of absence:  ", counts[["absence"]], "\n\n")

  cat("Mover edges (verdict changes with the scale or differs from the\n")
  cat("chosen-scale verdict):", n_movers, "of", n_edges, "\n")
  if(n_saturated > 0) {
    cat(n_saturated, "edge(s) selected or excluded at every scale considered\n")
  }
  cat(
    "Reliable grid points (importance ESS >=", x$ess_floor, "):",
    reliable, "of", nrow(x$grid), "\n\n"
  )
  cat("Use plot() for the inclusion-Bayes-factor-versus-scale curves,\n")
  cat("and $edges for the per-edge table.\n")
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


#' @title Plot a Prior Sensitivity Check
#'
#' @description
#' Draws the inclusion-Bayes-factor-versus-scale curves over the relative slab
#' scale (scale / chosen scale) on a log axis. The top panel shows one
#' \eqn{\log_{10}} inclusion-Bayes-factor curve per edge, with the undecided
#' band shaded, the chosen scale marked, mover edges at full opacity and stable
#' edges faint, and thin per-chain curves underneath for the Monte Carlo
#' fragility check. The bottom panel is a forest of scale-stability intervals
#' for the mover edges.
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
  rel = x$grid$relative
  thr = log10(x$evidence_threshold)
  edges = x$edges
  n_edges = nrow(edges)
  movers = which(edges$mover & !edges$saturated)

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(2, 1), mar = c(4.2, 4.2, 2.2, 1))

  # --- Top panel: log10 BF curves ---------------------------------------------
  finite_bf = x$log10_bf[is.finite(x$log10_bf)]
  ylim = if(length(finite_bf)) {
    range(c(finite_bf, -thr, thr))
  } else {
    c(-thr, thr) * 2
  }
  graphics::plot(NA, NA,
    xlim = range(rel), ylim = ylim, log = "x",
    xlab = "relative slab scale (scale / chosen scale)",
    ylab = expression(log[10] * " inclusion Bayes factor"),
    main = "Inclusion Bayes factor vs. slab scale"
  )
  graphics::rect(min(rel), -thr, max(rel), thr,
    col = grDevices::adjustcolor(pal[["undecided"]], 0.12), border = NA
  )
  graphics::abline(
    h = c(-thr, 0, thr), col = "grey70",
    lty = c(2, 1, 2)
  )
  graphics::abline(v = 1, col = "grey40", lty = 3)

  # Thin per-chain curves for the mover edges (fragility check).
  for(e in movers) {
    col_e = grDevices::adjustcolor(pal[[edges$chosen_scale_verdict[e]]], 0.25)
    for(ch in x$log10_bf_by_chain) {
      graphics::lines(rel, ch[, e], col = col_e, lwd = 0.5)
    }
  }
  # Pooled curves: stable faint, movers solid.
  for(e in seq_len(n_edges)) {
    if(edges$saturated[e]) next
    is_mover = edges$mover[e]
    col_e = grDevices::adjustcolor(
      pal[[edges$chosen_scale_verdict[e]]], if(is_mover) 0.9 else 0.2
    )
    graphics::lines(rel, x$log10_bf[, e],
      col = col_e,
      lwd = if(is_mover) 2 else 1
    )
  }
  n_saturated = sum(edges$saturated)
  if(n_saturated > 0) {
    graphics::mtext(
      paste0(
        n_saturated,
        " edge(s): presence or absence at every scale considered"
      ),
      side = 3, line = 0.2, cex = 0.75, adj = 1, col = "grey40"
    )
  }

  # --- Bottom panel: forest of stability intervals for movers -----------------
  if(length(movers) == 0L) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No mover edges: every verdict is scale-robust.",
      cex = 1
    )
    return(invisible(x))
  }
  ord = movers[order(edges$stability_upper[movers] -
    edges$stability_lower[movers])]
  show = utils::head(ord, max_labels)
  yy = seq_along(show)
  xr = range(c(edges$stability_lower[show], edges$stability_upper[show], 1),
    na.rm = TRUE
  )
  graphics::plot(NA, NA,
    xlim = xr, ylim = c(0.5, length(show) + 0.5),
    log = "x", yaxt = "n", xlab = "relative slab scale",
    ylab = "", main = "Scale-stability intervals (mover edges)"
  )
  graphics::abline(v = 1, col = "grey40", lty = 3)
  for(k in seq_along(show)) {
    e = show[k]
    graphics::segments(edges$stability_lower[e], yy[k],
      edges$stability_upper[e], yy[k],
      col = pal[[edges$chosen_scale_verdict[e]]], lwd = 3
    )
    graphics::points(1, yy[k],
      pch = 18,
      col = pal[[edges$chosen_scale_verdict[e]]]
    )
  }
  graphics::axis(2, at = yy, labels = edges$edge[show], las = 1, cex.axis = 0.7)
  invisible(x)
}
