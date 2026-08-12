# ==============================================================================
# Posterior centrality
# ==============================================================================
#
# A centrality is a function of the network, so evaluating it on every posterior
# draw of the network gives its posterior distribution directly. Nothing here is
# an estimate plus an error bar bolted on afterwards: the draws come from the
# same model-averaged posterior as the edge weights.
# ==============================================================================


#' @title Extract Posterior Centrality
#'
#' @description
#' Evaluates a node centrality on every posterior draw of the network, giving
#' the posterior distribution of each node's centrality.
#'
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()]) or
#'   `bgmCompare` (from [bgmCompare()]).
#' @param measure Character; the centrality to evaluate. Currently `"strength"`
#'   (default), the sum of the absolute weights of a node's edges.
#' @param group For a `bgmCompare` fit: a single group index, giving that
#'   group's centrality, or two indices, giving the difference in centrality
#'   between them (the first minus the second). Default `1`. Ignored for a
#'   `bgm` fit, which has one network.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix of class `bgms_centrality` with one row per
#'   posterior draw and one column per variable, carrying the `measure` as an
#'   attribute. Use [summary()] for posterior means, credible intervals, and
#'   the probability of being the most central node, and [plot()] for the
#'   ordered interval display.
#'
#' @details
#' The centrality is evaluated on the model-averaged pairwise draws, in which an
#' edge excluded at a given iteration contributes exactly zero. Structural
#' uncertainty therefore propagates into the centrality without any extra step:
#' a node whose edges are themselves uncertain gets a wide centrality
#' posterior, which is the honest summary. Conditioning on the included-only
#' draws instead would report each node's centrality in the models where its
#' edges happen to be present, which is a different and generally larger
#' quantity.
#'
#' Strength centrality sums the absolute edge weights, so positive and negative
#' associations both add to a node's total involvement rather than cancelling.
#'
#' For a [bgmCompare()] fit each group's network is rebuilt on every draw as
#' `baseline + (P %*% differences)`, with the fit's own contrast projection `P`,
#' rather than from posterior means. That is what carries the uncertainty
#' through: with `group = c(1, 2)` the credible interval is the interval of the
#' *difference* in a node's centrality, which answers whether the groups differ
#' in it directly, where two separately drawn intervals do not. A draw in which
#' a difference indicator is zero gives both groups the same edge weight and so
#' contributes exactly zero to the difference, which is why this is computed per
#' draw and not from summaries.
#'
#' Read a centrality difference with care, and as numbers rather than as a
#' picture ([summary()]; there is no plot method for it): strength sums
#' absolute weights, so a difference of zero can mean identical networks or
#' compensating edge differences, and its sign says nothing about which edges
#' moved. [verdicts()] is the per-difference evidence.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' strength = extract_centrality(fit)
#' summary(strength)
#' }
#'
#' @seealso [extract_pairwise_interactions()] for the draws this is computed
#'   from, [verdicts()] for which of those edges the data settle
#' @family extractors
#' @export
extract_centrality = function(bgms_object, measure = "strength", group = 1, ...) {
  UseMethod("extract_centrality")
}


# ------------------------------------------------------------------
# strength_from_pairwise
# ------------------------------------------------------------------
# Strength centrality of every node, on each row of a draws x pairs matrix of
# edge weights.
#
# The incidence comes from the edge order rather than from the column names: a
# variable name containing a hyphen would defeat splitting "A-B" back into its
# nodes.
#
# @param samples        draws x pairs matrix of edge weights.
# @param num_variables  Number of nodes.
#
# Returns: a draws x nodes matrix.
# ------------------------------------------------------------------
strength_from_pairwise = function(samples, num_variables) {
  pairs = which(upper.tri(matrix(0, num_variables, num_variables)), arr.ind = TRUE)
  pairs = pairs[order(pairs[, "row"], pairs[, "col"]), , drop = FALSE]
  stopifnot(nrow(pairs) == ncol(samples))

  absolute = abs(samples)
  vapply(
    seq_len(num_variables),
    function(v) {
      rowSums(absolute[, pairs[, 1] == v | pairs[, 2] == v, drop = FALSE])
    },
    numeric(nrow(samples))
  )
}

#' @inheritParams extract_centrality
#' @exportS3Method
#' @noRd
extract_centrality.bgms = function(bgms_object, measure = "strength", group = 1, ...) {
  measure = match.arg(measure, choices = "strength")

  samples = extract_pairwise_interactions(bgms_object)
  nodes = extract_arguments(bgms_object)$data_columnnames

  structure(
    strength_from_pairwise(samples, length(nodes)),
    dimnames = list(NULL, nodes),
    class = c("bgms_centrality", "matrix", "array"),
    measure = measure
  )
}

#' @inheritParams extract_centrality
#' @exportS3Method
#' @noRd
extract_centrality.bgmCompare = function(bgms_object, measure = "strength",
                                         group = 1, ...) {
  measure = match.arg(measure, choices = "strength")

  arguments = extract_arguments(bgms_object)
  nodes = arguments$data_columnnames
  num_variables = as.integer(arguments$num_variables)
  num_groups = as.integer(arguments$num_groups)
  projection = arguments$projection

  if(!is.numeric(group) || length(group) < 1L || length(group) > 2L ||
    anyNA(group) || any(group != as.integer(group)) ||
    any(group < 1L) || any(group > num_groups)) {
    stop(
      "Argument 'group' must be one group index, for that group's centrality, ",
      "or two, for the difference between them. This fit has ", num_groups,
      " groups."
    )
  }
  group = as.integer(group)
  if(length(group) == 2L && group[1] == group[2]) {
    stop("The two entries of 'group' must be different groups.")
  }

  # Per-draw group networks, not posterior means: the pairwise draws hold the
  # baseline block followed by one block per contrast, and a group's weights are
  # baseline + (P %*% differences) with the fit's own projection.
  draws = do.call(rbind, get_raw_samples(bgms_object)$pairwise)
  num_pairs = ncol(draws) / num_groups
  baseline = draws[, seq_len(num_pairs), drop = FALSE]
  differences = lapply(seq_len(num_groups - 1L), function(k) {
    draws[, k * num_pairs + seq_len(num_pairs), drop = FALSE]
  })

  weights_of = function(g) {
    out = baseline
    for(k in seq_along(differences)) {
      out = out + projection[g, k] * differences[[k]]
    }
    out
  }

  labels = compare_group_labels(arguments, num_groups)
  if(length(group) == 1L) {
    centrality = strength_from_pairwise(weights_of(group), num_variables)
    label = sprintf("%s %s centrality", group_tag(labels, group), measure)
  } else {
    centrality = strength_from_pairwise(weights_of(group[1]), num_variables) -
      strength_from_pairwise(weights_of(group[2]), num_variables)
    label = sprintf(
      "difference in %s centrality (%s - %s)",
      measure, group_tag(labels, group[1]), group_tag(labels, group[2])
    )
  }

  structure(
    centrality,
    dimnames = list(NULL, nodes),
    class = c("bgms_centrality", "matrix", "array"),
    measure = measure,
    label = label,
    group = group
  )
}


# ------------------------------------------------------------------
# centrality_label
# ------------------------------------------------------------------
# What the axis and the summary column call the quantity. bgm() has one
# network and one label; bgmCompare() names the group or the contrast.
# ------------------------------------------------------------------
centrality_label = function(x) {
  attr(x, "label") %||% sprintf("%s centrality", attr(x, "measure"))
}


# ------------------------------------------------------------------
# is_centrality_difference
# ------------------------------------------------------------------
# Whether the object holds a difference between two groups' centralities,
# in which case zero is the reference rather than the smallest value.
# ------------------------------------------------------------------
is_centrality_difference = function(x) {
  length(attr(x, "group")) == 2L
}


#' @title Summarize Posterior Centrality
#'
#' @description
#' Posterior mean, credible interval, and probability of being the most central
#' node, per variable.
#'
#' @param object An object of class `bgms_centrality`, from
#'   [extract_centrality()].
#' @param probs Numeric of length two; the credible-interval quantiles. Default
#'   `c(0.025, 0.975)`.
#' @param ... Ignored.
#'
#' @return A data frame with one row per variable, ordered by decreasing
#'   posterior mean, with columns `node`, `mean`, `lower`, `upper`, and
#'   `p_most_central`, the posterior probability that the node has the largest
#'   centrality of all nodes.
#'
#'   For a difference between two groups the last column is `p_positive`
#'   instead: the posterior probability that the node's centrality is higher in
#'   the first group than in the second. Which node is *most* central is not the
#'   question a difference answers.
#'
#' @details
#' A difference has a point mass at exactly zero, from the draws in which every
#' one of the node's difference indicators is excluded and the two groups share
#' the network. `p_positive` and its mirror therefore need not sum to one, and a
#' node can have a positive posterior mean with `p_positive` well below `0.5`:
#' the remaining mass is on no difference at all, which is the model averaging
#' reporting itself.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' summary(extract_centrality(fit))
#' }
#'
#' @seealso [extract_centrality()]
#' @family extractors
#' @export
summary.bgms_centrality = function(object, probs = c(0.025, 0.975), ...) {
  nodes = colnames(object)
  interval = apply(object, 2, stats::quantile, probs = probs)

  out = data.frame(
    node = nodes,
    mean = colMeans(object),
    lower = interval[1, ],
    upper = interval[2, ],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if(is_centrality_difference(object)) {
    # A difference has no "most central" node to win; the direct read is which
    # side of zero the node's difference falls on.
    out$p_positive = colMeans(object > 0)
  } else {
    # Which node is most central is itself uncertain, and the share of draws a
    # node wins is the direct read of that uncertainty. A near-tie between the
    # top nodes shows here as two middling probabilities, where a ranking of
    # posterior means would show a winner.
    top = factor(nodes[max.col(object, ties.method = "first")], levels = nodes)
    out$p_most_central = as.numeric(table(top)) / nrow(object)
  }
  out[order(out$mean, decreasing = TRUE), , drop = FALSE]
}


#' @title Plot Posterior Centrality
#'
#' @description
#' Draws each node's posterior mean centrality with its credible interval,
#' nodes ordered by mean.
#'
#' @param x An object of class `bgms_centrality`, from [extract_centrality()].
#' @param probs Numeric of length two; the credible-interval quantiles. Default
#'   `c(0.025, 0.975)`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly. Called for the side effect of drawing.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' plot(extract_centrality(fit))
#' }
#'
#' @seealso [extract_centrality()]
#' @family extractors
#' @export
plot.bgms_centrality = function(x, probs = c(0.025, 0.975), ...) {
  if(is_centrality_difference(x)) {
    stop(
      "A difference-centrality plot is not offered: strength centrality sums ",
      "the absolute weights of a node's edges, so a between-group difference ",
      "in it conflates which edges differ with how their signs cancel and has ",
      "no edge-level reading. summary() reports the difference as numbers."
    )
  }
  summ = summary(x, probs = probs)
  # Largest mean at the top of the panel.
  summ = summ[order(summ$mean), , drop = FALSE]
  n = nrow(summ)

  style = bgms_style()
  # The node names live in the left margin, so the margin is sized for the
  # names this fit actually has rather than for a width that happened to fit
  # the author's example.
  left = margin_lines_for(summ$node, cex = style$cex_axis, pad = 1.4)
  style = bgms_panel_par(mar = c(4.6, left, 4.0, 2.2))
  on.exit(graphics::par(style$old_par), add = TRUE)

  # Difference objects are rejected above, so the scale is always a plain
  # centrality scale: no zero anchoring and no zero reference line.
  x_axis = bgms_axis_range(c(summ$lower, summ$upper))

  graphics::plot(NA, NA,
    xlim = x_axis$lim, ylim = c(0.5, n + 0.5),
    axes = FALSE, xlab = "", ylab = "", main = ""
  )
  bgms_axis(1, x_axis$at, style = style)
  graphics::axis(2,
    at = seq_len(n), labels = summ$node, las = 1, tick = FALSE,
    line = -0.5, cex.axis = style$cex_axis, col.axis = style$ink
  )
  graphics::segments(summ$lower, seq_len(n), summ$upper, seq_len(n),
    col = grDevices::adjustcolor(style$accent, 0.45), lwd = 4.5, lend = 1
  )
  graphics::points(summ$mean, seq_len(n),
    pch = 16, col = style$accent, cex = 1.3
  )

  # The axis carries the quantity. What the dot and the bar are is in the Rd:
  # a figure does not need a sentence under it to be read.
  graphics::mtext(centrality_label(x),
    side = 1, line = 3.0, cex = style$cex_lab, col = style$ink
  )

  invisible(x)
}
