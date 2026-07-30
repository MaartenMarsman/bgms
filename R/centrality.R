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
#' @param bgms_object A fitted model object of class `bgms` (from [bgm()]).
#' @param measure Character; the centrality to evaluate. Currently `"strength"`
#'   (default), the sum of the absolute weights of a node's edges.
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
extract_centrality = function(bgms_object, measure = "strength", ...) {
  UseMethod("extract_centrality")
}

#' @inheritParams extract_centrality
#' @exportS3Method
#' @noRd
extract_centrality.bgms = function(bgms_object, measure = "strength", ...) {
  measure = match.arg(measure, choices = "strength")

  samples = extract_pairwise_interactions(bgms_object)
  nodes = extract_arguments(bgms_object)$data_columnnames
  num_variables = length(nodes)

  # Incidence from the edge order rather than from the column names: a variable
  # name containing a hyphen would defeat splitting "A-B" back into its nodes.
  pairs = which(upper.tri(matrix(0, num_variables, num_variables)), arr.ind = TRUE)
  pairs = pairs[order(pairs[, "row"], pairs[, "col"]), , drop = FALSE]
  stopifnot(nrow(pairs) == ncol(samples))

  absolute = abs(samples)
  strength = vapply(
    seq_len(num_variables),
    function(v) {
      rowSums(absolute[, pairs[, 1] == v | pairs[, 2] == v, drop = FALSE])
    },
    numeric(nrow(samples))
  )

  structure(
    strength,
    dimnames = list(NULL, nodes),
    class = c("bgms_centrality", "matrix", "array"),
    measure = measure
  )
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

  # Which node is most central is itself uncertain, and the share of draws a
  # node wins is the direct read of that uncertainty. A near-tie between the
  # top nodes shows here as two middling probabilities, where a ranking of
  # posterior means would show a winner.
  top = factor(nodes[max.col(object, ties.method = "first")], levels = nodes)

  out = data.frame(
    node = nodes,
    mean = colMeans(object),
    lower = interval[1, ],
    upper = interval[2, ],
    p_most_central = as.numeric(table(top)) / nrow(object),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
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
  summ = summary(x, probs = probs)
  # Largest mean at the top of the panel.
  summ = summ[order(summ$mean), , drop = FALSE]
  n = nrow(summ)

  ink = "grey25"
  muted = "grey55"
  accent = mover_palette()[1]

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  label_width = max(graphics::strwidth(summ$node, units = "inches", cex = 0.85))
  graphics::par(
    mar = c(4.1, 1.2 + 6 * label_width, 2.6, 1.6), mgp = c(2.4, 0.6, 0),
    col.axis = ink, col.lab = ink, col.main = ink
  )

  graphics::plot(NA, NA,
    xlim = range(c(summ$lower, summ$upper)), ylim = c(0.5, n + 0.5),
    axes = FALSE, xlab = sprintf(
      "%s centrality (posterior mean and %g%% credible interval)",
      attr(x, "measure"), 100 * diff(probs)
    ),
    ylab = "", main = ""
  )
  graphics::axis(1, col = muted, col.ticks = muted)
  graphics::axis(2,
    at = seq_len(n), labels = summ$node, las = 1, tick = FALSE,
    line = -0.5, cex.axis = 0.85
  )
  graphics::segments(summ$lower, seq_len(n), summ$upper, seq_len(n),
    col = grDevices::adjustcolor(accent, 0.45), lwd = 4, lend = 1
  )
  graphics::points(summ$mean, seq_len(n), pch = 16, col = accent, cex = 1.1)

  invisible(x)
}
