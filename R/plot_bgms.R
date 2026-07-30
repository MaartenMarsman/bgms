# ==============================================================================
# Plot methods for fitted models
# ==============================================================================
#
# The default picture is the reporting object, not an estimate: edges are drawn
# by what the data settle about them, so the network a reader sees and the
# verdict table they cite are the same statement. Network layout comes from
# qgraph, a suggested package; the wider set of network displays lives in
# easybgm, which is what JASP wraps.
# ==============================================================================


# ------------------------------------------------------------------
# verdict_edge_colors
# ------------------------------------------------------------------
# Edge colours for the verdict-encoded network: sign-carrying accents for edges
# the data place in the model, one recessive grey for the undecided. The two
# accents are the Okabe-Ito blue/vermillion pair, so the sign survives common
# forms of colour-vision deficiency; the verdict itself is carried by line type
# as well, never by colour alone.
#
# @param weight   Model-averaged edge weights.
# @param verdict  Verdict per edge.
#
# Returns: character vector of colours.
# ------------------------------------------------------------------
verdict_edge_colors = function(weight, verdict) {
  palette = mover_palette()
  out = rep("grey65", length(weight))
  present = verdict == "presence"
  out[present & weight >= 0] = palette[1]
  out[present & weight < 0] = palette[2]
  out
}


#' @title Plot a Fitted bgms Model
#'
#' @description
#' Draws the model-averaged network with edges encoded by what the data settle
#' about them, or one of the other standard displays.
#'
#' @param x A fitted model object of class `bgms`, from [bgm()].
#' @param type Character; which display to draw. `"network"` (default) is the
#'   verdict-encoded model-averaged network; `"centrality"` is the posterior
#'   strength centrality of [extract_centrality()].
#' @param evidence_threshold Numeric > 1; the inclusion Bayes factor separating
#'   evidence of presence from undecided, as in [verdicts()]. Default `10`.
#' @param layout Layout passed to [qgraph::qgraph()]. Default `"spring"`.
#' @param legend Logical; draw the edge legend. Default `TRUE`.
#' @param ... Passed to [qgraph::qgraph()] for `type = "network"`, and to
#'   [plot.bgms_centrality()] otherwise.
#'
#' @return `x`, invisibly. Called for the side effect of drawing.
#'
#' @details
#' Edges are drawn by verdict rather than by estimate. An edge with evidence of
#' presence is drawn solid, with its width scaled by the model-averaged weight
#' and its colour carrying the sign; an undecided edge is drawn as a thin
#' dotted grey line, because the data neither place it in the network nor rule
#' it out; an edge with evidence of absence is not drawn. The default picture
#' therefore says the same thing as [verdicts()] at the same threshold, instead
#' of showing one unqualified network of point estimates.
#'
#' Drawing the network requires the suggested package qgraph. The three-panel
#' evidence display, the structure plots, and the other network displays live
#' in the easybgm package, which builds on these fits.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' plot(fit)
#' plot(fit, type = "centrality")
#' }
#'
#' @seealso [verdicts()] for the table the picture encodes,
#'   [extract_centrality()], [plot_edge_posterior()] for one edge in detail
#' @family posterior-methods
#' @export
plot.bgms = function(x,
                     type = c("network", "centrality"),
                     evidence_threshold = 10,
                     layout = "spring",
                     legend = TRUE,
                     ...) {
  type = match.arg(type)
  if(type == "centrality") {
    plot(extract_centrality(x), ...)
    return(invisible(x))
  }

  if(!requireNamespace("qgraph", quietly = TRUE)) {
    stop(
      "Drawing the network needs the qgraph package, which is suggested rather ",
      "than required by bgms. Install it with install.packages(\"qgraph\"), or ",
      "use verdicts() for the same information as a table."
    )
  }

  edges = verdicts(x, evidence_threshold = evidence_threshold)
  weight = colMeans(extract_pairwise_interactions(x))
  variables = extract_arguments(x)$data_columnnames
  num_variables = length(variables)

  pairs = which(upper.tri(matrix(0, num_variables, num_variables)), arr.ind = TRUE)
  pairs = pairs[order(pairs[, "row"], pairs[, "col"]), , drop = FALSE]

  drawn = !is.na(edges$verdict) & edges$verdict != "absence"
  if(!any(drawn)) {
    stop(
      "No edge reaches evidence of presence or sits undecided at this ",
      "threshold, so there is nothing to draw. verdicts(fit) reports the ",
      "evidence for every edge."
    )
  }

  verdict = as.character(edges$verdict[drawn])
  arguments = list(
    input = cbind(pairs[drawn, , drop = FALSE], weight[drawn]),
    nNodes = num_variables,
    labels = variables,
    directed = FALSE,
    layout = layout,
    edge.color = verdict_edge_colors(weight[drawn], verdict),
    lty = ifelse(verdict == "presence", 1L, 3L),
    # qgraph fades an edge toward the background in proportion to its weight,
    # which would wash the verdict encoding out: a settled but weak edge would
    # come out as faint as an undecided one. Width already carries the weight.
    fade = FALSE,
    minimum = 0
  )
  user = list(...)
  do.call(qgraph::qgraph, utils::modifyList(arguments, user))

  if(isTRUE(legend)) {
    graphics::legend("bottomleft",
      legend = c("present, positive", "present, negative", "undecided"),
      col = c(mover_palette()[1:2], "grey65"),
      lty = c(1L, 1L, 3L), lwd = c(2.4, 2.4, 1.2),
      bty = "n", cex = 0.75, text.col = "grey25", xpd = NA
    )
  }
  invisible(x)
}


#' @title Plot the Posterior of One Edge Weight
#'
#' @description
#' Draws the spike-and-slab posterior of a single edge weight: a stem at zero
#' whose height is the posterior probability that the edge is absent, and the
#' slab of the weight given that it is present.
#'
#' @param bgms_object A fitted model object of class `bgms`, from [bgm()] with
#'   `edge_selection = TRUE`.
#' @param variable1,variable2 The two variables naming the edge. Either names
#'   or column positions.
#' @param evidence_threshold Numeric > 1; the threshold used for the verdict in
#'   the panel's title, as in [verdicts()]. Default `10`.
#' @param binwidth Width of the weight bin the slab's height is expressed per,
#'   so that the stem and the slab are on one scale. Default `0.01`.
#' @param ... Ignored.
#'
#' @return `bgms_object`, invisibly. Called for the side effect of drawing.
#'
#' @details
#' The two parts are commensurable, which is the point of the picture: the stem
#' is the probability of exactly zero, and the slab is drawn as probability per
#' weight bin, so the visible areas divide the posterior mass as the model
#' does. A decisive edge is a short stem under a well-separated slab; an
#' undecided edge splits its mass between them.
#'
#' The y axis carries no numbers. Its unit is a display choice, and the two
#' quantities a reader needs -- the probability of absence and the location of
#' the weight -- are printed and on the x axis respectively.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' plot_edge_posterior(fit, "intrusion", "dreams")
#' }
#'
#' @seealso [verdicts()], [extract_pairwise_interactions()]
#' @family posterior-methods
#' @export
plot_edge_posterior = function(bgms_object, variable1, variable2,
                               evidence_threshold = 10, binwidth = 0.01, ...) {
  check_evidence_threshold(evidence_threshold)
  samples = extract_pairwise_interactions(bgms_object)
  variables = extract_arguments(bgms_object)$data_columnnames

  first = resolve_variable(variable1, variables, "variable1")
  second = resolve_variable(variable2, variables, "variable2")
  if(first == second) {
    stop("Arguments 'variable1' and 'variable2' must name two different variables.")
  }
  if(first > second) {
    swap = first
    first = second
    second = swap
  }
  label = paste(variables[first], variables[second], sep = "-")
  if(!label %in% colnames(samples)) {
    stop("No edge between '", variables[first], "' and '", variables[second], "'.")
  }

  draws = samples[, label]
  p_spike = mean(draws == 0)
  slab = draws[draws != 0]
  if(length(slab) < 2L) {
    stop(
      "The edge '", label, "' was included in fewer than two retained draws, ",
      "so there is no slab to draw. verdicts(fit) reports its evidence."
    )
  }

  edges = verdicts(bgms_object, evidence_threshold = evidence_threshold)
  row = edges[edges$parameter == label, , drop = FALSE]
  bayes_factor = row$bf
  title = sprintf(
    "%s\n%s, BF = %s", label, as.character(row$verdict),
    if(is.finite(bayes_factor) && bayes_factor >= 1) {
      sprintf("%.1f", bayes_factor)
    } else if(is.finite(bayes_factor)) {
      sprintf("%.3f", bayes_factor)
    } else {
      format(bayes_factor)
    }
  )

  ink = "grey25"
  muted = "grey55"
  accent = mover_palette()[1]

  span = max(abs(slab))
  density = stats::density(slab, from = -1.1 * span, to = 1.1 * span, n = 512)
  height = density$y * (1 - p_spike) * binwidth

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mar = c(5.4, 1.6, 3.6, 1.6), mgp = c(2.3, 0.6, 0),
    col.axis = ink, col.lab = ink, col.main = ink
  )

  graphics::plot(NA, NA,
    xlim = range(c(density$x, 0)), ylim = c(0, 1.18 * max(c(height, p_spike))),
    axes = FALSE, xlab = "Edge weight", ylab = "", main = title, cex.main = 1
  )
  graphics::axis(1, col = muted, col.ticks = muted)

  # Fill without a border, then trace the top only: a bordered polygon would
  # close along y = 0 and draw a baseline across the whole weight range.
  graphics::polygon(
    c(density$x, rev(density$x)), c(height, rep(0, length(height))),
    col = grDevices::adjustcolor(accent, 0.25), border = NA
  )
  graphics::lines(density$x, height, col = accent, lwd = 1.6)
  graphics::segments(0, 0, 0, p_spike, col = ink, lwd = 2.6, lend = 1)
  graphics::points(0, p_spike, pch = 16, col = ink, cex = 1.2)
  graphics::text(0, p_spike,
    labels = if(p_spike < 0.01) "<.01" else sub("^0", "", sprintf("%.2f", p_spike)),
    pos = 3, offset = 0.5, cex = 0.85, col = ink
  )
  graphics::mtext("P(edge absent) at the stem; slab on the same scale",
    side = 1, line = 3.8, cex = 0.75, col = muted
  )

  invisible(bgms_object)
}


# ------------------------------------------------------------------
# resolve_variable
# ------------------------------------------------------------------
# Turn a variable name or column position into a column index.
#
# @param value      A name or a position.
# @param variables  The fit's variable names.
# @param argument   Argument name, for the error message.
#
# Returns: an integer column index.
# ------------------------------------------------------------------
resolve_variable = function(value, variables, argument) {
  if(is.character(value)) {
    index = match(value, variables)
    if(is.na(index)) {
      stop(
        "Argument '", argument, "' is '", value, "', which is not one of the ",
        "model's variables: ", paste(variables, collapse = ", "), "."
      )
    }
    return(index)
  }
  if(is.numeric(value) && length(value) == 1L &&
    value >= 1 && value <= length(variables) && value == as.integer(value)) {
    return(as.integer(value))
  }
  stop(
    "Argument '", argument, "' must be a variable name or a column position ",
    "between 1 and ", length(variables), "."
  )
}
