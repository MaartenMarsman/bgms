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


# ------------------------------------------------------------------
# verdict_network_input
# ------------------------------------------------------------------
# The weight matrix qgraph is given, with the per-edge colours and line types
# in the order it reads that matrix.
#
# A weighted edgelist would be the direct expression of "draw these edges", but
# qgraph derives the node set from the edgelist and its nNodes argument does not
# reliably override that, so a network that leaves a node out fails rather than
# drawing the node alone. A square matrix states the node set in its dimensions.
# qgraph then reads its non-zero upper triangle in column-major order, which is
# the order the colour and line-type vectors are put in here.
#
# @param weight         Weight per pair, in row-major upper-triangle order.
# @param verdict        Verdict per pair, in the same order.
# @param pairs          The row-major upper-triangle index.
# @param num_variables  Number of nodes.
#
# Returns: list(weights = matrix, edge.color, lty), or NULL when nothing is
# drawable, in which case the caller draws the nodes alone.
# ------------------------------------------------------------------
verdict_network_input = function(weight, verdict, pairs, num_variables) {
  drawn = !is.na(verdict) & verdict != "absence"
  weights = matrix(0, num_variables, num_variables)
  if(!any(drawn)) {
    return(list(weights = weights, edge.color = character(0), lty = integer(0)))
  }
  # An exactly zero weight would be dropped as a non-edge, which would silently
  # shift every later colour onto the wrong edge. Nudge it instead.
  value = weight[drawn]
  value[value == 0] = .Machine$double.eps
  index = pairs[drawn, , drop = FALSE]
  weights[index] = value
  weights[index[, 2:1, drop = FALSE]] = value

  order_read = which(upper.tri(weights) & weights != 0, arr.ind = TRUE)
  position = match(
    paste(order_read[, 1], order_read[, 2]),
    paste(index[, 1], index[, 2])
  )
  list(
    weights = weights,
    edge.color = verdict_edge_colors(value, verdict[drawn])[position],
    lty = ifelse(verdict[drawn][position] == "presence", 1L, 3L)
  )
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

  network = verdict_network_input(
    weight, as.character(edges$verdict), pairs, num_variables
  )
  arguments = list(
    input = network$weights,
    labels = variables,
    directed = FALSE,
    layout = layout,
    edge.color = network$edge.color,
    lty = network$lty,
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


# ------------------------------------------------------------------
# compare_difference_verdicts
# ------------------------------------------------------------------
# Split a bgmCompare verdict table into its two indicator families, in the
# orders the drawing needs them: pairwise differences in the row-major upper
# triangle, main-effect differences one per variable.
#
# @param x                   A fitted bgmCompare object.
# @param evidence_threshold  Threshold passed to verdicts().
#
# Returns: list(pairwise = verdict per pair, main = verdict per variable,
#   pairs = the row-major upper-triangle index, main_selected = whether the
#   main-effect difference indicators were ever updated).
# ------------------------------------------------------------------
compare_difference_verdicts = function(x, evidence_threshold) {
  table = verdicts(x, evidence_threshold = evidence_threshold)
  num_variables = as.integer(extract_arguments(x)$num_variables)
  index = compare_indicator_index(num_variables)
  is_main = index[, 1] == index[, 2]

  pairs = which(upper.tri(matrix(0, num_variables, num_variables)), arr.ind = TRUE)
  pairs = pairs[order(pairs[, "row"], pairs[, "col"]), , drop = FALSE]

  main = as.character(table$verdict[is_main])
  list(
    pairwise = as.character(table$verdict[!is_main]),
    main = main,
    pairs = pairs,
    # Without main_difference_selection those indicators are never updated and
    # carry no verdict, so the node channel has nothing to say.
    main_selected = !all(is.na(main))
  )
}


# ------------------------------------------------------------------
# main_difference_nodes
# ------------------------------------------------------------------
# Node shape and border colour encoding the main-effect difference verdicts.
# Shape carries the settled/unsettled distinction, so the encoding does not
# ride on colour alone; the border colour grades it.
#
# @param verdict        Verdict per variable, or all NA when never updated.
# @param main_selected  Whether the indicators were updated at all.
#
# Returns: list(shape, border_color).
# ------------------------------------------------------------------
main_difference_nodes = function(verdict, main_selected) {
  n = length(verdict)
  if(!main_selected) {
    return(list(shape = rep("circle", n), border_color = rep("grey55", n)))
  }
  shape = ifelse(verdict %in% "presence", "square", "circle")
  border = rep("grey80", n)
  border[verdict %in% "undecided"] = "grey55"
  border[verdict %in% "presence"] = mover_palette()[1]
  list(shape = shape, border_color = border)
}


#' @title Plot a Fitted bgmCompare Model
#'
#' @description
#' Draws the group differences the data settle, as one network, or the groups'
#' own networks beside it on a shared layout.
#'
#' @param x A fitted model object of class `bgmCompare`, from [bgmCompare()].
#' @param type Character; which display to draw. `"difference"` (default) is
#'   the verdict-encoded network of group differences; `"groups"` draws each
#'   group's own network and the difference panel on one shared layout;
#'   `"centrality"` is the posterior strength centrality of
#'   [extract_centrality()].
#' @param evidence_threshold Numeric > 1; the inclusion Bayes factor separating
#'   evidence of a difference from undecided, as in [verdicts()]. Default `10`.
#' @param group For `type = "centrality"`: passed to [extract_centrality()], so
#'   a single index gives that group's centrality and two give the difference.
#'   Default `1`.
#' @param layout Layout passed to [qgraph::qgraph()]. Default `"spring"`.
#' @param legend Logical; draw the legend. Default `TRUE`.
#' @param ... Passed to [qgraph::qgraph()], or to [plot.bgms_centrality()] for
#'   `type = "centrality"`.
#'
#' @return `x`, invisibly. Called for the side effect of drawing.
#'
#' @details
#' The default picture is about differences, because differences are what
#' [bgmCompare()] parameterizes. A pairwise difference with evidence of presence
#' is drawn solid, its width scaled by the size of the difference and its colour
#' carrying the sign; an undecided difference is a thin dotted grey line; a
#' difference the data rule out is not drawn. The picture therefore says the
#' same thing as `verdicts(fit)` at the same threshold.
#'
#' A network with no edges left is a result, not an error: "the groups do not
#' differ anywhere" is a common and correct finding, so the nodes are drawn on
#' their own and the subtitle says so.
#'
#' Main-effect differences are not edges. When `main_difference_selection =
#' TRUE` gave them their own indicators, their verdicts are carried on the
#' nodes: a square node with an accented border where the data settle a
#' main-effect difference, a circle with a grey border where they leave it
#' undecided, and a faint border where they rule it out. Shape carries the
#' settled/unsettled distinction so the encoding does not rest on colour alone.
#' Under the default `main_difference_selection = FALSE` those indicators are
#' never updated and have no verdict, so every node is drawn alike and the
#' subtitle says the channel is empty. `verdicts()` remains the place to read
#' main-effect differences precisely; the nodes are a summary of it.
#'
#' Drawing needs the suggested package qgraph.
#'
#' @examples
#' \donttest{
#' fit = bgmCompare(
#'   x = Wenchuan[, 1:5],
#'   group_indicator = rep(1:2, length.out = nrow(Wenchuan)),
#'   display_progress = "none"
#' )
#' plot(fit)
#' plot(fit, type = "groups")
#' }
#'
#' @seealso [verdicts()] for the table the picture encodes,
#'   [extract_centrality()], [prior_sensitivity_check()]
#' @family posterior-methods
#' @export
plot.bgmCompare = function(x,
                           type = c("difference", "groups", "centrality"),
                           evidence_threshold = 10,
                           group = 1,
                           layout = "spring",
                           legend = TRUE,
                           ...) {
  type = match.arg(type)
  if(type == "centrality") {
    plot(extract_centrality(x, group = group), ...)
    return(invisible(x))
  }

  if(!requireNamespace("qgraph", quietly = TRUE)) {
    stop(
      "Drawing the network needs the qgraph package, which is suggested rather ",
      "than required by bgms. Install it with install.packages(\"qgraph\"), or ",
      "use verdicts() for the same information as a table."
    )
  }

  arguments = extract_arguments(x)
  variables = arguments$data_columnnames
  num_variables = length(variables)
  num_groups = as.integer(arguments$num_groups)

  found = compare_difference_verdicts(x, evidence_threshold)
  differences = x@posterior_mean_pairwise_differences
  if(is.list(differences)) {
    stop(
      "Drawing group differences is implemented for two groups; this fit has ",
      num_groups, " groups and so ", num_groups - 1L, " contrasts, which are ",
      "not one network. Use verdicts() for the evidence and ",
      "extract_group_params() for each group's parameters."
    )
  }
  weight = differences[cbind(found$pairs[, 1], found$pairs[, 2])]

  if(type == "groups") {
    compare_group_panels(x, found, weight, variables, num_groups, layout, legend, ...)
    return(invisible(x))
  }

  nodes = main_difference_nodes(found$main, found$main_selected)
  draw_difference_network(
    weight, found$pairwise, found$pairs, variables, num_variables,
    nodes, layout, legend, found$main_selected, ...
  )
  invisible(x)
}


# ------------------------------------------------------------------
# draw_difference_network
# ------------------------------------------------------------------
# The difference panel: edges by difference verdict, nodes by main-effect
# difference verdict. Returns the layout qgraph used, so the group panels can
# share it.
# ------------------------------------------------------------------
draw_difference_network = function(weight, verdict, pairs, variables,
                                   num_variables, nodes, layout, legend,
                                   main_selected, ..., title = NULL) {
  drawn = !is.na(verdict) & verdict != "absence"
  # An empty difference network is the finding "the groups do not differ
  # anywhere", and a matrix of zeros draws the nodes alone.
  network = verdict_network_input(weight, verdict, pairs, num_variables)

  arguments = list(
    input = network$weights,
    labels = variables,
    directed = FALSE,
    layout = layout,
    edge.color = network$edge.color,
    lty = network$lty,
    shape = nodes$shape,
    border.color = nodes$border_color,
    fade = FALSE,
    minimum = 0,
    title = title,
    DoNotPlot = FALSE
  )
  result = do.call(qgraph::qgraph, utils::modifyList(arguments, list(...)))

  subtitle = character(0)
  if(!any(drawn)) {
    subtitle = c(subtitle, "no difference reaches presence or undecided")
  }
  if(!main_selected) {
    subtitle = c(subtitle, "main-effect differences not selected")
  }
  if(length(subtitle)) {
    # Along the top, left-aligned: qgraph leaves no bottom margin, and the
    # centre of the panel is where the nodes are.
    graphics::mtext(paste(subtitle, collapse = "; "),
      side = 3, line = -1, adj = 0, cex = 0.7, col = "grey55"
    )
  }
  if(isTRUE(legend)) {
    keys = c("difference, positive", "difference, negative", "undecided")
    colors = c(mover_palette()[1:2], "grey65")
    line = c(1L, 1L, 3L)
    if(main_selected) {
      keys = c(keys, "node: main-effect difference")
      colors = c(colors, mover_palette()[1])
      line = c(line, NA)
    }
    graphics::legend("bottomleft",
      legend = keys, col = colors,
      lty = line, lwd = c(rep(2.4, 2), 1.2, rep(NA, length(keys) - 3L)),
      pch = c(rep(NA, 3L), rep(22L, length(keys) - 3L)),
      pt.cex = 1.2, bty = "n", cex = 0.7, text.col = "grey25", xpd = NA
    )
  }
  invisible(result$layout)
}


# ------------------------------------------------------------------
# compare_group_panels
# ------------------------------------------------------------------
# Each group's own network and the difference panel, on one layout so a reader
# compares by position. The layout is computed once from the difference panel's
# node set and reused, which is what makes the three comparable by eye.
# ------------------------------------------------------------------
compare_group_panels = function(x, found, weight, variables, num_groups,
                                layout, legend, ...) {
  num_variables = length(variables)
  nodes = main_difference_nodes(found$main, found$main_selected)

  # Posterior-mean group networks, in the same row-major upper-triangle order
  # as the pairs.
  effects = extract_group_params(x)$pairwise_effects_groups
  group_weight = lapply(seq_len(num_groups), function(g) effects[, g])

  # One layout for every panel: qgraph computes it once with DoNotPlot, on the
  # union of the groups' networks, so a node sits in the same place throughout.
  union_weight = Reduce(pmax, lapply(group_weight, abs))
  union_matrix = matrix(0, num_variables, num_variables)
  union_matrix[found$pairs] = union_weight
  union_matrix[found$pairs[, 2:1, drop = FALSE]] = union_weight
  shared = qgraph::qgraph(
    input = union_matrix, labels = variables, directed = FALSE,
    layout = layout, DoNotPlot = TRUE
  )$layout

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1L, num_groups + 1L))

  for(g in seq_len(num_groups)) {
    weights = matrix(0, num_variables, num_variables)
    weights[found$pairs] = group_weight[[g]]
    weights[found$pairs[, 2:1, drop = FALSE]] = group_weight[[g]]
    panel = list(
      input = weights,
      labels = variables, directed = FALSE,
      layout = shared, fade = FALSE, minimum = 0,
      # qgraph's default green/red sign pair is the one colour-vision
      # deficiency most often collapses; the difference panel's Okabe-Ito pair
      # carries the sign here too, so all three panels read alike.
      posCol = mover_palette()[1], negCol = mover_palette()[2],
      shape = nodes$shape, border.color = nodes$border_color,
      title = sprintf("group %d", g)
    )
    do.call(qgraph::qgraph, utils::modifyList(panel, list(...)))
  }
  draw_difference_network(
    weight, found$pairwise, found$pairs, variables, num_variables,
    nodes, shared, legend, found$main_selected, ...,
    title = "difference"
  )
  invisible(NULL)
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
  title = edge_panel_title(label, as.character(row$verdict), row$log_bf)

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
# edge_panel_title
# ------------------------------------------------------------------
# Two-line title of an edge posterior panel: the edge, then its verdict and
# the evidence as a natural log Bayes factor. A saturated edge gets an
# inequality; printing exp() of a three-figure log Bayes factor would fill the
# title with a hundred digits none of which the run resolves.
#
# @param label    The edge label.
# @param verdict  The verdict, as a character string.
# @param log_bf   Natural log inclusion Bayes factor.
#
# Returns: a length-one character string.
# ------------------------------------------------------------------
edge_panel_title = function(label, verdict, log_bf) {
  sprintf("%s\n%s, log BF %s", label, verdict, format_log_bf(log_bf))
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
