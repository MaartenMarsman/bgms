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

  # The weights and the verdicts are both in the fit's raw indicator order, so
  # the pair positions have to be read off in that order too.
  pairs = indicator_pair_index(x, num_variables)

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
    main_pip = table$pip[is_main],
    pairs = pairs,
    # Without main_difference_selection those indicators are never updated and
    # carry no verdict, so the node channel has nothing to say.
    main_selected = !all(is.na(main))
  )
}


# ------------------------------------------------------------------
# main_difference_nodes
# ------------------------------------------------------------------
# Node-ring encoding of the main-effect difference evidence: each node wears a
# ring (qgraph's pie channel) filled to its difference indicator's posterior
# inclusion probability -- a full ring is 1, half a ring 0.5 -- coloured by the
# verdict. The fill fraction carries the number, so the encoding does not ride
# on colour alone. Without main_difference_selection there is no indicator and
# no ring.
#
# @param verdict        Verdict per variable, or all NA when never updated.
# @param pip            Posterior inclusion probability per variable.
# @param main_selected  Whether the indicators were updated at all.
#
# Returns: list(pie, pie_color), both NULL when no ring is drawn.
# ------------------------------------------------------------------
main_difference_nodes = function(verdict, pip, main_selected) {
  if(!main_selected) {
    return(list(pie = NULL, pie_color = NULL))
  }
  pie = pip
  pie[!is.finite(pie)] = 0
  color = rep("grey55", length(verdict))
  color[verdict %in% "presence"] = mover_palette()[1]
  color[verdict %in% "absence"] = "grey80"
  list(pie = pie, pie_color = color)
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
#' @param group For `type = "centrality"`: a single group index, passed to
#'   [extract_centrality()] for that group's centrality. A difference in
#'   centrality (two indices) can be extracted and summarized but has no plot;
#'   see [extract_centrality()] for the interpretation caveat. Default `1`.
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
#' TRUE` gave them their own indicators, their evidence is carried on the
#' nodes: each node wears a ring filled to its difference indicator's
#' posterior inclusion probability -- a full ring is probability 1, half a
#' ring 0.5 -- coloured by the verdict (accented for presence, grey for
#' undecided, faint for absence). The fill fraction carries the number, so
#' the encoding does not rest on colour alone. Under the default
#' `main_difference_selection = FALSE` those indicators do not exist, so no
#' ring is drawn and the subtitle names the setting. `verdicts()` remains the
#' place to read main-effect differences precisely; the rings are a summary
#' of it.
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
    if(length(group) == 2L) {
      stop(
        "A difference-centrality plot is not offered: strength centrality ",
        "sums the absolute weights of a node's edges, so a between-group ",
        "difference in it conflates which edges differ with how their signs ",
        "cancel and has no edge-level reading. Plot one group at a time ",
        "(group = 1), or inspect the difference as numbers with ",
        "summary(extract_centrality(fit, group = c(1, 2)))."
      )
    }
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

  nodes = main_difference_nodes(found$main, found$main_pip, found$main_selected)
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
    fade = FALSE,
    minimum = 0,
    title = title,
    DoNotPlot = FALSE
  )
  if(!is.null(nodes$pie)) {
    arguments$pie = nodes$pie
    arguments$pieColor = nodes$pie_color
  }
  # qgraph folds its mar argument into the coordinate range (the nodes live in
  # [-1, 1]), so a wide bottom margin opens a strip under the network that the
  # subtitle and the legend draw into without touching a node.
  arguments$mar = c(8, 3, 3, 3)
  result = do.call(qgraph::qgraph, utils::modifyList(arguments, list(...)))
  usr = graphics::par("usr")

  subtitle = character(0)
  if(!any(drawn)) {
    subtitle = c(subtitle, "no difference reaches presence or undecided")
  }
  if(!main_selected) {
    subtitle = c(subtitle,
      "main-effect differences not under selection (main_difference_selection = FALSE)"
    )
  }
  if(length(subtitle)) {
    graphics::text(
      usr[1] + 0.01 * diff(usr[1:2]), -1.22,
      paste(subtitle, collapse = "; "),
      adj = c(0, 1), cex = 0.7, col = "grey55", xpd = NA
    )
  }
  if(isTRUE(legend)) {
    keys = c("difference, positive", "difference, negative", "undecided")
    colors = c(mover_palette()[1:2], "grey65")
    line = c(1L, 1L, 3L)
    if(!is.null(nodes$pie)) {
      keys = c(keys, "node ring: P(main-effect difference); full ring = 1")
      colors = c(colors, mover_palette()[1])
      line = c(line, NA)
    }
    graphics::legend("bottomleft",
      legend = keys, col = colors,
      lty = line, lwd = c(rep(2.4, 2), 1.2, rep(NA, length(keys) - 3L)),
      pch = c(rep(NA, 3L), rep(21L, length(keys) - 3L)),
      pt.cex = 1.1, pt.lwd = 2, bty = "n", cex = 0.7, text.col = "grey25",
      xpd = NA
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
  nodes = main_difference_nodes(found$main, found$main_pip, found$main_selected)

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
      title = sprintf("group %d", g)
    )
    if(!is.null(nodes$pie)) {
      panel$pie = nodes$pie
      panel$pieColor = nodes$pie_color
    }
    # The difference panel opens a bottom strip for its subtitle and legend;
    # the group panels match it so the three networks share their extent.
    panel$mar = c(8, 3, 3, 3)
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
#' Draws one edge's weight the way JASP draws a parameter: the posterior
#' density against the prior it was updated from, with the evidence for the
#' edge as a filled probability wheel and the estimate printed beside it.
#'
#' @param bgms_object A fitted model object of class `bgms`, from [bgm()].
#' @param variable1,variable2 The two variables naming the edge. Either names
#'   or column positions.
#' @param evidence_threshold Numeric > 1; the threshold used for the verdict
#'   printed on the panel, as in [verdicts()]. Default `10`. Ignored for a fit
#'   without edge selection, which has no inclusion verdict to read.
#' @param binwidth `r lifecycle::badge("deprecated")` The panel no longer
#'   expresses the weight as probability per bin, so this has nothing to set;
#'   it is warned about and ignored.
#' @param ... Ignored.
#'
#' @return `bgms_object`, invisibly. Called for the side effect of drawing.
#'
#' @details
#' The panel shows two curves on one density scale. The solid accented curve is
#' the posterior of the edge weight; the dashed grey curve is the prior it was
#' updated from, computed in closed form from the fit's own `interaction_prior`
#' rather than assumed, so a panel drawn from a `cauchy_prior(scale = 2.5)` fit
#' and one drawn from the `normal_prior(scale = 1)` default do not look alike.
#'
#' What the two curves are, and how the evidence is read off them, follows the
#' model the fit actually used.
#'
#' \strong{With edge selection} (the default), the posterior of the weight is a
#' spike-and-slab: a point mass at zero and a continuous part. The panel draws
#' the continuous part \emph{conditional on the edge being included}, which is a
#' genuine density integrating to one, so the y axis carries numbers. The mass
#' at zero is not drawn as a stem competing with that density; it is carried by
#' the probability wheel, whose accented share is the posterior inclusion
#' probability and whose pale share is the probability that the edge is absent.
#' The evidence is printed as the natural log of the inclusion Bayes factor
#' from [extract_inclusion_bf()], which `bgms` estimates by Rao-Blackwellizing
#' the indicator draws. It is \emph{not} a ratio of densities at zero, so the
#' panel does not draw Savage-Dickey ordinates: dots at zero would assert an
#' estimator the package does not use.
#'
#' \strong{Without edge selection} there is no indicator and no point mass; the
#' posterior of the weight is continuous and the Savage-Dickey density ratio is
#' the licensed estimator of the inclusion Bayes factor. The panel then draws
#' the JASP figure exactly: both ordinates at zero are marked with grey dots,
#' and their ratio -- prior over posterior -- is the Bayes factor for the edge,
#' printed on the same natural-log scale. The wheel is filled by
#' \eqn{BF/(1 + BF)}, which is the posterior probability of the edge at equal
#' prior odds, so the wheel still shows a probability, and it is labelled with
#' JASP's `data|H1` / `data|H0` pair.
#'
#' The prior ordinate at zero is exact. The posterior ordinate is estimated
#' from the draws by a Gaussian kernel density with the Sheather-Jones
#' bandwidth (`stats::density(bw = "SJ")`), evaluated at zero; JASP's own
#' implementations use a logspline fit, which `bgms` does not adopt because it
#' would add a dependency for one number. The two estimators agree closely
#' where the posterior is smooth near zero and diverge where it is not, which
#' is the regime in which a Savage-Dickey Bayes factor is unreliable whatever
#' fits its density.
#'
#' On a block of continuous variables the slab sits on an entry of a precision
#' matrix, whose joint prior is that slab times a prior on the diagonal,
#' restricted to the positive-definite cone. The curve drawn is the slab -- the
#' density the sampler evaluates for that edge, and the one the prior
#' sensitivity machinery reweights on -- which is the marginal prior of the
#' entry only up to that restriction. The restriction pulls the true marginal
#' in, so on a continuous block the drawn prior is the wider of the two. On a
#' discrete block the parameters are unconstrained and the slab is the marginal
#' prior exactly.
#'
#' The median and the 95% credible interval printed at the top right are of the
#' same posterior the density shows -- conditional on inclusion under edge
#' selection, unconditional without it -- and the interval is repeated as a bar
#' above the curve.
#'
#' An edge the data rule out is still drawn. When no retained draw included it
#' there is no conditional posterior to show, so the panel draws the prior, the
#' wheel (nearly all pale), and the evidence, and says so in its caption:
#' decisive absence is a result, not a failure.
#'
#' The panel follows the package's plotting conventions (see the internal
#' `R/plot_style.R`): offset axes, no box, large type, no headline title. The
#' edge and its verdict are printed as annotations, where the rest of the
#' numbers are.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' plot_edge_posterior(fit, "intrusion", "dreams")
#' }
#'
#' @seealso [verdicts()], [extract_pairwise_interactions()],
#'   [extract_inclusion_bf()]
#' @family posterior-methods
#' @export
plot_edge_posterior = function(bgms_object, variable1, variable2,
                               evidence_threshold = 10,
                               binwidth = lifecycle::deprecated(), ...) {
  if(lifecycle::is_present(binwidth)) {
    lifecycle::deprecate_warn(
      "0.2.0.0", "plot_edge_posterior(binwidth =)",
      details = paste(
        "The panel draws the posterior of the weight as a density conditional",
        "on inclusion, not as probability per weight bin, so there is no bin",
        "to set a width for. The inclusion probability is now shown as a",
        "filled wheel. The argument is ignored."
      )
    )
  }
  check_evidence_threshold(evidence_threshold)
  if(inherits(bgms_object, "bgmCompare")) {
    stop(
      "plot_edge_posterior() draws an edge of a single network, and a ",
      "bgmCompare() fit parameterizes differences between networks rather ",
      "than one edge weight. Use plot(fit) for the difference network and ",
      "verdicts(fit) for the per-difference evidence."
    )
  }

  samples = extract_pairwise_interactions(bgms_object)
  arguments = extract_arguments(bgms_object)
  variables = arguments$data_columnnames

  first = resolve_variable(variable1, variables, "variable1")
  second = resolve_variable(variable2, variables, "variable2")
  if(first == second) {
    stop(
      "Arguments 'variable1' and 'variable2' both resolve to '",
      variables[first], "', but an edge joins two different variables. ",
      "Name the other end of the edge, for example plot_edge_posterior(fit, '",
      variables[first], "', '",
      variables[if(first == 1L) 2L else 1L], "')."
    )
  }
  label = edge_column_label(variables[first], variables[second],
    colnames(samples))
  if(is.null(label)) {
    stop("No edge between '", variables[first], "' and '", variables[second], "'.")
  }

  draws = samples[, label]
  prior = edge_slab_prior(bgms_object)

  if(isTRUE(arguments$edge_selection)) {
    evidence = edge_selection_evidence(bgms_object, label, evidence_threshold)
    panel = edge_panel_selection(label, draws, prior,
      evidence$pip, evidence$log_bf, evidence$verdict)
  } else {
    panel = edge_panel_savage_dickey(label, draws, prior)
  }
  draw_edge_panel(panel)

  invisible(bgms_object)
}


# ------------------------------------------------------------------
# edge_column_label
# ------------------------------------------------------------------
# The column name a pair of variables is stored under, in whichever order the
# fit stores it. A single-network fit names an edge in variable order, but a
# mixed fit lays its pairwise draws out by block -- discrete, continuous, cross
# -- and names a cross edge by its discrete end first, which need not be the
# earlier column. Looking for one orientation only would report "no edge"
# for an edge that exists.
#
# @param a, b     The two variable names.
# @param columns  The pairwise draws' column names.
#
# Returns: the matching column name, or NULL.
# ------------------------------------------------------------------
edge_column_label = function(a, b, columns) {
  forward = paste(a, b, sep = "-")
  if(forward %in% columns) {
    return(forward)
  }
  reverse = paste(b, a, sep = "-")
  if(reverse %in% columns) {
    return(reverse)
  }
  NULL
}


# ------------------------------------------------------------------
# edge_slab_prior
# ------------------------------------------------------------------
# The fit's own slab prior on one edge weight, as a function of the weight.
#
# The family and its scale are read from the spec rather than assumed: the
# default changed in 0.2.0 (from Cauchy to Normal), and a panel that hardcoded
# either would misdraw half the fits.
#
# No change of variable is needed. The sampler evaluates the slab on the
# association-scale parameter -- for a continuous block that is -K_ij/2 rather
# than the precision entry itself (ggm_model.cpp, mixed_mrf_gradient.cpp) --
# and extract_pairwise_interactions() reports the draws in exactly that frame
# for every model type. It returns, element for element, the theta that
# anchor_draws() reweights the sensitivity curve on, which is defined as the
# frame the slab prior applies to.
#
# @param bgms_object  The fit.
#
# Returns: list(density = function(w), quantile = function(p), family, scale),
#   or NULL when the fit carries no spec to read the prior from, or carries a
#   slab family with no closed form here.
# ------------------------------------------------------------------
edge_slab_prior = function(bgms_object) {
  spec = get_fit_spec(bgms_object)
  if(is.null(spec) || is.null(spec$prior$interaction_prior_type)) {
    return(NULL)
  }
  prior = spec$prior
  family = tolower(prior$interaction_prior_type)
  scale = prior$pairwise_scale
  alpha = prior$interaction_alpha
  beta = prior$interaction_beta

  density = switch(family,
    normal = function(w) stats::dnorm(w, 0, scale),
    cauchy = function(w) stats::dcauchy(w, 0, scale),
    # mu = logit(Y) with Y ~ Beta(alpha, beta), so the density of mu carries
    # the logistic Jacobian. Far out in either tail the Beta density overflows
    # while the Jacobian underflows; the product tends to zero, so that is what
    # is returned rather than the NaN the floating-point form produces.
    `beta-prime` = function(w) {
      p = stats::plogis(w)
      out = stats::dbeta(p, alpha, beta) * p * (1 - p)
      out[!is.finite(out)] = 0
      out
    },
    NULL
  )
  if(is.null(density)) {
    return(NULL)
  }

  quantile = switch(family,
    normal = function(p) stats::qnorm(p, 0, scale),
    cauchy = function(p) stats::qcauchy(p, 0, scale),
    `beta-prime` = function(p) stats::qlogis(stats::qbeta(p, alpha, beta))
  )

  list(density = density, quantile = quantile, family = family, scale = scale)
}


# ------------------------------------------------------------------
# conditional_density
# ------------------------------------------------------------------
# Kernel density of a set of draws on a grid covering them and zero. The
# Sheather-Jones bandwidth is the package's choice everywhere a density is
# read at a point; it can fail on a degenerate sample, in which case the
# default rule stands in rather than the call failing.
#
# @param draws  The draws.
# @param n      Grid size.
#
# Returns: list(x, y), or NULL when there is nothing to smooth.
# ------------------------------------------------------------------
conditional_density = function(draws, n = 512L) {
  if(length(draws) < 2L || !any(is.finite(draws))) {
    return(NULL)
  }
  span = diff(range(draws))
  if(!is.finite(span) || span <= 0) {
    return(NULL)
  }
  fit = tryCatch(
    stats::density(draws, bw = "SJ", n = n),
    error = function(e) NULL, warning = function(w) NULL
  )
  if(is.null(fit)) {
    fit = stats::density(draws, n = n)
  }
  list(x = fit$x, y = fit$y)
}


# ------------------------------------------------------------------
# edge_selection_evidence
# ------------------------------------------------------------------
# The verdicts() row of one edge: its inclusion probability, its evidence, and
# its verdict. A mixed fit names its indicators in block order, which need not
# be variable order, so the other orientation is tried before giving up.
#
# @param bgms_object        The fit.
# @param label              The edge label, as the pairwise draws name it.
# @param evidence_threshold The verdict threshold.
#
# Returns: list(pip, log_bf, verdict); all NA when the edge has no row.
# ------------------------------------------------------------------
edge_selection_evidence = function(bgms_object, label, evidence_threshold) {
  edges = verdicts(bgms_object, evidence_threshold = evidence_threshold)
  row = edges[edges$parameter == label, , drop = FALSE]
  if(!nrow(row)) {
    ends = strsplit(label, "-", fixed = TRUE)[[1]]
    if(length(ends) == 2L) {
      flipped = paste(ends[2], ends[1], sep = "-")
      row = edges[edges$parameter == flipped, , drop = FALSE]
    }
  }
  if(!nrow(row)) {
    return(list(pip = NA_real_, log_bf = NA_real_, verdict = NA_character_))
  }
  list(
    pip = row$pip[1], log_bf = row$log_bf[1],
    verdict = as.character(row$verdict[1])
  )
}


# ------------------------------------------------------------------
# edge_panel_selection
# ------------------------------------------------------------------
# The panel of an edge from a fit with edge selection: a conditional posterior,
# the prior, and the inclusion probability on the wheel. No Savage-Dickey
# ordinates -- the Bayes factor printed here is the Rao-Blackwellized indicator
# Bayes factor, not a ratio of densities at zero, and marking the ordinates
# would claim otherwise.
#
# @param label    The edge label.
# @param draws    The edge's pairwise draws.
# @param prior    From edge_slab_prior(), or NULL.
# @param pip      Posterior inclusion probability.
# @param log_bf   Natural log inclusion Bayes factor.
# @param verdict  The verdict, as a character string, or NA.
#
# Returns: a panel description for draw_edge_panel().
# ------------------------------------------------------------------
edge_panel_selection = function(label, draws, prior, pip, log_bf, verdict) {
  slab = draws[draws != 0]
  if(!is.finite(pip)) {
    pip = mean(draws != 0)
  }
  posterior = conditional_density(slab)

  style = bgms_style()
  list(
    label = label,
    subtitle = verdict_phrase(verdict),
    posterior = posterior,
    prior = prior,
    dots = NULL,
    wheel_prob = pip,
    wheel_labels = NULL,
    evidence = c(
      paste("PIP", format_probability(pip)),
      paste("log BF", format_log_bf(log_bf))
    ),
    estimate = if(is.null(posterior)) NULL else estimate_lines(slab),
    interval = if(is.null(posterior)) NULL else stats::quantile(
      slab, c(0.025, 0.975), names = FALSE
    ),
    caption = if(is.null(posterior)) {
      "No retained draw included this edge. Pale share of the wheel: P(absent)."
    } else {
      "Density: the weight given inclusion. Pale share of the wheel: P(absent)."
    },
    style = style
  )
}


# ------------------------------------------------------------------
# edge_panel_savage_dickey
# ------------------------------------------------------------------
# The panel of an edge from a fit without edge selection. There is no
# indicator, the posterior of the weight is continuous, and the Savage-Dickey
# density ratio at zero is the licensed estimator of the inclusion Bayes
# factor -- so this is the JASP figure, ordinate dots and all.
#
# @param label  The edge label.
# @param draws  The edge's pairwise draws.
# @param prior  From edge_slab_prior(), or NULL.
#
# Returns: a panel description for draw_edge_panel().
# ------------------------------------------------------------------
edge_panel_savage_dickey = function(label, draws, prior) {
  posterior = conditional_density(draws)

  prior_at_zero = if(is.null(prior)) NA_real_ else prior$density(0)
  posterior_at_zero = if(is.null(posterior)) {
    NA_real_
  } else {
    stats::approx(posterior$x, posterior$y, xout = 0, rule = 2)$y
  }
  # BF for the edge: prior ordinate over posterior ordinate at zero. A
  # posterior that has left zero entirely gives an infinite Bayes factor,
  # which format_log_bf() prints as the capped inequality.
  log_bf = log(prior_at_zero) - log(posterior_at_zero)
  bf = exp(min(log_bf, 700))
  wheel = bf / (1 + bf)
  if(!is.finite(wheel)) {
    wheel = if(is.finite(log_bf) && log_bf > 0) 1 else NA_real_
  }

  dots = if(is.finite(prior_at_zero) && is.finite(posterior_at_zero)) {
    list(x = c(0, 0), y = c(prior_at_zero, posterior_at_zero))
  } else {
    NULL
  }

  list(
    label = label,
    subtitle = "no edge selection",
    posterior = posterior,
    prior = prior,
    dots = dots,
    wheel_prob = wheel,
    wheel_labels = c("data|H1", "data|H0"),
    evidence = paste("log BF", format_log_bf(log_bf)),
    estimate = estimate_lines(draws),
    interval = stats::quantile(draws, c(0.025, 0.975), names = FALSE),
    caption = paste(
      "Savage-Dickey ratio at zero (grey dots).",
      "Accented share: P(edge | data), equal prior odds."
    ),
    style = bgms_style()
  )
}


# ------------------------------------------------------------------
# verdict_phrase
# ------------------------------------------------------------------
# The verdict as a panel says it. "presence" and "absence" are readings of the
# evidence and are named as such; "undecided" is not evidence of anything, so
# it stands alone rather than being wrapped in "evidence of".
#
# @param verdict  The verdict, as a character string, or NA.
#
# Returns: a length-one character string, or NULL for an unread verdict.
# ------------------------------------------------------------------
verdict_phrase = function(verdict) {
  if(is.null(verdict) || is.na(verdict)) {
    return(NULL)
  }
  switch(verdict,
    presence = "evidence of presence",
    absence = "evidence of absence",
    undecided = "undecided",
    verdict
  )
}


# ------------------------------------------------------------------
# estimate_lines
# ------------------------------------------------------------------
# The two estimate lines a JASP panel prints at its top right.
#
# @param draws  The draws the estimate summarizes.
#
# Returns: a length-two character vector.
# ------------------------------------------------------------------
estimate_lines = function(draws) {
  quantiles = stats::quantile(draws, c(0.025, 0.5, 0.975), names = FALSE)
  # A weight that rounds to nothing is nothing, not "-0.00": the sign of a
  # rounded-away quantity is not information the run established.
  shown = round(quantiles, 2)
  shown[shown == 0] = 0
  c(
    sprintf("median = %.2f", shown[2]),
    sprintf("95%% CI [%.2f, %.2f]", shown[1], shown[3])
  )
}


# ------------------------------------------------------------------
# edge_panel_window
# ------------------------------------------------------------------
# The weight range the panel is drawn over: the posterior's own support where
# there is one, and the prior's central 90% where there is not, always
# including zero, which is the point the whole figure is about.
#
# @param panel  The panel description.
#
# Returns: a length-two numeric.
# ------------------------------------------------------------------
edge_panel_window = function(panel) {
  if(!is.null(panel$posterior)) {
    return(range(c(panel$posterior$x, 0)))
  }
  if(!is.null(panel$prior)) {
    edge = panel$prior$quantile(0.95)
    return(c(-abs(edge), abs(edge)))
  }
  c(-1, 1)
}


# ------------------------------------------------------------------
# draw_edge_panel
# ------------------------------------------------------------------
# Render a panel description in the package style. The layout is the JASP
# prior-and-posterior figure's: curves in the plot region, everything a reader
# would otherwise look for in a title arranged across the top margin -- the
# edge and its verdict at the left, the estimate at the right, the evidence on
# and beside the wheel between them.
#
# @param panel  A panel description from edge_panel_selection() or
#   edge_panel_savage_dickey().
# ------------------------------------------------------------------
draw_edge_panel = function(panel) {
  style = bgms_panel_par()
  on.exit(graphics::par(style$old_par), add = TRUE)

  window = edge_panel_window(panel)
  x_axis = bgms_axis_range(window)
  grid = seq(x_axis$lim[1], x_axis$lim[2], length.out = 1024)

  prior_y = if(is.null(panel$prior)) NULL else panel$prior$density(grid)
  posterior_y = if(is.null(panel$posterior)) {
    NULL
  } else {
    stats::approx(panel$posterior$x, panel$posterior$y, xout = grid,
      yleft = 0, yright = 0)$y
  }
  # A fit old enough to carry no spec and an edge no draw included leave both
  # curves empty; max() of nothing warns, so the empty case is named instead.
  heights = c(posterior_y, prior_y, panel$dots$y)
  peak = if(length(heights)) max(heights, na.rm = TRUE) else NA_real_
  if(!is.finite(peak) || peak <= 0) peak = 1
  # The tick range starts at zero and the data region drops just below it, so
  # the density baseline sits above the x axis rather than on it.
  y_axis = bgms_axis_range(c(0, style$stretch * peak))

  graphics::plot(NA, NA,
    xlim = x_axis$lim, ylim = y_axis$lim,
    axes = FALSE, xlab = "", ylab = "", main = ""
  )
  bgms_axis(1, x_axis$at, style = style)
  bgms_axis(2, y_axis$at, style = style)
  graphics::mtext("Edge weight",
    side = 1, line = 2.9, cex = style$cex_lab, col = style$ink
  )
  graphics::mtext("Density",
    side = 2, line = 3.3, las = 0, cex = style$cex_lab, col = style$ink
  )

  if(!is.null(prior_y)) {
    graphics::lines(grid, prior_y,
      col = style$muted, lwd = style$lwd_curve, lty = 3
    )
  }
  if(!is.null(posterior_y)) {
    graphics::polygon(
      c(grid, rev(grid)), c(posterior_y, rep(0, length(posterior_y))),
      col = grDevices::adjustcolor(style$accent, style$fill_alpha), border = NA
    )
    graphics::lines(grid, posterior_y,
      col = style$accent, lwd = style$lwd_curve
    )
  }
  if(!is.null(panel$dots)) {
    graphics::points(panel$dots$x, panel$dots$y,
      pch = 21, col = style$ink, bg = "grey75", cex = 1.5, lwd = style$lwd_axis
    )
  }

  # The credible interval, as a bar above the curve it belongs to.
  if(!is.null(panel$interval) && all(is.finite(panel$interval))) {
    bar = 1.1 * peak
    graphics::arrows(
      panel$interval[1], bar, panel$interval[2], bar,
      angle = 90, code = 3, length = 0.05,
      col = style$ink, lwd = style$lwd_axis
    )
  }

  legend_labels = c(
    if(!is.null(posterior_y)) "Posterior",
    if(!is.null(prior_y)) "Prior"
  )
  if(length(legend_labels)) {
    # The compendium puts the legend on the side the mass is not; with the
    # posterior usually right of zero that is the left, and a left-heavy
    # posterior moves it rather than being drawn over.
    left_heavy = !is.null(posterior_y) &&
      sum(posterior_y[grid < mean(x_axis$lim)]) > sum(posterior_y) / 2
    graphics::legend(
      if(left_heavy) x_axis$lim[2] else x_axis$lim[1], max(y_axis$at),
      legend = legend_labels,
      lty = if(length(legend_labels) == 2L) c(1, 3) else if(is.null(posterior_y)) 3 else 1,
      col = if(length(legend_labels) == 2L) c(style$accent, style$muted) else {
        if(is.null(posterior_y)) style$muted else style$accent
      },
      lwd = style$lwd_curve, bty = "n", cex = style$cex_legend,
      xjust = if(left_heavy) 1 else 0, yjust = 1,
      x.intersp = 0.7, seg.len = 1.3, text.col = style$ink
    )
  }

  draw_edge_annotations(panel, x_axis, style)

  graphics::mtext(panel$caption,
    side = 1, line = 4.4, cex = style$cex_caption, col = style$muted
  )
  invisible(NULL)
}


# ------------------------------------------------------------------
# draw_edge_annotations
# ------------------------------------------------------------------
# The top-margin band: identity and verdict at the left, estimate at the right,
# and the probability wheel with the evidence beside it between them. Positions
# are figure fractions, so the band holds its shape while the axis scale moves
# under it -- the same device JASP's own code uses.
#
# @param panel   The panel description.
# @param x_axis  The x tick/limit list.
# @param style   The style list.
# ------------------------------------------------------------------
draw_edge_annotations = function(panel, x_axis, style) {
  old_xpd = graphics::par("xpd")
  graphics::par(xpd = NA)
  on.exit(graphics::par(xpd = old_xpd), add = TRUE)

  left = x_axis$lim[1]
  right = x_axis$lim[2]

  annotation_block(left, panel_y(0.975), panel$label,
    adj = 0, cex = style$cex_label, col = style$ink, style = style
  )
  if(!is.null(panel$subtitle)) {
    annotation_block(left, panel_y(0.920), panel$subtitle,
      adj = 0, cex = style$cex_annotation, col = style$muted, style = style
    )
  }
  if(!is.null(panel$estimate)) {
    annotation_block(right, panel_y(0.975), panel$estimate,
      adj = 1, cex = style$cex_annotation, col = style$ink, style = style
    )
  }

  # A labelled wheel needs room above and below it for its labels, so it hangs
  # lower in the band than a bare one.
  wheel_y = panel_y(if(is.null(panel$wheel_labels)) 0.815 else 0.795)
  wheel = probability_wheel(
    panel_x(0.20), wheel_y, panel$wheel_prob,
    labels = panel$wheel_labels, style = style
  )
  gap = 0.6 * graphics::strwidth("M", cex = style$cex_annotation)
  top = wheel_y + if(length(panel$evidence) == 2L) {
    0.9 * graphics::strheight("Mg", cex = style$cex_annotation)
  } else {
    0.35 * graphics::strheight("Mg", cex = style$cex_annotation)
  }
  annotation_block(wheel$x + wheel$dx + gap, top, panel$evidence,
    adj = 0, cex = style$cex_annotation, col = style$ink, style = style
  )
  invisible(NULL)
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
