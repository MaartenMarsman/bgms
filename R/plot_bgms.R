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
# network_unit
# ------------------------------------------------------------------
# The wording one network display uses for the thing it draws. A bgms fit
# draws edges; a bgmCompare fit draws differences between two groups' edges.
# Everything else about the two pictures -- the three-panel split, the
# encoding within each panel, the shared layout, the thresholds -- is one
# routine, so the wording is the only thing that has to be told apart, and it
# is told apart here rather than in two drawing routines that can drift.
#
# @param kind  "edge" or "difference".
#
# Returns: a list of strings.
# ------------------------------------------------------------------
network_unit = function(kind = c("edge", "difference")) {
  kind = match.arg(kind)
  if(kind == "edge") {
    return(list(
      kind = kind,
      weights_label = "Edge weights",
      panels = c(
        presence = "evidence of presence",
        absence = "evidence of absence",
        undecided = "undecided"
      )
    ))
  }
  list(
    kind = kind,
    weights_label = "Difference weights",
    panels = c(
      presence = "difference supported",
      absence = "difference ruled out",
      undecided = "undecided"
    )
  )
}


# ------------------------------------------------------------------
# threshold_rules
# ------------------------------------------------------------------
# The classification rule each panel stands for, in raw Bayes factors. A
# reader does not think in logs, so the number printed under a panel title is
# the Bayes factor itself even though every per-edge number the package prints
# elsewhere is a natural log. The rule is printed because the threshold is the
# caller's to choose: two figures of the same fit at different thresholds are
# different classifications, and nothing else on the figure would say so.
#
# @param evidence_threshold  The inclusion Bayes factor threshold.
#
# Returns: a named character vector, one rule per panel.
# ------------------------------------------------------------------
threshold_rules = function(evidence_threshold) {
  high = format_bayes_factor(evidence_threshold)
  low = format_bayes_factor(1 / evidence_threshold)
  c(
    presence = sprintf("BF > %s", high),
    absence = sprintf("BF < %s", low),
    undecided = sprintf("%s < BF < %s", low, high)
  )
}


# ------------------------------------------------------------------
# panel_edge_matrix
# ------------------------------------------------------------------
# The square matrix qgraph is given for one evidence panel.
#
# A weighted edgelist would be the direct expression of "draw these edges", but
# qgraph derives the node set from the edgelist and its nNodes argument does not
# reliably override that, so a network that leaves a node out fails rather than
# drawing the node alone. A square matrix states the node set in its dimensions,
# and an empty panel is a matrix of zeros, which draws the nodes on their own --
# which is what an all-decided fit's undecided panel should look like.
#
# @param value          Value per selected pair (a weight, or 1 for a panel
#                       drawn at uniform width).
# @param keep           Logical, which pairs this panel draws.
# @param pairs          The row-major upper-triangle index.
# @param num_variables  Number of nodes.
#
# Returns: a symmetric num_variables x num_variables matrix.
# ------------------------------------------------------------------
panel_edge_matrix = function(value, keep, pairs, num_variables) {
  out = matrix(0, num_variables, num_variables)
  if(!any(keep)) {
    return(out)
  }
  # An exactly zero weight would be dropped as a non-edge, which would silently
  # shift every later colour onto the wrong edge. Nudge it instead.
  drawn = value[keep]
  drawn[drawn == 0] = .Machine$double.eps
  index = pairs[keep, , drop = FALSE]
  out[index] = drawn
  out[index[, 2:1, drop = FALSE]] = drawn
  out
}


# ------------------------------------------------------------------
# matrix_edge_colors
# ------------------------------------------------------------------
# Sign-carrying colours for a panel's edges, in the order qgraph reads them.
#
# qgraph reads the non-zero upper triangle of its input in column-major order,
# which is not the row-major order the verdict table is in; a colour vector
# built in the wrong order lands every colour on the wrong edge. Reading the
# signs off the matrix itself is what makes the order right by construction.
#
# The two accents are the Okabe-Ito blue/vermillion pair, so the sign survives
# common forms of colour-vision deficiency.
#
# @param m  A panel matrix from panel_edge_matrix().
#
# Returns: a character vector, one colour per drawn edge.
# ------------------------------------------------------------------
matrix_edge_colors = function(m) {
  palette = mover_palette()
  order_read = which(upper.tri(m) & m != 0, arr.ind = TRUE)
  if(!nrow(order_read)) {
    return(character(0))
  }
  ifelse(m[order_read] >= 0, palette[1], palette[2])
}


# ------------------------------------------------------------------
# shared_network_layout
# ------------------------------------------------------------------
# One layout, computed from every pair, so panels drawn from subsets of those
# pairs match node for node. This is the whole point of the three-panel
# display: a reader compares the panels by position, and a node that moved
# between them would make that impossible.
#
# @param weight         Weight per pair.
# @param pairs          The row-major upper-triangle index.
# @param num_variables  Number of nodes.
# @param layout         Layout passed to qgraph.
# @param labels         Node labels.
#
# Returns: the layout matrix qgraph computed.
# ------------------------------------------------------------------
shared_network_layout = function(weight, pairs, num_variables, layout, labels) {
  full = panel_edge_matrix(abs(weight), rep(TRUE, length(weight)), pairs,
    num_variables)
  qgraph::qgraph(
    input = full, labels = labels, directed = FALSE,
    layout = layout, DoNotPlot = TRUE
  )$layout
}


# ------------------------------------------------------------------
# contrast_magnitude
# ------------------------------------------------------------------
# How much a pair differs anywhere, when it differs in more than one place.
# On more than two groups bgmCompare holds a posterior mean difference matrix
# per contrast, and a pair has K - 1 of them under the single indicator they
# share. The largest of them in absolute value is the summary used to lay the
# network out -- and only to lay it out. It is never drawn: the panels it feeds
# are unweighted, precisely because no one of these numbers is the pair's
# difference. Taking the largest rather than a mean keeps a pair that differs
# sharply in one contrast from being averaged into the middle of the picture.
#
# @param differences  List of per-contrast posterior mean difference matrices.
# @param pairs        The row-major upper-triangle index.
#
# Returns: one non-negative number per pair, in the order of `pairs`.
# ------------------------------------------------------------------
contrast_magnitude = function(differences, pairs) {
  index = cbind(pairs[, 1], pairs[, 2])
  per_contrast = vapply(differences, function(d) abs(d[index]),
    numeric(nrow(pairs)))
  apply(matrix(per_contrast, nrow = nrow(pairs)), 1L, max)
}


# ------------------------------------------------------------------
# network_panel_par
# ------------------------------------------------------------------
# The graphical parameters a network panel is drawn under, and the qgraph
# margin that goes with them. qgraph sets par("mar") itself, so a margin
# reserved before the call does not survive it; its own `mar` folds into the
# coordinate range instead, and the top strip that opens there is where the
# panel title is drawn. The bottom strip is narrow: it once held the sign key,
# and when the key went the space it was holding went with it, so the network
# fills the panel rather than sitting above a reserved blank.
#
# Every panel of a figure uses the same two, whatever it draws, because panels
# that differ in margin differ in scale and stop being comparable.
#
# These are constants and not measurements, which is a real limitation on a
# panel far from square. qgraph writes its coordinate range straight from
# `mar` while drawing nodes at a physical size, so the two cannot be brought
# into agreement from here: widening the margin pulls the layout in without
# shrinking a node. Three panels in one row of a 7x7 device therefore still
# crowd, and the display wants a wide one. Review finding F-115, routed to a
# later batch.
#
# Returns: list(mar = par margin, qgraph_mar = qgraph's margin).
# ------------------------------------------------------------------
network_panel_par = function() {
  list(mar = c(0.4, 0.4, 0.4, 0.4), qgraph_mar = c(3, 3, 6.5, 3))
}


# ------------------------------------------------------------------
# draw_network_panel
# ------------------------------------------------------------------
# One network panel: the drawing, then its title and the classification rule
# it stands for. The title names the panel and counts what is in it; it never
# states a conclusion, which is print()'s and verdicts()' business.
#
# @param m          The panel matrix.
# @param variables  Node labels.
# @param shared     The shared layout.
# @param title      Panel title, already composed.
# @param rule       Muted second line, or NULL.
# @param nodes      Node ring channel from main_difference_nodes(), or NULL.
# @param style      The style list.
# @param ...        Passed to qgraph::qgraph().
# ------------------------------------------------------------------
draw_network_panel = function(m, variables, shared, title, rule = NULL,
                              nodes = NULL, style = bgms_style(), ...) {
  geometry = network_panel_par()
  graphics::par(mar = geometry$mar)
  arguments = list(
    input = m,
    labels = variables,
    directed = FALSE,
    layout = shared,
    # qgraph fades an edge toward the background in proportion to its weight,
    # which would wash the encoding out: a supported but weak edge would come
    # out as faint as one drawn at uniform width. Width already carries what
    # there is to carry.
    fade = FALSE,
    minimum = 0,
    color = style$pale,
    border.color = style$muted,
    label.color = style$ink,
    label.scale.equal = TRUE,
    vsize = 16,
    mar = geometry$qgraph_mar
  )
  if(!is.null(nodes$pie)) {
    arguments$pie = nodes$pie
    arguments$pieColor = nodes$pie_color
  }
  result = do.call(qgraph::qgraph, utils::modifyList(arguments, list(...)))

  # The title sits in the strip qgraph's own margin opened above the nodes,
  # inset and spaced by its own type height so it holds its place on any
  # device. It is the largest type on the figure: it is what a reader looks at
  # first and reads across three panels rather than up close.
  usr = graphics::par("usr")
  left = usr[1] + 0.03 * diff(usr[1:2])
  # Sized by measurement rather than by a constant that happens to fit one
  # device: the title is drawn as large as the style asks for, or as large as
  # the panel has room for, whichever is smaller. A narrow device shrinks it
  # instead of clipping it.
  available = 0.94 * diff(usr[1:2])
  cex_title = style$cex_panel_title
  wide = graphics::strwidth(title, cex = cex_title)
  if(wide > available) cex_title = cex_title * available / wide
  height = graphics::strheight("Mg", cex = cex_title)

  graphics::text(left, usr[4] - 0.55 * height, title,
    adj = c(0, 1), cex = cex_title, font = 2, col = style$ink, xpd = NA
  )
  if(!is.null(rule)) {
    graphics::text(left, usr[4] - 1.62 * height, rule,
      adj = c(0, 1), cex = style$cex_annotation, col = style$muted, xpd = NA
    )
  }
  invisible(result)
}


# ------------------------------------------------------------------
# draw_evidence_panels
# ------------------------------------------------------------------
# The package's network display: three panels on one shared layout -- the
# pairs the data support, the pairs the data rule out, the pairs the data
# cannot decide.
#
# A single drawing has to make every pair either an edge or a blank, and a
# blank cannot say whether the data ruled the pair out or simply had too little
# to say. Because there is an inclusion Bayes factor for every pair, the choice
# does not have to be made: the two kinds of blank get a panel each.
#
# Only the supported panel is weighted. There, line width is the posterior mean
# association and colour its sign, because that is where the effect sizes live.
# The other two are drawn at uniform width, dashed and dotted, because for those
# pairs the classification is the result and a weight would suggest otherwise.
#
# The supported panel goes unweighted too when a pair has no single magnitude to
# draw -- a bgmCompare fit on more than two groups, where one indicator per pair
# is shared across K - 1 contrasts. The split is as well defined there as it
# ever is; only the width channel has nothing to carry, so it is dropped rather
# than filled with a summary the model never formed.
#
# Both plot methods call this; a bgms edge and a bgmCompare difference are
# split, laid out, encoded and titled by the same code.
#
# @param weight              Posterior mean per pair, row-major upper triangle.
#                            With weighted = FALSE it only orders the layout.
# @param verdict             Verdict per pair, in the same order.
# @param pairs               The row-major upper-triangle index.
# @param variables           Node labels.
# @param evidence_threshold  The threshold the verdicts were read at.
# @param unit                From network_unit().
# @param nodes               Node ring channel, or NULL.
# @param weighted            Whether the supported panel carries width and sign.
# @param layout, ...         As the plot methods take them.
#
# Returns: invisibly, the shared layout.
# ------------------------------------------------------------------
draw_evidence_panels = function(weight, verdict, pairs, variables,
                                evidence_threshold, unit, nodes = NULL,
                                weighted = TRUE, layout = "spring", ...) {
  style = bgms_style()
  num_variables = length(variables)
  shared = shared_network_layout(weight, pairs, num_variables, layout, variables)
  rules = threshold_rules(evidence_threshold)

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1L, 3L))

  classes = c("presence", "absence", "undecided")
  for(class in classes) {
    keep = !is.na(verdict) & verdict == class
    carries_weight = weighted && class == "presence"
    m = panel_edge_matrix(
      if(carries_weight) weight else rep(1, length(weight)), keep, pairs,
      num_variables
    )
    draw_network_panel(
      m, variables, shared,
      title = sprintf("%s: %d", unit$panels[[class]], sum(keep)),
      rule = rules[[class]],
      nodes = nodes, style = style,
      edge.color = switch(class,
        # Unweighted, the supported panel has no sign to colour by and takes
        # the plain ink; the classification is the whole of what it says.
        presence = if(carries_weight) matrix_edge_colors(m) else style$ink,
        # A ruled-out pair is a finding and takes the darker ink; an undecided
        # one is not, and recedes.
        absence = grDevices::adjustcolor(style$ink, 0.8),
        undecided = style$muted
      ),
      lty = switch(class, presence = 1L, absence = 2L, undecided = 3L),
      edge.width = if(carries_weight) NULL else 1.6,
      ...
    )
  }
  invisible(shared)
}


# ------------------------------------------------------------------
# draw_weight_network
# ------------------------------------------------------------------
# The display for a fit that was not run under selection: one panel, every
# pair drawn, width the posterior mean and colour its sign.
#
# There is no indicator and so no inclusion Bayes factor, which means there is
# no evidence to split the pairs by -- the analysis is an estimation one, and
# the honest picture is the estimate. The panel title says which channel the
# figure is drawn in, so a reader never has to work out whether a wide line
# means a large association or strong evidence for one.
#
# @param weight     Posterior mean per pair, row-major upper triangle.
# @param pairs      The row-major upper-triangle index.
# @param variables  Node labels.
# @param unit       From network_unit().
# @param layout, ... As the plot methods take them.
#
# Returns: invisibly, the layout qgraph used.
# ------------------------------------------------------------------
draw_weight_network = function(weight, pairs, variables, unit,
                               layout = "spring", ...) {
  style = bgms_style()
  num_variables = length(variables)
  m = panel_edge_matrix(weight, rep(TRUE, length(weight)), pairs, num_variables)

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  result = draw_network_panel(
    m, variables, layout,
    title = sprintf("%s: %d", unit$weights_label, sum(weight != 0)),
    rule = "posterior mean, every pair drawn",
    style = style,
    edge.color = matrix_edge_colors(m),
    lty = 1L,
    ...
  )
  invisible(result$layout)
}


#' @title Plot a Fitted bgms Model
#'
#' @description
#' Draws the network as three panels split by what the data settle about each
#' pair, or one of the other standard displays.
#'
#' @param x A fitted model object of class `bgms`, from [bgm()].
#' @param type Character; which display to draw. `"network"` (default) is the
#'   edge evidence plot; `"centrality"` is the posterior strength centrality of
#'   [extract_centrality()].
#' @param evidence_threshold Numeric > 1; the inclusion Bayes factor separating
#'   evidence of presence from undecided, as in [verdicts()]. Default `10`.
#' @param layout Layout passed to [qgraph::qgraph()]. Default `"spring"`.
#' @param ... Passed to [qgraph::qgraph()] for `type = "network"`, and to
#'   [plot.bgms_centrality()] otherwise.
#'
#' @return `x`, invisibly. Called for the side effect of drawing.
#'
#' @details
#' \strong{The edge evidence plot.} A single network drawing has to make every
#' pair either an edge or a blank, and a blank cannot say whether the data ruled
#' the pair out or simply had too little to say about it. `bgm()` returns an
#' inclusion Bayes factor for every pair, so that choice does not have to be
#' made: the network is drawn as three panels on one shared layout -- the pairs
#' the data support, the pairs the data rule out, and the pairs the data cannot
#' decide. Each panel is titled with what it holds and how many pairs are in it,
#' and with the rule that put them there, stated as a Bayes factor rather than
#' its logarithm. The classification is the one [verdicts()] reports at the same
#' `evidence_threshold`.
#'
#' Only the first panel is weighted. There, line width is the posterior mean
#' pairwise association and colour carries its sign -- blue for a positive
#' association, vermillion for a negative one, the Okabe-Ito pair, so the sign
#' survives common forms of colour-vision deficiency -- because that is where
#' the effect sizes are. The other two panels are drawn at uniform width, dashed
#' for evidence of absence and dotted for undecided: for those pairs the
#' classification is the result, and a width would suggest an effect size that
#' the data have either ruled out or not established.
#'
#' The figure carries no key. Each panel is titled with the evidence class it
#' holds and the rule that defines it, which is what a key would otherwise
#' repeat; the sign convention is documented here rather than reprinted on
#' every figure.
#'
#' The layout is computed once from every pair and reused, so a node sits in the
#' same place in all three panels and a reader compares them by position.
#'
#' A fit with no edge left in a panel is a result, not a failure. An all-absence
#' fit fills the second panel and leaves the first empty, which is the honest
#' picture of it.
#'
#' \strong{Without edge selection} there is no indicator and so no inclusion
#' Bayes factor, and nothing to split the pairs by: the analysis is an
#' estimation one. `plot()` then draws one panel with every pair on it, width
#' the posterior mean association and colour its sign, and titles it
#' `"Edge weights"` so that the channel the figure is drawn in is never in
#' doubt.
#'
#' Drawing the network requires the suggested package qgraph. The structure
#' plots and the other network displays live in the easybgm package, which
#' builds on these fits.
#'
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:5], display_progress = "none")
#' plot(fit)
#' plot(fit, type = "centrality")
#' }
#'
#' @seealso [verdicts()] for the table the panels encode,
#'   [extract_centrality()], [plot_edge_posterior()] for one edge in detail
#' @family posterior-methods
#' @export
plot.bgms = function(x,
                     type = c("network", "centrality"),
                     evidence_threshold = 10,
                     layout = "spring",
                     ...) {
  type = match.arg(type)
  if(type == "centrality") {
    plot(extract_centrality(x), ...)
    return(invisible(x))
  }
  require_qgraph()

  weight = colMeans(extract_pairwise_interactions(x))
  arguments = extract_arguments(x)
  variables = arguments$data_columnnames
  num_variables = length(variables)

  # The weights and the verdicts are both in the fit's raw indicator order, so
  # the pair positions have to be read off in that order too.
  pairs = indicator_pair_index(x, num_variables)
  unit = network_unit("edge")

  if(!isTRUE(arguments$edge_selection)) {
    draw_weight_network(weight, pairs, variables, unit, layout, ...)
    return(invisible(x))
  }

  edges = verdicts(x, evidence_threshold = evidence_threshold)
  draw_evidence_panels(
    weight = weight,
    verdict = as.character(edges$verdict),
    pairs = pairs,
    variables = variables,
    evidence_threshold = evidence_threshold,
    unit = unit,
    layout = layout,
    ...
  )
  invisible(x)
}


# ------------------------------------------------------------------
# require_qgraph
# ------------------------------------------------------------------
# qgraph is suggested rather than required, so every drawing entry point has
# to say so in the same words and point at the table that carries the same
# information.
# ------------------------------------------------------------------
require_qgraph = function() {
  if(!requireNamespace("qgraph", quietly = TRUE)) {
    stop(
      "Drawing the network needs the qgraph package, which is suggested rather ",
      "than required by bgms. Install it with install.packages(\"qgraph\"), or ",
      "use verdicts() for the same information as a table."
    )
  }
  invisible(TRUE)
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
    pairwise_log_bf = table$log_bf[!is_main],
    main = main,
    main_pip = table$pip[is_main],
    main_log_bf = table$log_bf[is_main],
    pairs = pairs,
    # Without main_difference_selection those indicators are never updated and
    # carry no verdict, so the node channel has nothing to say.
    main_selected = !all(is.na(main))
  )
}


# ------------------------------------------------------------------
# main_difference_nodes
# ------------------------------------------------------------------
# Node encoding of the main-effect difference evidence: each node wears a ring
# around its circle -- qgraph's pie channel -- filled to its difference
# indicator's posterior inclusion probability, a full ring being 1 and a half
# ring 0.5, coloured by the verdict. The fill fraction carries the number, so
# the encoding does not ride on colour alone. Without main_difference_selection
# there is no indicator and no ring.
#
# The ring is drawn around the node rather than beside it, which is where a
# reader of a qgraph network looks for a node's own quantity. The package's
# standalone probability wheel is the mark for a probability a panel prints;
# on a network the ring is the mark.
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
  color = rep(bgms_style()$muted, length(verdict))
  color[verdict %in% "presence"] = mover_palette()[1]
  color[verdict %in% "absence"] = bgms_style()$pale
  list(pie = pie, pie_color = color)
}


#' @title Plot a Fitted bgmCompare Model
#'
#' @description
#' Draws the group differences as three panels split by what the data settle
#' about each pair, or the groups' own networks on a shared layout.
#'
#' @param x A fitted model object of class `bgmCompare`, from [bgmCompare()].
#' @param type Character; which display to draw. `"difference"` (default) is
#'   the difference evidence plot; `"groups"` draws each group's own network on
#'   one shared layout; `"centrality"` is the posterior strength centrality of
#'   [extract_centrality()].
#' @param evidence_threshold Numeric > 1; the inclusion Bayes factor separating
#'   evidence of a difference from undecided, as in [verdicts()]. Default `10`.
#' @param group For `type = "centrality"`: a single group index, passed to
#'   [extract_centrality()] for that group's centrality. A difference in
#'   centrality (two indices) can be extracted and summarized but has no plot;
#'   see [extract_centrality()] for the interpretation caveat. Default `1`.
#' @param layout Layout passed to [qgraph::qgraph()]. Default `"spring"`.
#' @param max_panels For `type = "groups"`: how many group networks are drawn
#'   at once. A fit with more groups than this is drawn a page at a time.
#'   Default `3`.
#' @param page For `type = "groups"`: which page of `max_panels` networks to
#'   draw. Default `1`.
#' @param ... Passed to [qgraph::qgraph()], or to [plot.bgms_centrality()] for
#'   `type = "centrality"`.
#'
#' @return `x`, invisibly. Called for the side effect of drawing.
#'
#' @details
#' \strong{The difference evidence plot.} The default picture is about
#' differences, because differences are what [bgmCompare()] parameterizes, and
#' it is the display [plot.bgms()] uses, read for differences: three panels on
#' one shared layout -- the pairs whose difference the data support, the pairs
#' whose difference the data rule out, and the pairs the data cannot decide.
#' Each panel is titled with what it holds and how many pairs are in it, and
#' with the rule that put them there, stated as a Bayes factor rather than its
#' logarithm. The classification is the one `verdicts(fit)` reports at the same
#' `evidence_threshold`.
#'
#' Only the first panel is weighted: line width is the posterior mean difference
#' and colour carries its sign -- blue for a positive difference, vermillion
#' for a negative one, the Okabe-Ito pair, so the sign survives common forms of
#' colour-vision deficiency. Which group a positive difference favours follows
#' the contrast coding, which [extract_group_params()] reports per group. The
#' other two are drawn at uniform width, dashed for a difference the data rule
#' out and dotted for undecided, because for those pairs the classification is
#' the result.
#'
#' The figure carries no key: the panel titles name the evidence class and the
#' rule that defines it, and the sign convention is documented here.
#'
#' "The groups do not differ anywhere" is a common and correct finding, and it
#' is what a filled second panel and an empty first panel say.
#'
#' \strong{More than two groups.} [bgmCompare()] gives each pair a single
#' inclusion indicator shared across all `K - 1` contrasts, so the three-way
#' split is exactly as well defined for `K > 2` as it is for two groups and the
#' panels are the same. What changes is the width channel: a pair then has
#' `K - 1` posterior mean differences rather than one, and no single number is
#' "the" difference. The first panel is therefore drawn unweighted for `K > 2`
#' -- uniform width, plain ink, no sign colour -- because there the
#' classification is the whole of what the panel reports. Read the magnitudes
#' where they are per group: `plot(fit, type = "groups")` for the picture,
#' [extract_group_params()] for the numbers.
#'
#' The panels need the split to exist, and more than two groups \emph{without}
#' `difference_selection` is the one case where it does not. Selection off
#' means weights are the display -- as it does for a [bgm()] fit -- and beyond
#' two groups the weights of a pair are `K - 1` numbers rather than one, whose
#' honest weighted picture is the groups themselves. `plot()` therefore draws
#' the `type = "groups"` panels in that case, with the same layout, the same
#' `max_panels` paging and the same passthrough to [qgraph::qgraph()].
#' [extract_group_params()] remains the way to read the differences as numbers.
#'
#' \strong{Main-effect differences are not edges.} When
#' `main_difference_selection = TRUE` gave them their own indicators, their
#' evidence is carried on the nodes: each node wears a ring filled to its
#' difference indicator's posterior inclusion probability -- a full ring is
#' probability 1, half a ring 0.5 -- coloured by the verdict (accented for
#' presence, grey for undecided, faint for absence). The fill fraction carries
#' the number, so the encoding does not rest on colour alone. Under the default
#' `main_difference_selection = FALSE` those indicators do not exist and no ring
#' is drawn. `verdicts()` remains the place to read main-effect differences
#' precisely; the rings are a summary of it.
#'
#' \strong{Without difference selection} there is no indicator and so no
#' inclusion Bayes factor to split the pairs by. `plot()` then draws one panel
#' with every pair on it, width the posterior mean difference and colour its
#' sign, titled `"Difference weights"` so that the channel the figure is drawn
#' in is never in doubt. On more than two groups a pair has `K - 1` differences
#' and no one of them is the weight to draw, so `plot()` draws the groups'
#' own networks instead -- the `type = "groups"` display, reached without
#' asking for it, because with selection off that is what the weighted picture
#' of such a fit is.
#'
#' \strong{`type = "groups"`} draws each group's own posterior mean network on
#' the layout the difference display uses, so a node sits in the same place
#' throughout and a reader compares by position. Those panels are estimates,
#' not evidence, so every pair is drawn with its weight. A fit with more groups
#' than `max_panels` is paged rather than squeezed. The difference evidence is
#' the default display and is not repeated here.
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
#' @seealso [verdicts()] for the table the panels encode,
#'   [extract_centrality()], [prior_sensitivity_check()]
#' @family posterior-methods
#' @export
plot.bgmCompare = function(x,
                           type = c("difference", "groups", "centrality"),
                           evidence_threshold = 10,
                           group = 1,
                           layout = "spring",
                           max_panels = 3L,
                           page = 1L,
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
  require_qgraph()

  arguments = extract_arguments(x)
  variables = arguments$data_columnnames
  num_variables = length(variables)
  num_groups = as.integer(arguments$num_groups)
  unit = network_unit("difference")

  differences = x@posterior_mean_pairwise_differences
  pairs = which(upper.tri(matrix(0, num_variables, num_variables)),
    arr.ind = TRUE)
  pairs = pairs[order(pairs[, "row"], pairs[, "col"]), , drop = FALSE]

  if(type == "groups") {
    compare_group_panels(x, pairs, variables, num_groups, layout,
      max_panels, page, ...)
    return(invisible(x))
  }

  # Two groups make one contrast, so a pair has one difference and the
  # supported panel can carry it as width. More groups make K - 1 contrasts
  # under the one indicator the pair shares, so the split is unchanged and only
  # the width channel is dropped. Magnitudes are then read from type = "groups"
  # and extract_group_params().
  weighted = !is.list(differences)
  weight = if(weighted) {
    differences[cbind(pairs[, 1], pairs[, 2])]
  } else {
    contrast_magnitude(differences, pairs)
  }

  if(!isTRUE(arguments$difference_selection)) {
    # Selection off means weights are the display, exactly as it does for a
    # bgm() fit. On two groups the weights of a pair are one number and the
    # difference network is that display. On more than two they are K - 1
    # numbers, and the honest weighted picture of them is not a summary of the
    # differences but the groups themselves -- which is the display
    # `type = "groups"` already draws, on the same layout and with the same
    # paging. So this dispatches there rather than refusing.
    if(!weighted) {
      compare_group_panels(x, pairs, variables, num_groups, layout,
        max_panels, page, ...)
      return(invisible(x))
    }
    draw_weight_network(weight, pairs, variables, unit, layout, ...)
    return(invisible(x))
  }

  found = compare_difference_verdicts(x, evidence_threshold)
  nodes = main_difference_nodes(found$main, found$main_pip, found$main_selected)
  draw_evidence_panels(
    weight = weight,
    verdict = found$pairwise,
    pairs = found$pairs,
    variables = variables,
    evidence_threshold = evidence_threshold,
    unit = unit,
    nodes = nodes,
    weighted = weighted,
    layout = layout,
    ...
  )
  invisible(x)
}


# ------------------------------------------------------------------
# compare_group_panels
# ------------------------------------------------------------------
# Each group's own posterior mean network, on the layout the difference
# display uses, so a node sits in the same place throughout and a reader
# compares by position.
#
# These panels are estimates rather than evidence: bgmCompare's indicators are
# on the differences, not on either group's edges, so there is no per-group
# inclusion Bayes factor to split a group's own network by. Every pair is drawn
# with its weight, and the panel title says so.
#
# Three networks side by side is already at the limit of what a default device
# can show; a six-group fit crammed into one row is not a figure. Pages, on the
# calibration check's pattern, are what makes the display honest for K > 2.
#
# @param x           The fit.
# @param pairs       The row-major upper-triangle index.
# @param variables   Node labels.
# @param num_groups  Number of groups.
# @param layout      Layout passed to qgraph for the shared construction.
# @param max_panels  Networks per page.
# @param page        Which page.
# @param ...         Passed to qgraph::qgraph().
# ------------------------------------------------------------------
compare_group_panels = function(x, pairs, variables, num_groups, layout,
                                max_panels = 3L, page = 1L, ...) {
  check_positive_integer(max_panels, "max_panels")
  check_positive_integer(page, "page")
  max_panels = as.integer(max_panels)
  page = as.integer(page)

  num_pages = max(1L, ceiling(num_groups / max_panels))
  if(page > num_pages) {
    stop(
      "Argument 'page' is ", page, ", but ", num_groups, " group",
      if(num_groups == 1L) "" else "s", " at ", max_panels,
      " panels a page make ", num_pages, " page",
      if(num_pages == 1L) "" else "s", "."
    )
  }
  first = (page - 1L) * max_panels + 1L
  shown = seq.int(first, min(first + max_panels - 1L, num_groups))
  if(num_pages > 1L && isTRUE(getOption("bgms.verbose", TRUE))) {
    message(
      "Showing page ", page, " of ", num_pages, " (", num_groups,
      " groups). Draw the rest with page = ",
      paste(setdiff(seq_len(num_pages), page), collapse = ", "), "."
    )
  }

  style = bgms_style()
  num_variables = length(variables)
  labels = compare_group_labels(extract_arguments(x), num_groups)

  # Posterior-mean group networks, in the same row-major upper-triangle order
  # as the pairs.
  effects = extract_group_params(x)$pairwise_effects_groups
  group_weight = lapply(seq_len(num_groups), function(g) effects[, g])

  # One layout for every panel, computed on the union of the groups' networks
  # so that a node sits in the same place on every page as well as in every
  # panel of one page.
  union_weight = Reduce(pmax, lapply(group_weight, abs))
  shared = shared_network_layout(union_weight, pairs, num_variables, layout,
    variables)

  old_par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  # Every page gets the same number of panel slots, so the last page's networks
  # are drawn at the same size as the first page's and the two pages can be
  # compared. A page with fewer groups than slots leaves the rest blank.
  graphics::par(mfrow = c(1L, min(max_panels, num_groups)))

  for(g in shown) {
    m = panel_edge_matrix(group_weight[[g]], rep(TRUE, nrow(pairs)), pairs,
      num_variables)
    draw_network_panel(
      m, variables, shared,
      title = group_tag(labels, g),
      rule = "posterior mean network",
      style = style,
      edge.color = matrix_edge_colors(m),
      lty = 1L,
      ...
    )
  }
  invisible(shared)
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
#' @param evidence_threshold Numeric > 1; the threshold [verdicts()] is called
#'   at to read this edge's row. Default `10`. The panel prints no verdict, and
#'   neither the inclusion probability nor the Bayes factor depends on the
#'   threshold, so this does not change what is drawn. Ignored for a fit
#'   without edge selection, which has no inclusion row to read.
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
#' the JASP figure: both ordinates at zero are marked with grey dots, and their
#' ratio -- prior over posterior -- is the Bayes factor for the edge, printed on
#' the same natural-log scale. The wheel is filled by \eqn{BF/(1 + BF)}, which
#' is the posterior probability that the edge is there when the two hypotheses
#' are equally likely before seeing the data; the filled share is that
#' probability and the pale share its complement. JASP labels those two shares
#' `data|H1` and `data|H0`; this panel does not, because the notation cannot be
#' read without knowing the convention, and says it here instead.
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
#' wheel (nearly all pale) and the evidence, and the missing accented curve is
#' itself the statement: decisive absence is a result, not a failure.
#'
#' The panel follows the package's plotting conventions (see the internal
#' `R/plot_style.R`): offset axes, no box, large type, no headline title. It
#' carries no caption; what each mark means is stated above rather than
#' reprinted under every figure the panel draws. The edge is named as an
#' annotation, where the rest of the numbers are.
#'
#' No verdict word is printed. What the panel shows is the evidence -- the
#' wheel, the inclusion probability and the log Bayes factor -- and the reading
#' those license is the reader's to make at a threshold they choose;
#' [verdicts()] is where the package states verdicts, and it names the
#' threshold it used.
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
      evidence$pip, evidence$log_bf)
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
# The verdicts() row of one edge: its inclusion probability and its evidence.
# A mixed fit names its indicators in block order, which need not be variable
# order, so the other orientation is tried before giving up.
#
# Neither number depends on the threshold -- the threshold only decides which
# verdict word verdicts() attaches to them, and the panel prints no verdict --
# so it is passed on to verdicts() and has no further say here.
#
# @param bgms_object        The fit.
# @param label              The edge label, as the pairwise draws name it.
# @param evidence_threshold The threshold verdicts() is called at.
#
# Returns: list(pip, log_bf); both NA when the edge has no row.
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
    return(list(pip = NA_real_, log_bf = NA_real_))
  }
  list(pip = row$pip[1], log_bf = row$log_bf[1])
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
# The verdict is deliberately not shown. The panel's business is the evidence:
# the wheel, the inclusion probability and the log Bayes factor say where the
# edge stands, and a reader draws the conclusion those numbers license.
# verdicts() is where the package speaks in verdicts, at a threshold the reader
# chose; a word printed here would have asserted one silently.
#
# @param label    The edge label.
# @param draws    The edge's pairwise draws.
# @param prior    From edge_slab_prior(), or NULL.
# @param pip      Posterior inclusion probability.
# @param log_bf   Natural log inclusion Bayes factor.
#
# Returns: a panel description for draw_edge_panel().
# ------------------------------------------------------------------
edge_panel_selection = function(label, draws, prior, pip, log_bf) {
  slab = draws[draws != 0]
  if(!is.finite(pip)) {
    pip = mean(draws != 0)
  }
  posterior = conditional_density(slab)

  style = bgms_style()
  list(
    label = label,
    subtitle = NULL,
    posterior = posterior,
    prior = prior,
    dots = NULL,
    wheel_prob = pip,
    wheel_labels = NULL,
    evidence = c(
      format_inclusion(pip),
      paste("log BF", display_log_bf(log_bf))
    ),
    estimate = if(is.null(posterior)) NULL else estimate_lines(slab),
    interval = if(is.null(posterior)) NULL else stats::quantile(
      slab, c(0.025, 0.975), names = FALSE
    ),
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
    # The wheel carried JASP's data|H1 / data|H0 tags above and below it. They
    # are notation, not language: a reader who has not met the convention
    # cannot decode them, and they crowded the corner the panel's identity and
    # evidence already share. What they said is said in the Rd, where an
    # explanation can be as long as it needs to be without costing a figure
    # anything.
    wheel_labels = NULL,
    evidence = paste("log BF", display_log_bf(log_bf)),
    estimate = estimate_lines(draws),
    interval = stats::quantile(draws, c(0.025, 0.975), names = FALSE),
    style = bgms_style()
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
