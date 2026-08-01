# ==============================================================================
# The package's plotting law: the JASP / R Graph Compendium panel style
# ==============================================================================
#
# Every base-graphics panel bgms draws follows one set of conventions, borrowed
# from the figures of the R Graph Compendium (Wagenmakers & Gronau) and from
# JASP's own plots, whose working constants live in jaspGraphs. The rules, in
# the order they matter:
#
#   1. NO BIG TITLES. A JASP plot carries annotations, not a headline. What a
#      title would say -- the quantity, the evidence, the estimate -- is said by
#      the numbers printed in the top margin, where a reader looking for them
#      will find them next to the thing they describe. A panel may carry a
#      compact identifying label; it does not carry a `main =` banner.
#   2. OFFSET AXES. Each axis is drawn only across the range of its own tick
#      marks (a Tufte range frame; `geom_rangeframe()` in jaspGraphs), and the
#      data region is slightly larger than that range, so the two axes do not
#      meet in a corner and the frame does not close. `bgms_panel_par()` returns
#      the epsilon that opens that gap.
#   3. NO BOX. `bty = "n"`; the axes are the only frame.
#   4. LARGE FONTS. jaspGraphs works at a 17pt base with axis titles at 1.2x and
#      legends at 1.25x of it. The base-graphics equivalents are collected in
#      `bgms_style()` and are deliberately larger than R's defaults.
#   5. NEEDLESS INK OMITTED. Grey ink rather than black, no grid, no legend box,
#      fills at low alpha, one accent colour per panel taken from
#      `mover_palette()`.
#   6. PROBABILITIES ARE WHEELS. Any probability a panel reports is also shown
#      as a filled wheel (`probability_wheel()`), whose filled share is the
#      number. This is the same encoding the compare plots put on their node
#      rings, so one visual language carries probability everywhere; printed
#      text carries log Bayes factors, which are not probabilities and get no
#      wheel.
#
# These helpers are the reference implementation of that style. A new panel
# should reach for them rather than re-deriving margins and colours.
# ==============================================================================


# ------------------------------------------------------------------
# bgms_style
# ------------------------------------------------------------------
# The style constants, in one place, so that a panel never hardcodes them.
#
# Typography follows jaspGraphs' graph options (base font 17pt, axis titles at
# 1.2x, legends at 1.25x) expressed as base-graphics `cex` multipliers, and the
# compendium's own base-R figures, which run cex.axis = 1.2, cex.lab = 1.5 and
# lwd = 2 against R's default of 1.
#
# Colours are three greys and one accent: `ink` for anything a reader reads,
# `muted` for structure (axes, captions), `pale` for the unfilled share of a
# wheel, and `accent` for the single quantity the panel is about.
#
# Returns: a named list of style constants.
# ------------------------------------------------------------------
bgms_style = function() {
  list(
    # Colours
    ink = "grey25",
    muted = "grey55",
    pale = "grey88",
    accent = mover_palette()[1],
    fill_alpha = 0.22,

    # Typography (base-graphics cex multipliers)
    cex_axis = 1.2,
    cex_lab = 1.4,
    cex_label = 1.2,
    cex_annotation = 1.1,
    cex_legend = 1.15,
    cex_caption = 0.8,

    # Line weights
    lwd_curve = 2,
    lwd_axis = 1.2,
    lwd_wheel = 1.4,

    # Geometry
    mar = c(5.8, 5.0, 8.4, 2.2),
    mgp = c(3.1, 0.8, 0),
    # The data region overshoots the tick range by this fraction on each side,
    # which is what makes the axes read as offset.
    eps = 0.045,
    # Headroom above the tallest curve, for the interval bar.
    stretch = 1.25,
    # Radius of a probability wheel, in inches, so it is the same size on any
    # device.
    wheel_radius = 0.26
  )
}


# ------------------------------------------------------------------
# bgms_panel_par
# ------------------------------------------------------------------
# Open a panel in the package style: set the graphical parameters and hand
# back the constants the caller still needs (the axis epsilon, the colours,
# the type sizes).
#
# The caller is responsible for restoring `par`; this function returns the old
# settings for that purpose, and the idiom is
#
#   style = bgms_panel_par()
#   on.exit(graphics::par(style$old_par), add = TRUE)
#
# `mar` is the compendium's asymmetric one: a tall top margin, because that is
# where a JASP panel puts its annotations, and a wide left margin, because the
# y axis carries numbers and a label.
#
# @param mar    Margin override, in lines. Default from bgms_style().
# @param ...    Further arguments passed to graphics::par().
#
# Returns: the bgms_style() list with `old_par` added.
# ------------------------------------------------------------------
bgms_panel_par = function(mar = NULL, ...) {
  style = bgms_style()
  if(!is.null(mar)) {
    style$mar = mar
  }
  old_par = graphics::par(no.readonly = TRUE)
  graphics::par(
    mar = style$mar,
    mgp = style$mgp,
    las = 1,
    bty = "n",
    col.axis = style$ink,
    col.lab = style$ink,
    col.main = style$ink,
    cex.axis = style$cex_axis,
    cex.lab = style$cex_lab,
    ...
  )
  style$old_par = old_par
  style
}


# ------------------------------------------------------------------
# bgms_axis_range
# ------------------------------------------------------------------
# Tick positions and the slightly wider data range that offsets the axis from
# them. `pretty()` chooses the ticks, as jaspGraphs' break logic does; the data
# range is the tick range grown by the style's epsilon, so the axis line stops
# short of the corner.
#
# @param values  The values the axis has to cover.
# @param n       Target number of ticks.
# @param eps     Overshoot fraction; default from bgms_style().
#
# Returns: list(at = tick positions, lim = data range).
# ------------------------------------------------------------------
bgms_axis_range = function(values, n = 5, eps = bgms_style()$eps) {
  at = pretty(range(values, finite = TRUE), n = n)
  span = diff(range(at))
  if(!is.finite(span) || span <= 0) {
    span = max(abs(at), 1)
  }
  list(at = at, lim = range(at) + c(-1, 1) * eps * span)
}


# ------------------------------------------------------------------
# bgms_axis
# ------------------------------------------------------------------
# Draw one offset axis: the line spans the ticks and nothing more.
#
# @param side    1 or 2.
# @param at      Tick positions, from bgms_axis_range().
# @param labels  Tick labels; default formats to one decimal, as the
#                compendium's figures do.
# @param style   The style list.
# ------------------------------------------------------------------
bgms_axis = function(side, at, labels = NULL, style = bgms_style()) {
  if(is.null(labels)) {
    digits = if(all(abs(at - round(at, 1)) < 1e-8)) 1L else 2L
    labels = formatC(at, digits = digits, format = "f")
  }
  graphics::axis(
    side, at = at, labels = labels,
    col = style$muted, col.ticks = style$muted, col.axis = style$ink,
    lwd = style$lwd_axis, cex.axis = style$cex_axis
  )
}


# ------------------------------------------------------------------
# probability_wheel
# ------------------------------------------------------------------
# The package's one picture of a probability: a wheel whose filled share is the
# number. `prob` of the circle is drawn in the accent colour and the rest in a
# pale grey, so the reading survives a monochrome print and a reader who cannot
# separate the two hues -- the fraction, not the colour, is the message.
#
# The wedge is centred on the top of the circle, as JASP's own wheels are, so
# that a wheel near 0 and a wheel near 1 are distinguishable at a glance rather
# than by chasing a start angle.
#
# The circle is drawn from cos()/sin() and closed with polygon(), which is what
# the compendium does; the radius is given in inches and converted through the
# current user coordinates, so the wheel is round on any device and any axis
# scale, and can be placed in a margin under `xpd = TRUE`.
#
# @param x,y     Centre, in user coordinates.
# @param prob    Filled share, in [0, 1]. A non-finite value draws an empty
#                wheel rather than failing.
# @param radius  Radius in inches. Default from bgms_style().
# @param labels  Optional length-2 character vector printed above and below the
#                wheel, as JASP labels its data|H1 / data|H0 wheels.
# @param col     Fill colour of the `prob` share. Default the style accent.
# @param style   The style list.
#
# Returns: invisibly, list(x, y, dx, dy) -- the centre and the radius in user
#   units on each axis, so a caller can place text beside the wheel.
# ------------------------------------------------------------------
probability_wheel = function(x, y, prob, radius = NULL, labels = NULL,
                             col = NULL, style = bgms_style()) {
  if(is.null(radius)) radius = style$wheel_radius
  if(is.null(col)) col = style$accent
  if(!is.finite(prob)) prob = 0
  prob = min(max(prob, 0), 1)

  # Inches to user units, separately per axis, so the circle stays round.
  usr = graphics::par("usr")
  pin = graphics::par("pin")
  dx = radius * diff(usr[1:2]) / pin[1]
  dy = radius * diff(usr[3:4]) / pin[2]

  circle = seq(0, 2 * pi, length.out = 181)

  # Always draw the minority share as the wedge, on top of a full disc of the
  # majority colour. A wedge of nearly the whole circle closes on itself and
  # leaves a visible seam where its two radii meet; drawing the small share
  # instead keeps every wedge well under half a turn. The accented share still
  # sits at the top of the wheel either way, because the pale complement of a
  # top-centred wedge is a bottom-centred one.
  minority = min(prob, 1 - prob)
  wedge_at_top = prob <= 0.5
  graphics::polygon(
    x + dx * cos(circle), y + dy * sin(circle),
    col = if(wedge_at_top) style$pale else col, border = NA
  )
  if(minority > 1e-9) {
    centre = if(wedge_at_top) pi / 2 else -pi / 2
    half = pi * minority
    arc = seq(centre - half, centre + half, length.out = 181)
    graphics::polygon(
      c(x, x + dx * cos(arc)), c(y, y + dy * sin(arc)),
      col = if(wedge_at_top) col else style$pale, border = NA
    )
  }
  graphics::lines(
    x + dx * cos(circle), y + dy * sin(circle),
    col = style$muted, lwd = style$lwd_wheel
  )

  if(!is.null(labels)) {
    gap = 0.55 * graphics::strheight("M", cex = style$cex_annotation)
    graphics::text(x, y + dy + gap, labels[1],
      cex = style$cex_annotation, col = style$ink, adj = c(0.5, 0)
    )
    graphics::text(x, y - dy - gap, labels[2],
      cex = style$cex_annotation, col = style$ink, adj = c(0.5, 1)
    )
  }

  invisible(list(x = x, y = y, dx = dx, dy = dy))
}


# ------------------------------------------------------------------
# annotation_block
# ------------------------------------------------------------------
# A stack of annotation lines placed at a corner of the panel: the median, the
# credible interval, the evidence. This is where a JASP panel puts what a title
# would otherwise say, so the placement rule is worth having in one function --
# lines are spaced by their own height rather than by a guessed constant, and
# the block grows downward from its anchor whatever the device size.
#
# @param x,y      Anchor, in user coordinates. `y` is the top of the first line.
# @param lines    Character vector, one element per line.
# @param adj      Horizontal alignment: 0 anchors the left edge, 1 the right.
# @param cex      Type size; default the style's annotation size.
# @param col      Colour; default the style ink.
# @param spacing  Line pitch, as a multiple of line height.
# @param style    The style list.
#
# Returns: invisibly, the y coordinate of the last baseline drawn.
# ------------------------------------------------------------------
annotation_block = function(x, y, lines, adj = 0, cex = NULL, col = NULL,
                            spacing = 1.5, style = bgms_style()) {
  if(is.null(cex)) cex = style$cex_annotation
  if(is.null(col)) col = style$ink
  lines = lines[!is.na(lines)]
  if(!length(lines)) {
    return(invisible(y))
  }
  pitch = spacing * graphics::strheight("Mg", cex = cex)
  at = y - pitch * (seq_along(lines) - 1)
  for(k in seq_along(lines)) {
    graphics::text(x, at[k], lines[k], adj = c(adj, 1), cex = cex, col = col)
  }
  invisible(at[length(at)])
}


# ------------------------------------------------------------------
# panel_y / panel_x
# ------------------------------------------------------------------
# Figure-relative placement in user coordinates. The annotation band of a JASP
# panel sits in the top margin, whose extent in user units depends on the
# device; addressing it as a fraction of the figure region keeps the layout
# fixed while the axis scale moves under it. JASP's own code does the same
# through grconvertY(..., "ndc", "user").
#
# @param frac  Fraction of the figure region, 0 at the bottom/left.
#
# Returns: the corresponding user coordinate.
# ------------------------------------------------------------------
panel_y = function(frac) {
  graphics::grconvertY(frac, from = "nfc", to = "user")
}

panel_x = function(frac) {
  graphics::grconvertX(frac, from = "nfc", to = "user")
}


# ------------------------------------------------------------------
# format_probability
# ------------------------------------------------------------------
# A probability in the form a JASP panel prints it: two decimals, no leading
# zero, and an inequality at the ends rather than a rounded ".00" or "1.00"
# that would claim the run resolved a probability it did not.
#
# The relation is part of the returned string, as it is in format_log_bf(), so
# that a caller writes paste("PIP", format_probability(p)) and gets either
# "PIP = .79" or "PIP > .99" without composing the operator itself.
#
# @param p  A probability.
#
# Returns: a length-one character string.
# ------------------------------------------------------------------
format_probability = function(p) {
  if(!is.finite(p)) {
    return("= NA")
  }
  if(p >= 0.995) {
    return("> .99")
  }
  if(p <= 0.005) {
    return("< .01")
  }
  sprintf("= %s", sub("^0", "", sprintf("%.2f", p)))
}
