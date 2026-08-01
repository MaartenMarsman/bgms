# ==============================================================================
# The package's plotting style module (R/plot_style.R)
# ==============================================================================
#
# These helpers are the reference implementation the other panels will be
# restyled onto, so what is tested here is the contract they offer rather than
# the appearance they produce: the constants a caller reads, the geometry of
# an offset axis, the arithmetic of a wheel, and the placement rule of an
# annotation block.
# ==============================================================================

test_that("bgms_panel_par sets the style and hands back the old settings", {
  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  before = graphics::par(no.readonly = TRUE)
  style = bgms_panel_par()
  on.exit(graphics::par(style$old_par), add = TRUE, after = FALSE)

  # No box, axis labels running horizontally, and type larger than R's default
  # of 1 -- the three conventions a reader would notice first.
  expect_equal(graphics::par("bty"), "n")
  expect_equal(graphics::par("las"), 1L)
  expect_gt(graphics::par("cex.axis"), 1)
  expect_equal(graphics::par("mar"), style$mar)

  # The caller can restore exactly what was there.
  expect_equal(style$old_par$bty, before$bty)
  expect_equal(style$old_par$mar, before$mar)
})

test_that("bgms_panel_par takes a margin override", {
  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  style = bgms_panel_par(mar = c(1, 2, 3, 4))
  on.exit(graphics::par(style$old_par), add = TRUE, after = FALSE)
  expect_equal(graphics::par("mar"), c(1, 2, 3, 4))
  expect_equal(style$mar, c(1, 2, 3, 4))
})

test_that("an offset axis stops short of the data range on both sides", {
  axis = bgms_axis_range(c(-0.3, 0.9))

  # The ticks are pretty numbers covering the values...
  expect_true(min(axis$at) <= -0.3)
  expect_true(max(axis$at) >= 0.9)
  # ...and the data region overshoots them, which is what opens the gap at the
  # corner. An axis drawn at these ticks cannot close into a box.
  expect_lt(axis$lim[1], min(axis$at))
  expect_gt(axis$lim[2], max(axis$at))

  eps = bgms_style()$eps
  span = diff(range(axis$at))
  expect_equal(axis$lim, range(axis$at) + c(-1, 1) * eps * span)
})

test_that("a density axis starts its ticks at zero and drops below it", {
  # The lowest tick is the density baseline; the data region reaches under it,
  # which is what lifts the curve off the x axis line.
  axis = bgms_axis_range(c(0, 12))
  expect_equal(min(axis$at), 0)
  expect_lt(axis$lim[1], 0)
  expect_gt(axis$lim[2], 12)
})

test_that("a constant value still yields a drawable axis", {
  axis = bgms_axis_range(c(2, 2))
  expect_true(all(is.finite(axis$lim)))
  expect_lt(axis$lim[1], axis$lim[2])
})

test_that("the probability wheel fills the share it is given", {
  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)

  wheel = probability_wheel(0.5, 0.5, 0.25, radius = 0.2)
  expect_equal(wheel$x, 0.5)
  expect_equal(wheel$y, 0.5)
  # The radius is given in inches and converted per axis, so the wheel is round
  # on the device rather than in user units.
  expect_gt(wheel$dx, 0)
  expect_gt(wheel$dy, 0)

  # The ends are drawn rather than refused, and a probability that was never
  # estimated draws an empty wheel instead of failing.
  expect_silent(probability_wheel(0.5, 0.5, 0, radius = 0.1))
  expect_silent(probability_wheel(0.5, 0.5, 1, radius = 0.1))
  expect_silent(probability_wheel(0.5, 0.5, NaN, radius = 0.1))
  expect_silent(probability_wheel(0.5, 0.5, NA_real_, radius = 0.1))
  # Out-of-range input is clamped, not wrapped around the circle.
  expect_silent(probability_wheel(0.5, 0.5, 1.4, radius = 0.1))
  expect_silent(probability_wheel(0.5, 0.5, -0.2, radius = 0.1))
})

test_that("an annotation block stacks downward from its anchor", {
  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)

  last = annotation_block(0.1, 0.9, c("median = 0.31", "95% CI [0.12, 0.48]"))
  expect_lt(last, 0.9)

  # One line ends where it starts, and nothing to say draws nothing.
  expect_equal(annotation_block(0.1, 0.9, "one"), 0.9)
  expect_equal(annotation_block(0.1, 0.9, character(0)), 0.9)
  expect_equal(annotation_block(0.1, 0.9, NA_character_), 0.9)
})

test_that("a probability prints with its relation, as a log Bayes factor does", {
  expect_equal(format_probability(0.786), "= .79")
  expect_equal(format_probability(0.5), "= .50")
  # The ends get an inequality: a run that never saw the edge absent has not
  # established that the probability is exactly one.
  expect_equal(format_probability(1), "> .99")
  expect_equal(format_probability(0.999), "> .99")
  expect_equal(format_probability(0), "< .01")
  expect_equal(format_probability(0.001), "< .01")
  expect_equal(format_probability(NaN), "= NA")

  # The caller composes "PIP" + this, exactly as it composes "log BF" +
  # format_log_bf(), so the two evidence lines of a panel read alike.
  expect_equal(paste("PIP", format_probability(0.786)), "PIP = .79")
  expect_equal(paste("PIP", format_probability(1)), "PIP > .99")
})
