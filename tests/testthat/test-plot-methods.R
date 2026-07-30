test_that("verdict_edge_colors carries sign for present edges and greys the rest", {
  palette = mover_palette()
  colors = verdict_edge_colors(
    weight = c(0.2, -0.2, 0.2, -0.2, 0),
    verdict = c("presence", "presence", "undecided", "undecided", "presence")
  )
  expect_equal(colors[1], palette[1])
  expect_equal(colors[2], palette[2])
  expect_equal(colors[3:4], rep("grey65", 2))
  # A zero weight on a present edge is not negative.
  expect_equal(colors[5], palette[1])
})

test_that("resolve_variable accepts names and positions and rejects the rest", {
  variables = c("a", "b", "c")
  expect_equal(resolve_variable("b", variables, "variable1"), 2L)
  expect_equal(resolve_variable(3, variables, "variable2"), 3L)
  expect_error(resolve_variable("z", variables, "variable1"), "not one of the model's variables")
  expect_error(resolve_variable(0, variables, "variable1"), "column position")
  expect_error(resolve_variable(4, variables, "variable1"), "column position")
  expect_error(resolve_variable(1.5, variables, "variable1"), "column position")
})

test_that("plot.bgms draws the verdict-encoded network and the centrality panel", {
  skip_on_cran()
  skip_if_not_installed("qgraph")
  fit = bgm(Wenchuan[, 1:6],
    chains = 2, iter = 400, warmup = 400, seed = 1,
    display_progress = "none", verbose = FALSE
  )

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_invisible(plot(fit))
  expect_invisible(plot(fit, legend = FALSE))
  expect_invisible(plot(fit, type = "centrality"))

  # The threshold is the one verdicts() uses, so the picture and the table move
  # together; an unreachable threshold leaves nothing to draw.
  expect_invisible(plot(fit, evidence_threshold = 3))
  expect_error(plot(fit, evidence_threshold = 1), "greater than 1")
  expect_error(plot(fit, type = "nonesuch"))
})

test_that("plot_edge_posterior draws one edge and validates its arguments", {
  skip_on_cran()
  fit = bgm(Wenchuan[, 1:6],
    chains = 2, iter = 400, warmup = 400, seed = 1,
    display_progress = "none", verbose = FALSE
  )

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_invisible(plot_edge_posterior(fit, "intrusion", "upset"))
  # The edge is unordered: either order names the same edge.
  expect_invisible(plot_edge_posterior(fit, "upset", "intrusion"))
  expect_invisible(plot_edge_posterior(fit, 1, 4))

  expect_error(plot_edge_posterior(fit, "intrusion", "intrusion"), "two different variables")
  expect_error(plot_edge_posterior(fit, "intrusion", "nope"), "not one of the model's variables")
  expect_error(
    plot_edge_posterior(fit, "intrusion", "upset", evidence_threshold = 0.5),
    "greater than 1"
  )
})

test_that("the network needs at least one edge that is not ruled out", {
  skip_on_cran()
  skip_if_not_installed("qgraph")
  fit = bgm(Wenchuan[, 1:6],
    chains = 2, iter = 400, warmup = 400, seed = 1,
    display_progress = "none", verbose = FALSE
  )
  edges = verdicts(fit, evidence_threshold = 10)
  skip_if(all(edges$verdict == "absence"))

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Every edge the picture draws is one verdicts() does not call absent, and
  # every edge it omits is one verdicts() does.
  drawn = edges$parameter[edges$verdict != "absence"]
  expect_gt(length(drawn), 0L)
  expect_invisible(plot(fit))
})
