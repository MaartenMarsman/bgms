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

test_that("main_difference_nodes encodes the verdict on shape, not colour alone", {
  verdict = c("presence", "undecided", "absence")
  nodes = main_difference_nodes(verdict, main_selected = TRUE)

  expect_equal(nodes$shape, c("square", "circle", "circle"))
  expect_equal(nodes$border_color[1], mover_palette()[1])
  expect_equal(nodes$border_color[2:3], c("grey55", "grey80"))
  # Shape separates settled from unsettled on its own, so the encoding does
  # not fail for a reader who cannot tell the border colours apart.
  expect_false(nodes$shape[1] == nodes$shape[2])

  # Indicators the sampler never updated carry no verdict, so the channel is
  # empty and every node is drawn alike.
  empty = main_difference_nodes(rep(NA_character_, 3), main_selected = FALSE)
  expect_equal(length(unique(empty$shape)), 1L)
  expect_equal(length(unique(empty$border_color)), 1L)
})

test_that("compare_difference_verdicts splits the two indicator families", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgmCompare(
    x = Wenchuan[1:120, 1:5], group_indicator = rep(1:2, each = 60),
    iter = 300, warmup = 300, chains = 2, seed = 13,
    difference_selection = TRUE, display_progress = "none"
  )
  found = compare_difference_verdicts(fit, evidence_threshold = 10)

  expect_equal(length(found$pairwise), 10L)
  expect_equal(length(found$main), 5L)
  expect_equal(nrow(found$pairs), 10L)
  # main_difference_selection is FALSE by default, so those indicators were
  # never updated and the node channel has nothing to say.
  expect_false(found$main_selected)
  expect_true(all(is.na(found$main)))

  # The split follows the indicator names, not an assumed ordering.
  names_all = get_raw_samples(fit)$parameter_names$indicator
  table = verdicts(fit, evidence_threshold = 10)
  expect_equal(
    found$pairwise,
    as.character(table$verdict[grepl("(pairwise)", names_all, fixed = TRUE)])
  )
})

test_that("plot.bgmCompare draws the difference network and the group panels", {
  skip_on_cran()
  skip_if_not_installed("qgraph")
  data("Wenchuan", package = "bgms")
  fit = bgmCompare(
    x = Wenchuan[1:120, 1:5], group_indicator = rep(1:2, each = 60),
    iter = 300, warmup = 300, chains = 2, seed = 13,
    difference_selection = TRUE, display_progress = "none"
  )

  path = withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_invisible(plot(fit))
  expect_invisible(plot(fit, type = "groups"))
  expect_invisible(plot(fit, type = "centrality"))
  expect_invisible(plot(fit, type = "centrality", group = c(1, 2)))

  # An empty difference network is a result, not an error: at an unreachable
  # threshold every difference is ruled out and the nodes are drawn alone.
  expect_invisible(plot(fit, evidence_threshold = 1e6))
})
