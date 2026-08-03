# The joint precision-graph specification realizes pi(Gamma) * Z(Gamma) as the
# graph marginal, so the edge-inclusion prior a fit actually applies is not the
# nominal edge prior. These tests pin the advisory notice: which cells of the
# spec matrix it fires on, which wording it uses, and that it is silenced by
# verbose = FALSE.

notice = function(spec = "joint", model_type = "ggm", edge_selection = TRUE,
                  edge_prior = "Bernoulli", num_continuous = 5,
                  verbose = TRUE) {
  # setup.R silences bgms.verbose for the suite; the notice is advisory and
  # follows that flag, so the firing cells have to re-enable it.
  withr::local_options(bgms.verbose = verbose)
  zratio_joint_realized_prior_notice(
    precision_graph_prior = spec, model_type = model_type,
    edge_selection = edge_selection, edge_prior = edge_prior,
    num_continuous = num_continuous
  )
}

test_that("the realized-prior notice fires for every edge prior under joint", {
  for(ep in c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block")) {
    expect_message(notice(edge_prior = ep), "realized edge-inclusion prior")
  }
})

test_that("the notice separates a learned from a fixed inclusion probability", {
  expect_message(notice(edge_prior = "Bernoulli"), "nothing absorbs the tilt")
  for(ep in c("Beta-Bernoulli", "Stochastic-Block")) {
    expect_message(notice(edge_prior = ep), "hyperparameter update is corrected")
  }
})

test_that("the notice is silent outside the tilted-block cells", {
  # Hierarchical targets the nominal edge prior, so there is nothing to report.
  expect_no_message(notice(spec = "hierarchical"))
  # No edge selection: one graph, no reweighting across graphs.
  expect_no_message(notice(edge_selection = FALSE))
  # No continuous precision block to tilt (ordinal MRF), and a mixed model with
  # a single continuous variable has no continuous-continuous edge.
  expect_no_message(notice(model_type = "omrf"))
  expect_no_message(notice(model_type = "mixed_mrf", num_continuous = 1))
  expect_message(notice(model_type = "mixed_mrf", num_continuous = 2))
})

test_that("the notice follows the advisory verbose flag", {
  expect_no_message(notice(verbose = FALSE))
})

test_that("bgm wires the realized-prior notice into the joint spec", {
  skip_on_cran()
  set.seed(3)
  y = matrix(rnorm(30 * 5), 30, 5)
  fit_args = list(
    x = y, variable_type = "continuous", iter = 10, warmup = 10,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 2),
    update_method = "gibbs", chains = 1, cores = 1, seed = 1,
    display_progress = "none", verbose = TRUE
  )
  # verbose = TRUE is what makes the notice fire, and it also un-silences the
  # rest of the fit's advisory output (the correction build). Collect every
  # message and read the notice out of it, so the others neither satisfy the
  # expectation nor escape to the console.
  fit_messages = function(extra = list()) {
    msgs = character(0)
    withCallingHandlers(
      do.call(bgm, c(fit_args, extra)),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    msgs
  }
  fires = function(msgs) any(grepl("realized edge-inclusion prior", msgs, fixed = TRUE))

  # The notice is about the joint specification, which is no longer the default
  # (F-010), so the firing case names it.
  expect_true(fires(fit_messages(list(precision_graph_prior = "joint"))))
  expect_false(fires(fit_messages(list(precision_graph_prior = "hierarchical"))))
  # ... and the default, now hierarchical, does not fire it either.
  expect_false(fires(fit_messages()))
})
