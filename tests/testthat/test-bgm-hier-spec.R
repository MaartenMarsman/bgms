# Tests for bgm(precision_graph_prior = "hierarchical"): eligibility validation
# and the end-to-end fit with the Z-ratio engine and trust gauge attached.
#
# Local tier keeps the eligibility errors, one gamma-shape fit smoke, and the
# gauge-off / joint-default wiring. The multi-method Cauchy sweep, the trust-
# gauge attachment, and the mixed-data breadth are calibration/breadth checks
# for the (settled) hierarchical engine and run in the BGMS_RUN_SLOW_TESTS tier.

skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run the hierarchical engine sweep"
  )
}

hier_test_data = function(q = 10, n = 40, seed = 4) {
  set.seed(seed)
  matrix(rnorm(n * q), n, q)
}

ordinal_test_data = function(q = 4, n = 50, seed = 6) {
  set.seed(seed)
  matrix(sample(0:3, n * q, replace = TRUE), n, q)
}

test_that("hierarchical rejects a slab it cannot serve", {
  # The one configuration where the choice means something and is unsupported:
  # a continuous precision block under edge selection, with a beta-prime slab.
  expect_error(
    bgm(
      x = hier_test_data(), variable_type = "continuous",
      interaction_prior = beta_prime_prior(),
      precision_graph_prior = "hierarchical",
      update_method = "gibbs", display_progress = "none", verbose = FALSE
    ),
    "normal or Cauchy"
  )
})

test_that("hierarchical is accepted where the choice is vacuous", {
  skip_on_cran()
  # These three cells used to error. They carry no between-model move, so the
  # two specifications coincide and the argument is tolerated rather than
  # rejected. The fixed-graph cell is accepted silently -- there is nothing to
  # report -- and the two without a precision block say so once.
  expect_no_error(
    bgm(
      x = hier_test_data(q = 6), variable_type = "continuous",
      iter = 50, warmup = 50,
      interaction_prior = normal_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 1, rate = 2),
      edge_selection = FALSE,
      precision_graph_prior = "hierarchical",
      update_method = "gibbs", chains = 1, cores = 1, seed = 2,
      display_progress = "none", verbose = FALSE
    )
  )
  expect_message(
    bgm(
      x = ordinal_test_data(), iter = 50, warmup = 50,
      interaction_prior = normal_prior(scale = 0.5),
      precision_graph_prior = "hierarchical",
      update_method = "adaptive-metropolis", chains = 1, cores = 1, seed = 2,
      display_progress = "none", verbose = TRUE
    ),
    "has no effect"
  )
  X = cbind(ordinal_test_data(q = 3), rnorm(50))
  colnames(X) = paste0("V", seq_len(4))
  expect_message(
    bgm(
      x = X, variable_type = c(rep("ordinal", 3), "continuous"),
      iter = 50, warmup = 50,
      interaction_prior = normal_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 1, rate = 2),
      precision_graph_prior = "hierarchical",
      update_method = "adaptive-metropolis", chains = 1, cores = 1, seed = 2,
      display_progress = "none", verbose = TRUE
    ),
    "has no effect"
  )
})

test_that("vacuity is settled before the slab", {
  skip_on_cran()
  # These two rows pin the check ordering, which is what a refactor would most
  # plausibly break. Without a precision block there is no slab of the
  # precision prior for the beta-prime error to be about, so the vacuity
  # message is the honest report; with one, the error stands.
  expect_message(
    bgm(
      x = ordinal_test_data(), iter = 50, warmup = 50,
      interaction_prior = beta_prime_prior(),
      precision_graph_prior = "hierarchical",
      update_method = "adaptive-metropolis", chains = 1, cores = 1, seed = 2,
      display_progress = "none", verbose = TRUE
    ),
    "has no effect"
  )
  expect_error(
    bgm(
      x = hier_test_data(q = 6), variable_type = "continuous",
      iter = 50, warmup = 50,
      interaction_prior = beta_prime_prior(),
      precision_graph_prior = "hierarchical",
      update_method = "gibbs", chains = 1, cores = 1, seed = 2,
      display_progress = "none", verbose = FALSE
    ),
    "normal or Cauchy"
  )
})

test_that("a fixed graph runs the joint path exactly", {
  skip_on_cran()
  # With no between-model move the joint path applies no correction either
  # (ggm_edge_prior_correction() returns NULL for edge_selection = FALSE), so
  # the two specifications reduce to the same sampler by construction. This
  # test confirms that reduction; it is not what establishes it.
  fit_args = list(
    x = hier_test_data(q = 6), variable_type = "continuous",
    iter = 100, warmup = 100,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 2),
    edge_selection = FALSE,
    update_method = "gibbs", chains = 1, cores = 1, seed = 21,
    display_progress = "none", verbose = TRUE
  )
  msgs = character(0)
  fit_h = withCallingHandlers(
    do.call(bgm, c(fit_args, list(precision_graph_prior = "hierarchical"))),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  # The surface build is skipped, not merely quiet: there are no edge moves to
  # approximate, so the fit must not pay for it.
  expect_false(any(grepl("Building normalizing-constant", msgs, fixed = TRUE)))
  expect_length(msgs, 0L)
  # The trust gauge is skipped for the same reason -- nothing to audit.
  expect_null(fit_h@zratio_diag)
  # The request is still recorded as the user made it.
  expect_equal(fit_h@arguments$precision_graph_prior, "hierarchical")

  fit_j = suppressMessages(
    do.call(bgm, c(fit_args, list(precision_graph_prior = "joint")))
  )
  expect_identical(fit_h@raw_samples$pairwise, fit_j@raw_samples$pairwise)
  expect_identical(fit_h@raw_samples$main, fit_j@raw_samples$main)
})

test_that("a vacuous hierarchical fit round-trips through a refit", {
  skip_on_cran()
  skip_unless_slow()
  # An ordinal fit records no precision_graph_prior in $arguments -- the
  # argument has no referent there, so the message is the record. The refit
  # paths read the stored spec rather than $arguments, so that absence must
  # not reach them; this runs the one refit entry point ordinal models have.
  fit = suppressMessages(bgm(
    x = ordinal_test_data(), iter = 200, warmup = 200,
    interaction_prior = normal_prior(scale = 1),
    precision_graph_prior = "hierarchical",
    update_method = "adaptive-metropolis", chains = 1, cores = 1, seed = 4,
    display_progress = "none", verbose = FALSE
  ))
  expect_null(fit@arguments$precision_graph_prior)
  expect_no_error(
    suppressMessages(prior_sensitivity_check(
      fit, anchors = c(0.63, 1.6), iter = 100, warmup = 100,
      cores = 1, seed = 2
    ))
  )
})

test_that("the vacuity message follows bgms.verbose", {
  skip_on_cran()
  withr::local_options(bgms.verbose = FALSE)
  expect_no_message(
    bgm(
      x = ordinal_test_data(), iter = 50, warmup = 50,
      interaction_prior = normal_prior(scale = 0.5),
      precision_graph_prior = "hierarchical",
      update_method = "adaptive-metropolis", chains = 1, cores = 1, seed = 2,
      display_progress = "none", verbose = FALSE
    )
  )
})

test_that("the hierarchical spec accepts a gamma-shape diagonal", {
  skip_on_cran()
  # The trust gauge is opt-in (off by default); enable it to check attachment.
  withr::local_options(bgms.zratio_gauge_sweeps = 2L)
  Y = hier_test_data(q = 8)
  fit = bgm(
    x = Y, variable_type = "continuous",
    iter = 100, warmup = 150,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 2, rate = 2),
    precision_graph_prior = "hierarchical",
    update_method = "gibbs", chains = 1, cores = 1, seed = 7,
    display_progress = "none", verbose = FALSE
  )
  s = summary(fit)
  expect_true(all(is.finite(s$pairwise$mean)))
  expect_false(is.null(fit@zratio_diag))
})

test_that("the hierarchical spec runs at the smallest analyses", {
  skip_on_cran()
  # A bipartite bridge needs 2 + 2 nodes, so at 2 or 3 variables the bipartite
  # anchor grid is empty and no surface is built. That used to abort the fit
  # ("replacement has 1 row, data has 0"); it now routes to the additive path,
  # which is exact at these block sizes. From 4 variables up both families
  # anchor and the surface is built as usual.
  withr::local_options(
    bgms.zratio_surface_cache = FALSE,
    bgms.correction_table_cache = FALSE
  )
  for(q in 2:5) {
    outcome = tryCatch(
      {
        fit = bgm(
          x = hier_test_data(q = q, n = 60), variable_type = "continuous",
          iter = 100, warmup = 150,
          interaction_prior = normal_prior(scale = 0.5),
          precision_scale_prior = gamma_prior(shape = 1, rate = 2),
          precision_graph_prior = "hierarchical",
          update_method = "gibbs", chains = 1, cores = 1, seed = 7,
          display_progress = "none", verbose = FALSE
        )
        if(all(is.finite(summary(fit)$pairwise$mean))) "fitted" else "non-finite"
      },
      error = function(e) conditionMessage(e)
    )
    expect_equal(outcome, "fitted", info = sprintf("q = %d", q))
  }
  expect_true(bgms:::zratio_anchor_grids_empty(2L))
  expect_true(bgms:::zratio_anchor_grids_empty(3L))
  for(q in c(4L, 5L)) {
    surf = bgms:::zratio_build_surfaces(
      bgms:::zratio_constants(0.5 * log(q), 2),
      max_size = q, cores = 1L
    )
    expect_false(is.null(surf$bip), label = sprintf("bip surface at q = %d", q))
  }
})

test_that("the hierarchical spec accepts a Cauchy slab on every update method", {
  skip_on_cran()
  skip_unless_slow()
  withr::local_options(bgms.zratio_gauge_sweeps = 2L)
  Y = hier_test_data(q = 8)
  for(method in c("nuts", "adaptive-metropolis", "gibbs")) {
    fit = bgm(
      x = Y, variable_type = "continuous",
      iter = 100, warmup = 150,
      interaction_prior = cauchy_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 1, rate = 2),
      precision_graph_prior = "hierarchical",
      update_method = method, chains = 1, cores = 1, seed = 7,
      display_progress = "none", verbose = FALSE
    )
    s = summary(fit)
    expect_true(all(is.finite(s$pairwise$mean)), info = method)
    expect_true(all(s$indicator$mean >= 0 & s$indicator$mean <= 1),
      info = method
    )
    expect_equal(fit@arguments$precision_graph_prior, "hierarchical", info = method)
    expect_false(is.null(fit@zratio_diag), info = method)
  }
})

test_that("bgm fits the hierarchical spec and attaches the trust gauge", {
  skip_on_cran()
  skip_unless_slow()
  withr::local_options(bgms.zratio_gauge_sweeps = 2L)
  Y = hier_test_data(q = 12)
  fit = bgm(
    x = Y, variable_type = "continuous",
    iter = 150, warmup = 250,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 2),
    edge_prior = beta_bernoulli_prior(2, 4),
    precision_graph_prior = "hierarchical",
    update_method = "gibbs", chains = 2, cores = 2, seed = 11,
    display_progress = "none", verbose = FALSE
  )
  zd = fit@zratio_diag
  expect_false(is.null(zd))
  expect_equal(nrow(zd$per_chain), 2L)
  expect_true(all(is.finite(zd$per_chain$flip_rate)))
  expect_true(all(zd$per_chain$n_ent >= 0))
  expect_false(zd$flagged)
  expect_equal(fit@arguments$precision_graph_prior, "hierarchical")
  # The joint-path hyperparameter correction must not run on this path;
  # inclusion-parameter samples come from the clean conjugate draw.
  expect_equal(length(fit@inclusion_parameter_samples), 2L)
})

test_that("the trust gauge runs by default on the hierarchical path", {
  skip_on_cran()
  # Default (no bgms.zratio_gauge_sweeps option): the post-sampling diagnostic
  # runs and its summary is attached.
  Y = hier_test_data(q = 8)
  fit_args = list(
    x = Y, variable_type = "continuous",
    iter = 100, warmup = 150,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 2),
    precision_graph_prior = "hierarchical",
    update_method = "gibbs", chains = 1, cores = 1, seed = 7,
    display_progress = "none", verbose = FALSE
  )
  fit = do.call(bgm, fit_args)
  expect_equal(fit@arguments$precision_graph_prior, "hierarchical")
  expect_false(is.null(fit@zratio_diag))
  expect_equal(nrow(fit@zratio_diag$per_chain), 1L)

  # The option remains the off switch.
  withr::local_options(bgms.zratio_gauge_sweeps = 0L)
  expect_null(do.call(bgm, fit_args)@zratio_diag)
})

test_that("a continuous fit defaults to the hierarchical spec (F-010)", {
  skip_on_cran()
  # Re-derived from the flipped default: an unnamed precision_graph_prior on
  # continuous data now resolves to "hierarchical", which means the fit records
  # that value and the trust gauge runs, where the joint default recorded
  # "joint" and left zratio_diag NULL. The joint path is still one argument
  # away, and the second half asserts it is unchanged.
  Y = hier_test_data(q = 6)
  fit_args = list(
    x = Y, variable_type = "continuous",
    iter = 100, warmup = 150,
    update_method = "gibbs", chains = 1, cores = 1, seed = 3,
    display_progress = "none", verbose = FALSE
  )
  fit = do.call(bgm, fit_args)
  expect_equal(fit@arguments$precision_graph_prior, "hierarchical")
  expect_false(is.null(fit@zratio_diag))

  joint = do.call(bgm, c(fit_args, list(precision_graph_prior = "joint")))
  expect_equal(joint@arguments$precision_graph_prior, "joint")
  expect_null(joint@zratio_diag)
})

test_that("the vacuous-spec notice reports a request, not the default", {
  # An ordinal model has no continuous precision block, so the argument has no
  # referent. Since F-010 the value arrives on every fit, so the advisory has
  # to distinguish a user's request from an inherited default -- otherwise the
  # flagship ordinal fit gains a message about an argument nobody named.
  withr::local_options(bgms.verbose = TRUE)
  expect_message(
    zratio_vacuous_spec_notice(has_precision_block = FALSE, explicit = TRUE),
    "no effect for this model"
  )
  expect_no_message(
    zratio_vacuous_spec_notice(has_precision_block = FALSE, explicit = FALSE)
  )
  expect_no_message(
    zratio_vacuous_spec_notice(has_precision_block = TRUE, explicit = TRUE)
  )
})

test_that("bgm keeps the ordinal default fit silent about the spec", {
  skip_on_cran()
  x = ordinal_test_data(q = 4)
  msgs = function(extra = list()) {
    out = character(0)
    withCallingHandlers(
      do.call(bgm, c(list(
        x = x, iter = 50, warmup = 300, chains = 1, cores = 1,
        display_progress = "none", verbose = TRUE
      ), extra)),
      message = function(m) {
        out <<- c(out, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    out
  }
  vacuous = function(m) any(grepl("no effect for this model", m, fixed = TRUE))

  expect_false(vacuous(msgs()))
  expect_true(vacuous(msgs(list(precision_graph_prior = "hierarchical"))))
})

test_that("mixed data supports the hierarchical spec on the continuous block", {
  skip_on_cran()
  skip_unless_slow()
  withr::local_options(bgms.zratio_gauge_sweeps = 2L)
  set.seed(9)
  n = 60
  X = cbind(
    matrix(sample(0:2, n * 3, replace = TRUE), n, 3),
    matrix(rnorm(n * 8), n, 8)
  )
  colnames(X) = paste0("V", seq_len(11))
  vt = c(rep("ordinal", 3), rep("continuous", 8))

  # One continuous variable: no precision block to normalize, so the argument
  # is vacuous and accepted with the has-no-effect message.
  expect_message(
    bgm(
      x = X[, 1:4], variable_type = vt[1:4],
      iter = 50, warmup = 50,
      interaction_prior = normal_prior(scale = 0.5),
      precision_scale_prior = gamma_prior(shape = 1, rate = 2),
      precision_graph_prior = "hierarchical",
      update_method = "adaptive-metropolis", chains = 1, cores = 1, seed = 2,
      display_progress = "none", verbose = TRUE
    ),
    "has no effect"
  )

  fit = bgm(
    x = X, variable_type = vt,
    iter = 120, warmup = 200,
    interaction_prior = normal_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 2),
    edge_prior = beta_bernoulli_prior(2, 4),
    precision_graph_prior = "hierarchical",
    update_method = "adaptive-metropolis", chains = 1, cores = 1, seed = 5,
    display_progress = "none", verbose = FALSE
  )
  zd = fit@zratio_diag
  expect_false(is.null(zd))
  expect_true(is.finite(zd$per_chain$flip_rate))
  expect_false(zd$flagged)
  expect_equal(fit@arguments$precision_graph_prior, "hierarchical")
})
