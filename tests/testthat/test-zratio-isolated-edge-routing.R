# Isolated-edge routing past the surface's validated shape range.
#
# Beyond .zratio_surface_shape_hi the mediating correction is switched off and
# every edge is served log(psi0), the exact ratio for an edge with no mediating
# structure. The alternative -- the additive saddle -- returns zero on
# common-neighbour blocks and discards the whole ratio, so this is a strict
# improvement, bounded by the measured mediation (2.84e-04 nats; see
# zratio_mediation_off).
#
# CI-only on CRAN's clock (skip_heavy_guard_on_cran below; see helper-tiers.R).
#
# The flag crosses four layers: the R spec, the C++ spec reader, the engine's
# hot path, and the chain runner's counter export. A flag set in one and ignored
# in another is exactly how the earlier deploy gate failed, so the probes here
# go through the DEPLOYED entry points -- zratio_test_spec_eval takes the same
# spec list the samplers take and builds the engine through the same reader, and
# the fit probe takes a real chain -- rather than through a component API.

skip_heavy_guard_on_cran()

DELTA = 0.5 * log(12)

# One constants build per cell, shared across the file (~5 s each otherwise).
zc_hi = suppressWarnings(bgms:::zratio_constants(DELTA, 2, alpha = 12))
zc_one = bgms:::zratio_constants(DELTA, 2, alpha = 1)

# CN block: `size` common neighbours of edge (1, 2), complete among themselves.
cn_block = function(size) {
  q = size + 2L
  G = matrix(0L, q, q)
  for(v in 3:q) G[1, v] = G[v, 1] = G[2, v] = G[v, 2] = 1L
  for(a in 3:(q - 1)) for(b in (a + 1):q) G[a, b] = G[b, a] = 1L
  G
}

# Bipartite bridge block: `size` nodes split between the exclusive neighbour
# sets, complete across the split.
bip_block = function(size) {
  q = size + 2L
  na_ = size %/% 2L
  A = 2L + seq_len(na_)
  B = setdiff(3:q, A)
  G = matrix(0L, q, q)
  for(v in A) G[1, v] = G[v, 1] = 1L
  for(v in B) G[2, v] = G[v, 2] = 1L
  for(a in A) for(b in B) G[a, b] = G[b, a] = 1L
  G
}

test_that("the deployed route serves the isolated-edge value past the shape range", {
  # THE deployment probe. It goes through the spec the sampler is handed, so it
  # fails if R sets mediation_off and any layer below drops it -- which is the
  # failure this routing was most at risk of.
  spec = bgms:::zratio_spec_list(zc_hi, gauge_sweeps = 0L)
  spec = bgms:::zratio_attach_surface(spec, zc_hi, size = 40, cores = 1L)
  expect_true(spec$mediation_off)
  expect_null(spec$surface)

  iso = log(zc_hi$psi0)
  for(fam in c("cn", "bip")) {
    for(size in c(4, 12, 30, 60)) {
      G = if(fam == "cn") cn_block(size) else bip_block(size)
      r = bgms:::zratio_test_spec_eval(spec, G, matrix(c(1L, 2L), 1, 2))
      # Equality with the isolated-edge value, not mere difference from the
      # additive one: a route that merely differs is not the route that was
      # designed, and on a CN block the additive value happens to be 0.
      expect_identical(
        r$log_zratio[1], iso,
        label = sprintf("deployed value, %s block of %d", fam, size)
      )
      expect_true(r$mediation_off)
      expect_false(r$has_surface)
      expect_equal(r$n_isolated, 1)
      expect_equal(r$n_pred, 0)
      expect_equal(r$n_add, 0)
    }
  }
})

test_that("the routing is the flag, not the shape: forcing it off restores the additive path", {
  # Disconfirming control. If the engine re-derived the policy from alpha (the
  # defect the deploy gate carried), clearing the flag would change nothing.
  spec = bgms:::zratio_spec_list(zc_hi, gauge_sweeps = 0L)
  spec$mediation_off = FALSE
  G = cn_block(12)
  r = bgms:::zratio_test_spec_eval(spec, G, matrix(c(1L, 2L), 1, 2))
  expect_false(r$mediation_off)
  expect_equal(r$n_isolated, 0)
  expect_equal(r$n_add, 1)
  expect_false(isTRUE(all.equal(r$log_zratio[1], log(zc_hi$psi0))))
})

test_that("shapes inside the validated range keep the surface route", {
  # Pinning: the routing must not reach any shape the surface still serves.
  for(alpha in c(0.5, 1, 2, 5, 10)) {
    expect_false(
      bgms:::zratio_mediation_off(list(alpha = alpha, eta = 2)),
      label = sprintf("mediation_off at shape %g", alpha)
    )
  }
  # Below the range the additive path is unchanged: no surface, no routing.
  expect_false(bgms:::zratio_mediation_off(list(alpha = 0.25, eta = 2)))
  expect_true(bgms:::zratio_mediation_off(list(alpha = 10.5, eta = 2)))
})

test_that("the alpha = 1 spec route is bit-identical to the surface path", {
  # A synthetic surface (no anchor build) attached to a real alpha = 1 spec:
  # the spec route must reproduce the surface entry point exactly, so adding
  # the flag left the deployed alpha = 1 numerics untouched.
  mkfam = function(seed) {
    set.seed(seed)
    list(
      c1 = round(rnorm(9) * 0.1, 4),
      c2 = c(-4, round(rnorm(8) * 0.02, 4)),
      size_lo = 3, size_hi = 40, dens_lo = 0.1, dens_hi = 1.0,
      l1_lo = -10, l1_hi = 10, l2_lo = -10, l2_hi = 10, size_min = 3
    )
  }
  surf = list(cn = mkfam(11), bip = mkfam(22))
  spec = bgms:::zratio_spec_list(zc_one, gauge_sweeps = 0L)
  spec$surface = surf
  expect_false(spec$mediation_off)

  for(size in c(6, 20)) {
    G = cn_block(size)
    r_spec = bgms:::zratio_test_spec_eval(spec, G, matrix(c(1L, 2L), 1, 2))
    r_comp = bgms:::zratio_test_surface_eval(
      G, 1, 2, zc_one$addc, zc_one$tg, zc_one$ihat, zc_one$ghat, zc_one$wt,
      zc_one$psi0, surf, zc_one$delta, zc_one$eta, FALSE, 1
    )
    expect_identical(r_spec$log_zratio[1], r_comp$log_zratio)
    expect_true(r_spec$has_surface)
    expect_false(r_spec$mediation_off)
    expect_equal(r_spec$n_isolated, 0)
  }
})

test_that("the eta > 2 notice fires past the measured range and nowhere inside it", {
  # zratio_surface_fence_message reads only (alpha, eta), so this needs no
  # constants build. Collect every message the call emits, so an assertion that
  # a wording is ABSENT is made against the whole output rather than against
  # whichever message testthat happened to catch first.
  said = function(alpha, eta) {
    out = character(0)
    withCallingHandlers(
      bgms:::zratio_surface_fence_message(list(alpha = alpha, eta = eta)),
      message = function(m) {
        out <<- c(out, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    paste(out, collapse = " ")
  }
  covered = "measured bound does not cover"

  expect_match(said(12, 3), "isolated-edge ratio")
  expect_match(said(12, 3), covered)
  expect_match(said(12, 2), "isolated-edge ratio")
  expect_false(grepl(covered, said(12, 2)))
  expect_false(grepl(covered, said(12, 1)))
  # Inside the range a build failure still reports itself as a build failure,
  # and below the range the additive wording is unchanged and carries no bound.
  expect_match(said(2, 5), "surface build failed")
  expect_match(said(0.25, 5), "additive path")
  expect_false(grepl(covered, said(0.25, 5)))
  expect_false(grepl("isolated-edge", said(0.25, 5)))
})

test_that("the post-fit notice fires on the counters, and scopes its bound by eta", {
  # Counter-driven, so it reports what the chains did rather than what the spec
  # intended. The spec-build message is the other half and fires before the
  # chains launch; this one reaches a fit run with verbose = FALSE.
  said = function(chains, eta) {
    out = character(0)
    withCallingHandlers(
      bgms:::zratio_isolated_route_notice(chains, eta),
      message = function(m) {
        out <<- c(out, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    paste(out, collapse = " ")
  }
  routed = list(list(zratio = list(counters = list(
    n_isolated = 1200, n_pred = 0, n_add = 0
  ))))
  covered = "measured bound does not cover"

  expect_match(said(routed, 2), "isolated-edge normalizer ratio")
  expect_match(said(routed, 2), "1,200 edge evaluations")
  expect_match(said(routed, 2), "0.00028 nats")
  expect_false(grepl(covered, said(routed, 2)))
  expect_match(said(routed, 4), covered)
  expect_false(grepl("0.00028 nats", said(routed, 4)))
  # A non-finite rate cannot be shown to be covered, so it is not claimed to be.
  expect_match(said(routed, NA_real_), covered)

  # Silent where the route was not taken, and on chain output that predates the
  # counter: a missing tally is not a tally of zero to be reported either way.
  not_routed = list(list(zratio = list(counters = list(
    n_isolated = 0, n_pred = 5000, n_add = 0
  ))))
  legacy = list(list(zratio = list(counters = list(n_pred = 5000))))
  expect_identical(said(not_routed, 2), "")
  expect_identical(said(legacy, 2), "")
  expect_identical(said(list(list()), 2), "")
})

test_that("the collapse counter counts exactly the guard hits, and nothing else", {
  # The additive kernel's zero-collapse, in the cell it is actually deployed in
  # (a Gamma diagonal shape below .zratio_surface_shape_lo). The contract is the
  # COUNTER, not the discarded value: asserting that log_zratio is 0 there would
  # pin a known defect as the specification, which is how the deploy-gate bug
  # survived review. So the value is only ever used as the reference the counter
  # is checked against, never as an expectation of its own.
  zc = suppressWarnings(bgms:::zratio_constants(DELTA, 1, alpha = 0.25))
  spec = bgms:::zratio_attach_surface(
    bgms:::zratio_spec_list(zc, gauge_sweeps = 0L), zc, size = 40, cores = 1L
  )
  # The cell must actually be the additive one, or this tests nothing.
  expect_false(isTRUE(spec$mediation_off))
  expect_null(spec$surface)

  collapsed = 0L
  for(size in 3:30) {
    r = bgms:::zratio_test_spec_eval(spec, cn_block(size), matrix(c(1L, 2L), 1, 2))
    hit = r$log_zratio[1] == 0
    collapsed = collapsed + as.integer(hit)
    expect_equal(
      r$n_collapsed, as.integer(hit),
      label = sprintf("collapse tally at CN block %d", size)
    )
    expect_equal(unname(r$n_add), 1)
    expect_equal(unname(r$n_isolated), 0)
    expect_equal(unname(r$n_pred), 0)
  }
  # The sweep must straddle the boundary, or the equality above is vacuous:
  # a run of all-zero tallies would pass it just as well.
  expect_gt(collapsed, 0)
  expect_lt(collapsed, 28)

  # The reported size is the common-neighbour count, the scale the documented
  # boundary is quoted in.
  big = bgms:::zratio_test_spec_eval(spec, cn_block(30), matrix(c(1L, 2L), 1, 2))
  expect_equal(unname(big$max_collapse_size), 30)
})

test_that("the surface and isolated-edge routes cannot contaminate the tally", {
  # The collapse belongs to the additive kernel alone. A cell served by the
  # surface, and a cell served by the isolated-edge route, must both leave the
  # counter at zero however large the block.
  mkfam = function(seed) {
    set.seed(seed)
    list(
      c1 = round(rnorm(9) * 0.1, 4),
      c2 = c(-4, round(rnorm(8) * 0.02, 4)),
      size_lo = 3, size_hi = 40, dens_lo = 0.1, dens_hi = 1.0,
      l1_lo = -10, l1_hi = 10, l2_lo = -10, l2_hi = 10, size_min = 3
    )
  }
  surfaced = bgms:::zratio_spec_list(zc_one, gauge_sweeps = 0L)
  surfaced$surface = list(cn = mkfam(11), bip = mkfam(22))
  r = bgms:::zratio_test_spec_eval(surfaced, cn_block(30), matrix(c(1L, 2L), 1, 2))
  expect_equal(unname(r$n_collapsed), 0)
  expect_equal(unname(r$max_collapse_size), 0)

  routed = bgms:::zratio_spec_list(zc_hi, gauge_sweeps = 0L)
  r2 = bgms:::zratio_test_spec_eval(routed, cn_block(30), matrix(c(1L, 2L), 1, 2))
  expect_true(r2$mediation_off)
  expect_equal(unname(r2$n_collapsed), 0)
})

test_that("the collapse notice fires on the counters and is silent otherwise", {
  said = function(chains) {
    out = character(0)
    withCallingHandlers(
      bgms:::zratio_collapse_notice(chains),
      message = function(m) {
        out <<- c(out, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    paste(out, collapse = " ")
  }
  hit = list(list(zratio = list(counters = list(
    n_collapsed = 1400, max_collapse_size = 26,
    n_collapsed_retained = 900, max_collapse_size_retained = 21, n_add = 5000
  ))))
  txt = said(hit)
  expect_match(txt, "900 edge evaluations in the retained sweeps")
  expect_match(txt, "21 variables")
  expect_match(txt, "1,400 such evaluations over the whole run")
  expect_match(txt, "discarded rather than approximated")
  expect_match(txt, "12 and 32 common-neighbour variables")
  expect_match(txt, "summarize_zratio_gauge")

  # Warmup-only: the sampler starts from a complete graph, so an oversized
  # block early on says nothing about the posterior. Reporting that as though
  # it described the stored draws would be a false alarm on every fit past the
  # boundary, which is why the tally is split by phase at the engine.
  warm = list(list(zratio = list(counters = list(
    n_collapsed = 1400, max_collapse_size = 26,
    n_collapsed_retained = 0, max_collapse_size_retained = 0, n_add = 5000
  ))))
  wtxt = said(warm)
  expect_match(wtxt, "during warmup only")
  expect_match(wtxt, "stored draws are unaffected")
  expect_match(wtxt, "discarded rather than approximated")
  expect_false(grepl("retained sweeps the hierarchical", wtxt))

  # Silent where nothing collapsed, and on chain output that predates the
  # counter -- a missing tally is not a tally of zero to be reported.
  expect_identical(said(list(list(zratio = list(counters = list(
    n_collapsed = 0, n_add = 5000
  ))))), "")
  expect_identical(said(list(list(zratio = list(counters = list(n_add = 5000))))), "")
  expect_identical(said(list(list())), "")
})

test_that("a fit past the shape range reports the isolated-edge counter", {
  # The chain-runner layer: the counter is filled in the model and named in
  # chain_runner, and nothing else spans the two. A real chain is the only
  # thing that crosses that boundary.
  set.seed(11)
  q = 8L
  y = matrix(rnorm(60 * q), 60, q)
  spec = suppressWarnings(bgms:::bgm_spec(
    x = y, model_type = "ggm", variable_type = "continuous",
    interaction_prior_type = "normal", pairwise_scale = 0.5,
    scale_prior_type = "gamma", scale_shape = 12, scale_rate = NA_real_,
    scale_eta = 2, precision_graph_prior = "hierarchical",
    update_method = "gibbs", iter = 60, warmup = 60, chains = 1, cores = 1,
    seed = 3, display_progress = "none", verbose = FALSE
  ))
  raw = suppressWarnings(bgms:::run_sampler(spec))
  ct = raw[[1]]$zratio$counters
  expect_true("n_isolated" %in% names(ct))
  expect_gt(ct[["n_isolated"]], 0)
  # Every evaluation took the route; neither of the other two branches ran.
  expect_equal(unname(ct[["n_pred"]]), 0)
  expect_equal(unname(ct[["n_add"]]), 0)

  # ... and the counter a real chain produced drives the notice, which closes
  # the loop from the spec through the sampler to the user-facing output. The
  # fit above ran with verbose = FALSE, so the spec-build message never fired
  # and this is the only thing that would tell its user anything.
  out = character(0)
  withCallingHandlers(
    bgms:::zratio_isolated_route_notice(raw, eta = 2),
    message = function(m) {
      out <<- c(out, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_match(paste(out, collapse = " "), "isolated-edge normalizer ratio")
})
