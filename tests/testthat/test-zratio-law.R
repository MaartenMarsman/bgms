# Tests for the mu-first CPA analytic CN law (zratio_law_moments): a dormant
# deterministic anchor engine kept as large-q insurance and NOT wired into the
# default build (the Monte-Carlo block oracle is the anchor source; see
# zratio_law.h). Because the engine is dormant, the whole file runs in the
# nightly tier (T1, BGMS_RUN_SLOW_TESTS), except the all-MC oracle cell
# comparison, which is weekly certification (T2, BGMS_RUN_CERTIFICATION).
# Two acceptance channels:
#   (a) porting fidelity  -- reproduces the companion R eval_mu_law on a fixture
#       of certified cells to machine precision (both port identical numerics);
#   (a') physical accuracy -- certified law cells match the all-MC block oracle
#       within Monte-Carlo noise;
#   plus gate parity (uncertified cells are rejected, not silently dressed) and
#   determinism.

fixture_path = testthat::test_path("fixtures", "zratio_law_reference.rds")

test_that("law reproduces the companion eval_mu_law reference (fidelity)", {
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to exercise the dormant CN law engine"
  )
  skip_if_not(file.exists(fixture_path), "zratio_law_reference.rds not generated")
  fx = readRDS(fixture_path)
  for(cell in fx) {
    a = cell$anchors
    for(r in seq_len(nrow(a))) {
      res = zratio_law_moments(cell$eta, cell$delta, a$D[r], a$n[r], 3000L)
      ref_cert = is.finite(a$S1[r])
      # Gate parity: the port certifies exactly the cells the reference did.
      expect_equal(res$certified, ref_cert,
        label = sprintf("cert eta=%g n=%d dens=%.1f", cell$eta, a$n[r], a$dens[r])
      )
      if(ref_cert) {
        # Same algorithm from the same cold start: agreement is FP-reordering
        # only (tighter than 1e-6 relative on both moments).
        expect_equal(res$S1, a$S1[r],
          tolerance = 1e-6,
          label = sprintf("S1 eta=%g n=%d dens=%.1f", cell$eta, a$n[r], a$dens[r])
        )
        expect_equal(res$S2, a$S2[r],
          tolerance = 1e-6,
          label = sprintf("S2 eta=%g n=%d dens=%.1f", cell$eta, a$n[r], a$dens[r])
        )
      }
    }
  }
})

test_that("certified law cells match the all-MC oracle within noise", {
  skip_unless_certification()
  eta = 2
  delta = 0.5 * log(50)
  zc = zratio_constants(delta, eta, alpha = 1, slab = "normal")
  cells = data.frame(n = c(8L, 12L, 16L), dens = c(1.0, 1.0, 0.9))
  for(r in seq_len(nrow(cells))) {
    n = cells$n[r]
    dens = cells$dens[r]
    law = zratio_law_moments(eta, delta, dens * (n - 1), n, 3000L)
    expect_true(law$certified)
    # Average the oracle over a few seeds to suppress single-chain noise; the
    # law sits < 1% from the MC mean on both moments (measured ~0.4%/0.5%).
    s1 = numeric(4)
    s2 = numeric(4)
    for(k in 1:4) {
      mc = zratio_anchor_cn(n, dens, zc, 5000L, 1000L, 7000L + k)
      s1[k] = mc$S1
      s2[k] = mc$S2
    }
    expect_equal(law$S1, mean(s1),
      tolerance = 0.03,
      label = sprintf("law vs MC S1 n=%d", n)
    )
    expect_equal(law$S2, mean(s2),
      tolerance = 0.03,
      label = sprintf("law vs MC S2 n=%d", n)
    )
  }
})

test_that("the gate rejects an unconverged (hard-corner) solve", {
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to exercise the dormant CN law engine"
  )
  # eta = 1 sparse frontier: cold start does not reach the fixed point within
  # budget, so the solve self-reports a large psi_gap and the cell is dropped
  # (moments NA), routing the caller to the MC fallback rather than dressing a
  # non-converged mu.
  res = zratio_law_moments(1, 0.5 * log(50), 0.8 * (16 - 1), 16, 3000L)
  expect_false(res$certified)
  expect_true(is.na(res$S1) && is.na(res$S2))
  expect_true(abs(res$psi_gap) > 1e-4)
})

test_that("the law solve is deterministic", {
  skip_if(
    !identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    "Set BGMS_RUN_SLOW_TESTS=true to exercise the dormant CN law engine"
  )
  eta = 2
  delta = 0.5 * log(50)
  a = zratio_law_moments(eta, delta, 12, 10, 3000L)
  b = zratio_law_moments(eta, delta, 12, 10, 3000L)
  expect_identical(a$S1, b$S1)
  expect_identical(a$S2, b$S2)
  expect_identical(a$psi_gap, b$psi_gap)
})
