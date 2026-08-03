# Test-suite tiers.
#
# The suite runs at three cadences. Which tier a block belongs to is decided by
# the CLASS OF BUG it catches, not by how long it takes
# (dev/plans/backlog/2026-07-29_test-suite-tiering-audit_AUDIT.md).
#
#   T0  every run (local devtools::test(), and every push/PR)
#       Contracts, validation, wiring, the product surface, and the tight
#       numerical unit guards that catch a shifted constant in one run.
#       No gate. Budget ~90 s.
#
#       CRAN runs a SUBSET of T0: the heavy internal numerical guard files
#       are excluded there via skip_heavy_guard_on_cran() below (maintainer
#       ruling, 2026-08-04). CRAN's Windows incoming pretest measured T0 at
#       506 s (~6x this machine) and the whole check at 13 min against
#       CRAN's 10-min ceiling. The guards catch regressions WE introduce,
#       not platform breakage, and run on every push anyway; what CRAN
#       checking is for -- "does the package work on this platform" -- stays:
#       the product-surface smokes (test-bgm.R, test-bgmCompare.R, methods,
#       extractors, validation) still run on CRAN in full.
#
#   T1  nightly heartbeat (daily, develop)   BGMS_RUN_SLOW_TESTS=true
#       "Does the settled math still hold tonight": graph-law and prior-chain
#       identities, single surface-vs-gold cells, the RB saturation boundary,
#       concordance smokes. Budget <= 60 min on the 2-core runner. Gated by
#       each file's own skip_unless_slow()-style helper.
#
#   T2  weekly certification (Sunday, develop)  BGMS_RUN_CERTIFICATION=true
#       The Monte-Carlo machinery: the trust gauge's harm controls (moved
#       here from the T1 charter by maintainer ruling, 2026-08-03, F-103 —
#       their statistic is a per-fit Monte-Carlo quantity with real
#       seed-to-seed dispersion, so they are distributional assertions
#       whatever their runtime), SBC suites, full parameter-recovery
#       sweeps, full NUTS-vs-MH condition grids, n = 2e6 MC channels, and the
#       refit cross-validations. Gated by skip_unless_certification() below.
#
# The weekly workflow sets BOTH variables, so a T2 run also carries T1. The
# nightly workflow sets BGMS_RUN_SLOW_TESTS only, so T2 blocks skip there and
# say so.

# File-level gate for the heavy internal numerical guard files (see the T0
# charter note above for the ruling and the rationale). Called at the top
# level of: test-rb-inclusion-probabilities.R, test-prior-interface.R,
# test-zratio-isolated-edge-routing.R, test-zratio-engine.R.
skip_heavy_guard_on_cran = function() {
  skip_on_cran()
}

skip_unless_certification = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_CERTIFICATION"), "true"),
    message = paste0(
      "Weekly certification tier (T2): set BGMS_RUN_CERTIFICATION=true to run. ",
      "BGMS_RUN_SLOW_TESTS=true alone is the nightly heartbeat (T1) and does ",
      "not enable this block."
    )
  )
}
