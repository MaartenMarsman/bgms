# TEMPORARY investigation probe (not for merge).
#
# Intermittent Windows CI failures of
# 'test-regressions-2.R:355' ("serial and parallel dispatch produce identical
# draws") reproduced on immediate retry. This probe repeats the exact
# comparison many times in one job to measure the mismatch rate and localise
# the first divergence, so we can tell a rare flake from a near-deterministic
# regression. Windows-only; a no-op elsewhere.

test_that("serial/parallel draw-identity mismatch rate (windows probe)", {
  if (Sys.info()[["sysname"]] != "Windows") {
    skip("windows-only determinism stress probe")
  }

  reps = 25L
  mism = 0L
  first_detail = ""

  for (r in seq_len(reps)) {
    set.seed(11)
    p = 4
    x = matrix(sample(0:2, 120 * p, replace = TRUE), ncol = p)
    colnames(x) = paste0("V", seq_len(p))

    serial = bgm(
      x, iter = 100, warmup = 100, chains = 2, cores = 1, seed = 7,
      display_progress = "none", verbose = FALSE
    )
    parallel = bgm(
      x, iter = 100, warmup = 100, chains = 2, cores = 2, seed = 7,
      display_progress = "none", verbose = FALSE
    )

    same_pw = identical(serial$raw_samples$pairwise, parallel$raw_samples$pairwise)
    same_main = identical(serial$raw_samples$main, parallel$raw_samples$main)
    if (!same_pw || !same_main) {
      mism = mism + 1L
      if (!nzchar(first_detail)) {
        which_chain = NA_integer_
        for (c in seq_along(serial$raw_samples$pairwise)) {
          if (!identical(serial$raw_samples$pairwise[[c]],
                         parallel$raw_samples$pairwise[[c]])) {
            which_chain = c
            break
          }
        }
        first_detail = sprintf(
          "first at rep %d (pairwise_ok=%s main_ok=%s, chain %s)",
          r, same_pw, same_main, which_chain
        )
      }
    }
  }

  expect_equal(
    mism, 0L,
    info = sprintf("[[PROBE]] mismatches %d/%d; %s", mism, reps, first_detail)
  )
})
