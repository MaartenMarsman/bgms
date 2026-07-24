# TEMPORARY diagnostic (not for merge).
#
# Question: under RcppParallel 6.0.0 (oneTBB 2022) on Windows, the second chain
# of a parallel run diverges bit-for-bit from the serial run. Is that a
# CORRECTNESS bug (chain samples the wrong posterior) or just NON-REPRODUCIBILITY
# (a different but equally valid trajectory)?
#
# Method: compare POSTERIOR SUMMARIES, not raw draws. If serial and parallel
# recover the same posterior means within Monte Carlo error -- even though the
# raw draws differ -- each chain is a valid sample and the issue is purely
# reproducibility. Reports the numbers via a deliberate diagnostic failure so
# they surface in the CI log. Windows-only; a no-op elsewhere.

test_that("serial vs parallel posterior agreement (windows correctness probe)", {
  if (Sys.info()[["sysname"]] != "Windows") {
    skip("windows-only correctness probe")
  }

  set.seed(11)
  p <- 4
  x <- matrix(sample(0:2, 120 * p, replace = TRUE), ncol = p)
  colnames(x) <- paste0("V", seq_len(p))

  it <- 4000L
  wu <- 2000L
  s <- bgm(x, iter = it, warmup = wu, chains = 2, cores = 1, seed = 7,
           display_progress = "none", verbose = FALSE)
  q <- bgm(x, iter = it, warmup = wu, chains = 2, cores = 2, seed = 7,
           display_progress = "none", verbose = FALSE)

  fmt <- function(comp) {
    ls <- s$raw_samples[[comp]]
    lq <- q$raw_samples[[comp]]
    # Pooled posterior mean across both chains, and posterior SD for scale.
    pooled_s <- colMeans(do.call(rbind, ls))
    pooled_q <- colMeans(do.call(rbind, lq))
    post_sd <- apply(do.call(rbind, ls), 2, stats::sd)
    # Chain 2 is the one that differs bit-for-bit under the regression.
    c2s <- colMeans(ls[[2]])
    c2q <- colMeans(lq[[2]])
    # Chain 1 vs chain 2 within the serial run: how far apart do two *valid*
    # independently-seeded chains land? A same-target reference scale.
    c1s <- colMeans(ls[[1]])
    sprintf(paste0("%s: d_pooled=%.4f d_chain2(serial-vs-par)=%.4f ",
                   "ref_chain1-vs-chain2=%.4f post_sd~%.3f rough_MCSE~%.4f"),
            comp,
            max(abs(pooled_s - pooled_q)),
            max(abs(c2s - c2q)),
            max(abs(c1s - c2s)),
            max(post_sd),
            max(post_sd) / sqrt(it))
  }

  bw1 <- identical(s$raw_samples$pairwise[[1]], q$raw_samples$pairwise[[1]])
  bw2 <- identical(s$raw_samples$pairwise[[2]], q$raw_samples$pairwise[[2]])

  fail(paste0(
    "[[CORRECTNESS PROBE — deliberate diagnostic fail]] ",
    "bitwise chain1_ident=", bw1, " chain2_ident=", bw2, " | ",
    fmt("pairwise"), " || ", fmt("main")
  ))
})
