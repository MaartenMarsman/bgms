# The anchor grid is filtered to sizes at or below the build's size cap, so a
# cap that lands between two tiers used to train the hull at the lower tier and
# extrapolate everything above it (a cap of 40 kept the size-36 tier, dropped
# the size-42 tier, and left blocks of 37-38 outside the hull). These tests pin
# the cap tier that closes that gap, without paying for an anchor build.

# Largest mediating block an analysis on `cap` variables can present: the block
# is the edge's common-neighbour or bridge set, which excludes both endpoints.
reachable = function(cap) cap - 2L

test_that("the anchor grids cover every reachable block size at any cap", {
  for(cap in 4:44) {
    g = bgms:::zratio_anchor_grids(cap)
    expect_gte(max(g$cn$n), reachable(cap))
    expect_lte(max(g$cn$n), cap)
    expect_gte(max(g$bip$n), reachable(cap))
    expect_lte(max(g$bip$n), cap)
  }
})

test_that("the size-40 cap trains a hull that covers its own blocks", {
  # The regression cell: without the cap tier this grid stopped at 36.
  g = bgms:::zratio_anchor_grids(40L)
  expect_equal(max(g$cn$n), 40)
  expect_equal(sum(g$cn$n == 40), 4L) # two densities, two replicates
  expect_equal(unique(g$cn$sweeps[g$cn$n == 40]), 800L)
})

test_that("a cap that already sits on a tier gets no duplicate tier", {
  # Size 22 is a bipartite grid tier and 44 is within two of the size-42 CN
  # tier, so neither cap needs a top-up.
  expect_equal(max(bgms:::zratio_anchor_grids(44L)$cn$n), 42)
  expect_equal(sum(bgms:::zratio_anchor_grids(44L)$cn$n == 44), 0L)
  bip22 = bgms:::zratio_anchor_grids(22L)$bip
  expect_equal(max(bip22$n), 22)
  expect_equal(sum(bip22$n == 22), 8L) # four densities, two replicates
})

test_that("zratio_cap_tier leaves a grid alone when a tier is close enough", {
  jobs = expand.grid(n = c(4, 10, 18), d = c(0.8, 0.9))
  expect_identical(bgms:::zratio_cap_tier(jobs, 20L, c(0.8, 0.9)), jobs)
  expect_identical(bgms:::zratio_cap_tier(jobs, 18L, c(0.8, 0.9)), jobs)
  topped = bgms:::zratio_cap_tier(jobs, 21L, c(0.8, 0.9))
  expect_equal(max(topped$n), 21)
  expect_equal(nrow(topped), nrow(jobs) + 2L)
})
