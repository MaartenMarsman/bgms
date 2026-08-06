# --------------------------------------------------------------------------- #
# The two guards zratio_constants() raises on the way to a constant set.
#
# They answer different questions and are gated differently:
#
#  - the certified-range warning is about the diagonal SHAPE. The fixed
#    quadrature grids were scored at shapes {2, 10, 12, 15, 20}; outside
#    [.zratio_constants_shape_lo, .zratio_constants_shape_hi] nothing is
#    measured. alpha = 1 is the exponential diagonal the grids are built for,
#    so only a non-unit shape can leave the band.
#
#  - the truncation warning is about the WIDTH of the c-grid, which is fixed at
#    cmax = 42 for every cell. Whether the pair integrals have decayed by then
#    depends on (delta, eta) as much as on alpha, and small eta widens the
#    diagonal at any shape. It is the only detector that the grid is wide
#    enough, and it used to sit inside the shape gate -- so the default
#    exponential path, which is what nearly every fit runs, was never checked.
# --------------------------------------------------------------------------- #

# zratio_constants serves a session cache keyed on the cell and returns before
# either guard on a hit, so a test that wants to observe a warning has to ask
# for a cold cell.
clear_zratio_constants_cache = function() {
  cache = bgms:::zratio_constants_cache
  rm(list = ls(cache, all.names = TRUE), envir = cache)
}

test_that("the truncation guard fires at the default shape alpha = 1", {
  clear_zratio_constants_cache()

  # A genuine violating cell, not an injected one: at eta = 0.05 the diagonal
  # is wide enough that the pair integrals still carry ~12% of their peak at
  # the end of the c-grid, five orders of magnitude past the 1e-6 threshold.
  # Both channels are checked, so both are pinned here.
  pair = bgms:::zratio_pair_integrals(delta = 1, sigma = 1, beta = 0.05)
  tail_g = pair$gv[length(pair$gv)] / pair$gv[1]
  tail_i = pair$ispike(max(pair$cg)) / pair$ispike(0)
  expect_gt(tail_g, 1e-6)
  expect_gt(tail_i, 1e-6)

  expect_warning(
    bgms:::zratio_constants(delta = 1, eta = 0.05, alpha = 1),
    "retain visible mass at the grid edge"
  )
})

test_that("the shape warning stays gated on a non-exponential diagonal", {
  clear_zratio_constants_cache()

  # alpha = 1 with a well-decayed grid: neither guard has anything to say.
  expect_no_warning(bgms:::zratio_constants(delta = 1, eta = 1, alpha = 1))

  clear_zratio_constants_cache()
  # The truncating cell above is alpha = 1, so it must NOT collect the
  # certified-range warning -- splitting the two guards must not have merged
  # their conditions.
  msgs = character(0)
  withCallingHandlers(
    bgms:::zratio_constants(delta = 1, eta = 0.05, alpha = 1),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("retain visible mass", msgs, fixed = TRUE)))
  expect_false(any(grepl("certified range", msgs, fixed = TRUE)))
})

test_that("a shape outside the certified range still warns about the shape", {
  skip_on_cran()
  clear_zratio_constants_cache()

  hi = bgms:::.zratio_constants_shape_hi
  msgs = character(0)
  withCallingHandlers(
    bgms:::zratio_constants(delta = 1, eta = 1, alpha = hi + 5),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("certified range", msgs, fixed = TRUE)))

  # Inside the range a non-unit shape is silent on both counts.
  clear_zratio_constants_cache()
  expect_no_warning(bgms:::zratio_constants(delta = 1, eta = 1, alpha = 3))
})
