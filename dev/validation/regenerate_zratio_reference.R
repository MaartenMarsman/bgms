# Regenerator for tests/testthat/fixtures/zratio_reference.rds, the bit-level
# reference the Z-ratio engine parity tests run against.
#
# The fixture holds four cells, and they are NOT interchangeable. The fit-time
# constant builder `zratio_constants(delta, eta)` works in the standardized
# cell, so it reproduces a stored cell only where sigma = 1 and eta = sigma *
# beta = beta; cells 2 and 3 carry sigma = 2 and were built through a different
# parameterization. Regenerating all four uniformly through the sigma = 1
# mapping silently corrupts those two (measured 94% and 97% shifts on the
# attempt that prompted this script). They are therefore frozen here, and the
# freeze is asserted byte-for-byte rather than trusted.
#
# `ihat`/`ghat` are also not the only vintage-carrying fields: `evals` stores
# engine log-ratios computed against them, so rebuilding the grids without the
# evaluations leaves the fixture internally inconsistent. Every dependent
# `evals` entry is rebuilt from the live engine in the same pass.
#
# Usage, from the repo root:
#   Rscript dev/validation/regenerate_zratio_reference.R --check   # verify only
#   Rscript dev/validation/regenerate_zratio_reference.R --write   # rewrite it
#
# --check recomputes everything and reports the deltas without touching the
# fixture, so the committed reference can be confirmed reproducible at any time.

suppressMessages(devtools::load_all(".", quiet = TRUE))

fixture_path = "tests/testthat/fixtures/zratio_reference.rds"

# Per-cell parameters, encoded rather than read back from the fixture: the
# assertion below is what catches a reordered or re-parameterized cell before
# anything is overwritten.
cell_params = data.frame(
  cell  = 1:4,
  delta = c(0, 0.5, 1.151292546, 1),
  sigma = c(1, 2, 2, 1),
  beta  = c(0.5, 0.5, 0.5, 1)
)

# The three engine variants differ only in how many trailing constants they are
# handed. The OLS correction slots are inert -- the engine reads the first six
# -- so `direct` and `clamp` reproduce `base`, which is what the parity test
# asserts. They are rebuilt anyway so a regenerated cell carries one vintage
# throughout; nothing reads them.
variants = c("base", "direct", "clamp")

TOL_STATIC = 1e-12  # addc/psi0: must not move beyond floating-point noise
tol_report = function(x) formatC(x, format = "g", digits = 4)

args = commandArgs(trailingOnly = TRUE)
do_write = "--write" %in% args
if(!do_write && !("--check" %in% args)) {
  stop("pass --check or --write", call. = FALSE)
}

fx_old = readRDS(fixture_path)
stopifnot(length(fx_old) == nrow(cell_params))

# Guard 1: the encoded parameters must describe the fixture in hand.
for(i in cell_params$cell) {
  cl = fx_old[[i]]
  ok = isTRUE(all.equal(cl$delta, cell_params$delta[i], tolerance = 1e-9)) &&
    isTRUE(all.equal(cl$sigma, cell_params$sigma[i], tolerance = 1e-12)) &&
    isTRUE(all.equal(cl$beta, cell_params$beta[i], tolerance = 1e-12))
  if(!ok) {
    stop(
      "cell ", i, " does not match the encoded parameters: fixture has ",
      sprintf("delta=%.10g sigma=%.10g beta=%.10g", cl$delta, cl$sigma, cl$beta),
      call. = FALSE
    )
  }
}

# Guard 2: which cells the standardized builder is allowed to touch is derived
# from sigma, not hardcoded by index, so the rule travels with the data.
regen = cell_params$cell[cell_params$sigma == 1]
frozen = setdiff(cell_params$cell, regen)
stopifnot(identical(regen, c(1L, 4L)), identical(frozen, c(2L, 3L)))
cat("regenerating cells", paste(regen, collapse = ", "),
    "| frozen cells", paste(frozen, collapse = ", "), "\n\n")

fx_new = fx_old
for(i in regen) {
  cl = fx_old[[i]]
  zc = bgms:::zratio_constants(cl$delta, cl$sigma * cl$beta)

  # Guard 3: the saddle nodes and weights are grid geometry, not integrals.
  # The c-grid append holds the node spacing, so these must be untouched; if
  # they move, the change is not the one this script is written for.
  stopifnot(identical(zc$tg, cl$tg), identical(zc$wt, cl$wt))

  # Guard 4: addc and psi0 do not read the saddle tail, so they must stay put
  # to floating-point noise. A real move here means the correction reached
  # further than the tail and the blast radius needs re-deriving.
  d_addc = max(abs(zc$addc - cl$addc6) / pmax(abs(cl$addc6), 1e-300))
  d_psi0 = abs(zc$psi0 - cl$psi0) / abs(cl$psi0)
  if(max(d_addc, d_psi0) > TOL_STATIC) {
    stop("cell ", i, ": addc/psi0 moved by ", tol_report(max(d_addc, d_psi0)),
         ", beyond the ", tol_report(TOL_STATIC), " noise budget", call. = FALSE)
  }

  d_ihat = max(abs(zc$ihat - cl$ihat))
  d_ghat = max(abs(zc$ghat - cl$ghat))
  cat(sprintf("cell %d (delta=%.6g eta=%g)\n", i, cl$delta, cl$sigma * cl$beta))
  cat(sprintf("  addc/psi0 rel move  %s (kept as stored)\n",
              tol_report(max(d_addc, d_psi0))))
  cat(sprintf("  ihat max abs move   %s\n", tol_report(d_ihat)))
  cat(sprintf("  ghat max abs move   %s  (%s relative to ghat[1] = %.6g)\n",
              tol_report(d_ghat), tol_report(d_ghat / abs(cl$ghat[1])),
              cl$ghat[1]))

  # addc6/psi0 are kept as stored: they are the same number to floating point,
  # and the evaluations below are rebuilt against the stored ones so the cell
  # stays self-consistent.
  cl$ihat = zc$ihat
  cl$ghat = zc$ghat

  addc13 = c(cl$addc6, cl$fc_coef, 1)
  addc23 = c(addc13, cl$hull)
  worst = 0
  for(gn in names(cl$graphs)) {
    g = cl$graphs[[gn]]
    eu = which(upper.tri(g), arr.ind = TRUE)
    edges = cbind(eu[, 1], eu[, 2])
    for(variant in variants) {
      addc = switch(variant, base = cl$addc6, direct = addc13, clamp = addc23)
      res = bgms:::zratio_test_eval(
        g, edges, addc, cl$tg, cl$ihat, cl$ghat, cl$wt, cl$psi0
      )
      key = paste(gn, variant, sep = "_")
      if(identical(variant, "base")) {
        worst = max(worst, max(abs(as.numeric(res$log_zratio) -
                                     as.numeric(cl$evals[[key]]))))
      }
      cl$evals[[key]] = as.numeric(res$log_zratio)
    }
  }
  cat(sprintf("  log Z-ratio max move %s  (over %d edges, all graphs)\n\n",
              tol_report(worst), length(cl$evals[[1]])))
  fx_new[[i]] = cl
}

# Guard 5: the frozen cells are byte-identical, checked on the serialized
# objects rather than by inspection.
frozen_ok = vapply(frozen, function(i) {
  identical(serialize(fx_old[[i]], NULL), serialize(fx_new[[i]], NULL))
}, logical(1))
if(!all(frozen_ok)) {
  stop("frozen cells changed: ", paste(frozen[!frozen_ok], collapse = ", "),
       call. = FALSE)
}
cat("frozen cells", paste(frozen, collapse = ", "), "byte-identical: OK\n")

if(!do_write) {
  cat("\n--check: fixture not written\n")
  quit(save = "no")
}

saveRDS(fx_new, fixture_path, version = 2)

# Guard 6: the freeze must survive the round trip through saveRDS, not just the
# in-memory edit.
fx_rt = readRDS(fixture_path)
rt_ok = vapply(frozen, function(i) {
  identical(serialize(fx_old[[i]], NULL), serialize(fx_rt[[i]], NULL))
}, logical(1))
if(!all(rt_ok)) {
  stop("frozen cells changed on the round trip: ",
       paste(frozen[!rt_ok], collapse = ", "), call. = FALSE)
}
cat("written to", fixture_path, "and re-read: frozen cells still identical\n")
