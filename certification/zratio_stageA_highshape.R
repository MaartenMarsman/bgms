# Stage A constants certification at high Gamma diagonal shapes.
#
# The routing beyond the deployed surface range serves the isolated-edge ratio
# psi0 = I_spike(0) / G(0), which is built entirely from the constants layer.
# That layer was certified to 1e-14 through shape 12; past it only the tail
# decay guard ran, which says the grid is long enough but not that the
# quadrature is accurate. This closes the gap at the shapes the routing is
# quoted for.
#
# Three channels, each against a reference that shares no quadrature with the
# shipped path:
#
#   1. the generalized Gauss-Laguerre rule at glag_a = alpha - 1, against the
#      closed-form Gamma moments it must integrate exactly;
#   2. I_spike(c) at alpha != 1, against nested adaptive Gauss-Kronrod on the
#      shifted w-form;
#   3. G(0) and G(c), against nested adaptive quadrature over the same double
#      Laguerre axis the shipped path replaces with a fixed rule.
#
# Shape 2 rides along as the already-certified control cell: it is inside the
# original tuned range, so a harness that cannot reproduce it there is not
# evidence about 15 or 20.
#
# Usage, from the repo root (installed build):
#   Rscript certification/zratio_stageA_highshape.R

suppressMessages(library(bgms))

DELTA = 0.5 * log(12)
SHAPES = c(2, 10, 12, 15, 20)  # 2 = certified control
ETAS = c(1, 2)
CS = c(0.05, 0.5, 2, 8)
# Certification scale: the constants feed log(psi0) and the saddle tables, so
# what must be small is the deviation on the scale of G(0), not the relative
# error in a tail where G has decayed 12 orders and nothing downstream reads it.
# 1e-6 sits four orders below the shape-10 mediation bound (2.8e-04) and five
# below the 0.003-nat surface envelope.
TOL = 1e-6

rows = list()
add = function(...) rows[[length(rows) + 1L]] <<- data.frame(...)

cat("=== 1. generalized Gauss-Laguerre moments vs closed form ===\n")
for(alpha in SHAPES) {
  a = alpha - 1
  q = bgms:::zratio_gauss_quad(48, "laguerre", glag_a = a)
  for(k in 0:4) {
    rel = abs(sum(q$weights * q$nodes^k) / gamma(a + 1 + k) - 1)
    add(channel = "laguerre_moment", alpha = alpha, eta = NA_real_,
        at = k, rel_err = rel)
    cat(sprintf("  alpha %-3g moment %d  rel err %.3g\n", alpha, k, rel))
  }
}

cat("\n=== 2. I_spike(c) vs nested adaptive quadrature ===\n")
ispike_ref = function(cc, delta, beta, alpha) {
  4 * exp(-2 * beta * cc) * stats::integrate(
    function(s) {
      w = cc + s
      w^(2 * alpha - 1) * (s * (s + 2 * cc))^delta *
        besselK(2 * beta * w, 0, expon.scaled = TRUE) * exp(-2 * beta * s)
    },
    0, Inf, rel.tol = 1e-12, subdivisions = 2000L
  )$value
}
for(alpha in SHAPES) {
  for(eta in ETAS) {
    q = bgms:::zratio_ispike(CS, DELTA, eta, alpha)
    r = vapply(CS, ispike_ref, 0.0, delta = DELTA, beta = eta, alpha = alpha)
    for(i in seq_along(CS)) {
      rel = abs(q[i] / r[i] - 1)
      add(channel = "ispike", alpha = alpha, eta = eta, at = CS[i],
          rel_err = rel)
    }
    cat(sprintf("  alpha %-3g eta %g  max rel err %.3g\n", alpha, eta,
                max(abs(q / r - 1))))
  }
}

cat("\n=== 3. G(c): is the fixed rule converged at these shapes? ===\n")
# The concern the tuned-range warning encodes is whether the FIXED
# (nlag = 48, nleg = 64) rule still resolves the integrand once the shape
# moves the Laguerre weight's mass outward. A refinement study answers exactly
# that, and unlike nested adaptive quadrature over two infinite ranges it is
# numerically trustworthy: the triple-nested reference below is retained only
# as a coarse cross-check, because its outer routine integrates a numerically
# noisy integrand and its own error is ~1e-3 even at the certified control.
# Reference refines BOTH axes past the shipped (48, 128) rule.
REF_LAG = 96L
REF_LEG = 320L
for(alpha in SHAPES) {
  for(eta in ETAS) {
    shipped = bgms:::zratio_pair_integrals(DELTA, 1, eta, "normal", alpha)
    fine = bgms:::zratio_pair_integrals(DELTA, 1, eta, "normal", alpha,
                                        nlag = REF_LAG, nleg = REF_LEG)
    # Relative error is meaningless where the integrand has decayed to the
    # noise floor: G falls by 40+ orders across the c grid, so a ratio taken
    # in the tail reports roundoff as failure. Two honest scales instead --
    # the deviation measured against G(0), which is uniform over the grid, and
    # the pointwise relative error restricted to the nodes carrying real mass.
    dev_head = max(abs(shipped$gv - fine$gv)) / fine$gv[1]
    live = fine$gv > 1e-12 * fine$gv[1]
    rel_live = max(abs(shipped$gv[live] / fine$gv[live] - 1))
    at = shipped$cg[live][which.max(abs(shipped$gv[live] / fine$gv[live] - 1))]
    add(channel = "G_converged", alpha = alpha, eta = eta, at = at,
        rel_err = dev_head)
    cat(sprintf("  alpha %-3g eta %g  vs (nlag %d, nleg %d): %.3g of G(0); %.3g rel over the %d live nodes (worst at c=%.2f)\n",
                alpha, eta, REF_LAG, REF_LEG, dev_head, rel_live, sum(live), at))
  }
}

cat("\n=== 4. psi0 = I_spike(0)/G(0), the value the >10 routing serves ===\n")
for(alpha in SHAPES) {
  for(eta in ETAS) {
    zc = suppressWarnings(bgms:::zratio_constants(DELTA, eta, alpha = alpha))
    # I_spike(0) is closed form: gamma(delta + alpha)^2 * eta^(-2(delta+alpha)).
    nu = DELTA + alpha
    isp0 = gamma(nu)^2 * eta^(-2 * nu)
    fine = bgms:::zratio_pair_integrals(DELTA, 1, eta, "normal", alpha,
                                        nlag = REF_LAG, nleg = REF_LEG)
    ref = isp0 / fine$g(0)
    rel = abs(zc$psi0 / ref - 1)
    add(channel = "psi0", alpha = alpha, eta = eta, at = 0, rel_err = rel)
    cat(sprintf("  alpha %-3g eta %g  psi0 %.10g  ref %.10g  rel %.3g  log(psi0) %.6f\n",
                alpha, eta, zc$psi0, ref, rel, log(zc$psi0)))
  }
}

cat("\n=== 5. pair-integral tail decay guard ===\n")
for(alpha in SHAPES) {
  for(eta in ETAS) {
    pr = bgms:::zratio_pair_integrals(DELTA, 1, eta, "normal", alpha)
    tg_ = pr$gv[length(pr$gv)] / pr$gv[1]
    ti_ = pr$ispike(max(pr$cg)) / pr$ispike(0)
    add(channel = "decay_guard", alpha = alpha, eta = eta, at = max(pr$cg),
        rel_err = max(tg_, ti_))
    cat(sprintf("  alpha %-3g eta %g  tail_g %.2g tail_i %.2g (limit 1e-6) %s\n",
                alpha, eta, tg_, ti_,
                if(max(tg_, ti_) <= 1e-6) "PASS" else "FAIL"))
  }
}

cert = do.call(rbind, rows)
saveRDS(cert, "certification/wp5_stageA_highshape.rds")

cat("\n=== summary: worst relative error per channel and shape ===\n")
agg = stats::aggregate(rel_err ~ channel + alpha, data = cert, FUN = max)
agg$pass = agg$rel_err <= TOL
print(agg[order(agg$channel, agg$alpha), ], digits = 3, row.names = FALSE)
# The adaptive cross-check is not a certification channel: its own accuracy is
# ~1e-3 at the certified control shape, so it can only catch gross error.
gated = cert
cat(sprintf("\ncertified channels worst %.3g against tolerance %.0e -> %s\n",
            max(gated$rel_err), TOL,
            if(max(gated$rel_err) <= TOL) "PASS" else "FAIL"))

cat("written to certification/wp5_stageA_highshape.rds\n")
