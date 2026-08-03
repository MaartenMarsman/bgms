# Deployed-route certification for the isolated-edge routing above shape 10.
#
# The routing switches the mediating correction off and serves log(psi0). The
# claim attached to it is an irrelevance bound, not an accuracy figure: the
# whole mediated correction is small enough at these shapes that any method
# returning the isolated-edge value is wrong by at most that. Certifying it
# therefore needs two different things measured, and one of them is a route,
# not a number.
#
#   Arm 1  DEPLOYMENT. The value a fit receives equals log(psi0) exactly. This
#          goes through zratio_test_spec_eval, which builds the engine from the
#          same spec list the samplers consume, through the same reader. A
#          component entry point cannot certify this: the previous deploy gate
#          was wrong for an entire validation programme precisely because every
#          score went through a component API while fits took another branch.
#   Arm 2  MAGNITUDE. The deployed value against block-Gibbs gold, alongside
#          what the additive kernel would have returned on the same block. The
#          deployed error IS the mediation, so this measures the bound rather
#          than assuming it.
#   Arm 3  ROUTE COUNTERS. n_isolated accounts for every evaluation and neither
#          the surface nor the additive branch ran.
#
# eta 2 is the bound-setting end (mediation grows with eta, measured ~3.3x
# from eta 1 to eta 2 at shape 10), so arm 2 runs there.
#
# Usage, from the repo root (installed build):
#   Rscript certification/zratio_isolated_route_cert.R --cores=12

suppressMessages(library(bgms))
source("certification/zratio_gold_bank.R")

opt = function(name, default) {
  a = grep(paste0("^--", name, "="), commandArgs(TRUE), value = TRUE)
  if(length(a) == 0L) default else sub(paste0("^--", name, "="), "", a[1])
}
cores = as.integer(opt("cores", "12"))

DELTA = 0.5 * log(12)
SHAPES = c(12, 15, 20)
options(bgms.zratio_surface_cache = FALSE, bgms.correction_table_cache = FALSE,
        bgms.verbose = FALSE)

zc_of = function(alpha, eta) {
  suppressWarnings(bgms:::zratio_constants(DELTA, eta, alpha = alpha))
}
spec_of = function(zc) {
  # The spec a fit hands the sampler, built by the same two functions the
  # sampler paths call. attach_surface returns the spec unchanged past the
  # shape range, which is itself part of what is being certified.
  bgms:::zratio_attach_surface(
    bgms:::zratio_spec_list(zc, gauge_sweeps = 0L),
    zc, size = 80, cores = 1L
  )
}

cat("=== arm 1: the deployed value is the isolated-edge value ===\n")
id_rows = list()
for(alpha in SHAPES) {
  for(eta in c(1, 2)) {
    zc = zc_of(alpha, eta)
    spec = spec_of(zc)
    stopifnot(isTRUE(spec$mediation_off), is.null(spec$surface))
    iso = log(zc$psi0)
    for(fam in c("cn", "bip")) {
      for(size in c(10, 24, 42, 100)) {
        G = gold_block(fam, size)
        r = bgms:::zratio_test_spec_eval(spec, G, matrix(c(1L, 2L), 1, 2))
        id_rows[[length(id_rows) + 1L]] = data.frame(
          alpha = alpha, eta = eta, family = fam, size = size,
          deployed = r$log_zratio[1], isolated = iso,
          abs_diff = abs(r$log_zratio[1] - iso),
          n_isolated = r$n_isolated, n_pred = r$n_pred, n_add = r$n_add
        )
      }
    }
    cat(sprintf("  shape %-3g eta %g  max |deployed - log(psi0)| over 8 blocks: %.3g\n",
                alpha, eta,
                max(vapply(id_rows[(length(id_rows) - 7L):length(id_rows)],
                           function(d) d$abs_diff, 0.0))))
  }
}
ident = do.call(rbind, id_rows)
cat(sprintf("\nidentity worst %.3g over %d cells -> %s\n", max(ident$abs_diff),
            nrow(ident), if(max(ident$abs_diff) == 0) "PASS (exact)" else "FAIL"))

cat("\n=== arm 3: route counters ===\n")
route_ok = all(ident$n_isolated == 1) && all(ident$n_pred == 0) &&
  all(ident$n_add == 0)
cat(sprintf("  every evaluation isolated, none predicted, none additive: %s\n",
            if(route_ok) "PASS" else "FAIL"))

cat("\n=== arm 2: the deployed value against block-Gibbs gold (eta 2) ===\n")
# Shape 12 at k = 42 and k = 100 is already banked; 15 and 20 are new cells and
# are computed and banked here at the same budget.
GOLD_SWEEPS = 4000L
GOLD_BURN = 500L
GOLD_SEEDS = c(4L, 8L, 12L)
SIZES = list("12" = c(42, 100), "15" = 42, "20" = 42)
spot_rows = list()
for(alpha in SHAPES) {
  zc = zc_of(alpha, 2)
  spec = spec_of(zc)
  cell = list(delta = DELTA, eta = 2, alpha = alpha, slab = "normal")
  iso = log(zc$psi0)
  for(fam in c("cn", "bip")) {
    g = gold_bank(zc, cell, fam, SIZES[[as.character(alpha)]], GOLD_SWEEPS,
                  GOLD_BURN, GOLD_SEEDS, verbose = TRUE, cores = cores)
    for(r in seq_len(nrow(g))) {
      size = g$size[r]
      G = gold_block(fam, size)
      dep = bgms:::zratio_test_spec_eval(
        spec, G, matrix(c(1L, 2L), 1, 2)
      )$log_zratio[1]
      add = bgms:::zratio_test_eval(
        G, matrix(c(1L, 2L), 1, 2), zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt,
        zc$psi0
      )$log_zratio[1]
      spot_rows[[length(spot_rows) + 1L]] = data.frame(
        alpha = alpha, family = fam, size = size,
        gold = g$mean[r], gold_sd = g$sd[r], isolated = iso, deployed = dep,
        err_deployed = abs(dep - g$mean[r]), additive = add,
        err_additive = abs(add - g$mean[r])
      )
    }
  }
  utils::flush.console()
}
spots = do.call(rbind, spot_rows)
print(spots[, c("alpha", "family", "size", "gold", "gold_sd", "deployed",
                "err_deployed", "additive", "err_additive")],
      digits = 4, row.names = FALSE)

BOUND = 2.84e-4
cat(sprintf("\nworst deployed error %.3g against the recorded shape-10 bound %.3g -> %s\n",
            max(spots$err_deployed), BOUND,
            if(max(spots$err_deployed) <= BOUND) "inside" else "OUTSIDE"))
cat(sprintf("worst additive error  %.3g (%.0fx the deployed one)\n",
            max(spots$err_additive),
            max(spots$err_additive) / max(spots$err_deployed)))
# Does the reference move? An irrelevance bound is only readable if the gold it
# is measured against is resolved at that scale.
cat(sprintf("gold sd worst %.3g; deployed error / gold sd worst %.1f\n",
            max(spots$gold_sd), max(spots$err_deployed / spots$gold_sd)))

saveRDS(list(identity = ident, spots = spots),
        "certification/zratio_isolated_route_cert.rds")
cat("\nwritten to certification/zratio_isolated_route_cert.rds\n")
