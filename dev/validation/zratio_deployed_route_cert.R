# Deployed-route certification.
#
# Every accuracy figure in this validation program was scored on `logR`, the
# surface path. A fit takes `log_zratio`, the hot path. Those were not the same
# function: the engine carried its own alpha == 1 test left over from the
# alpha = 1-only migration, so a surface built and attached at a non-unit shape
# was ignored and the additive kernel served instead. Component-level scoring
# cannot certify deployment, so this certifies the route.
#
# Two arms:
#
#   1. Identity. `log_zratio` must equal `logR` to 1e-12 across shapes,
#      families and block sizes. This is a consistency assertion, not a
#      re-validation -- the surface is already scored against gold.
#   2. One end-to-end gold spot per shape, through `log_zratio`, because
#      "the composition is what I think it is" deserves one measured
#      confirmation per shape rather than an argument.
#
# Usage, from the repo root (installed build):
#   Rscript dev/validation/zratio_deployed_route_cert.R --cores=12

suppressMessages(library(bgms))
source("dev/validation/zratio_gold_bank.R")

opt = function(name, default) {
  a = grep(paste0("^--", name, "="), commandArgs(TRUE), value = TRUE)
  if(length(a) == 0L) default else sub(paste0("^--", name, "="), "", a[1])
}
cores = as.integer(opt("cores", "12"))

DELTA = 0.5 * log(12)
SHAPES = c(0.5, 1, 2, 3, 5, 10)
options(bgms.zratio_surface_cache = FALSE, bgms.correction_table_cache = FALSE,
        bgms.verbose = FALSE)

cat("=== arm 1: deployed route == surface path ===\n")
id_rows = list()
for(alpha in SHAPES) {
  zc = suppressWarnings(bgms:::zratio_constants(DELTA, 2, alpha = alpha))
  surf = bgms:::zratio_build_surfaces(zc, max_size = 44L, cores = cores)
  if(is.null(surf)) stop("no surface at shape ", alpha, call. = FALSE)
  for(fam in c("cn", "bip")) {
    for(size in c(10, 24, 42)) {
      G = gold_block(fam, size)
      r = bgms:::zratio_test_surface_eval(
        G, 1, 2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
        surf, zc$delta, zc$eta, FALSE, alpha
      )
      id_rows[[length(id_rows) + 1L]] = data.frame(
        alpha = alpha, family = fam, size = size,
        log_zratio = r$log_zratio, logR = r$logR,
        abs_diff = abs(r$log_zratio - r$logR)
      )
    }
  }
  cat(sprintf("  shape %-4g max |log_zratio - logR| over 6 blocks: %.3g\n",
              alpha,
              max(vapply(id_rows[(length(id_rows) - 5L):length(id_rows)],
                         function(d) d$abs_diff, 0.0))))
}
ident = do.call(rbind, id_rows)
worst_id = max(ident$abs_diff)
cat(sprintf("\nidentity worst %.3g over %d cells -> %s\n", worst_id, nrow(ident),
            if(worst_id <= 1e-12) "PASS" else "FAIL"))

cat("\n=== arm 2: end-to-end gold spot per shape, through log_zratio ===\n")
# eta 2, CN, size 42 -- the family and size where the additive kernel it
# replaces is worst, scored at the budget already banked for each shape.
gold_budget = function(alpha) {
  if(abs(alpha - 1) < 1e-12) list(s = 2000L, b = 500L, sd = c(11L, 23L, 37L, 51L, 67L))
  else if(alpha %in% c(0.5, 2)) list(s = 4000L, b = 500L, sd = c(11L, 23L, 37L))
  else list(s = 4000L, b = 500L, sd = c(4L, 8L, 12L))
}
spot_rows = list()
for(alpha in SHAPES) {
  zc = suppressWarnings(bgms:::zratio_constants(DELTA, 2, alpha = alpha))
  surf = bgms:::zratio_build_surfaces(zc, max_size = 80L, cores = cores)
  cell = list(delta = DELTA, eta = 2, alpha = alpha, slab = "normal")
  gb = gold_budget(alpha)
  g = gold_bank(zc, cell, "cn", 42, gb$s, gb$b, gb$sd, verbose = FALSE,
                cores = cores)
  G = gold_block("cn", 42)
  r = bgms:::zratio_test_surface_eval(
    G, 1, 2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    surf, zc$delta, zc$eta, FALSE, alpha
  )
  add = bgms:::zratio_test_eval(
    G, matrix(c(1, 2), 1, 2), zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0
  )$log_zratio[1]
  spot_rows[[length(spot_rows) + 1L]] = data.frame(
    alpha = alpha, gold = g$mean, gold_sd = g$sd, deployed = r$log_zratio,
    err_deployed = abs(r$log_zratio - g$mean), err_additive = abs(add - g$mean)
  )
  cat(sprintf("  shape %-4g gold %.6f (sd %.1g)  deployed %.6f  err %.3g   (additive err %.3g)\n",
              alpha, g$mean, g$sd, r$log_zratio, abs(r$log_zratio - g$mean),
              abs(add - g$mean)))
  utils::flush.console()
}
spots = do.call(rbind, spot_rows)

saveRDS(list(identity = ident, spots = spots),
        "dev/validation/wp5_deployed_route_cert.rds")
cat("\n=== summary ===\n")
print(spots, digits = 3, row.names = FALSE)
cat(sprintf("\nidentity worst %.3g (limit 1e-12); deployed gold error worst %.3g\n",
            worst_id, max(spots$err_deployed)))
cat("written to dev/validation/wp5_deployed_route_cert.rds\n")
