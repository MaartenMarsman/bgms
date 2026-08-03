# Mediation magnitude per Gamma diagonal shape, and the shape-10 irrelevance
# bound.
#
# Two claim types run in this validation program and must not be blended:
#
#   Accuracy claim   -- the surface tracks block-Gibbs gold to within X, and the
#                       reference resolves the comparison. Scored elsewhere.
#   Irrelevance claim -- the whole mediated correction is bounded by Y over the
#                       scored band, so EVERY method's error is at most Y by
#                       construction, whichever approximation runs.
#
# The second is not a weaker version of the first. A Gamma diagonal at a large
# shape concentrates the precision diagonal, the matrix goes diagonally
# dominant, and mediation dies; past that point the deployment guarantee rests
# on the correction being too small to feel rather than on the surface being
# right about it.
#
# Mediation is measured against the isolated-edge ratio log(psi0), which is what
# the engine returns for an edge with no mediating structure at all, so
#
#   mediation(block) = gold(block) - log(psi0)
#
# is the entire quantity any correction exists to supply. The bound Y is the max
# of |mediation| over the scored band, not one block.
#
# Blocks are the canonical single-component ones at density 1, which is where
# mediation is largest (every mediating path present); a sparser spot check
# confirms the direction rather than assuming it.
#
# Usage, from the repo root (installed build):
#   Rscript certification/zratio_mediation_bound.R

suppressMessages(library(bgms))
source("certification/zratio_gold_bank.R")

DELTA = 0.5 * log(12)
SHAPES = c(2, 3, 5, 10)
ETAS = c(1, 2)
SIZES = list(cn = c(20, 26, 32, 38, 42, 100), bip = c(24, 30, 36, 42, 100))
SEEDS = list("2" = c(11L, 23L, 37L), "3" = c(4L, 8L, 12L),
             "5" = c(4L, 8L, 12L), "10" = c(4L, 8L, 12L))

bank = readRDS(gold_bank_path())
key_mean = function(fam, eta, alpha, n) {
  ks = vapply(SEEDS[[as.character(alpha)]], function(sd_) {
    gold_bank_key(list(delta = DELTA, eta = eta, alpha = alpha,
                       slab = "normal"), fam, n, 4000L, 500L, sd_)
  }, character(1))
  v = vapply(ks, function(k) if(is.null(bank[[k]])) NA_real_ else bank[[k]],
             numeric(1))
  if(anyNA(v)) NA_real_ else mean(v)
}

rows = list()
for(alpha in SHAPES) {
  for(eta in ETAS) {
    zc = suppressWarnings(bgms:::zratio_constants(DELTA, eta, alpha = alpha))
    iso = log(zc$psi0)
    for(fam in names(SIZES)) {
      for(n in SIZES[[fam]]) {
        gm = key_mean(fam, eta, alpha, n)
        rows[[length(rows) + 1L]] = data.frame(
          alpha = alpha, eta = eta, family = fam, size = n,
          isolated = iso, gold = gm, mediation = gm - iso
        )
      }
    }
  }
}
med = do.call(rbind, rows)
med = med[!is.na(med$gold), ]

cat("=== mediation magnitude by shape (max |gold - isolated edge| over band) ===\n")
by_shape = do.call(rbind, lapply(split(med, med$alpha), function(d) data.frame(
  alpha = d$alpha[1], n_cells = nrow(d),
  min_med = min(d$mediation), max_med = max(d$mediation),
  max_abs = max(abs(d$mediation))
)))
rownames(by_shape) = NULL
print(by_shape, digits = 3)

cat("\n=== by shape and eta (eta dependence) ===\n")
by_eta = do.call(rbind, lapply(split(med, list(med$alpha, med$eta), drop = TRUE),
                               function(d) data.frame(
  alpha = d$alpha[1], eta = d$eta[1], n_cells = nrow(d),
  max_abs_mediation = max(abs(d$mediation)),
  at_size = d$size[which.max(abs(d$mediation))],
  family = d$family[which.max(abs(d$mediation))]
)))
rownames(by_eta) = NULL
by_eta = by_eta[order(by_eta$alpha, by_eta$eta), ]
print(by_eta, digits = 3)

cat("\n=== shape-10 detail (the irrelevance bound) ===\n")
d10 = med[med$alpha == 10, ]
print(d10[order(d10$eta, d10$family, d10$size), ], digits = 4, row.names = FALSE)
cat(sprintf("\nY = max |mediation| over the shape-10 band = %.3g nats (%d cells)\n",
            max(abs(d10$mediation)), nrow(d10)))

# Density spot check: mediation should fall as the mediating block thins, so the
# density-1 band is the worst case rather than an arbitrary slice.
cat("\n=== density spot check at shape 10 (does density 1 maximize mediation?) ===\n")
thin_cn = function(size, frac, seed) {
  G = gold_block("cn", size)
  q = size + 2L
  set.seed(seed)
  for(a in 3:(q - 1)) for(b in (a + 1):q) {
    if(G[a, b] == 1L && stats::runif(1) > frac) G[a, b] = G[b, a] = 0L
  }
  G
}
for(eta in ETAS) {
  zc = suppressWarnings(bgms:::zratio_constants(DELTA, eta, alpha = 10))
  iso = log(zc$psi0)
  for(frac in c(1, 0.5)) {
    G = if(frac == 1) gold_block("cn", 42) else thin_cn(42, frac, 31L)
    g = bgms:::zratio_test_gold_moments(
      G, 1, 2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
      zc$delta, zc$eta, 4000L, 500L, 4L, FALSE, 10
    )
    cat(sprintf("  eta %g cn k=42 clique density %.1f -> mediation %.3g\n",
                eta, frac, g$logR - iso))
  }
}

saveRDS(med, "certification/wp5_mediation.rds")
cat("\nwritten to certification/wp5_mediation.rds\n")
