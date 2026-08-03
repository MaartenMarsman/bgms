# Where the additive Z-ratio kernel collapses to zero, and why.
#
# The additive kernel returns exactly 0 on common-neighbour mediating blocks --
# it discards the whole log-ratio rather than approximating it coarsely. That
# was first read as a large-shape effect, then measured at the exponential shape
# too, so the "progressive in shape" framing was retracted. This maps the actual
# dependence (family, block size, density, eta, shape) and tests one mechanism.
#
# The additive path (ZRatioEngine::log_zratio) accumulates
#
#   S1 = ncn * addc[1] + cne * addc[3] + bre * addc[5]
#   S2 = ncn * addc[2] + cne * addc[4] + bre * addc[6]      (R indexing)
#
# and hands (S1, S2) to saddle_ratio, whose first line is
#
#   if (s1 <= 0 || s2 <= 0) return 1.0;
#
# i.e. a guard returning log-ratio 0 exactly. The CN-CN edge constants are
# differences -- addc[3] = e2_1 - 2 w1, addc[4] = e2_2 - 2 w2, the clique-2
# moment net of its two endpoint nodes -- so they can be negative, while cne
# grows as k^2 against ncn ~ k. HYPOTHESIS: on a CN block S1 crosses zero at a
# size the density sets, the guard fires, and the kernel returns 0. On a
# bipartite block only the bridge channel contributes (ncn = cne = 0) and no
# subtraction occurs, so nothing crosses.
#
# The script does not assume this. It measures the deployed additive value, the
# predicted (S1, S2), and checks the guard against both, so the mechanism either
# reproduces the map exactly or it is reported as unexplained.
#
# Usage, from the repo root (installed build):
#   Rscript certification/zratio_additive_collapse_map.R

suppressMessages(library(bgms))

DELTA = 0.5 * log(12)
# 0.1 and 0.25 are the shapes that matter for exposure: they are BELOW
# .zratio_surface_shape_lo, which is where the additive kernel is what a fit
# actually gets. Everything from 0.5 up is served by the surface or, above 10,
# by the isolated-edge route, so those rows characterize the kernel rather than
# the deployed behaviour. Both slab families are covered because both reach the
# additive path the same way.
SHAPES = c(0.1, 0.25, 0.5, 1, 2, 3, 5, 10, 12, 20)
ETAS = c(1, 2)
SLABS = c("normal", "cauchy")
SIZES = c(3, 4, 5, 6, 8, 10, 12, 16, 20, 26, 32, 42, 60)
DENS = c(0.35, 0.5, 0.7, 1.0)
# Shapes at or above this are not served by the additive kernel in any deployed
# cell; kept in the map to characterize the mechanism, excluded from the
# exposure statement.
FENCE_LO = 0.5

# CN block: `size` common neighbours of the edge (1, 2), with a `dens` share of
# the CN-CN edges present (deterministic thinning, so the map is reproducible).
# Every CN node keeps both endpoint edges, so ncn = size at any density.
cn_block = function(size, dens, seed) {
  q = size + 2L
  G = matrix(0L, q, q)
  for(v in 3:q) G[1, v] = G[v, 1] = G[2, v] = G[v, 2] = 1L
  pairs = if(size >= 2) t(utils::combn(3:q, 2)) else matrix(0L, 0, 2)
  keep = if(nrow(pairs) == 0) integer(0) else {
    set.seed(seed)
    sort(sample.int(nrow(pairs), round(dens * nrow(pairs))))
  }
  for(r in keep) {
    G[pairs[r, 1], pairs[r, 2]] = G[pairs[r, 2], pairs[r, 1]] = 1L
  }
  G
}

# Bipartite bridge block: `size` nodes split between the exclusive neighbour
# sets of the two endpoints, with a `dens` share of the A-B bridges present.
# Isolated nodes drop out of the block, which the count function handles.
bip_block = function(size, dens, seed) {
  q = size + 2L
  na_ = size %/% 2L
  A = 2L + seq_len(na_)
  B = setdiff(3:q, A)
  G = matrix(0L, q, q)
  for(v in A) G[1, v] = G[v, 1] = 1L
  for(v in B) G[2, v] = G[v, 2] = 1L
  pairs = expand.grid(a = A, b = B)
  set.seed(seed)
  keep = sort(sample.int(nrow(pairs), round(dens * nrow(pairs))))
  for(r in keep) {
    G[pairs$a[r], pairs$b[r]] = G[pairs$b[r], pairs$a[r]] = 1L
  }
  G
}

# The engine's own mediating-block counts, recomputed in R from the definition
# in extract_block_: common neighbours of (1, 2), edges among them, and bridge
# edges between the exclusive neighbour sets. Cross-checked below against the
# deployed value through the saddle map, which fails loudly if these are wrong.
block_counts = function(G) {
  q = nrow(G)
  rest = 3:q
  cn = rest[G[1, rest] == 1L & G[2, rest] == 1L]
  ei = rest[G[1, rest] == 1L & G[2, rest] == 0L]
  ej = rest[G[2, rest] == 1L & G[1, rest] == 0L]
  cne = if(length(cn) < 2) 0L else sum(G[cn, cn]) / 2L
  bre = if(length(ei) == 0 || length(ej) == 0) 0L else sum(G[ei, ej])
  list(ncn = length(cn), cne = as.integer(cne), bre = as.integer(bre))
}

cat("=== 1. the additive constants, by cell ===\n")
cat("addc = (w1, w2, ce1, ce2, cb1, cb2): CN node, CN-CN edge, bridge.\n")
cat("ce = clique-2 moment net of its two endpoint nodes, hence a difference.\n\n")
const_rows = list()
zcs = list()
for(alpha in SHAPES) {
  for(eta in ETAS) {
    for(slab in SLABS) {
      key = paste(alpha, eta, slab)
      zc = suppressWarnings(
        bgms:::zratio_constants(DELTA, eta, alpha = alpha, slab = slab)
      )
      zcs[[key]] = zc
      const_rows[[length(const_rows) + 1L]] = data.frame(
        alpha = alpha, eta = eta, slab = slab,
        w1 = zc$addc[1], ce1 = zc$addc[3], cb1 = zc$addc[5],
        w2 = zc$addc[2], ce2 = zc$addc[4], cb2 = zc$addc[6]
      )
    }
  }
}
consts = do.call(rbind, const_rows)
print(consts, digits = 4, row.names = FALSE)
cat(sprintf("\nce1 < 0 in %d of %d cells; cb1 < 0 in %d.\n",
            sum(consts$ce1 < 0), nrow(consts), sum(consts$cb1 < 0)))

cat("\n=== 2. the map: deployed additive value over (family, size, density) ===\n")
rows = list()
for(alpha in SHAPES) {
  for(eta in ETAS) {
    for(slab in SLABS) {
      zc = zcs[[paste(alpha, eta, slab)]]
      for(fam in c("cn", "bip")) {
        for(k in SIZES) {
          for(d in DENS) {
            G = if(fam == "cn") cn_block(k, d, 101L) else bip_block(k, d, 202L)
            ct = block_counts(G)
            s1 = ct$ncn * zc$addc[1] + ct$cne * zc$addc[3] + ct$bre * zc$addc[5]
            s2 = ct$ncn * zc$addc[2] + ct$cne * zc$addc[4] + ct$bre * zc$addc[6]
            lz = bgms:::zratio_test_eval(
              G, matrix(c(1L, 2L), 1, 2), zc$addc, zc$tg, zc$ihat, zc$ghat,
              zc$wt, zc$psi0
            )$log_zratio[1]
            rows[[length(rows) + 1L]] = data.frame(
              alpha = alpha, eta = eta, slab = slab, family = fam, size = k,
              dens = d, ncn = ct$ncn, cne = ct$cne, bre = ct$bre,
              s1 = s1, s2 = s2, log_zratio = lz
            )
          }
        }
      }
    }
  }
}
map = do.call(rbind, rows)
map$collapsed = map$log_zratio == 0
map$guard = map$s1 <= 0   # s2 is floored to 1e-3 before the saddle, s1 is not

cat(sprintf("\ncollapsed (log_zratio exactly 0) in %d of %d cells\n",
            sum(map$collapsed), nrow(map)))
cat("\nby family:\n")
print(table(family = map$family, collapsed = map$collapsed))
cat("\nby slab (CN only):\n")
cn_all = map[map$family == "cn", ]
print(table(slab = cn_all$slab, collapsed = cn_all$collapsed))
cat("\nby shape (CN only):\n")
print(table(shape = cn_all$alpha, collapsed = cn_all$collapsed))
cat("\nby eta (CN only):\n")
print(table(eta = cn_all$eta, collapsed = cn_all$collapsed))

cat("\n--- EXPOSURE: the cells the additive kernel actually serves ---\n")
cat("Shapes below", FENCE_LO, "only; at or above it the surface or the\n")
cat("isolated-edge route serves instead.\n")
dep = map[map$alpha < FENCE_LO, ]
print(table(slab = dep$slab, family = dep$family, collapsed = dep$collapsed))
cat(sprintf("collapsed in %d of %d deployed-cell rows\n",
            sum(dep$collapsed), nrow(dep)))

cat("\n=== 3. mechanism: does the s1 <= 0 guard reproduce the map exactly? ===\n")
agree = map$collapsed == map$guard
cat(sprintf("guard predicts the collapse in %d of %d cells (%d disagreements)\n",
            sum(agree), nrow(map), sum(!agree)))
if(any(!agree)) {
  cat("disagreeing cells:\n")
  print(head(map[!agree, ], 20), digits = 4, row.names = FALSE)
}
# Second leg: where the guard does NOT fire, the deployed value must equal the
# saddle map at the predicted (S1, S2). That certifies the counts and the
# accumulation, so a passing first leg is not an accident of two wrong things.
live = map[!map$guard, ]
pred = vapply(seq_len(nrow(live)), function(r) {
  zc = zcs[[paste(live$alpha[r], live$eta[r], live$slab[r])]]
  # The additive branch replaces a non-positive S2 by kS2Floor = 1e-3 before
  # the saddle; it does not raise a small positive one.
  s2 = if(live$s2[r] <= 0) 1e-3 else live$s2[r]
  log(bgms:::zratio_test_saddle(
    live$s1[r], s2, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0
  ))
}, 0.0)
cat(sprintf("non-collapsed cells: max |deployed - saddle(S1, S2)| = %.3g over %d cells\n",
            max(abs(pred - live$log_zratio)), nrow(live)))

cat("\n=== 4. the collapse boundary in (size, density), CN family ===\n")
cat("smallest CN size that collapses, per (shape, eta, density);",
    "NA = no size in the grid collapses\n")
cn = map[map$family == "cn", ]
bnd = do.call(rbind, lapply(
  split(cn, list(cn$alpha, cn$eta, cn$slab, cn$dens), drop = TRUE),
  function(d) {
    hit = d$size[d$collapsed]
    data.frame(alpha = d$alpha[1], eta = d$eta[1], slab = d$slab[1],
               dens = d$dens[1],
               first_collapse = if(length(hit) == 0) NA_integer_ else min(hit))
  }))
rownames(bnd) = NULL
bnd = bnd[order(bnd$alpha, bnd$eta, bnd$slab, bnd$dens), ]
print(bnd, row.names = FALSE)

cat("\n=== 5. closed-form boundary against the measured one ===\n")
# S1 = k w1 + cne ce1 with cne = dens * k (k - 1) / 2, so with ce1 < 0 the sign
# flips at k = 1 + 2 w1 / (dens * |ce1|). Measured against the grid's own first
# collapsing size, which is coarse from above by construction.
bnd$k_star = NA_real_
for(r in seq_len(nrow(bnd))) {
  zc = zcs[[paste(bnd$alpha[r], bnd$eta[r], bnd$slab[r])]]
  ce1 = zc$addc[3]
  bnd$k_star[r] = if(ce1 >= 0) Inf else {
    1 + 2 * zc$addc[1] / (bnd$dens[r] * abs(ce1))
  }
}
print(bnd, digits = 4, row.names = FALSE)

cat("\n=== 6. the boundary at density 1, scanned at every size ===\n")
# The grid above is coarse, so a k_star between two grid points reads as a
# mismatch. At density 1 the edge count is exact (k(k-1)/2, no thinning), so a
# size-by-size scan turns the comparison into an equality check.
fine_rows = list()
for(alpha in SHAPES) {
  for(eta in ETAS) {
    for(slab in SLABS) {
      zc = zcs[[paste(alpha, eta, slab)]]
      collapsed_k = integer(0)
      for(k in 3:80) {
        G = cn_block(k, 1.0, 101L)
        lz = bgms:::zratio_test_eval(
          G, matrix(c(1L, 2L), 1, 2), zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt,
          zc$psi0
        )$log_zratio[1]
        if(lz == 0) collapsed_k = c(collapsed_k, k)
      }
      ce1 = zc$addc[3]
      ks = if(ce1 >= 0) Inf else 1 + 2 * zc$addc[1] / abs(ce1)
      fine_rows[[length(fine_rows) + 1L]] = data.frame(
        alpha = alpha, eta = eta, slab = slab,
        first_collapse = if(length(collapsed_k) == 0) NA_integer_ else min(collapsed_k),
        predicted = if(is.finite(ks)) ceiling(ks) else NA_integer_,
        contiguous = length(collapsed_k) == 0 ||
          identical(collapsed_k, min(collapsed_k):80L)
      )
    }
  }
}
fine = do.call(rbind, fine_rows)
# A predicted boundary past the top of the scan predicts "no collapse here",
# which is what an NA in the scan column means; scoring that as a mismatch
# would be scoring the scan range, not the formula.
SCAN_HI = 80L
fine$agree = ifelse(
  is.na(fine$predicted) | fine$predicted > SCAN_HI,
  is.na(fine$first_collapse),
  !is.na(fine$first_collapse) & fine$first_collapse == fine$predicted
)
print(fine, row.names = FALSE)
cat(sprintf("\nclosed-form boundary matches the scan in %d of %d cells (predictions past %d count as 'no collapse in range'); collapse is an upper tail in %d of %d\n",
            sum(fine$agree), nrow(fine), SCAN_HI, sum(fine$contiguous),
            nrow(fine)))

saveRDS(list(constants = consts, map = map, boundary = bnd, fine = fine),
        "certification/zratio_additive_collapse.rds")
cat("\nwritten to certification/zratio_additive_collapse.rds\n")
