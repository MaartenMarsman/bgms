# Is the additive zero-collapse reachable by a real fit?
#
# The collapse map fixed delta = 0.5*log(12). But bgm() resolves delta from the
# model dimension (0.5*log(p)), and the largest common-neighbour block an edge
# can have is p - 2. So the deployed question is not "what is k* at one delta"
# but "does k*(delta(p)) ever fall below p - 2". MEASURED here.
suppressMessages(devtools::load_all(".", quiet = TRUE))
QS = c(6, 8, 10, 12, 16, 20, 30, 50, 80, 120)
rows = list()
for(q in QS) {
  delta = 0.5 * log(q)
  for(alpha in c(0.1, 0.25)) {
    for(eta in c(1, 2)) {
      zc = suppressWarnings(bgms:::zratio_constants(delta, eta, alpha = alpha))
      w1 = zc$addc[1]
      ce1 = zc$addc[3]
      kstar = if(ce1 >= 0) Inf else 1 + 2 * w1 / abs(ce1)
      rows[[length(rows) + 1L]] = data.frame(
        q = q, delta = delta, alpha = alpha, eta = eta, ce1 = ce1,
        k_star = kstar, k_max = q - 2, reachable = (q - 2) >= kstar
      )
    }
  }
}
d = do.call(rbind, rows)
print(d, digits = 4, row.names = FALSE)
cat(sprintf("\nreachable at default delta in %d of %d (q, shape, eta) cells\n",
            sum(d$reachable), nrow(d)))
saveRDS(d, "certification/zratio_collapse_reachability.rds")
