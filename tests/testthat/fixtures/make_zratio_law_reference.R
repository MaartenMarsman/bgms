# Generates tests/testthat/fixtures/zratio_law_reference.rds: companion R
# eval_mu_law (via law_mu_cpp) moments at a handful of (eta, delta, n, dens)
# CN cells. The bgms C++ port (zratio_law_moments) must reproduce these to
# porting tolerance. The generator needs the (non-public) companion ggm_paper
# sources; the committed .rds is the frozen reference the test reads, so this
# script only needs re-running if the companion law changes. Run from the bgms
# package root, pointing BGMS_GGM_PAPER at the ggm_paper checkout:
#   BGMS_GGM_PAPER=/path/to/ggm_paper \
#     Rscript tests/testthat/fixtures/make_zratio_law_reference.R

COMPANION <- Sys.getenv("BGMS_GGM_PAPER")
if(!nzchar(COMPANION)) {
  stop("Set BGMS_GGM_PAPER to the ggm_paper checkout (the companion sources).")
}
# Resolve against the (existing) package root before the setwd below:
# normalizePath on a not-yet-existing file returns the path unresolved, which
# would make OUT relative to the companion checkout when regenerating.
OUT <- file.path(
  normalizePath("."),
  "tests", "testthat", "fixtures", "zratio_law_reference.rds"
)

setwd(COMPANION)
Sys.setenv(SURFB_FUNCS_ONLY = "1", S1_FUNCS_ONLY = "1", OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
})
source("R/scripts/surface-b.R")   # re_param2 (SURFB_FUNCS_ONLY: no driver)
source("R/scripts/s1-correction.R")   # law stack + cpp ports
s1_setup()
LAW_BUDGET <<- 3000L               # well-converged reference

q <- 50
cells <- list(
  list(eta = 2, delta = 0.5 * log(q)),
  list(eta = 1, delta = 0.5 * log(q))
)
grid <- expand.grid(
  n = c(6L, 10L, 16L, 22L), dens = c(0.8, 1.0), KEEP.OUT.ATTRS = FALSE
)

fixture <- list()
for (cell in cells) {
  re_param2(cell$eta, cell$delta)
  rows <- list()
  for (r in seq_len(nrow(grid))) {
    n <- grid$n[r]
    dens <- grid$dens[r]
    D <- dens * (n - 1)
    res <- law_mu_cpp(D = D, n = n)
    rows[[length(rows) + 1]] <- data.frame(
      n = n, dens = dens, D = D,
      S1 = res$out$S1d, S2 = res$out$S2d, psi_gap = res$out$psi_gap
    )
  }
  fixture[[length(fixture) + 1]] <- list(
    eta = cell$eta, delta = cell$delta,
    anchors = do.call(rbind, rows)
  )
}
saveRDS(fixture, OUT)
cat("wrote", OUT, "\n")
for (fc in fixture) {
  cat(sprintf("--- eta=%g delta=%.4f ---\n", fc$eta, fc$delta))
  print(fc$anchors)
}
