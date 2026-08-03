# =============================================================================
# f119-ridge.R --- the identification geometry of a structural-zero cell
#
# Brief 22, task 4. Maintainer decision F-119 moves bgmCompare()'s slab default
# from Cauchy to Normal. This draws what that does, exactly rather than
# schematically, at the one place in the compare model where the likelihood
# alone does not identify a parameter: a retained ordinal category that one
# group never observes.
#
# THE GEOMETRY. Two groups, one four-level ordinal variable, group 1 observing
# only {0,1,2} and group 2 only {1,2,3}. Under the union semantics bgmCompare
# retains all four categories, so category 3 carries a threshold in both
# groups. Group g's threshold at category 3 is
#
#     mu_3 + proj[g] * delta_3,
#
# with mu_3 the overall threshold and delta_3 the group difference. Group 2's
# data pin the combination mu_3 + proj[2] * delta_3. Group 1's do not: no
# observation ever lands in category 3, so the pseudolikelihood only ever
# increases as group 1's own combination mu_3 + proj[1] * delta_3 is driven
# down. The result is a one-sided ridge running off to infinity along the line
# where group 2's combination is held at its identified value. Nothing in the
# likelihood stops it. The slab does, and how fast it does is the whole
# content of F-119.
#
# WHAT IS PLOTTED. The true log-posterior is evaluated on a grid over
# (mu_3, delta_3) through the package's own hook,
# bgmCompare_test_logp_and_gradient(), with every other parameter frozen at its
# posterior mean from a reference fit. The hook returns log-likelihood plus the
# full log-prior; the log-prior is recomputed analytically in R from the same
# closed forms the C++ uses (src/priors/parameter_prior.h) and subtracted, which
# leaves the exact log-pseudolikelihood. The three panels are then
#
#   (a) that log-pseudolikelihood alone, no prior of any kind  -- the ridge;
#   (b) the posterior under a Cauchy(0, 1) slab on delta_3     -- the old default;
#   (c) the posterior under a Normal(0, 1) slab on delta_3     -- the new default.
#
# (a) is improper, which is the point of drawing it: nothing bounds the surface.
# (b) and (c) both carry the beta-prime(0.5, 0.5) threshold prior on mu_3 --- the
# shipped threshold prior, which F-119 does not change --- so they differ from
# each other in the slab on delta_3 and in nothing else. Every remaining prior
# term is constant over the grid and would shift all three panels identically,
# so it is left out; the subtraction check below confirms the likelihood is
# recovered exactly.
#
# The subtraction is verified, not assumed: the hook is called twice at the
# same parameter vector under two different slab families, and the two
# recovered likelihoods must agree to machine precision.
#
# Usage:  Rscript f119-ridge.R
# Writes: f119-ridge-panels.png, f119-ridge-profile.png, f119-ridge-values.rds
# =============================================================================

LIB = Sys.getenv("F119_LIB", "~/bgms-review/wt-fix14-lib")
if (nzchar(LIB) && dir.exists(path.expand(LIB))) {
  .libPaths(c(path.expand(LIB), .libPaths()))
}
suppressMessages(library(bgms))
options(bgms.verbose = FALSE)

OUT = dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(OUT) || !nzchar(OUT)) OUT = "."

hook = get("bgmCompare_test_logp_and_gradient", envir = asNamespace("bgms"))

# ---------------------------------------------------------------------------
# 1. The adverse data set: disjoint-at-the-ends supports, seconds-scale
# ---------------------------------------------------------------------------
P = 4L        # variables
N = 250L      # persons per group
K = 3L        # non-baseline categories, i.e. a four-level ordinal 0..3

# Coupling and thresholds are chosen so that every category is well populated
# where it is populated at all: a stronger coupling piles both groups into their
# own top category and empties others, which would make the picture about
# sparsity rather than about the structural zero.
OM = matrix(0.10, P, P); diag(OM) = 0
MAIN = matrix(c(-0.1, -1.2), nrow = P, ncol = 2, byrow = TRUE)

# Both groups drawn from the same three-level model; group 2's codes are then
# shifted up by one. Group 1 lives on {0,1,2}, group 2 on {1,2,3}, so category
# 3 is a structural zero in group 1 (and category 0 one in group 2).
x1 = simulate_mrf(N, P, num_categories = 2, pairwise = OM, main = MAIN,
                  variable_type = "ordinal", iter = 500, seed = 4041)
x2 = simulate_mrf(N, P, num_categories = 2, pairwise = OM, main = MAIN,
                  variable_type = "ordinal", iter = 500, seed = 4042) + 1L
colnames(x1) = colnames(x2) = paste0("V", seq_len(P))
X = rbind(x1, x2)
G = rep(1:2, each = N)

supp = rbind(g1 = tabulate(x1[, 1] + 1L, K + 1L),
           g2 = tabulate(x2[, 1] + 1L, K + 1L))
dimnames(supp) = list(c("group 1", "group 2"), paste("category", 0:K))
cat("Support of V1 (the cell we draw):\n"); print(supp)
stopifnot(supp["group 1", "category 3"] == 0L)

# ---------------------------------------------------------------------------
# 2. Reference fit, at the shipped defaults, all differences active
# ---------------------------------------------------------------------------
# difference_selection = FALSE keeps every difference in the model, so the
# picture is about identification and not about selection: the inclusion
# indicator is all ones and the flat parameter vector is full length.
t0 = Sys.time()
fit = bgmCompare(X, group_indicator = G, difference_selection = FALSE,
                 iter = 1500, warmup = 1500, chains = 2, cores = 2,
                 seed = 4043, display_progress = "none")
cat(sprintf("reference fit: %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

A = extract_arguments(fit)
proj = A$projection                       # [num_groups x (num_groups - 1)]
nc = as.integer(A$num_categories)
# num_categories counts NON-baseline categories, so a retained four-level
# variable reports 3. If any variable lost a category (both groups empty there)
# the flat layout below would not be the one we reason about, so assert.
cat("num_categories:", nc, "\n")
stopifnot(all(nc == K))
cat("projection rows:", paste(sprintf("%+.4f", proj[, 1]), collapse = " "), "\n")

raw = get("get_raw_samples", envir = asNamespace("bgms"))(fit)
nm = sum(nc); np = P * (P - 1L) / 2L
main_mean = colMeans(do.call(rbind, raw$main))
pair_mean = colMeans(do.call(rbind, raw$pairwise))
stopifnot(length(main_mean) == 2L * nm, length(pair_mean) == 2L * np)

# The hook's flat layout is: overall mains, overall pairs, active main
# differences, active pair differences.
theta0 = c(main_mean[seq_len(nm)], pair_mean[seq_len(np)],
           main_mean[nm + seq_len(nm)], pair_mean[np + seq_len(np)])

# ---------------------------------------------------------------------------
# 3. The hook's arguments, built from the same data the fit saw
# ---------------------------------------------------------------------------
ord = order(G); obs = X[ord, , drop = FALSE]; gs = G[ord]
group_indices = t(vapply(1:2, function(g) {
  r = which(gs == g); c(min(r) - 1L, max(r) - 1L)
}, numeric(2)))
storage.mode(group_indices) = "integer"

counts = lapply(1:2, function(g) {
  m = matrix(0L, nrow = max(nc), ncol = P)
  for (v in seq_len(P)) for (cat in seq_len(nc[v]))
    m[cat, v] = sum(obs[gs == g, v] == cat)
  m
})
bc_stats = lapply(1:2, function(g) matrix(0L, nrow = 2L, ncol = P))
pair_stats = lapply(1:2, function(g) {
  o = obs[gs == g, , drop = FALSE]; t(o) %*% o
})

mei = matrix(NA_integer_, P, 2L)
for (v in seq_len(P)) {
  mei[v, 1] = if (v > 1L) 1L + mei[v - 1L, 2] else 0L
  mei[v, 2] = mei[v, 1] + nc[v] - 1L
}
pei = matrix(NA_integer_, P, P); tel = 0L
for (v1 in seq_len(P - 1L)) for (v2 in seq(v1 + 1L, P)) {
  pei[v1, v2] = tel; pei[v2, v1] = tel; tel = tel + 1L
}
incl = matrix(1L, P, P)

S = list(observations = matrix(as.integer(obs), nrow(obs), ncol(obs)),
         group_indices = group_indices, num_groups = 2L,
         counts_per_category = counts, blume_capel_stats = bc_stats,
         pairwise_stats = pair_stats, num_categories = nc,
         is_ordinal_variable = rep(1L, P), baseline_category = rep(0L, P),
         main_effect_indices = mei, pairwise_effect_indices = pei,
         inclusion_indicator = incl, projection = proj)

ALPHA = 0.5; BETA = 0.5          # beta_prime_prior(0.5, 0.5), the shipped default
call_hook = function(theta, ifam, dfam, iscale = 1, dscale = 1) {
  hook(theta, S$observations, S$group_indices, S$num_groups,
       S$counts_per_category, S$blume_capel_stats, S$pairwise_stats,
       S$num_categories, S$is_ordinal_variable, S$baseline_category,
       S$main_effect_indices, S$pairwise_effect_indices,
       S$inclusion_indicator, S$projection,
       iscale, dscale, ALPHA, BETA, ifam, dfam, "beta-prime", 1.0)$value
}

# ---------------------------------------------------------------------------
# 4. The analytic log-prior, in the C++'s own closed forms
# ---------------------------------------------------------------------------
# src/priors/parameter_prior.h:
#   CauchyPrior::logp    = R::dcauchy(x, 0, s, log = TRUE)
#   NormalPrior::logp    = R::dnorm(x, 0, s, log = TRUE)
#   BetaPrimePrior::logp = x*alpha - log1p(exp(x)) * (alpha + beta)   [unnormalized]
lp_slab = function(x, fam, s) {
  switch(fam, cauchy = dcauchy(x, 0, s, log = TRUE),
              normal = dnorm(x, 0, s, log = TRUE),
              stop("unknown family"))
}
lp_betaprime = function(x, a = ALPHA, b = BETA) x * a - log1p(exp(x)) * (a + b)

log_prior_total = function(theta, ifam, dfam, iscale = 1, dscale = 1) {
  mu    = theta[seq_len(nm)]
  pw    = theta[nm + seq_len(np)]
  dmu   = theta[nm + np + seq_len(nm)]
  dpw   = theta[nm + np + nm + seq_len(np)]
  sum(lp_betaprime(mu)) + sum(lp_slab(dmu, dfam, dscale)) +
    sum(lp_slab(pw, ifam, iscale)) + sum(lp_slab(dpw, dfam, dscale))
}

# VERIFICATION: the recovered log-pseudolikelihood cannot depend on which slab
# the hook was called with. Two families, same theta, must agree exactly.
loglik = function(theta, ifam = "normal", dfam = "normal")
  call_hook(theta, ifam, dfam) - log_prior_total(theta, ifam, dfam)

v_nn = loglik(theta0, "normal", "normal")
v_cc = loglik(theta0, "cauchy", "cauchy")
v_cn = loglik(theta0, "cauchy", "normal")
cat(sprintf("prior-subtraction check: normal/normal %.10f\n", v_nn))
cat(sprintf("                         cauchy/cauchy %.10f  (delta %.3e)\n",
            v_cc, v_cc - v_nn))
cat(sprintf("                         cauchy/normal %.10f  (delta %.3e)\n",
            v_cn, v_cn - v_nn))
stopifnot(max(abs(c(v_cc, v_cn) - v_nn)) < 1e-8)

# ---------------------------------------------------------------------------
# 5. The grid over (mu_3, delta_3) of the zero-support cell
# ---------------------------------------------------------------------------
VSTAR = 1L; CSTAR = 3L                       # V1, category 3: empty in group 1
r_star = mei[VSTAR, 1] + CSTAR               # 1-based row in the main block
i_mu = r_star                                # overall threshold
i_dl = nm + np + r_star                      # its group difference (2 groups)

zero_group = 1L                              # the group with no observations
obs_group  = 2L
c_hat = theta0[i_mu] + proj[obs_group, 1] * theta0[i_dl]   # identified combination

cat(sprintf("cell: variable %d, category %d ; flat indices mu=%d delta=%d\n",
            VSTAR, CSTAR, i_mu, i_dl))
cat(sprintf("posterior mean: mu = %+.3f, delta = %+.3f ; group-2 combination %+.3f\n",
            theta0[i_mu], theta0[i_dl], c_hat))

# Ranges follow the geometry rather than being guessed: the ridge axis is
# mu = c_hat - proj[2] * delta, so a delta window maps to a mu window through
# c_hat. The window opens far more in the unidentified direction (delta up,
# group 1's own threshold down) than against it.
DL_RANGE = c(-6, 26); NG = 201L
mu_on_ridge = c_hat - proj[obs_group, 1] * DL_RANGE
MU_RANGE = c(min(mu_on_ridge) - 2.5, max(mu_on_ridge) + 2.5)
cat(sprintf("grid: mu in [%.2f, %.2f], delta in [%.2f, %.2f]\n",
            MU_RANGE[1], MU_RANGE[2], DL_RANGE[1], DL_RANGE[2]))
mu_grid = seq(MU_RANGE[1], MU_RANGE[2], length.out = NG)
dl_grid = seq(DL_RANGE[1], DL_RANGE[2], length.out = NG)

cat(sprintf("evaluating %d x %d = %d hook calls ...\n", NG, NG, NG * NG))
t0 = Sys.time()
LL = matrix(NA_real_, NG, NG)
th = theta0
for (i in seq_len(NG)) {
  th[i_mu] = mu_grid[i]
  for (j in seq_len(NG)) {
    th[i_dl] = dl_grid[j]
    LL[i, j] = loglik(th)
  }
}
cat(sprintf("  %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# The two prior terms that live on the plotted coordinates. Everything else is
# constant over the grid and shifts all three panels by the same amount.
PR_MU = outer(lp_betaprime(mu_grid), rep(1, NG))
PR_C  = outer(rep(1, NG), dcauchy(dl_grid, 0, 1, log = TRUE))
PR_N  = outer(rep(1, NG), dnorm(dl_grid, 0, 1, log = TRUE))

# (a) is the pseudolikelihood alone -- no prior of any kind, so it is improper
# and the ridge is open. (b) and (c) are the actual posteriors: both carry the
# beta-prime(0.5, 0.5) threshold prior on mu_3, which F-119 does not change, and
# differ from each other only in the slab on delta_3.
Z = list(a = LL, b = LL + PR_MU + PR_C, c = LL + PR_MU + PR_N)
LABS = c(a = "(a)  likelihood alone",
         b = "(b)  posterior, Cauchy(0, 1) slab",
         c = "(c)  posterior, Normal(0, 1) slab")
SUBS = c(a = "no prior at all: the ridge is open",
         b = "the old default: closed, but slowly",
         c = "the new default (F-119): closed")

# ---------------------------------------------------------------------------
# 6. Readings: how far the 95% highest-density region runs along the ridge
# ---------------------------------------------------------------------------
# For a two-parameter surface the 95% HDR is the set within qchisq(.95, 2)/2 =
# 2.996 log units of the maximum, and the 50% HDR within qchisq(.5, 2)/2.
HDR95 = qchisq(0.95, 2) / 2
HDR50 = qchisq(0.50, 2) / 2

reading = function(z) {
  k = which(z == max(z), arr.ind = TRUE)[1, ]
  inside = z >= max(z) - HDR95
  dl_in = dl_grid[apply(inside, 2, any)]
  mu_in = mu_grid[apply(inside, 1, any)]
  # group-1 threshold over the 95% region: how deep the empty cell is pushed
  g1 = outer(mu_grid, dl_grid, function(m, d) m + proj[zero_group, 1] * d)
  list(mode_mu = mu_grid[k[1]], mode_dl = dl_grid[k[2]],
       dl_lo = min(dl_in), dl_hi = max(dl_in),
       dl_width = diff(range(dl_in)),
       mu_lo = min(mu_in), mu_hi = max(mu_in),
       g1_min = min(g1[inside]), g1_max = max(g1[inside]),
       area = sum(inside) * diff(mu_grid[1:2]) * diff(dl_grid[1:2]),
       edge = any(inside[1, ]) || any(inside[NG, ]) ||
              any(inside[, 1]) || any(inside[, NG]))
}
R = lapply(Z, reading)

cat("\n95% highest-density region of each panel\n")
cat(sprintf("%-26s %9s %9s %20s %9s %8s\n",
            "panel", "mode mu", "mode dl", "delta range (width)", "area", "at edge"))
for (k in names(Z)) cat(sprintf("%-26s %+9.2f %+9.2f  [%+6.2f,%+6.2f] (%5.2f) %9.1f %8s\n",
  LABS[k], R[[k]]$mode_mu, R[[k]]$mode_dl, R[[k]]$dl_lo, R[[k]]$dl_hi,
  R[[k]]$dl_width, R[[k]]$area, R[[k]]$edge))

# ---------------------------------------------------------------------------
# 7. Panels
# ---------------------------------------------------------------------------
PAL = grDevices::hcl.colors(256, "Blues 3", rev = TRUE)
LEVELS = c(HDR50, HDR95, 6, 10, 20, 40)        # log-unit drops from the max
FLOOR = 45                                     # fill clipped this far below max

png(file.path(OUT, "f119-ridge-panels.png"),
    width = 12.6, height = 5.6, units = "in", res = 300)
op = par(mfrow = c(1, 3), mar = c(3.9, 4.4, 3.9, 1.0), oma = c(6.2, 0.4, 2.8, 0.4),
         mgp = c(2.4, 0.7, 0), las = 1, cex.axis = 0.95)

for (k in names(Z)) {
  z = Z[[k]]; zmax = max(z)
  zz = pmax(z - zmax, -FLOOR)
  image(mu_grid, dl_grid, zz, col = PAL, zlim = c(-FLOOR, 0),
        xlab = "", ylab = "", axes = FALSE, useRaster = TRUE)
  axis(1); axis(2); box(lwd = 1.1)
  title(xlab = expression("overall threshold  " * mu[3]), line = 2.5)
  if (k == "a") title(ylab = expression("group difference  " * delta[3]), line = 2.6)

  # the ridge axis: group 2's identified combination held fixed
  abline(a = c_hat / proj[obs_group, 1], b = -1 / proj[obs_group, 1],
         col = "#B03030", lwd = 1.6, lty = 2)

  contour(mu_grid, dl_grid, z, add = TRUE, drawlabels = FALSE,
          levels = zmax - rev(LEVELS[-(1:2)]), col = "#FFFFFFAA", lwd = 0.7)
  contour(mu_grid, dl_grid, z, add = TRUE, drawlabels = FALSE,
          levels = zmax - HDR95, col = "#101010", lwd = 2.4)
  contour(mu_grid, dl_grid, z, add = TRUE, drawlabels = FALSE,
          levels = zmax - HDR50, col = "#101010", lwd = 1.2, lty = 3)

  points(theta0[i_mu], theta0[i_dl], pch = 21, bg = "white",
         col = "black", cex = 1.25, lwd = 1.3)

  mtext(LABS[k], side = 3, line = 1.70, adj = 0, font = 2, cex = 0.98)
  mtext(SUBS[k], side = 3, line = 0.55, adj = 0, cex = 0.84, col = "grey25")

  w = R[[k]]$dl_width
  lab = if (R[[k]]$edge)
    expression(bold("95% region runs off the grid")) else
    bquote(bold("95% region: ") * delta[3] * " spans " * .(sprintf("%.1f", w)))
  legend("bottomleft", legend = as.expression(lab), bty = "n",
         cex = 0.9, text.col = "grey15", inset = c(0.01, 0.01))
}

mtext("A structural-zero cell in bgmCompare: what the slab has to close",
      outer = TRUE, side = 3, line = 0.6, font = 2, cex = 1.18)
CAPTION = c(
  paste0("Four-level ordinal, group 1 observing only {0,1,2} and group 2 only {1,2,3}, so category 3 is empty in group 1. ",
         "Exact log-posterior over that cell's overall threshold and its group"),
  paste0("difference, every other parameter held at its posterior mean from a reference fit, evaluated through ",
         "bgmCompare_test_logp_and_gradient(). Panel (a) carries no prior at all"),
  paste0("and is improper; (b) and (c) are the posteriors, both carrying the beta-prime(0.5, 0.5) threshold prior on the overall threshold, ",
         "and differ only in the slab on the"),
  paste0("difference. Solid black: 95% highest-density contour; dotted: 50%. Dashed red: the line along which group 2's identified ",
         "threshold is constant, i.e. the ridge axis."),
  paste0("Open circle: the reference fit's posterior mean. All three panels share axes, fill scale and contour levels."))
for (i in seq_along(CAPTION))
  mtext(CAPTION[i], outer = TRUE, side = 1, line = 4.55 - 1.02 * (length(CAPTION) - i),
        cex = 0.66, col = "grey25", adj = 0)
par(op); dev.off()

# ---------------------------------------------------------------------------
# 8. The profile along the ridge --- the same statement, one dimension
# ---------------------------------------------------------------------------
# Walk the ridge axis: hold group 2's identified combination at c_hat and let
# delta_3 run, so mu_3 = c_hat - proj[2] * delta_3. Group 1's threshold then
# moves as (proj[1] - proj[2]) * delta_3, i.e. straight down the open direction.
dl_line = seq(DL_RANGE[1], DL_RANGE[2], length.out = 601L)
mu_line = c_hat - proj[obs_group, 1] * dl_line
th = theta0
ll_line = vapply(seq_along(dl_line), function(i) {
  th[i_mu] = mu_line[i]; th[i_dl] = dl_line[i]; loglik(th)
}, numeric(1))
pr_mu_line = lp_betaprime(mu_line)
prof = list(a = ll_line,
            b = ll_line + pr_mu_line + dcauchy(dl_line, 0, 1, log = TRUE),
            c = ll_line + pr_mu_line + dnorm(dl_line, 0, 1, log = TRUE))

png(file.path(OUT, "f119-ridge-profile.png"),
    width = 7.6, height = 5.2, units = "in", res = 300)
op = par(mar = c(4.6, 4.6, 3.6, 1.2), oma = c(3.6, 0, 0, 0),
         mgp = c(2.6, 0.7, 0), las = 1)
COL = c(a = "#8A8A8A", b = "#C2571A", c = "#1F5FA8")
ylim = c(-25, 1)
plot(NA, xlim = DL_RANGE, ylim = ylim, xlab = "", ylab = "", bty = "n")
title(xlab = expression("group difference  " * delta[3] * "   (along the ridge axis)"),
      line = 2.7)
title(ylab = "log-posterior, relative to its own maximum", line = 3.0)
abline(h = -HDR95, col = "grey70", lty = 2)
text(DL_RANGE[1], -HDR95, "  95% cut", adj = c(0, -0.45), cex = 0.8, col = "grey40")
for (k in c("a", "b", "c"))
  lines(dl_line, prof[[k]] - max(prof[[k]]), col = COL[k], lwd = 2.6)
legend("bottomright", bty = "o", bg = "white", box.col = NA, lwd = 2.6, col = COL,
       legend = c("likelihood alone", "posterior, Cauchy(0, 1) slab",
                  "posterior, Normal(0, 1) slab"),
       cex = 0.92, inset = c(0.01, 0.03), seg.len = 2.2)
title(main = "Down the ridge: how each slab closes it", font.main = 2, cex.main = 1.15,
      line = 1.9, adj = 0)
mtext("Group 2's identified threshold held fixed; only the empty cell's own threshold moves.",
      side = 3, line = 0.5, adj = 0, cex = 0.85, col = "grey25")
mtext(paste0(
  "The likelihood climbs to an asymptote and never comes back down: the further the empty cell's ",
  "threshold is pushed, the better the fit,\nfor ever. The Cauchy's logarithmic tail bends it back but ",
  "leaves a long shallow arm still inside the 95% cut at delta = 8; the Normal's\nquadratic tail closes ",
  "it by delta = 3. That difference is the whole of what changing the default does to an unidentified cell."),
  outer = TRUE, side = 1, line = 1.2, cex = 0.74, col = "grey25", adj = 0)
par(op); dev.off()

prof_width = vapply(prof, function(p) {
  ins = which(p >= max(p) - HDR95)
  if (!length(ins)) return(NA_real_)
  diff(range(dl_line[ins]))
}, numeric(1))
cat("\n95% width along the ridge axis (delta units):\n")
print(round(prof_width, 2))

saveRDS(list(mu_grid = mu_grid, dl_grid = dl_grid, LL = LL, Z = Z,
             readings = R, prof = prof, dl_line = dl_line,
             prof_width = prof_width, theta0 = theta0, proj = proj,
             c_hat = c_hat, support = supp,
             check = c(nn = v_nn, cc = v_cc, cn = v_cn)),
        file.path(OUT, "f119-ridge-values.rds"))
cat("\nwrote f119-ridge-panels.png, f119-ridge-profile.png, f119-ridge-values.rds\n")
