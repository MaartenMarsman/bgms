# Fit-time constants for the hierarchical-spec per-edge Z-ratio estimator.
#
# The hierarchical prior specification p(K | Gamma) = rho_Gamma(K) / Z(Gamma)
# carries the ratio J = Z(Gamma-)/Z(Gamma+) in every between-edge move. The
# deterministic estimator evaluates J from three local neighbourhood counts
# through a two-moment saddle over pair integrals of the tilted prior. The
# objects built here are its fit-time constants:
#   - pair integrals I_spike(c) (closed-form Bessel) and G(c) (quadrature),
#   - their cosine-transform tables on a t-grid (the saddle grid),
#   - the per-channel moment constants addc[0..5] for the additive counts
#     (common-neighbour node, CN-CN edge, bridge edge),
#   - psi0, the isolated-edge ratio I_spike(0)/G(0).
# Conventions (bare scale): K_ii ~ Exp(beta), slab K_ij ~ N(0, sigma^2),
# tilt |K|^delta. In bgms parameter units: sigma = 2 * pairwise_scale and
# beta = scale_rate / 2 (priors act on K/2). Reference implementation:
# SV/Z don-validation (sd_marginal_helpers.R, ks_validation_grid.R).

# Golub-Welsch Gauss quadrature nodes/weights. kind: "laguerre" (weight
# e^{-x} on (0, Inf)), "hermite" (weight e^{-x^2} on (-Inf, Inf)),
# "legendre" (weight 1 on (-1, 1)).
zratio_gauss_quad = function(n, kind) {
  k = seq_len(n - 1)
  if(kind == "laguerre") {
    a = 2 * seq_len(n) - 1
    b = k
    mu0 = 1
  } else if(kind == "hermite") {
    a = rep(0, n)
    b = sqrt(k / 2)
    mu0 = sqrt(pi)
  } else if(kind == "legendre") {
    a = rep(0, n)
    b = k / sqrt(4 * k^2 - 1)
    mu0 = 2
  } else {
    stop("unknown quadrature kind: ", kind)
  }
  J = diag(a)
  for(i in k) {
    J[i, i + 1] = b[i]
    J[i + 1, i] = b[i]
  }
  ev = eigen(J, symmetric = TRUE)
  ord = order(ev$values)
  list(
    nodes = ev$values[ord],
    weights = mu0 * ev$vectors[1, ord]^2
  )
}

# Closed-form spike pair integral
#   I_spike(c) = 2 Gamma(nu) beta^{-nu} |c|^nu K_nu(2 beta |c|), nu = delta + 1,
# with limit Gamma(nu)^2 beta^{-2 nu} at c = 0.
zratio_ispike = function(c_val, delta, beta) {
  nu = delta + 1
  ac = abs(c_val)
  out = numeric(length(ac))
  small = ac < 1e-8
  out[small] = gamma(nu)^2 * beta^(-2 * nu)
  if(any(!small)) {
    a = ac[!small]
    out[!small] = 2 * gamma(nu) * beta^(-nu) * a^nu * besselK(2 * beta * a, nu)
  }
  out
}

# Slab pair integral G(c) tabulated on [0, cmax] with linear interpolation.
# Inner S12 integral by the sin substitution (removes the boundary
# singularity of (S11 S22 - S12^2)^delta); outer (S11, S22) by
# Gauss-Laguerre with weight e^{-beta S}.
zratio_pair_integrals = function(
  delta, sigma, beta,
  cmax = 18, ngrid = 121, nlag = 48, nleg = 64
) {
  gl = zratio_gauss_quad(nlag, "laguerre")
  xq = gl$nodes / beta
  wq = gl$weights / beta
  lg = zratio_gauss_quad(nleg, "legendre")
  th = lg$nodes * (pi / 2)
  wt = lg$weights * (pi / 2)
  cth = cos(th)
  sth = sin(th)
  s1 = rep(xq, each = nlag)
  s2 = rep(xq, times = nlag)
  w_pair = rep(wq, each = nlag) * rep(wq, times = nlag)
  b = sqrt(s1 * s2)
  g_one = function(c_val) {
    cm = outer(b, sth)
    dn = dnorm(cm + c_val, 0, sigma)
    inner = (b^(2 * delta + 1)) * as.numeric(dn %*% (wt * cth^(2 * delta + 1)))
    sum(w_pair * inner)
  }
  cg = seq(0, cmax, length.out = ngrid)
  gv = vapply(cg, g_one, 0.0)
  list(
    ispike = function(c_val) zratio_ispike(c_val, delta, beta),
    g = approxfun(cg, gv, rule = 2),
    cg = cg,
    gv = gv
  )
}

# Cosine-transform saddle grid: Ihat(t) = 2 int_0^Cmax cos(t c) I_spike(c) dc
# and Ghat(t) likewise for G, on nt t-points with trapezoid end-weights.
zratio_saddle_grid = function(
  pair, Cmax = 40, nc = 8001L, Tmax = 160, nt = 801L
) {
  cg = seq(0, Cmax, length.out = nc)
  dc = cg[2] - cg[1]
  is_v = pair$ispike(cg)
  gv = pair$g(cg)
  wc = rep(dc, nc)
  wc[1] = wc[nc] = dc / 2
  tg = seq(0, Tmax, length.out = nt)
  ih = gh = numeric(nt)
  for(ss in seq(1, nt, by = 100L)) {
    ix = ss:min(ss + 99L, nt)
    m = cos(outer(tg[ix], cg))
    ih[ix] = 2 * as.numeric(m %*% (wc * is_v))
    gh[ix] = 2 * as.numeric(m %*% (wc * gv))
  }
  wt = rep(tg[2] - tg[1], nt)
  wt[1] = wt[nt] = wt[2] / 2
  list(tg = tg, ihat = ih, ghat = gh, wt = wt)
}

# Two-moment constants for the common-neighbour node channel: weighted
# moments of the single-node resolvent under the tilted diagonal prior,
#   w_k = sigma^{4k} int x^{delta+1} (x + t2)^{-(2k+1)} e^{-beta x} dx / I(1).
zratio_node_channel = function(delta, sigma, beta) {
  t2 = 2 * beta * sigma^2
  ip = function(p) {
    integrate(
      function(x) x^(delta + 1) * (x + t2)^(-p) * exp(-beta * x),
      0, Inf,
      rel.tol = 1e-10
    )$value
  }
  i1 = ip(1)
  c(sigma^4 * ip(3) / i1, sigma^8 * ip(5) / i1)
}

# Excess two-moment constants for the CN-CN edge channel: moments of the
# connected 2-clique block minus twice the single-node constants. Estimated
# by a seeded within-block Gibbs run (matches the reference implementation
# draw for draw at the same seed).
zratio_clique2_moments = function(
  delta, sigma, beta, n_mc = 20000, burn = 60, seed = 7
) {
  t2 = 2 * beta * sigma^2
  set.seed(seed)
  k_mat = diag(rexp(2, beta) + 2, 2)
  s2i = 1 / sigma^2
  p1 = p2 = w = numeric(n_mc)
  for(s in 1:(burn + n_mc)) {
    for(i in 1:2) {
      rest = setdiff(1:2, i)
      c_inv = solve(k_mat[rest, rest, drop = FALSE])
      m = 2 * beta * c_inv
      diag(m) = diag(m) + s2i
      r = chol(m)
      bvec = backsolve(r, rnorm(1))
      k_mat[i, rest] = bvec
      k_mat[rest, i] = bvec
      k_mat[i, i] = rgamma(1, delta + 1, beta) +
        as.numeric(t(bvec) %*% c_inv %*% bvec)
    }
    if(s > burn) {
      k = s - burn
      a2 = k_mat + t2 * diag(2)
      ri = solve(a2)
      w[k] = det(k_mat) / det(a2)
      p1[k] = sigma^4 * sum(diag(ri %*% ri))
      p2[k] = sigma^8 * sum(diag(ri %*% ri %*% ri %*% ri))
    }
  }
  ok = is.finite(w) & is.finite(p1) & is.finite(p2)
  w = w[ok]
  sw = sum(w)
  c(sum(w * p1[ok]) / sw, sum(w * p2[ok]) / sw)
}

# Two-moment constants for the bridge channel (edge from Si\Sj to Sj\Si):
# 2-node quadrature, Gauss-Laguerre on the diagonals x Gauss-Hermite on the
# coupling.
zratio_bridge_channel = function(delta, sigma, beta, nlag = 64, nher = 80) {
  t2 = 2 * beta * sigma^2
  gl = zratio_gauss_quad(nlag, "laguerre")
  xa = gl$nodes / beta
  wa = gl$weights
  gh = zratio_gauss_quad(nher, "hermite")
  zk = gh$nodes
  wh = gh$weights / sqrt(pi)
  kab = sigma * sqrt(2) * zk
  den = n1 = n2 = 0
  for(ia in 1:nlag) {
    kaa = xa[ia]
    for(ib in 1:nlag) {
      kbb = xa[ib]
      pw = wa[ia] * wa[ib]
      d = kaa * kbb - kab^2
      pd = d > 0
      if(!any(pd)) next
      da = d[pd] + t2 * kbb
      db = d[pd] + t2 * kaa
      w = (d[pd]^delta) * d[pd] / sqrt(da * db)
      base = pw * wh[pd]
      den = den + sum(base * w)
      n1 = n1 + sum(base * w * sigma^4 * kab[pd]^2 / (da * db))
      n2 = n2 + sum(base * w * sigma^8 * kab[pd]^4 / (da * db)^2)
    }
  }
  c(n1 / den, n2 / den)
}

# Session cache for zratio_constants: the constant set is deterministic per
# (delta, sigma, beta) cell, and one fit resolves the same cell more than
# once (sampler dispatch and diagnostics assembly).
zratio_constants_cache = new.env(parent = emptyenv())

# Full fit-time constant set for one (delta, sigma, beta) cell:
# addc[1..6] (R indexing) = (w1, w2, ce1, ce2, cb1, cb2), the saddle grid,
# and psi0 = I_spike(0)/G(0). The OLS-correction slots (7..13) and the hull
# box (14..23) are absent here; the warm-up calibrator appends them.
# The clique-2 channel draws seeded Monte Carlo samples, so the caller's RNG
# state is saved and restored. Results are served from a session cache keyed
# on the cell.
zratio_constants = function(delta, sigma, beta) {
  key = paste(
    format(delta, digits = 17), format(sigma, digits = 17),
    format(beta, digits = 17),
    sep = "_"
  )
  cached = zratio_constants_cache[[key]]
  if(!is.null(cached)) {
    return(cached)
  }
  has_seed = exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if(has_seed) {
    old_seed = get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(
      if(exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(list = ".Random.seed", envir = globalenv())
      },
      add = TRUE
    )
  }
  pair = zratio_pair_integrals(delta, sigma, beta)
  grid = zratio_saddle_grid(pair)
  w12 = zratio_node_channel(delta, sigma, beta)
  e2 = zratio_clique2_moments(delta, sigma, beta)
  cb = zratio_bridge_channel(delta, sigma, beta)
  addc = c(w12[1], w12[2], e2[1] - 2 * w12[1], e2[2] - 2 * w12[2], cb[1], cb[2])
  out = list(
    delta = delta,
    sigma = sigma,
    beta = beta,
    addc = addc,
    tg = grid$tg,
    ihat = grid$ihat,
    ghat = grid$ghat,
    wt = grid$wt,
    psi0 = pair$ispike(0) / pair$g(0)
  )
  assign(key, out, envir = zratio_constants_cache)
  out
}
