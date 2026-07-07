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
# Conventions: K_ii ~ Exp(beta), slab K_ij ~ N(0, sigma^2), tilt
# |K|^delta. The between-graph ratio Z(Gamma-)/Z(Gamma+) is invariant under
# the diagonal congruence Theta = A K A (scale standardization of the
# normalizer), so it depends on (delta, eta) alone, where in bgms parameter
# units eta = (2 * pairwise_scale) * (scale_rate / 2) = pairwise_scale *
# scale_rate (priors act on K/2). The fixed quadrature grids below (cmax,
# Cmax, Tmax, Laguerre ranges) are sized for the sigma = 1 frame, so every
# consumer builds in the standardized cell (delta, sigma = 1, beta = eta)
# via zratio_cell_constants(). Reference implementation: SV/Z don-validation
# (sd_marginal_helpers.R, ks_validation_grid.R).

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

# Cauchy-slab leg mixture in resolvent form: with omega ~ IG(1/2, 1/2)
# (equivalently omega = 1/(2t), t ~ Gamma(1/2, 1)),
#   E[omega^k (A + B omega)^-(k + 1/2)]
#     = sqrt(2/pi) sum_i w_i (2 t_i A + B)^-(k + 1/2)
# exactly, on the plain Gauss-Laguerre grid (t_i, w_i). The transformed
# integrand is analytic in t (the omega-side form carries a sqrt(t) factor
# that defeats polynomial rules), and every leg term in the channel
# moments has this (k, k + 1/2) power pairing.
zratio_omega_mixture = function(n) {
  q = zratio_gauss_quad(n, "laguerre")
  list(t = q$nodes, w = sqrt(2 / pi) * q$weights)
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
# Gauss-Laguerre with weight e^{-beta S}. The slab density multiplies at
# the raw entry (the tilt sees the Schur-shifted entry), so the Cauchy
# variant only swaps the density factor.
zratio_pair_integrals = function(
  delta, sigma, beta, slab = "normal",
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
  slab_dens = if(identical(slab, "cauchy")) {
    function(x) dcauchy(x, 0, sigma)
  } else {
    function(x) dnorm(x, 0, sigma)
  }
  g_one = function(c_val) {
    cm = outer(b, sth)
    dn = slab_dens(cm + c_val)
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
# Cauchy slab: the two legs from the node to the toggled endpoints carry
# independent mixture weights, K_leg | omega ~ N(0, sigma^2 omega), so the
# resolvent splits per leg, (x + t2)^{-(2k+1)} ->
# (x + t2 w_a)^{-(2k+1)/2} (x + t2 w_b)^{-(2k+1)/2} with prefactor
# (w_a w_b)^k, and each leg mixes on the omega grid.
zratio_node_channel = function(delta, sigma, beta, slab = "normal") {
  t2 = 2 * beta * sigma^2
  if(identical(slab, "cauchy")) {
    mx = zratio_omega_mixture(48)
    gl = zratio_gauss_quad(96, "laguerre")
    x = gl$nodes / beta
    cx = (gl$weights / beta) * x^(delta + 1)
    tx = 2 * outer(x, mx$t)
    mix = function(p) as.numeric((tx + t2)^(-p) %*% mx$w)
    den = sum(cx * mix(0.5)^2)
    return(c(
      sigma^4 * sum(cx * mix(1.5)^2) / den,
      sigma^8 * sum(cx * mix(2.5)^2) / den
    ))
  }
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
# draw for draw at the same seed). Cauchy slab: the block coupling runs
# omega-augmented (conjugate IG(1, 1/2 + k^2/(2 sigma^2)) refresh after
# each draw), the four legs draw fresh prior weights sqrt(omega) = 1/|z|
# per kept sweep, and the moments use the leg-dressed block recipe
#   Mi = (I + t2 Wi R Wi)^{-1},  Wt = sqrt(det Mi det Mj),
#   P = (Wi Mi Wi) R (Wj Mj Wj) R,  p1 = sigma^4 tr P, p2 = sigma^8 tr P^2,
# which reduces to the a2-resolvent form at omega = 1.
zratio_clique2_moments = function(
  delta, sigma, beta, slab = "normal", n_mc = 20000, burn = 60, seed = 7
) {
  t2 = 2 * beta * sigma^2
  cauchy = identical(slab, "cauchy")
  set.seed(seed)
  k_mat = diag(rexp(2, beta) + 2, 2)
  s2i = 1 / sigma^2
  om_e = 1
  p1 = p2 = w = numeric(n_mc)
  for(s in 1:(burn + n_mc)) {
    for(i in 1:2) {
      rest = setdiff(1:2, i)
      c_inv = solve(k_mat[rest, rest, drop = FALSE])
      m = 2 * beta * c_inv
      diag(m) = diag(m) + if(cauchy) 1 / (sigma^2 * om_e) else s2i
      r = chol(m)
      bvec = backsolve(r, rnorm(1))
      k_mat[i, rest] = bvec
      k_mat[rest, i] = bvec
      k_mat[i, i] = rgamma(1, delta + 1, beta) +
        as.numeric(t(bvec) %*% c_inv %*% bvec)
      if(cauchy) {
        om_e = (0.5 + bvec^2 / (2 * sigma^2)) / rexp(1)
      }
    }
    if(s > burn) {
      k = s - burn
      if(cauchy) {
        wsi = 1 / abs(rnorm(2))
        wsj = 1 / abs(rnorm(2))
        r_blk = solve(k_mat)
        mi = solve(diag(2) + t2 * (r_blk * outer(wsi, wsi)))
        mj = solve(diag(2) + t2 * (r_blk * outer(wsj, wsj)))
        p_mat = (mi * outer(wsi, wsi)) %*% r_blk %*%
          (mj * outer(wsj, wsj)) %*% r_blk
        w[k] = sqrt(det(mi) * det(mj))
        p1[k] = sigma^4 * sum(diag(p_mat))
        p2[k] = sigma^8 * sum(p_mat * t(p_mat))
      } else {
        a2 = k_mat + t2 * diag(2)
        ri = solve(a2)
        w[k] = det(k_mat) / det(a2)
        p1[k] = sigma^4 * sum(diag(ri %*% ri))
        p2[k] = sigma^8 * sum(diag(ri %*% ri %*% ri %*% ri))
      }
    }
  }
  ok = is.finite(w) & is.finite(p1) & is.finite(p2)
  w = w[ok]
  sw = sum(w)
  c(sum(w * p1[ok]) / sw, sum(w * p2[ok]) / sw)
}

# Two-moment constants for the bridge channel (edge from Si\Sj to Sj\Si):
# 2-node quadrature, Gauss-Laguerre on the diagonals x Gauss-Hermite on the
# coupling. Cauchy slab: the coupling integrates on its exact PD support by
# the sin substitution u = sqrt(kaa kbb) sin(theta) against the Cauchy
# density, and the two legs mix independently on the omega grid; at fixed
# theta the leg resolvents factor, Da = d + t2 w_a kbb, Db = d + t2 w_b kaa,
# with per-moment prefactors (w_a w_b)^k folded into the leg mixtures.
zratio_bridge_channel = function(delta, sigma, beta, slab = "normal",
                                 nlag = 64, nher = 80) {
  t2 = 2 * beta * sigma^2
  gl = zratio_gauss_quad(nlag, "laguerre")
  xa = gl$nodes / beta
  wa = gl$weights
  if(identical(slab, "cauchy")) {
    lgq = zratio_gauss_quad(64, "legendre")
    th = lgq$nodes * (pi / 2)
    wth = lgq$weights * (pi / 2)
    sth = sin(th)
    cth = cos(th)
    mx = zratio_omega_mixture(32)
    den = n1 = n2 = 0
    for(ia in 1:nlag) {
      kaa = xa[ia]
      for(ib in 1:nlag) {
        kbb = xa[ib]
        pw = wa[ia] * wa[ib]
        bmax = sqrt(kaa * kbb)
        u = bmax * sth
        d = bmax^2 * cth^2
        wu = pw * wth * dcauchy(u, 0, sigma) * bmax * cth
        da = 2 * outer(d, mx$t) + t2 * kbb
        db = 2 * outer(d, mx$t) + t2 * kaa
        fa0 = as.numeric(da^(-0.5) %*% mx$w)
        fb0 = as.numeric(db^(-0.5) %*% mx$w)
        fa1 = as.numeric(da^(-1.5) %*% mx$w)
        fb1 = as.numeric(db^(-1.5) %*% mx$w)
        fa2 = as.numeric(da^(-2.5) %*% mx$w)
        fb2 = as.numeric(db^(-2.5) %*% mx$w)
        base = wu * d^delta * d
        den = den + sum(base * fa0 * fb0)
        n1 = n1 + sum(base * sigma^4 * u^2 * fa1 * fb1)
        n2 = n2 + sum(base * sigma^8 * u^4 * fa2 * fb2)
      }
    }
    return(c(n1 / den, n2 / den))
  }
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

# Standardized cell for one fit's Z-ratio constants. The between-graph
# ratio depends on (delta, eta) only, so the constants take just those two:
# eta is the user-specified standardized rate when the scale prior carries
# one, else pairwise_scale * scale_rate (the same number up to rounding).
zratio_cell_constants = function(delta, pairwise_scale, scale_rate,
                                 scale_eta = NA_real_, slab = "normal") {
  eta = if(is.finite(scale_eta)) scale_eta else pairwise_scale * scale_rate
  zratio_constants(delta, eta, slab = slab)
}

# Session cache for zratio_constants: the constant set is deterministic per
# (delta, eta, slab) cell, and one fit resolves the same cell more than once
# (sampler dispatch and diagnostics assembly).
zratio_constants_cache = new.env(parent = emptyenv())

# Full fit-time constant set for one (delta, eta, slab) cell. The estimator
# is defined in the standardized frame, so the slab has unit scale and the
# diagonal rate is eta: the internal channel builders are evaluated at
# sigma = 1, beta = eta (the frame the quadrature grids are sized for), and
# no other scale is representable. addc[1..6] (R indexing) =
# (w1, w2, ce1, ce2, cb1, cb2), plus the saddle grid and
# psi0 = I_spike(0)/G(0). The OLS-correction slots (7..13) and the hull box
# (14..23) are absent here; the warm-up calibrator appends them. The
# clique-2 channel draws seeded Monte Carlo samples, so the caller's RNG
# state is saved and restored. Results are served from a session cache keyed
# on the cell.
zratio_constants = function(delta, eta, slab = "normal") {
  slab = match.arg(slab, c("normal", "cauchy"))
  sigma = 1
  beta = eta
  key = paste(
    format(delta, digits = 17), format(eta, digits = 17), slab,
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
  pair = zratio_pair_integrals(delta, sigma, beta, slab)
  grid = zratio_saddle_grid(pair)
  w12 = zratio_node_channel(delta, sigma, beta, slab)
  # The omega-augmented chain mixes slower than the Normal one, so the
  # Cauchy cell runs longer; the Normal cell keeps its draw-for-draw
  # reference length.
  e2 = zratio_clique2_moments(
    delta, sigma, beta, slab,
    n_mc = if(identical(slab, "cauchy")) 60000 else 20000
  )
  cb = zratio_bridge_channel(delta, sigma, beta, slab)
  addc = c(w12[1], w12[2], e2[1] - 2 * w12[1], e2[2] - 2 * w12[2], cb[1], cb[2])
  out = list(
    delta = delta,
    eta = eta,
    slab = slab,
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
