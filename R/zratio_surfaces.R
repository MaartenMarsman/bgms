# Option-B absolute-moment surfaces for the hierarchical-spec Z-ratio engine.
#
# One theta-independent smooth surface per family (common-neighbour cluster and
# bipartite bridge) predicts the absolute per-component log-S1 and log-S2 from
# (size, density), fit once per analysis at the deployment (eta, delta) to
# block-Gibbs anchors on synthetic components. The engine then decomposes each
# mediating block into disjoint components and sums the per-component surface
# moments (C++ log_zratio surface branch), replacing the online ridge-OLS
# correction. eta is a build parameter, not a switch: the surface is built at
# the analysis's own eta and used at every eta. Deployment is fenced to the
# validated alpha = 1 Normal-slab family (build_surfaces_allmc returns NULL
# otherwise, so the engine keeps the additive path).
#
# Anchors are drawn from the sampler's own kernel via the C++ bare-component
# oracle (zratio_block_oracle_moments). Short chains suffice: the low-order fit
# denoises the Monte-Carlo noise across anchors.

# Random connected graph on nn nodes with ne edges: a random spanning tree
# (uniform attachment over a permutation) plus extra edges filled at random.
zratio_rand_conn_graph = function(nn, ne, seed) {
  set.seed(seed)
  adj = matrix(0L, nn, nn)
  perm = sample(nn)
  for(i in seq_len(nn)[-1]) {
    j = perm[if(i == 2) 1 else sample(i - 1, 1)]
    adj[perm[i], j] = adj[j, perm[i]] = 1L
  }
  extra = ne - (nn - 1)
  if(extra > 0) {
    avail = which(upper.tri(adj) & adj == 0)
    pick = if(length(avail) == 1) avail else sample(avail, min(extra, length(avail)))
    adj[pick] = 1L
    adj = pmax(adj, t(adj))
  }
  adj
}

# Random connected bipartite graph (na_ x nb_) with ne edges; returns the full
# (na_ + nb_) adjacency and the side membership vector.
zratio_rand_conn_bip = function(na_, nb_, ne, seed) {
  set.seed(seed)
  bip = matrix(0L, na_, nb_)
  conn_a = 1
  conn_b = 1
  bip[1, 1] = 1L
  rest = c(if(na_ > 1) paste0("a", 2:na_), if(nb_ > 1) paste0("b", 2:nb_))
  for(nd in sample(rest)) {
    if(startsWith(nd, "a")) {
      i = as.integer(substring(nd, 2))
      bip[i, sample(conn_b, 1)] = 1L
      conn_a = c(conn_a, i)
    } else {
      j = as.integer(substring(nd, 2))
      bip[sample(conn_a, 1), j] = 1L
      conn_b = c(conn_b, j)
    }
  }
  extra = ne - sum(bip)
  if(extra > 0) {
    avail = which(bip == 0)
    pick = if(length(avail) == 1) avail else sample(avail, min(extra, length(avail)))
    bip[pick] = 1L
  }
  adj = rbind(
    cbind(matrix(0L, na_, na_), bip),
    cbind(t(bip), matrix(0L, nb_, nb_))
  )
  list(adj = adj, aside = c(rep(TRUE, na_), rep(FALSE, nb_)))
}

# One CN-cluster anchor: a random connected graph at (size, density), all nodes
# common neighbours (adjacent to both endpoints), moments from the oracle.
zratio_anchor_cn = function(n, dens, zc, sweeps, burn, seed) {
  e = max(n - 1, round(dens * choose(n, 2)))
  adj = zratio_rand_conn_graph(n, e, seed)
  all_rows = 0:(n - 1)
  r = zratio_block_oracle_moments(
    adj, all_rows, all_rows, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    zc$delta, zc$eta, as.integer(sweeps), as.integer(burn), seed
  )
  if(!isTRUE(r$ok)) return(NULL)
  data.frame(size = n, dens = e / choose(n, 2), S1 = r$S1, S2 = r$S2)
}

# One bipartite-bridge anchor (near-balanced canonical split).
zratio_anchor_bip = function(n, dens, zc, sweeps, burn, seed) {
  na_ = max(2, floor(n / 2))
  nb_ = n - na_
  if(nb_ < 2) return(NULL)
  e = min(max(round(dens * na_ * nb_), n - 1), na_ * nb_)
  bp = zratio_rand_conn_bip(na_, nb_, e, seed)
  ee = sum(bp$adj) / 2
  si = which(bp$aside) - 1
  sj = which(!bp$aside) - 1
  r = zratio_block_oracle_moments(
    bp$adj, si, sj, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    zc$delta, zc$eta, as.integer(sweeps), as.integer(burn), seed
  )
  if(!isTRUE(r$ok)) return(NULL)
  data.frame(size = n, dens = ee / (na_ * nb_), S1 = r$S1, S2 = r$S2)
}

# Fit one family's log-S1 / log-S2 surface over the explicit raw monomial basis
# [1, L, L^2, d, d^2, L*d, L^2*d, L*d^2, L^2*d^2], L = log(size). Explicit
# columns fix the 9-coefficient order the C++ deploy evaluates and let a
# rank-deficient design (few unique sizes/densities at tiny q) drop terms
# gracefully: aliased coefficients come back NA and are zeroed. Returns the two
# coefficient vectors, the (size, density) hull, the log-moment ranges (the
# deploy clamps predictions to these +/- 0.1), and size_min (below which the
# engine uses the additive per-component moment). NULL if no usable anchors.
zratio_fit_surface_family = function(anchors) {
  train = anchors[is.finite(anchors$S1) & anchors$S1 > 0 &
                  is.finite(anchors$S2) & anchors$S2 > 0, ]
  if(nrow(train) == 0) return(NULL)
  ll = log(train$size)
  dd = train$dens
  df = data.frame(
    y1 = log(train$S1), y2 = log(train$S2),
    m1 = ll, m2 = ll^2, m3 = dd, m4 = dd^2,
    m5 = ll * dd, m6 = ll^2 * dd, m7 = ll * dd^2, m8 = ll^2 * dd^2
  )
  rhs = "m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8"
  fit1 = stats::lm(stats::as.formula(paste("y1 ~", rhs)), data = df)
  fit2 = stats::lm(stats::as.formula(paste("y2 ~", rhs)), data = df)
  coef9 = function(fit) {
    v = unname(stats::coef(fit))   # [Intercept, m1..m8] in the deploy order
    v[!is.finite(v)] = 0
    v
  }
  list(
    c1 = coef9(fit1), c2 = coef9(fit2),
    size_lo = min(train$size), size_hi = max(train$size),
    dens_lo = min(train$dens), dens_hi = max(train$dens),
    l1_lo = min(df$y1), l1_hi = max(df$y1),
    l2_lo = min(df$y2), l2_hi = max(df$y2),
    size_min = min(train$size)
  )
}

# Build both family surfaces (CN + bipartite) for one analysis. `zc` is the
# fit-time cell (zratio_cell_constants); the surface is fenced to the validated
# alpha = 1 Normal-slab cell (NULL otherwise, so the engine keeps additive).
# max_size caps the anchor sizes at the reachable giant; components larger than
# the trained hull clamp to its edge at deploy. Sparse, law-informed placement
# with short chains (~800-1000 sweeps); the fit denoises. Returns
# list(cn = <family>, bip = <family>) or NULL.
build_surfaces_allmc = function(zc, max_size = 44L, cores = 1L,
                                seed0 = 700000L) {
  if(!identical(zc$slab, "normal") || abs(zc$alpha - 1) > 1e-12) return(NULL)
  cap = as.integer(max_size)

  cn_jobs = rbind(
    expand.grid(n = c(4, 6, 8, 10, 12, 15, 18, 22, 26, 30), d = c(0.7, 0.8, 0.9, 1.0)),
    expand.grid(n = c(3, 4),                                 d = c(0.5, 0.7, 0.85, 1.0)),
    expand.grid(n = c(4, 6, 8, 10, 12, 15),                  d = c(0.35, 0.5)),
    expand.grid(n = c(36, 42),                               d = c(0.8, 0.9))
  )
  cn_jobs = cn_jobs[cn_jobs$n <= cap, , drop = FALSE]
  cn_jobs = rbind(cn_jobs, cn_jobs)                          # 2 reps
  cn_jobs$sweeps = ifelse(cn_jobs$n >= 32, 800L, 1000L)

  bip_jobs = expand.grid(
    n = c(4, 6, 8, 10, 12, 14, 16, 18, 20, 22), d = c(0.55, 0.7, 0.85, 1.0)
  )
  bip_jobs = bip_jobs[bip_jobs$n <= cap, , drop = FALSE]
  bip_jobs = rbind(bip_jobs, bip_jobs)
  bip_jobs$sweeps = 1000L

  cn_rows = do.call(rbind, parallel::mclapply(seq_len(nrow(cn_jobs)), function(k) {
    zratio_anchor_cn(cn_jobs$n[k], cn_jobs$d[k], zc, cn_jobs$sweeps[k], 200L, seed0 + k)
  }, mc.cores = cores))
  bip_rows = do.call(rbind, parallel::mclapply(seq_len(nrow(bip_jobs)), function(k) {
    zratio_anchor_bip(bip_jobs$n[k], bip_jobs$d[k], zc, bip_jobs$sweeps[k], 200L, seed0 + 100000L + k)
  }, mc.cores = cores))

  if(is.null(cn_rows) || is.null(bip_rows)) return(NULL)
  cn = zratio_fit_surface_family(cn_rows)
  bip = zratio_fit_surface_family(bip_rows)
  if(is.null(cn) || is.null(bip)) return(NULL)
  list(cn = cn, bip = bip)
}

# Parallelism for the one-time surface build. Anchors self-seed per job, so the
# result is independent of the core count; only the wall time changes. Default
# is serial (portable, no fork surprises inside a user's own parallel context);
# raise it with options(bgms.zratio_surface_cores = <n>).
zratio_surface_build_cores = function() {
  cores = suppressWarnings(as.integer(getOption("bgms.zratio_surface_cores", 1L)))
  if(length(cores) != 1L || is.na(cores) || cores < 1L) cores = 1L
  cores
}
