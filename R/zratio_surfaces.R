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
# supported Gamma-shape range (Normal or Cauchy slab); a shape outside
# [.zratio_surface_shape_lo, .zratio_surface_shape_hi] returns NULL from
# zratio_build_surfaces, so the engine keeps the additive path.
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
# common neighbours (adjacent to both endpoints), moments from the oracle. The
# cell's Gamma shape is passed to the oracle: it selects the independence-
# Metropolis row update, without which the anchors would be drawn at the
# exponential shape while the constants carry the cell's own.
zratio_anchor_cn = function(n, dens, zc, sweeps, burn, seed) {
  e = max(n - 1, round(dens * choose(n, 2)))
  adj = zratio_rand_conn_graph(n, e, seed)
  all_rows = 0:(n - 1)
  r = zratio_block_oracle_moments(
    adj, all_rows, all_rows, zc$addc, zc$tg, zc$ihat, zc$ghat, zc$wt, zc$psi0,
    zc$delta, zc$eta, as.integer(sweeps), as.integer(burn), seed,
    slab_cauchy = identical(zc$slab, "cauchy"), alpha = zc$alpha
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
    zc$delta, zc$eta, as.integer(sweeps), as.integer(burn), seed,
    slab_cauchy = identical(zc$slab, "cauchy"), alpha = zc$alpha
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
# fit-time cell (zratio_cell_constants). Built for the Normal slab and for the
# Cauchy slab (a scale-mixture of normals the oracle draws directly, so the same
# anchor machinery covers it). Fenced to the supported Gamma-shape range: a
# shape outside it returns NULL and the engine keeps the additive path.
# max_size caps the anchor sizes at the reachable giant; components larger than
# the trained hull clamp to its edge at deploy. The anchor grids come from
# zratio_anchor_grids. Returns list(cn = <family>, bip = <family>) or NULL.
# Session cache for the one-time surface build. The surface depends only on the
# fit cell (delta, eta, alpha, slab) and the size cap, never on the data, so it
# is built once per cell and reused. Backed by disk (same directory and toggle
# convention as the edge-selection correction tables) so the build is also
# skipped across sessions.
.zratio_surface_cache = new.env(parent = emptyenv())

zratio_surface_cache_key = function(zc, max_size, seed0) {
  # The package version is part of the key: a release that changes the anchor
  # grids, the basis, or the fit must not be served a surface cached by an
  # earlier version (the one-time rebuild per cell is seconds). The vN tag
  # carries the same guarantee within a version, and is bumped whenever the
  # grids or the fit change during development.
  sprintf(
    "zratio_surf_v2_%s_delta%.8g_eta%.8g_alpha%.8g_%s_ms%d_sd%d",
    as.character(utils::packageVersion("bgms")),
    as.numeric(zc$delta), as.numeric(zc$eta), as.numeric(zc$alpha),
    as.character(zc$slab), as.integer(max_size), as.integer(seed0)
  )
}

# ------------------------------------------------------------------------------
# zratio_cap_tier
# ------------------------------------------------------------------------------
# Adds an anchor tier at the size cap when the surviving grid stops short of it.
# The grid filter keeps tiers at or below the cap, so a cap that lands between
# two tiers trains the hull at the lower one and every larger block is
# extrapolated: at 40 variables the cap of 40 dropped the size-42 tier and left
# a hull of 36. A tier within two sizes of the cap needs no top-up, because a
# mediating block excludes the edge's own two endpoints and so never exceeds
# cap - 2.
#
# @param jobs  Filtered anchor grid, columns n (size) and d (density).
# @param cap   Size cap for this build.
# @param dens  Densities to place at the cap tier.
#
# Returns: The grid, with a cap tier appended when one is needed.
# ------------------------------------------------------------------------------
zratio_cap_tier = function(jobs, cap, dens) {
  if(nrow(jobs) == 0L || cap - max(jobs$n) <= 2) return(jobs)
  rbind(jobs, expand.grid(n = as.numeric(cap), d = dens))
}

# Block-Gibbs sweeps per anchor, by anchor size. Anchor cost scales ~ n^3 *
# sweeps, so the large tiers run shorter chains; the fit denoises across
# anchors, and the measured hull accuracy at the top tiers is unaffected.
zratio_anchor_sweeps = function(n) {
  ifelse(n >= 46, 600L, ifelse(n >= 32, 800L, 1000L))
}

# Gamma-shape deployment range for the absolute-moment surface. Outside it the
# engine keeps the additive path. The surface is validated at three shapes --
# 0.5, 1, and 2, each scored against block-Gibbs gold at eta 1 and 2 on both
# component families -- and deploys on the interval they span; the interior is
# interpolated, not measured. The upper end is set by the anchor oracle, not by
# the surface: at a non-unit shape the oracle's row update is an
# independence-Metropolis step whose acceptance falls away from shape 1
# (measured at 20 nodes, density 0.9: 80% at shape 2, 70% at 2.5, 60% at 3,
# 35% at 4, 1-12% at 5), and with it the anchor Monte-Carlo error rises out of
# reach of any affordable sweep budget (at shape 5, 50-500x the shape-1 anchor
# error and not restored by 94x the sweeps).
#
# This pair is the single owner of the deployment policy: the C++ gate trusts
# whether a surface was attached and does not re-derive the range (see
# ZRatioEngine::log_zratio). The range now runs to 10. It is accuracy-validated
# against block-Gibbs gold at shapes 0.5, 1, 2, 3 and 5 -- every interior cell
# inside the 0.003-nat envelope -- and carries a different guarantee at 10,
# where the whole mediated correction is bounded by 2.8e-04 nats over the
# scored band at eta <= 2, so any method returning the isolated-edge value is
# wrong by at most that. The interior of the range is interpolated, not
# measured, at both ends.
.zratio_surface_shape_lo = 0.5
.zratio_surface_shape_hi = 10

# Largest diagonal rate the isolated-edge routing below is measured at. Every
# gate in this program was scored at eta 1 and 2, and mediation grows with eta
# (measured at shape 10: eta 2 exceeds eta 1 by ~3.3x, so the bound is set at
# the eta-2 end and does not transfer upward by argument). Past this the same
# route deploys -- the additive alternative is no better there and one rule is
# better than two -- but the bound is not claimed, and the user is told so.
.zratio_mediation_off_eta_hi = 2

# ------------------------------------------------------------------------------
# zratio_mediation_off
# ------------------------------------------------------------------------------
# TRUE when the fit deploys the isolated-edge ratio log(psi0) with the mediating
# correction switched off: the route past the top of the surface's validated
# shape range.
#
# The alternative there is the additive saddle, which is not a coarser version
# of the correction but a broken one on common-neighbour mediating blocks: it
# returns essentially zero and discards the entire log-ratio (MEASURED at
# shape 10, eta 2, k = 20..42: additive error 0.0478 against a gold of 0.0478,
# 16x outside the 0.003-nat envelope). Serving log(psi0) instead is exact for an
# unmediated edge and drops only the mediation, and past shape 10 the Gamma
# diagonal concentrates the precision diagonal until mediation dies:
#
#   MEASURED max |gold - log(psi0)|, eta 2, block-Gibbs gold
#     shape 10   2.84e-04   (22 cells: both families, both eta, sizes 20-100)
#     shape 12   2.18e-04   (k = 42 and k = 100 at the measured shape-12 bump)
#     shape 15   1.37e-04   (k = 42)
#     shape 20   6.22e-05   (k = 42)
#
# so the whole error of this route is at most 2.84e-04 nats, two orders below
# the 0.003-nat envelope every accuracy claim in this program lives in. The
# decay is NOT pointwise monotone in shape (shape 12 sits above shape 10), so
# the claim rests on the measured band maximum and not on a monotonicity
# argument. The constants feeding psi0 are certified to 4.3e-08 through shape 20
# (nleg = 128; R/zratio_tables.R).
#
# One rule for both component families: mediation is negligible for both at
# these shapes, and two fallback routes would be interface surface with nothing
# to buy. Below .zratio_surface_shape_lo the additive path is unchanged -- the
# mediation bound is a large-shape phenomenon and does not apply there.
zratio_mediation_off = function(zc) {
  isTRUE(zc$alpha > .zratio_surface_shape_hi)
}

# Sweep multiplier restoring the shape-1 anchor Monte-Carlo error at a non-unit
# shape, resolved by matching measured across-seed anchor spread rather than by
# the 1/acceptance heuristic (which understates the cost, since a rejected row
# repeats the previous state and leaves autocorrelation behind).
zratio_anchor_shape_multiplier = function(alpha) {
  if(abs(alpha - 1) < 1e-12) 1L else if(alpha < 1) 4L else 2L
}

# ------------------------------------------------------------------------------
# zratio_anchor_grids
# ------------------------------------------------------------------------------
# The (size, density, sweeps) anchor grids for both families at one size cap.
# Sparse, law-informed placement: dense high-size tiers where the components the
# engine meets live, a low-density tail at small sizes, and two replicates
# throughout (the low-order fit denoises the Monte-Carlo noise across anchors).
#
# @param cap  Largest anchor size for this build.
#
# Returns: list(cn = <grid>, bip = <grid>), columns n, d, sweeps.
# ------------------------------------------------------------------------------
zratio_anchor_grids = function(cap) {
  cn = rbind(
    expand.grid(n = c(4, 6, 8, 10, 12, 15, 18, 22, 26, 30), d = c(0.7, 0.8, 0.9, 1.0)),
    expand.grid(n = c(3, 4),                                 d = c(0.5, 0.7, 0.85, 1.0)),
    expand.grid(n = c(4, 6, 8, 10, 12, 15),                  d = c(0.35, 0.5)),
    expand.grid(n = c(36, 42),                               d = c(0.8, 0.9)),
    expand.grid(n = c(52, 64, 80),                           d = c(0.8, 0.9))
  )
  cn = zratio_cap_tier(cn[cn$n <= cap, , drop = FALSE], cap, dens = c(0.8, 0.9))
  cn = rbind(cn, cn)                                         # 2 reps

  bip = rbind(
    expand.grid(n = c(4, 6, 8, 10, 12, 14, 16, 18, 20, 22), d = c(0.55, 0.7, 0.85, 1.0)),
    expand.grid(n = c(30, 38, 52, 64, 80),                  d = c(0.7, 1.0))
  )
  bip = zratio_cap_tier(bip[bip$n <= cap, , drop = FALSE], cap, dens = c(0.7, 1.0))
  bip = rbind(bip, bip)

  # One sweeps rule for both families keeps the top tier affordable wherever
  # the cap places it.
  cn$sweeps = zratio_anchor_sweeps(cn$n)
  bip$sweeps = zratio_anchor_sweeps(bip$n)

  list(cn = cn, bip = bip)
}

# Trained size-hull cap for the anchor build. Components larger than this are
# extended along the surface's boundary slope at deploy, which is accurate but
# unanchored, so the cap sets where measured accuracy ends: at 80 the reachable
# giant of an ordinary large fit sits inside the hull, and the one-time build
# stays at 24 s serial (2 s on four cores at the small caps a modest fit uses).
# The build default and every sampler call site size through this one value, so
# the cap cannot drift between them.
.zratio_surface_size_cap = 80L

zratio_build_surfaces = function(zc, max_size = .zratio_surface_size_cap,
                                 cores = 1L, seed0 = 700000L) {
  if(zc$alpha < .zratio_surface_shape_lo ||
    zc$alpha > .zratio_surface_shape_hi) {
    return(NULL)
  }
  cap = as.integer(max_size)

  # Get-or-build: the build is data-independent, so a repeat fit of the same
  # cell returns the cached surface (session memory first, then disk) instead of
  # re-running the anchor sweeps. Disable with
  # options(bgms.zratio_surface_cache = FALSE); it also follows
  # options(bgms.correction_table_cache) so one switch governs both.
  use_cache = isTRUE(getOption(
    "bgms.zratio_surface_cache",
    correction_cache_enabled()
  ))
  cache_file = NULL
  if(use_cache) {
    key = zratio_surface_cache_key(zc, cap, seed0)
    hit = get0(key, envir = .zratio_surface_cache, inherits = FALSE)
    if(!is.null(hit)) return(hit)
    cache_dir = correction_cache_dir()
    cache_file = file.path(cache_dir, paste0(key, ".rds"))
    if(file.exists(cache_file)) {
      surf = tryCatch(readRDS(cache_file), error = function(e) NULL)
      if(!is.null(surf)) {
        assign(key, surf, envir = .zratio_surface_cache)
        return(surf)
      }
    }
  }

  # The anchors self-seed per job (zratio_rand_conn_graph/bip call set.seed),
  # which runs in-process at the default cores = 1 and would otherwise leave the
  # caller's .Random.seed advanced after every fit. Save and restore it, matching
  # zratio_constants (R/zratio_tables.R).
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

  grids = zratio_anchor_grids(cap)
  # A non-unit shape samples its anchors through an independence-Metropolis
  # step, so the same nominal budget buys fewer effective sweeps.
  shape_mult = zratio_anchor_shape_multiplier(zc$alpha)
  cn_jobs = grids$cn
  bip_jobs = grids$bip
  cn_jobs$sweeps = as.integer(cn_jobs$sweeps * shape_mult)
  bip_jobs$sweeps = as.integer(bip_jobs$sweeps * shape_mult)

  # One scheduling pool over both families, heaviest job first with dynamic
  # assignment (mc.preschedule = FALSE): anchor cost scales ~ n^3 * sweeps and
  # spans ~1000x across the grid, so static per-worker chunks leave cores idle
  # behind the giants, and a separate bip pass cannot fill the CN tail. Each
  # job keeps the seed it had in its own family grid and the anchor rows are
  # reassembled in grid order, so the fitted surfaces are identical for any
  # core count and any schedule.
  cn_jobs$fam = "cn"
  cn_jobs$seed = seed0 + seq_len(nrow(cn_jobs))
  bip_jobs$fam = "bip"
  bip_jobs$seed = seed0 + 100000L + seq_len(nrow(bip_jobs))
  jobs = rbind(cn_jobs, bip_jobs)
  ord = order(-(as.numeric(jobs$n)^3 * jobs$sweeps))
  job_fun = function(k) {
    if(jobs$fam[k] == "cn") {
      zratio_anchor_cn(jobs$n[k], jobs$d[k], zc, jobs$sweeps[k], 200L, jobs$seed[k])
    } else {
      zratio_anchor_bip(jobs$n[k], jobs$d[k], zc, jobs$sweeps[k], 200L, jobs$seed[k])
    }
  }
  # Fork where available; on Windows (no fork) use a socket cluster instead,
  # but only for builds big enough to repay the ~1-2s worker launch -- below
  # the size-16 grid tier the whole serial build is cheaper than the cluster
  # start. Both schedules are dynamic and jobs self-seed, so every path
  # returns identical surfaces. options(bgms.zratio_surface_psock) forces the
  # cluster branch on or off (tests; forking-hostile unix hosts).
  psock = getOption("bgms.zratio_surface_psock", NULL)
  use_psock = if(is.null(psock)) {
    .Platform$OS.type != "unix" && cores > 1L && max(jobs$n) >= 16
  } else {
    isTRUE(psock) && cores > 1L
  }
  # The cache missed, so the surfaces are about to be built for real; announce
  # the one-time cost so the pre-chain pause is not silent. n_workers is the
  # parallelism actually used, which drops to 1 on the Windows serial tier.
  n_workers = if(use_psock || .Platform$OS.type == "unix") cores else 1L
  verbose = isTRUE(getOption("bgms.verbose", TRUE))
  if(verbose) {
    message(
      "Building normalizing-constant corrections for the hierarchical ",
      "precision prior (", cap, " variables, ", n_workers,
      if(n_workers == 1L) " core)." else " cores)."
    )
  }
  t0 = proc.time()[["elapsed"]]
  res = vector("list", nrow(jobs))
  if(use_psock) {
    cl = parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    res[ord] = parallel::parLapplyLB(cl, ord, job_fun)
  } else {
    mc = if(.Platform$OS.type == "unix") cores else 1L
    res[ord] = parallel::mclapply(
      ord, job_fun,
      mc.cores = mc, mc.preschedule = FALSE
    )
  }
  cn_rows = do.call(rbind, res[jobs$fam == "cn"])
  bip_rows = do.call(rbind, res[jobs$fam == "bip"])

  if(is.null(cn_rows) || is.null(bip_rows)) return(NULL)
  cn = zratio_fit_surface_family(cn_rows)
  bip = zratio_fit_surface_family(bip_rows)
  if(is.null(cn) || is.null(bip)) return(NULL)
  surf = list(cn = cn, bip = bip)
  if(verbose) {
    secs = proc.time()[["elapsed"]] - t0
    message(
      "Correction build complete (",
      if(secs < 1) "< 1s" else paste0(round(secs), "s"), ")."
    )
  }
  if(use_cache) {
    assign(key, surf, envir = .zratio_surface_cache)
    tryCatch({
      dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
      saveRDS(surf, cache_file)
    }, error = function(e) NULL)
  }
  surf
}

# Parallelism for the one-time surface build. Anchors self-seed per job, so the
# result is independent of the core count; only the wall time changes. Defaults
# to the fit's own `cores`: the build runs before the chains launch, so those
# cores are idle during exactly this window. options(bgms.zratio_surface_cores)
# overrides. Unix parallelizes by forking; Windows by a socket cluster on large
# builds (small ones run serially there -- see zratio_build_surfaces).
zratio_surface_build_cores = function(fit_cores = 1L) {
  fallback = suppressWarnings(as.integer(fit_cores))
  if(length(fallback) != 1L || is.na(fallback) || fallback < 1L) fallback = 1L
  cores = suppressWarnings(as.integer(
    getOption("bgms.zratio_surface_cores", fallback)
  ))
  if(length(cores) != 1L || is.na(cores) || cores < 1L) cores = 1L
  cores
}

# Message the route taken when no surface is attached. Three cases, and they
# are different claims, so they get different wordings: a shape past the top of
# the validated range (isolated-edge routing, bounded), a shape below it
# (additive path, unchanged), or a failed build inside the range (which must not
# downgrade the fit silently). Shared by every sampler call site so the wordings
# cannot drift apart.
zratio_surface_fence_message = function(zc) {
  if(zratio_mediation_off(zc)) {
    message(
      "z-ratio: precision shape alpha = ", format(zc$alpha),
      " is past the absolute-moment surface's validated range (up to shape ",
      format(.zratio_surface_shape_hi),
      "), so the mediating correction is switched off and every edge gets the ",
      "isolated-edge ratio. At these shapes the Gamma diagonal concentrates ",
      "the precision diagonal and the whole mediated correction is at most ",
      "0.00028 nats, measured against a block-Gibbs reference at shapes 12, ",
      "15 and 20; that is the entire error of this route, and it is two ",
      "orders below the 0.003 nats the surface is claimed to within inside ",
      "its range."
    )
    if(zc$eta > .zratio_mediation_off_eta_hi) {
      message(
        "z-ratio: that bound was measured at diagonal rates eta up to ",
        format(.zratio_mediation_off_eta_hi), " and this fit runs at eta = ",
        format(zc$eta),
        ". Mediation grows with eta, so the same route deploys but the ",
        "measured bound does not cover this cell."
      )
    }
  } else if(zc$alpha < .zratio_surface_shape_lo) {
    message(
      "z-ratio: precision shape alpha = ", format(zc$alpha),
      " -> additive path (coarser correction). The absolute-moment surface is ",
      "scored against a block-Gibbs reference at shapes 0.5, 1, 2, 3 and 5, ",
      "and deploys on the range those points span up to shape ",
      format(.zratio_surface_shape_hi),
      "; the interior of that range is interpolated, not measured. Below it ",
      "the surface is unscored, and the additive path that serves instead is ",
      "measurably coarse on common-neighbour mediating blocks."
    )
  } else {
    message(
      "z-ratio: the absolute-moment surface build failed -> additive ",
      "path (coarser correction; the trust gauge quantifies the impact in ",
      "fit$zratio_diag)."
    )
  }
}

# ------------------------------------------------------------------------------
# zratio_spec_list
# ------------------------------------------------------------------------------
# The `zratio` spec one fit hands the sampler: the cell's constants, its
# (delta, eta, alpha, slab) identity, the trust-gauge sweep count, and the
# routing R resolved. The C++ side reads it in exactly one place
# (zratio_engine_from_spec, src/models/ggm/zratio_engine.h), which the GGM
# sampler, the mixed sampler, and the deployed-route test entry all share.
#
# `mediation_off` travels with the constants rather than being re-derived
# downstream: R owns the deployment policy, the engine deploys what it is
# handed. The one time that rule was broken -- the engine keeping its own copy
# of the shape range -- R widened the range, C++ did not, and every non-unit
# shape silently took the wrong branch.
#
# @param zc            Cell constants from zratio_cell_constants.
# @param gauge_sweeps  In-chain trust-gauge assessment sweeps (0 = off).
#
# Returns: the spec list, before any surface is attached.
# ------------------------------------------------------------------------------
zratio_spec_list = function(zc, gauge_sweeps) {
  list(
    addc = zc$addc, tg = zc$tg, ihat = zc$ihat, ghat = zc$ghat,
    wt = zc$wt, psi0 = zc$psi0,
    delta = zc$delta, eta = zc$eta, alpha = zc$alpha, slab = zc$slab,
    gauge_sweeps = gauge_sweeps,
    mediation_off = zratio_mediation_off(zc)
  )
}

# Build the Option-B surface for cell `zc`, sized on `size` variables, and, on a
# successful build, attach it to the `zratio` spec; otherwise keep whatever
# route the spec already carries -- isolated-edge past the validated shape
# range, additive elsewhere -- and message the reason when `verbose`.
# Centralizes the size cap, cores policy, and fence message the GGM, mixed, and
# prior sampler paths share. Returns the (possibly surface-carrying) `zratio`
# list.
zratio_attach_surface = function(zratio, zc, size, cores, verbose = FALSE) {
  surf = zratio_build_surfaces(
    zc,
    max_size = min(size, .zratio_surface_size_cap),
    cores = cores
  )
  if(!is.null(surf)) {
    zratio$surface = surf
  } else if(isTRUE(verbose)) {
    zratio_surface_fence_message(zc)
  }
  zratio
}

# Number of in-chain trust-gauge assessment sweeps. The gauge is a post-sampling
# diagnostic (chain_runner.cpp) and runs by default; options(bgms.zratio_gauge_
# sweeps = 0L) is the off switch. Its cost is fixed per chain (two sweeps, each
# referencing a capped number of edge moves), so it does not scale with iter:
# nothing on a sparse posterior, where no mediating block is non-trivial and the
# ratio is exact, and seconds on a dense large-q one. The prior sampler wires its
# own flag (sample_ggm_prior); this governs the deployed hierarchical path.
zratio_gauge_sweeps = function() {
  n = suppressWarnings(as.integer(getOption("bgms.zratio_gauge_sweeps", 2L)))
  if(length(n) != 1L || is.na(n) || n < 0L) n = 0L
  n
}
