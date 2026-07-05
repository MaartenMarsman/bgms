# ==============================================================================
# Z-ratio diagnostics (hierarchical GGM prior specification)
# ==============================================================================
#
# Post-sampling alarm suite for the per-edge Z-ratio kernel: measurement-only
# oracle audits on visited graphs (targeted, random, additive-zone channels),
# calibration-stream drift checks, regime flags, and the per-chain verdict.
# Reference protocol: SV/z_graph_prior notes/bgms-warmup-alarms-spec.md Sec 4-5
# and the measured operating characteristics in w5c3.
# ==============================================================================


# ------------------------------------------------------------------------------
# zratio_eta
# ------------------------------------------------------------------------------
# Standardized prior-scale regime index: eta = sigma * beta in the bare frame
# (sigma = 1, beta = eta after standardization). Drives the verdict threshold.
#
# @param zratio_spec  Z-ratio spec list with bare-scale sigma and beta.
#
# Returns: Numeric scalar eta.
# ------------------------------------------------------------------------------
zratio_eta = function(zratio_spec) {
  zratio_spec$sigma * zratio_spec$beta
}

# ------------------------------------------------------------------------------
# zratio_tau
# ------------------------------------------------------------------------------
# Verdict threshold per regime: 0.01 at eta = 1, 0.02 at eta = 2, 0.04 at
# eta >= 3 (healthy boundary cells read up to ~0.027 at eta = 3).
#
# @param eta  Regime index from zratio_eta().
#
# Returns: Numeric scalar threshold on the audit error.
# ------------------------------------------------------------------------------
zratio_tau = function(eta) {
  if(eta >= 3) 0.04 else if(eta >= 2) 0.02 else 0.01
}

# ------------------------------------------------------------------------------
# zratio_drift_flag
# ------------------------------------------------------------------------------
# End-of-warmup drift check on a selection-warm stream trace: compare the
# first-half and second-half means against twice the binomial half-width of
# a single graph's density around its expectation.
#
# @param trace    Numeric vector (density or theta over the selection-warm
#   warmup iterations).
# @param n_pairs  Number of candidate edges p(p-1)/2.
#
# Returns: Logical scalar; NA when the trace is shorter than 40 points.
# ------------------------------------------------------------------------------
zratio_drift_flag = function(trace, n_pairs) {
  n = length(trace)
  if(n < 40) {
    return(NA)
  }
  first_half = trace[seq_len(floor(n / 2))]
  second_half = trace[seq.int(floor(n / 2) + 1L, n)]
  center = mean(trace)
  half_width = 2 * sqrt(max(center * (1 - center), 1e-12) / n_pairs)
  abs(mean(second_half) - mean(first_half)) > half_width
}

# ------------------------------------------------------------------------------
# zratio_indicator_graph
# ------------------------------------------------------------------------------
# Rebuild the p x p adjacency from one column of the vectorized indicator
# samples (upper triangle including the diagonal, row-major).
#
# @param column     Integer vector of length p(p+1)/2.
# @param num_nodes  Number of nodes p.
#
# Returns: Integer p x p adjacency matrix with unit diagonal.
# ------------------------------------------------------------------------------
zratio_indicator_graph = function(column, num_nodes) {
  G = matrix(0L, num_nodes, num_nodes)
  e = 1L
  for(i in seq_len(num_nodes)) {
    for(j in i:num_nodes) {
      G[i, j] = column[e]
      G[j, i] = column[e]
      e = e + 1L
    }
  }
  diag(G) = 1L
  G
}

# ------------------------------------------------------------------------------
# zratio_audit_chain
# ------------------------------------------------------------------------------
# Alarm channels for one chain: scan visited graphs for mediating-block
# descriptors, pick audit candidates per targeting rule, and score each
# against the block-Gibbs local oracle.
#
# Channels (spec Sec 4): targeted picks by anchor distance, leverage, and
# prediction magnitude (fit-bearing runs); random visitation-weighted picks
# (the gate without a fit); additive-zone picks (maxbd < 2 blocks scored as
# |0 - oracle|, random plus highest-m). A1a out-of-hull fraction, m32
# closure reach, and the visited-density band are regime context.
#
# @param chain        One chain output holding indicator_samples and zratio.
# @param zratio_spec  Z-ratio spec list (tables, bare-scale constants).
# @param num_nodes    Number of nodes p.
# @param n_graphs     Visited graphs to scan.
# @param top_k        Picks per targeted rule.
# @param rand_k       Random picks (entangled and additive-zone channels).
# @param audit_sweep  Oracle block-Gibbs sweeps per audited block.
# @param seed         Seed for graph sampling, random picks, and the oracle.
#
# Returns: List with the channel maxima/medians, regime fractions, pick
#   count, and the per-pick audit table.
# ------------------------------------------------------------------------------
zratio_audit_chain = function(chain, zratio_spec, num_nodes, n_graphs, top_k,
                              rand_k, audit_sweep, seed) {
  zr = chain$zratio
  addc = as.numeric(zr$addc)
  have_fit = length(addc) >= 23 && addc[13] > 0.5
  inds = chain$indicator_samples
  n_iter = ncol(inds)
  set.seed(seed)
  cols = sort(sample.int(n_iter, min(n_graphs, n_iter)))

  scan = vector("list", length(cols))
  graphs = vector("list", length(cols))
  for(g in seq_along(cols)) {
    G = zratio_indicator_graph(inds[, cols[g]], num_nodes)
    graphs[[g]] = G
    rows = zratio_scan_graph(
      G, addc, zratio_spec$tg, zratio_spec$ihat, zratio_spec$ghat,
      zratio_spec$wt, zratio_spec$psi0
    )
    if(nrow(rows) > 0) {
      scan[[g]] = cbind(graph = g, rows)
    }
  }
  scan = do.call(rbind, scan)
  empty = list(
    aud_targeted_max = NA_real_, aud_rand_max = NA_real_,
    aud_addz_max = NA_real_, aud_targeted_med = NA_real_,
    a1a_out_frac = NA_real_, m32_frac = NA_real_, dens_out_frac = NA_real_,
    n_audit = 0L, picks = NULL
  )
  if(is.null(scan) || nrow(scan) == 0) {
    return(empty)
  }
  colnames(scan) = c(
    "graph", "i", "j", "ncn", "cne", "bre", "maxbd", "m", "dens", "pred",
    "clamped"
  )

  entangled = scan[scan[, "maxbd"] >= 2, , drop = FALSE]
  additive_zone = scan[scan[, "maxbd"] < 2 & scan[, "m"] >= 2, , drop = FALSE]

  # Regime context: out-of-hull fraction (trigger only), closure reach, and
  # the visited densities against the calibration stream's band.
  a1a = if(have_fit && nrow(entangled) > 0) {
    mean(entangled[, "clamped"])
  } else {
    NA_real_
  }
  m32 = if(nrow(entangled) > 0) mean(entangled[, "m"] > 32) else NA_real_
  dens_out = NA_real_
  warm_dens = as.numeric(zr$warmup_density)
  if(length(warm_dens) > 0) {
    n_pairs = num_nodes * (num_nodes - 1) / 2
    visited = vapply(
      graphs, function(G) mean(G[upper.tri(G)]), numeric(1)
    )
    center = mean(warm_dens)
    half_width = 2 * sqrt(max(center * (1 - center), 1e-12) / n_pairs)
    dens_out = mean(
      visited < min(warm_dens) - half_width |
        visited > max(warm_dens) + half_width
    )
  }

  # Targeting rules over the entangled blocks.
  picks = NULL
  pick_rows = function(rows, rule) {
    if(is.null(rows) || nrow(rows) == 0) {
      return(NULL)
    }
    data.frame(
      graph = rows[, "graph"], i = rows[, "i"], j = rows[, "j"], rule = rule
    )
  }
  if(nrow(entangled) > 0) {
    if(have_fit) {
      ax = zr$anchors_x
      feats = entangled[, c("bre", "m", "cne", "maxbd", "dens"), drop = FALSE]
      if(!is.null(ax) && nrow(ax) >= 2) {
        anchor_feats = ax[, -1, drop = FALSE]
        sds = apply(anchor_feats, 2, stats::sd)
        sds[!is.finite(sds) | sds < 1e-9] = 1
        scaled_anchors = sweep(anchor_feats, 2, sds, "/")
        dq = apply(sweep(feats, 2, sds, "/"), 1, function(z) {
          min(sqrt(rowSums(sweep(scaled_anchors, 2, z, "-")^2)))
        })
        xtx_inv = solve(crossprod(ax) + diag(1e-8, ncol(ax)))
        design = cbind(1, feats)
        lev = rowSums((design %*% xtx_inv) * design)
        picks = rbind(
          picks,
          pick_rows(
            entangled[order(-dq)[seq_len(min(top_k, nrow(entangled)))], ,
              drop = FALSE
            ],
            "dist"
          ),
          pick_rows(
            entangled[order(-lev)[seq_len(min(top_k, nrow(entangled)))], ,
              drop = FALSE
            ],
            "lev"
          ),
          pick_rows(
            entangled[
              order(-abs(entangled[, "pred"]))[
                seq_len(min(top_k, nrow(entangled)))
              ], ,
              drop = FALSE
            ],
            "pred"
          )
        )
      }
    }
    rand_idx = sample.int(nrow(entangled), min(rand_k, nrow(entangled)))
    picks = rbind(
      picks, pick_rows(entangled[rand_idx, , drop = FALSE], "rand")
    )
  }
  if(nrow(additive_zone) > 0) {
    addz_rand = sample.int(
      nrow(additive_zone), min(rand_k, nrow(additive_zone))
    )
    addz_top = order(-additive_zone[, "m"])[
      seq_len(min(top_k, nrow(additive_zone)))
    ]
    picks = rbind(
      picks,
      pick_rows(
        additive_zone[unique(c(addz_rand, addz_top)), , drop = FALSE], "addz"
      )
    )
  }
  if(is.null(picks) || nrow(picks) == 0) {
    return(empty)
  }

  # One audit call per visited graph covers all its picked edges; identical
  # (graph, edge) picks from different rules share one oracle run.
  key = paste(picks$graph, picks$i, picks$j)
  first = !duplicated(key)
  err = rep(NA_real_, nrow(picks))
  for(g in unique(picks$graph[first])) {
    sel = first & picks$graph == g
    edges = matrix(as.integer(c(picks$i[sel], picks$j[sel])), ncol = 2)
    audit = zratio_audit_edges(
      graphs[[g]], edges,
      addc, zratio_spec$tg, zratio_spec$ihat, zratio_spec$ghat,
      zratio_spec$wt, zratio_spec$psi0,
      zratio_spec$delta, zratio_spec$sigma, zratio_spec$beta,
      as.integer(audit_sweep), 30L, as.integer(seed + g)
    )
    err[sel] = audit$err
  }
  err[!first] = err[which(first)[match(key[!first], key[first])]]

  channel_max = function(rules) {
    e = err[picks$rule %in% rules]
    e = e[is.finite(e)]
    if(length(e)) max(e) else NA_real_
  }
  targeted = err[picks$rule %in% c("dist", "lev", "pred")]
  targeted = targeted[is.finite(targeted)]
  list(
    aud_targeted_max = channel_max(c("dist", "lev", "pred")),
    aud_rand_max = channel_max("rand"),
    aud_addz_max = channel_max("addz"),
    aud_targeted_med = if(length(targeted)) {
      stats::median(targeted)
    } else {
      NA_real_
    },
    a1a_out_frac = a1a, m32_frac = m32, dens_out_frac = dens_out,
    n_audit = sum(first),
    picks = cbind(picks, err = err)
  )
}

#' @title Summarize Z-Ratio Kernel Diagnostics
#'
#' @description Alarm suite for the hierarchical prior specification's
#' per-edge Z-ratio kernel. For each chain it audits the frozen kernel on
#' graphs the chain visited — targeted picks (anchor distance, leverage,
#' prediction magnitude), random picks, and additive-zone picks, each scored
#' as the absolute gap between the deployed correction and a measurement-only
#' block-Gibbs oracle — checks the calibration stream for end-of-warmup
#' drift, and reports regime context (out-of-hull fraction, closure reach,
#' visited-density band). The verdict flags a chain when the maximum of its
#' targeted (or, without a fit, random) and additive-zone audit errors
#' exceeds the regime threshold: 0.01 at \eqn{\eta = 1}, 0.02 at
#' \eqn{\eta = 2}, 0.04 at \eqn{\eta \ge 3}, with
#' \eqn{\eta = \sigma \beta} in the bare prior scale. A quiet verdict bounds
#' the pointwise approximation error on visited graphs; coherent sub-margin
#' bias that accumulates through inclusion-probability feedback (prior-only
#' chains at large \eqn{p}) is outside its reach.
#'
#' @param chains List of per-chain sampler outputs, each holding
#'   \code{indicator_samples} and the \code{zratio} diagnostics block
#'   (constant block, anchors, counters, warmup traces) that
#'   \code{sample_ggm} attaches when the hierarchical specification is
#'   active.
#' @param zratio_spec The Z-ratio specification list used for the run:
#'   quadrature tables \code{tg}, \code{ihat}, \code{ghat}, \code{wt},
#'   isolated-edge ratio \code{psi0}, and the bare-scale constants
#'   \code{delta}, \code{sigma}, \code{beta}.
#' @param num_nodes Integer: number of nodes \eqn{p}.
#' @param n_graphs Integer: visited graphs scanned per chain (default 12).
#' @param top_k Integer: picks per targeted rule and additive-zone top-m
#'   picks (default 6).
#' @param rand_k Integer: random picks per channel (default 6).
#' @param audit_sweep Integer: block-Gibbs sweeps per audited block
#'   (default 600).
#' @param seed Integer seed for graph sampling, random picks, and the audit
#'   oracle (default 1).
#' @param verbose Logical: print detected issues (default \code{TRUE}).
#'   Quiet chains print nothing.
#'
#' @return An invisible named list:
#'   \describe{
#'     \item{\code{per_chain}}{Data frame, one row per chain: engine
#'       counters (\code{n_oracle}, \code{n_anchors}, \code{n_pred},
#'       \code{n_add}, \code{n_clamp}, \code{frozen}), audit channel maxima
#'       (\code{aud_targeted_max}, \code{aud_rand_max},
#'       \code{aud_addz_max}), the verdict gate and \code{flagged}, drift
#'       flags (\code{drift_density}, \code{drift_theta}), and regime
#'       context (\code{a1a_out_frac}, \code{m32_frac},
#'       \code{dens_out_frac}).}
#'     \item{\code{eta}}{Regime index \eqn{\sigma \beta}.}
#'     \item{\code{tau}}{Verdict threshold at this \code{eta}.}
#'     \item{\code{verdict_flagged}}{Logical: any chain flagged.}
#'     \item{\code{calibration_incomplete}}{Logical: any chain's
#'       density or theta stream still drifting at the end of warmup.}
#'     \item{\code{audits}}{List of per-chain audit pick tables
#'       (graph, edge, rule, error).}
#'   }
#'
#' @examples
#' \donttest{
#' draws = sample_ggm_prior(
#'   p = 8, n_samples = 100, n_warmup = 200,
#'   interaction_prior = normal_prior(scale = 0.5),
#'   precision_scale_prior = gamma_prior(shape = 1, rate = 2),
#'   spec = "hierarchical", calibration_window = 50,
#'   verbose = FALSE
#' )
#' draws$zratio_diagnostics$per_chain
#' }
#'
#' @seealso \code{\link{sample_ggm_prior}}
#' @family diagnostics
#' @export
summarize_zratio_diagnostics = function(
  chains,
  zratio_spec,
  num_nodes,
  n_graphs = 12,
  top_k = 6,
  rand_k = 6,
  audit_sweep = 600,
  seed = 1,
  verbose = TRUE
) {
  chains = Filter(function(chain) !is.null(chain$zratio), chains)
  if(length(chains) == 0) {
    stop(
      "No Z-ratio diagnostics found in the chain outputs. They are only ",
      "recorded when the hierarchical prior specification is active."
    )
  }
  eta = zratio_eta(zratio_spec)
  tau = zratio_tau(eta)
  n_pairs = num_nodes * (num_nodes - 1) / 2

  rows = vector("list", length(chains))
  audits = vector("list", length(chains))
  for(c_idx in seq_along(chains)) {
    chain = chains[[c_idx]]
    counters = chain$zratio$counters
    audit = zratio_audit_chain(
      chain, zratio_spec, num_nodes, n_graphs, top_k, rand_k, audit_sweep,
      seed + 1000L * c_idx
    )
    audits[c_idx] = list(audit$picks)
    gate = if(!is.na(audit$aud_targeted_max)) {
      audit$aud_targeted_max
    } else {
      audit$aud_rand_max
    }
    verdict_stat = max(
      c(gate, audit$aud_addz_max)[is.finite(c(gate, audit$aud_addz_max))],
      -Inf
    )
    rows[[c_idx]] = data.frame(
      chain = c_idx,
      n_oracle = counters[["n_oracle"]],
      n_anchors = counters[["n_anchors"]],
      n_pred = counters[["n_pred"]],
      n_add = counters[["n_add"]],
      n_clamp = counters[["n_clamp"]],
      frozen = counters[["frozen"]] > 0,
      aud_targeted_max = audit$aud_targeted_max,
      aud_rand_max = audit$aud_rand_max,
      aud_addz_max = audit$aud_addz_max,
      gate = if(is.finite(verdict_stat)) verdict_stat else NA_real_,
      flagged = is.finite(verdict_stat) && verdict_stat > tau,
      drift_density = zratio_drift_flag(
        as.numeric(chain$zratio$warmup_density), n_pairs
      ),
      drift_theta = zratio_drift_flag(
        as.numeric(chain$zratio$warmup_theta), n_pairs
      ),
      a1a_out_frac = audit$a1a_out_frac,
      m32_frac = audit$m32_frac,
      dens_out_frac = audit$dens_out_frac,
      n_audit = audit$n_audit
    )
  }
  per_chain = do.call(rbind, rows)
  drifting = per_chain$drift_density | per_chain$drift_theta
  calibration_incomplete = any(drifting, na.rm = TRUE)
  verdict_flagged = any(per_chain$flagged)

  if(verbose && isTRUE(getOption("bgms.verbose", TRUE))) {
    issues = character(0)
    if(verdict_flagged) {
      flagged_chains = per_chain$chain[per_chain$flagged]
      issues = c(issues, sprintf(
        "Audit verdict: audit error %.4f > tau %.2f (eta = %g) in chain%s %s - the frozen Z-ratio kernel exceeds its error margin on visited graphs; edge-inclusion results may be biased",
        max(per_chain$gate[per_chain$flagged]), tau, eta,
        if(length(flagged_chains) > 1) "s" else "",
        paste(flagged_chains, collapse = ", ")
      ))
    }
    if(calibration_incomplete) {
      drift_chains = per_chain$chain[which(drifting)]
      issues = c(issues, sprintf(
        "Calibration incomplete: graph density or theta still drifting at the end of warmup in chain%s %s - increase warmup or the calibration window",
        if(length(drift_chains) > 1) "s" else "",
        paste(drift_chains, collapse = ", ")
      ))
    }
    m32_hot = !is.na(per_chain$m32_frac) & per_chain$m32_frac > 0.5
    if(any(m32_hot)) {
      issues = c(issues, sprintf(
        "Regime flag: over half the coupled-bridge blocks exceed m = 32 in chain%s %s - the two-moment closure degrades at this density; interpret with the additive-zone audit",
        if(sum(m32_hot) > 1) "s" else "",
        paste(per_chain$chain[m32_hot], collapse = ", ")
      ))
    }
    if(length(issues) > 0) {
      cat("Z-ratio issues:\n")
      for(issue in issues) {
        cat("  -", issue, "\n")
      }
    }
  }

  invisible(list(
    per_chain = per_chain,
    eta = eta,
    tau = tau,
    verdict_flagged = verdict_flagged,
    calibration_incomplete = calibration_incomplete,
    audits = audits
  ))
}
