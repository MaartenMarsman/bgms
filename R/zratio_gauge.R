#' @title Summarize the Hierarchical Prior Trust Gauge
#'
#' @description Reports the per-chain trust gauge for the hierarchical graph
#' prior. Under that prior the sampler decides each edge with a fast
#' approximation; in a set of assessment sweeps after sampling the gauge redoes
#' a subset of each chain's own edge decisions with the exact calculation and
#' records two statistics per chain:
#' \describe{
#'   \item{\code{flip_rate}}{The fraction of add/remove decisions that would
#'     come out differently under the exact calculation. A chain is flagged on
#'     this channel when \code{flip_rate} exceeds the tolerance by more than
#'     the exact reference's own Monte-Carlo noise. This channel detects error
#'     that changed decisions the chain actually made; it is insensitive to a
#'     small coherent error at chains whose decisions are far from their
#'     accept/reject boundaries.}
#'   \item{\code{harm_pred}}{The projected distortion of the mean posterior
#'     inclusion probability under the measured approximation error, a
#'     first-order (linear-response) quantity targeted at coherent error:
#'     \code{harm_pred = |mean(m_e s_e)| * A}, where \code{s_e} is the signed
#'     log-ratio error of audited edge \code{e} against the exact reference,
#'     \code{m_e = p_e (1 - p_e)} is that edge's inclusion sensitivity, and
#'     \code{A = 1 / (1 - g)} is the inclusion-probability feedback
#'     amplification with linearized gain
#'     \code{g = E m / (theta (1 - theta) (a + b + E))},
#'     \code{m = mean(p_e (1 - p_e))} over all edges, under a Beta-Bernoulli
#'     edge prior (for a fixed inclusion probability \code{g = 0}, so
#'     \code{A = 1}). When the per-pair audit stream is unavailable the
#'     unweighted form \code{|se_mean| * m * A} is used. A chain is flagged
#'     on this channel when \code{harm_pred} exceeds the tolerance and the
#'     weighted error is resolved above twice its standard error (per-edge
#'     cluster-robust spread plus the reference Monte-Carlo noise). This
#'     channel detects coherent error whose equilibrium effect exceeds the
#'     tolerance even when no individual decision visibly flips; it targets
#'     the mean-inclusion shift and does not bound edge-specific distortions
#'     that cancel in the mean. It is computed for Bernoulli and
#'     Beta-Bernoulli edge priors; under other priors it is \code{NA}.}
#' }
#'
#' The gauge reads the \code{zratio$gauge} block that the sampler attaches
#' under the hierarchical prior; it does not re-scan the stored draws.
#'
#' @details
#' The gauge audits only the edge moves whose mediating block is non-trivial
#' (two or more variables). On a sparse posterior no block reaches that size,
#' the normalizer ratio is exact, and the gauge reports nothing because there
#' is nothing to audit: silence there is exactness, not blindness. The cost
#' follows the same rule -- it is zero where the approximation is exact and
#' grows with the size of the blocks the chain actually visits.
#'
#' The audit is a sample. Each sweep references a capped number of edge moves,
#' so on a dense large graph a few tens of moves stand in for tens of
#' thousands of non-trivial ones (\code{n_ref} against \code{n_ent}, with the
#' remainder counted in \code{n_capped}). That sample resolves coherent error
#' -- error with a consistent sign across edges, which is what
#' \code{harm_pred} projects onto the inclusion-probability scale and what
#' shifts a recovered network. It does not resolve rare edge-specific
#' failures: \code{flip_rate} is the flip rate among audited decisions, not a
#' per-edge guarantee over all of them.
#'
#' @param chains List of per-chain sampler outputs, each carrying a
#'   \code{zratio$gauge} block (\code{flip_rate}, \code{noise_floor},
#'   \code{se_mean}, \code{se_sd}, \code{se_mcse}, \code{n_ent}, \code{n_ref},
#'   \code{n_capped}).
#' @param threshold Numeric flag threshold on \code{flip_rate} (default
#'   \code{0.01}).
#' @param verbose Logical: message flagged chains (default \code{TRUE}).
#' @param harm_inputs Optional list enabling the \code{harm_pred} channel:
#'   \code{pip} (a list with one numeric vector of posterior edge-inclusion
#'   probabilities per chain) and, for a Beta-Bernoulli edge prior, its shape
#'   parameters \code{a} and \code{b} (\code{NULL} for a fixed inclusion
#'   probability). When \code{NULL} (default) the \code{harm_pred} columns are
#'   \code{NA} and only the \code{flip_rate} channel flags.
#' @param harm_threshold Numeric flag threshold on \code{harm_pred}, in
#'   inclusion-probability units (default \code{0.01}).
#'
#' @return An invisible named list:
#'   \describe{
#'     \item{\code{per_chain}}{Data frame, one row per chain: \code{flip_rate},
#'       its \code{flag}, the signed mean and spread of the log-ratio error
#'       (\code{se_mean}, \code{se_sd}), the reference-noise component
#'       \code{se_mcse} and combined standard error \code{se_se} of
#'       \code{se_mean}, the reference \code{noise_floor}, the pair counts
#'       (\code{n_ent} non-trivial seen, \code{n_ref} referenced,
#'       \code{n_capped} cap hits), the mediating-block sizes the audit
#'       covered (\code{block_lo}, \code{block_hi}), and the harm channel
#'       (\code{amplification},
#'       \code{kappa} = the predicted mean-inclusion shift per nat of coherent
#'       error, \code{harm_pred}, \code{harm_flag}).}
#'     \item{\code{threshold}}{The flag tolerance on \code{flip_rate}.}
#'     \item{\code{harm_threshold}}{The flag tolerance on \code{harm_pred}.}
#'     \item{\code{flagged}}{Logical: any chain flagged on either channel.}
#'   }
#'
#' @section Fits that route around the correction entirely:
#'   Above the Gamma diagonal shape range the correction is scored on, the
#'   mediating correction is switched off and every edge is served the
#'   isolated-edge normalizer ratio, which is exact for an edge with no
#'   mediating structure. The per-chain \code{counters} vector records this as
#'   \code{n_isolated}, and a fit that took the route says so in a note
#'   reporting the measured bound on what it leaves out. The gauge still runs
#'   there and still measures that residual directly, so a flag on such a fit is
#'   as meaningful as on any other.
#'
#' @examples
#' \donttest{
#' draws = sample_ggm_prior(
#'   p = 8, n_samples = 100, n_warmup = 200,
#'   interaction_prior = normal_prior(scale = 0.5),
#'   precision_scale_prior = gamma_prior(shape = 1, rate = 2),
#'   spec = "hierarchical", verbose = FALSE
#' )
#' draws$zratio_diagnostics$per_chain
#' }
#'
#' @seealso \code{\link{sample_ggm_prior}}
#' @family diagnostics
#' @export
summarize_zratio_gauge = function(chains, threshold = 0.01, verbose = TRUE,
                                  harm_inputs = NULL, harm_threshold = 0.01) {
  # Keep the original chain indices: harm_inputs$pip is positional over ALL
  # chains, so a chain without gauge output (e.g. an interrupt during another
  # chain's sweeps) must not shift the pip alignment of the chains after it.
  keep = which(vapply(
    chains,
    function(ch) !is.null(ch$zratio) && !is.null(ch$zratio$gauge),
    logical(1)
  ))
  if(length(keep) == 0) {
    stop(
      "No Z-ratio trust-gauge output found in the chain outputs. It is ",
      "recorded only when the hierarchical prior specification is active and ",
      "the gauge is enabled."
    )
  }
  rows = lapply(keep, function(c_idx) {
    g = chains[[c_idx]]$zratio$gauge
    flip = as.numeric(g$flip_rate)
    floor = as.numeric(g$noise_floor)
    se_mean = as.numeric(g$se_mean)
    se_sd = as.numeric(g$se_sd)
    se_mcse = if(is.null(g$se_mcse)) NA_real_ else as.numeric(g$se_mcse)
    n_ref = as.integer(g$n_ref)
    se_se = if(is.finite(se_mcse) && is.finite(se_sd) && n_ref > 0) {
      sqrt(se_mcse^2 + se_sd^2 / n_ref)
    } else {
      NA_real_
    }

    amplification = NA_real_
    kappa = NA_real_
    harm_pred = NA_real_
    harm_flag = FALSE
    if(!is.null(harm_inputs) && length(harm_inputs$pip) >= c_idx &&
      !is.null(harm_inputs$pip[[c_idx]])) {
      pip = as.numeric(harm_inputs$pip[[c_idx]])
      m_bar = mean(pip * (1 - pip))
      n_edges = length(pip)
      gain = if(!is.null(harm_inputs$a) && !is.null(harm_inputs$b)) {
        a = as.numeric(harm_inputs$a)
        b = as.numeric(harm_inputs$b)
        theta_hat = (a + sum(pip)) / (a + b + n_edges)
        n_edges * m_bar / (theta_hat * (1 - theta_hat) * (a + b + n_edges))
      } else {
        0
      }
      amplification = 1 / (1 - min(gain, 0.98))
      kappa = m_bar * amplification

      # Sensitivity-weighted first-order predictor from the per-pair audit
      # stream (edge e with error s_e contributes m_e s_e); falls back to
      # the unweighted |se_mean| * m_bar form for outputs without the stream.
      pair_i = g$pair_i
      if(!is.null(pair_i) && length(pair_i) > 0 && n_ref > 0) {
        q = round((1 + sqrt(1 + 8 * n_edges)) / 2)
        i0 = pmin(as.integer(pair_i), as.integer(g$pair_j))
        j0 = pmax(as.integer(pair_i), as.integer(g$pair_j))
        # %/% binds tighter than *, so the row-offset product needs the
        # parentheses: idx = (i0 * (2q - i0 - 1)) %/% 2 + (j0 - i0).
        idx = (i0 * (2L * q - i0 - 1L)) %/% 2L + (j0 - i0)
        m_rec = pip[idx] * (1 - pip[idx])
        x = m_rec * as.numeric(g$pair_se)
        num = abs(mean(x))
        # Cluster-robust spread over unique audited edges plus the
        # reference-noise contribution of each record.
        n_rec = length(x)
        cl = split(seq_len(n_rec), idx)
        cr = sum(vapply(cl, function(ii) {
          (sum(x[ii]) - length(ii) * mean(x))^2
        }, numeric(1))) / n_rec^2
        noise2 = sum((m_rec * as.numeric(g$pair_mcse))^2) / n_rec^2
        se_num = sqrt(cr + noise2)
        harm_pred = num * amplification
        resolved = is.finite(se_num) && num > 2 * se_num
      } else {
        harm_pred = abs(se_mean) * m_bar * amplification
        resolved = is.finite(se_se) && abs(se_mean) > 2 * se_se
      }
      harm_flag = isTRUE(resolved && harm_pred > harm_threshold)
    }

    # Mediating-block sizes the audit actually covered; the remediation text
    # quotes them so a flag can be read against the trained hull.
    pair_m = if(is.null(g$pair_m)) integer(0) else as.integer(g$pair_m)
    data.frame(
      chain = c_idx,
      flip_rate = flip,
      flag = isTRUE(flip > threshold + floor),
      se_mean = se_mean,
      se_sd = se_sd,
      se_mcse = se_mcse,
      se_se = se_se,
      noise_floor = floor,
      n_ent = as.integer(g$n_ent),
      n_ref = n_ref,
      n_capped = as.integer(g$n_capped),
      block_lo = if(length(pair_m)) min(pair_m) else NA_integer_,
      block_hi = if(length(pair_m)) max(pair_m) else NA_integer_,
      amplification = amplification,
      kappa = kappa,
      harm_pred = harm_pred,
      harm_flag = harm_flag
    )
  })
  per_chain = do.call(rbind, rows)
  flagged = any(per_chain$flag) || any(per_chain$harm_flag)

  # Print flagged chains under a single header, matching the NUTS-issues block
  # (cat/stdout, one bullet per chain), then one shared remediation ladder. The
  # vignette pointer is emitted once by the output builder as a shared footer,
  # not here. The ladder never proposes an automatic switch: the joint
  # specification is a different model, not a more accurate version of this
  # one, so it is the last rung and is labelled as such.
  if(verbose && flagged && isTRUE(getOption("bgms.verbose", TRUE))) {
    audit = function(pc) {
      blocks = if(is.na(pc$block_lo)) {
        ""
      } else {
        sprintf("; mediating blocks %d-%d variables", pc$block_lo, pc$block_hi)
      }
      sprintf("audited %d of %d non-trivial edge moves%s",
        pc$n_ref, pc$n_ent, blocks
      )
    }
    cat("Graph-prior approximation issues:\n")
    for(i in seq_len(nrow(per_chain))) {
      pc = per_chain[i, ]
      if(pc$flag) {
        cat(sprintf(
          paste0(
            "  - Chain %d: %.1f%% of edge-toggle decisions in the ",
            "approximate chain differ from the exact reference (%s).\n"
          ),
          pc$chain, 100 * pc$flip_rate, audit(pc)
        ))
      }
      if(pc$harm_flag) {
        cat(sprintf(
          paste0(
            "  - Chain %d: the approximation shifts the inclusion ",
            "probabilities by an estimated %.2f (%s).\n"
          ),
          pc$chain, pc$harm_pred, audit(pc)
        ))
      }
    }
    cat(
      "  Raise options(bgms.zratio_gauge_sweeps) and refit to audit more edge\n",
      "  moves and resolve whether the signal is real. If it persists,\n",
      "  precision_graph_prior = \"joint\" avoids the approximation, but it\n",
      "  targets a different model: its graph marginal is the edge prior\n",
      "  reweighted by the per-graph normalizer, not the edge prior itself.\n",
      sep = ""
    )
  }

  invisible(list(
    per_chain = per_chain, threshold = threshold,
    harm_threshold = harm_threshold,
    flagged = flagged
  ))
}

# TRUE if any chain carries recorded trust-gauge output, i.e. the in-chain gauge
# actually ran. Lets callers skip the summary (which errors on empty input) when
# the gauge is disabled (options(bgms.zratio_gauge_sweeps = 0)).
zratio_gauge_present = function(chains) {
  isTRUE(any(vapply(
    chains,
    function(ch) !is.null(ch$zratio) && !is.null(ch$zratio$gauge),
    logical(1)
  )))
}

# One per-chain Z-ratio counter, read defensively: chain output that predates a
# counter reports it as absent rather than as zero, and a notice must not turn
# a missing tally into a claim.
zratio_counter = function(ct, nm) {
  if(is.null(ct) || !(nm %in% names(ct))) {
    return(0)
  }
  v = suppressWarnings(as.numeric(ct[[nm]]))
  if(length(v) != 1L || is.na(v)) 0 else v
}

# The per-chain counter blocks of a fit, dropping chains that carry none.
zratio_counter_blocks = function(chains) {
  Filter(
    function(ct) !is.null(ct),
    lapply(chains, function(ch) ch$zratio$counters)
  )
}

# ------------------------------------------------------------------------------
# zratio_isolated_route_notice
# ------------------------------------------------------------------------------
# One graceful, per-fit notice when the hierarchical prior served the
# isolated-edge ratio because the Gamma diagonal shape sits past the surface's
# validated range (zratio_mediation_off). The engine tallies `n_isolated` per
# chain, so this fires on what the fit actually did rather than on what the spec
# intended -- the spec-build message states the policy before the chains launch,
# and this states the outcome after they finish.
#
# It carries the bound, because the bound is the whole reason the route is
# acceptable: the value served is exact for an edge with no mediating structure,
# so the entire error is the mediation dropped, measured at no more than
# 2.8e-04 nats at diagonal rates up to .zratio_mediation_off_eta_hi. Past that
# rate the same route deploys and the notice says plainly that the measurement
# does not cover the cell, which is the only honest thing to say about it.
#
# @param chains  Raw per-chain sampler output.
# @param eta     The fit's standardized diagonal rate.
#
# Returns invisible(TRUE) when a notice was emitted.
# ------------------------------------------------------------------------------
zratio_isolated_route_notice = function(chains, eta) {
  counters = zratio_counter_blocks(chains)
  if(length(counters) == 0) {
    return(invisible(FALSE))
  }
  n_iso = sum(vapply(counters, zratio_counter, numeric(1), "n_isolated"))
  if(n_iso <= 0) {
    return(invisible(FALSE))
  }
  covered = is.finite(eta) && eta <= .zratio_mediation_off_eta_hi
  message(sprintf(
    paste0(
      "Note: the precision diagonal's Gamma shape is past the range the ",
      "hierarchical prior's edge correction is scored on, so all %s edge ",
      "evaluations in this fit used the isolated-edge normalizer ratio with ",
      "the mediating correction switched off. That value is exact for an edge ",
      "with no mediating structure, so the whole error is the mediating ",
      "correction it leaves out. %s"
    ),
    format(n_iso, big.mark = ",", scientific = FALSE),
    if(covered) {
      paste0(
        "At this diagonal rate that correction was measured against a ",
        "block-Gibbs reference at no more than 0.00028 nats, two orders below ",
        "the 0.003 nats the correction is held to inside its range."
      )
    } else {
      sprintf(
        paste0(
          "That correction grows with the diagonal rate and was measured only ",
          "at rates up to %s; this fit runs at %s, so the measured bound does ",
          "not cover it."
        ),
        format(.zratio_mediation_off_eta_hi), format(eta)
      )
    }
  ))
  invisible(TRUE)
}

# One graceful, per-fit notice when the hierarchical prior's fast edge correction
# was extrapolated beyond its validated block-size range. Mediating blocks larger
# than the trained surface hull are clamped at deploy (dense regions of large
# graphs); the C++ engine tallies how often per chain. This sums the tally and,
# if any block exceeded the hull, emits a single summary message. Independent of
# the trust gauge, so the signal reaches the user even with the gauge off.
zratio_extrapolation_notice = function(chains) {
  counters = zratio_counter_blocks(chains)
  if(length(counters) == 0) {
    return(invisible(FALSE))
  }
  n_extrap = sum(vapply(counters, zratio_counter, numeric(1), "n_extrap"))
  if(n_extrap <= 0) {
    return(invisible(FALSE))
  }
  n_pred = sum(vapply(counters, zratio_counter, numeric(1), "n_pred"))
  # The retained share is the one that describes the posterior: the sampler
  # initializes from a complete graph, so warmup alone puts every mediating
  # block past the hull for the first sweeps. Warmup is reported in brackets so
  # a warmup-only transient reads as what it is.
  n_extrap_ret = sum(vapply(counters, zratio_counter, numeric(1), "n_extrap_retained"))
  n_pred_ret = sum(vapply(counters, zratio_counter, numeric(1), "n_pred_retained"))
  max_ret = max(vapply(counters, zratio_counter, numeric(1), "max_extrap_size_retained"))
  max_all = max(vapply(counters, zratio_counter, numeric(1), "max_extrap_size"))
  pct_ret = if(n_pred_ret > 0) 100 * n_extrap_ret / n_pred_ret else 0
  pct_warm = if(n_pred - n_pred_ret > 0) {
    100 * (n_extrap - n_extrap_ret) / (n_pred - n_pred_ret)
  } else {
    NA_real_
  }
  # Chain output from before the phase split carries no retained counters, so
  # no phase can be claimed for it: report the whole-run share instead of
  # reassuring the reader that the stored draws are clean.
  has_split = any(vapply(
    counters,
    function(ct) all(c("n_pred_retained", "n_extrap_retained") %in% names(ct)),
    logical(1)
  ))
  if(!has_split) {
    message(sprintf(
      paste0(
        "Note: %.1f%% of the hierarchical prior's edge-correction evaluations ",
        "used a mediating block beyond its anchored size range (largest %d ",
        "variables). The correction is extended along the surface's own ",
        "boundary slope there. The trust gauge reports the measured ",
        "sensitivity per chain in fit$zratio_diag."
      ),
      if(n_pred > 0) 100 * n_extrap / n_pred else NA_real_, as.integer(max_all)
    ))
    return(invisible(TRUE))
  }
  if(n_extrap_ret <= 0) {
    message(sprintf(
      paste0(
        "Note: the hierarchical prior's edge correction was extended beyond ",
        "its anchored size range during warmup only (%.1f%% of warmup ",
        "evaluations, largest block %d variables); no retained sweep did. The ",
        "sampler starts from a complete graph, so this is the initial ",
        "transient and the stored draws are unaffected."
      ),
      pct_warm, as.integer(max_all)
    ))
    return(invisible(TRUE))
  }
  message(sprintf(
    paste0(
      "Note: %.1f%% of the hierarchical prior's edge-correction evaluations in ",
      "the retained sweeps used a mediating block beyond its anchored size ",
      "range (largest %d variables; warmup %.1f%%). The correction is extended ",
      "along the surface's own boundary slope there, which was measured against ",
      "a block-Gibbs reference at blocks of 90 to 150 variables at a median of ",
      "0.0006 nats and at most 0.0011 for common-neighbour blocks, and a median ",
      "of 0.0043 and at most 0.0060 for bipartite ones, against about 0.003 ",
      "inside the anchored range. Sparse graphs never ",
      "reach this. The trust gauge reports the measured sensitivity per chain ",
      "in fit$zratio_diag."
    ),
    pct_ret, as.integer(max_ret), pct_warm
  ))
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# zratio_harm_inputs
# ------------------------------------------------------------------------------
# Assembles the harm_inputs argument of summarize_zratio_gauge from per-chain
# edge-inclusion probabilities and the edge-prior identity. The harm channel is
# defined for the Bernoulli (fixed inclusion probability, no feedback) and
# Beta-Bernoulli (sampled inclusion probability) edge priors; for any other
# prior the channel is disabled (NULL return).
#
# @param pip         List with one numeric vector of posterior edge-inclusion
#                    probabilities per chain.
# @param edge_prior  Character edge-prior name (e.g. "Bernoulli",
#                    "Beta-Bernoulli").
# @param a,b         Numeric Beta-Bernoulli shape parameters; ignored for other
#                    priors.
#
# Returns: A list for harm_inputs, or NULL when the edge prior is not covered.
# ------------------------------------------------------------------------------
zratio_harm_inputs = function(pip, edge_prior, a = NULL, b = NULL) {
  if(identical(edge_prior, "Bernoulli")) {
    list(pip = pip, a = NULL, b = NULL)
  } else if(identical(edge_prior, "Beta-Bernoulli")) {
    list(pip = pip, a = a, b = b)
  } else {
    NULL
  }
}
