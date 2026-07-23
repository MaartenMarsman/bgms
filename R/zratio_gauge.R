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
#'       \code{n_capped} cap hits), and the harm channel (\code{amplification},
#'       \code{kappa} = the predicted mean-inclusion shift per nat of coherent
#'       error, \code{harm_pred}, \code{harm_flag}).}
#'     \item{\code{threshold}}{The flag tolerance on \code{flip_rate}.}
#'     \item{\code{harm_threshold}}{The flag tolerance on \code{harm_pred}.}
#'     \item{\code{flagged}}{Logical: any chain flagged on either channel.}
#'   }
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
  chains = Filter(
    function(ch) !is.null(ch$zratio) && !is.null(ch$zratio$gauge), chains
  )
  if(length(chains) == 0) {
    stop(
      "No Z-ratio trust-gauge output found in the chain outputs. It is ",
      "recorded only when the hierarchical prior specification is active and ",
      "the gauge is enabled."
    )
  }
  rows = lapply(seq_along(chains), function(c_idx) {
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
        idx = i0 * (2L * q - i0 - 1L) %/% 2L + (j0 - i0)
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
      amplification = amplification,
      kappa = kappa,
      harm_pred = harm_pred,
      harm_flag = harm_flag
    )
  })
  per_chain = do.call(rbind, rows)
  flagged = any(per_chain$flag) || any(per_chain$harm_flag)

  # Print flagged chains under a single header, matching the NUTS-issues block
  # (cat/stdout, one bullet per chain). The vignette pointer is emitted once by
  # the output builder as a shared footer, not here.
  if(verbose && flagged && isTRUE(getOption("bgms.verbose", TRUE))) {
    cat("Graph-prior approximation issues:\n")
    for(i in seq_len(nrow(per_chain))) {
      pc = per_chain[i, ]
      if(pc$flag) {
        cat(sprintf(
          paste0(
            "  - Chain %d: %.1f%% of edge-toggle decisions in the ",
            "approximate chain differ from the exact reference - the ",
            "z-ratio surface may be inaccurate here; consider ",
            "precision_graph_prior = \"joint\"\n"
          ),
          pc$chain, 100 * pc$flip_rate
        ))
      }
      if(pc$harm_flag) {
        cat(sprintf(
          paste0(
            "  - Chain %d: the approximation biases the inclusion ",
            "probabilities by an estimated %.2f - the z-ratio surface may ",
            "be inaccurate here; consider precision_graph_prior = \"joint\"\n"
          ),
          pc$chain, pc$harm_pred
        ))
      }
    }
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

# One graceful, per-fit notice when the hierarchical prior's fast edge correction
# was extrapolated beyond its validated block-size range. Mediating blocks larger
# than the trained surface hull are clamped at deploy (dense regions of large
# graphs); the C++ engine tallies how often per chain. This sums the tally and,
# if any block exceeded the hull, emits a single summary message. Independent of
# the trust gauge, so the signal reaches the user even with the gauge off.
zratio_extrapolation_notice = function(chains) {
  get_counter = function(ct, nm) {
    if(is.null(ct) || !(nm %in% names(ct))) {
      return(0)
    }
    v = suppressWarnings(as.numeric(ct[[nm]]))
    if(length(v) != 1L || is.na(v)) 0 else v
  }
  counters = Filter(
    function(ct) !is.null(ct),
    lapply(chains, function(ch) ch$zratio$counters)
  )
  if(length(counters) == 0) {
    return(invisible(FALSE))
  }
  n_extrap = sum(vapply(counters, get_counter, numeric(1), "n_extrap"))
  if(n_extrap <= 0) {
    return(invisible(FALSE))
  }
  n_pred = sum(vapply(counters, get_counter, numeric(1), "n_pred"))
  max_size = max(vapply(counters, get_counter, numeric(1), "max_extrap_size"))
  pct = if(n_pred > 0) 100 * n_extrap / n_pred else NA_real_
  message(sprintf(
    paste0(
      "Note: %.1f%% of the hierarchical prior's edge-correction evaluations ",
      "used a mediating block beyond its validated size range (largest %d ",
      "variables). The correction was extrapolated there, which can slightly ",
      "reduce edge-selection accuracy in dense regions of large graphs; sparse ",
      "graphs are unaffected. To assess the sensitivity, enable the trust gauge ",
      "with options(bgms.zratio_gauge_sweeps = 2L)."
    ),
    pct, as.integer(max_size)
  ))
  invisible(TRUE)
}

#' @title Harm-Channel Inputs for the Trust Gauge
#'
#' @description Assembles the \code{harm_inputs} argument of
#' \code{\link{summarize_zratio_gauge}} from per-chain edge-inclusion
#' probabilities and the edge-prior identity. The harm channel is defined for
#' the Bernoulli (fixed inclusion probability, no feedback) and Beta-Bernoulli
#' (sampled inclusion probability) edge priors; for any other prior the
#' channel is disabled (\code{NULL} return).
#'
#' @param pip List with one numeric vector of posterior edge-inclusion
#'   probabilities per chain.
#' @param edge_prior Character edge-prior name (e.g. \code{"Bernoulli"},
#'   \code{"Beta-Bernoulli"}).
#' @param a,b Numeric Beta-Bernoulli shape parameters; ignored for other
#'   priors.
#'
#' @return A list for \code{harm_inputs}, or \code{NULL} when the edge prior
#'   is not covered.
#'
#' @keywords internal
zratio_harm_inputs = function(pip, edge_prior, a = NULL, b = NULL) {
  if(identical(edge_prior, "Bernoulli")) {
    list(pip = pip, a = NULL, b = NULL)
  } else if(identical(edge_prior, "Beta-Bernoulli")) {
    list(pip = pip, a = a, b = b)
  } else {
    NULL
  }
}
