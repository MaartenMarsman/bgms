#' @title Summarize the Hierarchical Prior Trust Gauge
#'
#' @description Reports the per-chain trust gauge for the hierarchical graph
#' prior. Under that prior the sampler decides each edge with a fast
#' approximation; during sampling the gauge redoes a subset of each chain's own
#' edge decisions with the exact calculation and records two statistics per
#' chain:
#' \describe{
#'   \item{\code{flip_rate}}{The fraction of add/remove decisions that would
#'     come out differently under the exact calculation. A chain is flagged on
#'     this channel when \code{flip_rate} exceeds the tolerance by more than
#'     the exact reference's own Monte-Carlo noise. This channel detects error
#'     that changed decisions the chain actually made; it is insensitive to a
#'     small coherent error at chains whose decisions are far from their
#'     accept/reject boundaries.}
#'   \item{\code{harm_pred}}{The projected inclusion-probability distortion
#'     from the approximation error: \code{|se_mean| * m * A}, where
#'     \code{se_mean} is the signed mean log-ratio error against the exact
#'     reference, \code{m} is the chain's mean per-edge sensitivity
#'     \code{mean(p_e (1 - p_e))}, and \code{A = 1 / (1 - g)} is the
#'     inclusion-probability feedback amplification with linearized gain
#'     \code{g = E m / (theta (1 - theta) (a + b + E))} under a Beta-Bernoulli
#'     edge prior (for a fixed inclusion probability \code{g = 0}, so
#'     \code{A = 1}). A chain is flagged on this channel when
#'     \code{harm_pred} exceeds the tolerance and \code{se_mean} is resolved
#'     above its own standard error, \code{|se_mean| > 2 se_se} with
#'     \code{se_se = sqrt(se_mcse^2 + se_sd^2 / n_ref)}. This channel detects
#'     coherent error whose equilibrium effect exceeds the tolerance even when
#'     no individual decision visibly flips. It is computed for Bernoulli and
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
#'       \code{harm_pred}, \code{harm_flag}).}
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
      harm_pred = abs(se_mean) * m_bar * amplification
      resolved = is.finite(se_se) && abs(se_mean) > 2 * se_se
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
            "approximate chain differ from the exact reference - ",
            "increase calibration_window\n"
          ),
          pc$chain, 100 * pc$flip_rate
        ))
      }
      if(pc$harm_flag) {
        cat(sprintf(
          paste0(
            "  - Chain %d: the approximation error projects to a %.3f ",
            "distortion of the inclusion probabilities (tolerance %.2f) - ",
            "increase calibration_window; for evidence-free (prior-only) ",
            "use, prefer the joint specification\n"
          ),
          pc$chain, pc$harm_pred, harm_threshold
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
