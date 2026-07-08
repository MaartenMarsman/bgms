#' @title Summarize the In-Chain Z-Ratio Trust Gauge
#'
#' @description Reports the per-chain trust gauge for the hierarchical prior
#' specification's per-edge Z-ratio kernel. During sampling the gauge compares,
#' on each chain's own edge moves, the deployed approximate normalizing-constant
#' ratio against a block-local exact Monte-Carlo reference and records
#' \code{D} = the fraction of accept/reject decisions that would flip under the
#' exact ratio. A chain is flagged when \code{D} exceeds the tolerance by more
#' than the reference's own Monte-Carlo noise floor. The gauge flags; it does
#' not fix. \code{D} is in direct impact units, so the same tolerance applies at
#' any number of variables and any graph prior.
#'
#' The gauge reads the in-chain \code{zratio$gauge} block that the sampler
#' attaches under the hierarchical specification; it does not re-scan the stored
#' draws.
#'
#' @param chains List of per-chain sampler outputs, each carrying a
#'   \code{zratio$gauge} block (\code{D}, \code{noise_floor}, \code{se_mean},
#'   \code{se_sd}, \code{n_ent}, \code{n_ref}, \code{n_capped}).
#' @param threshold Numeric flag threshold on \code{D} (default \code{0.01}).
#' @param verbose Logical: message flagged chains (default \code{TRUE}).
#'
#' @return An invisible named list:
#'   \describe{
#'     \item{\code{per_chain}}{Data frame, one row per chain: \code{D}, the
#'       \code{flag}, the signed mean and spread of the log-ratio error
#'       (\code{se_mean}, \code{se_sd}), the reference \code{noise_floor}, and
#'       the pair counts (\code{n_ent} non-trivial seen, \code{n_ref}
#'       referenced, \code{n_capped} cap hits).}
#'     \item{\code{threshold}}{The flag tolerance on \code{D}.}
#'     \item{\code{flagged}}{Logical: any chain flagged.}
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
summarize_zratio_gauge = function(chains, threshold = 0.01, verbose = TRUE) {
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
    d_val = as.numeric(g$D)
    floor = as.numeric(g$noise_floor)
    data.frame(
      chain = c_idx,
      D = d_val,
      flag = isTRUE(d_val > threshold + floor),
      se_mean = as.numeric(g$se_mean),
      se_sd = as.numeric(g$se_sd),
      noise_floor = floor,
      n_ent = as.integer(g$n_ent),
      n_ref = as.integer(g$n_ref),
      n_capped = as.integer(g$n_capped)
    )
  })
  per_chain = do.call(rbind, rows)
  flagged = any(per_chain$flag)

  if(verbose) {
    for(i in seq_len(nrow(per_chain))) {
      pc = per_chain[i, ]
      if(pc$flag) {
        # Placeholder wording; the user-facing message is designed separately.
        message(sprintf(
          paste0("Chain %d: AMBER -- %.1f%% of edge moves would flip under ",
                 "the exact ratio (tolerance %.1f%%)."),
          pc$chain, 100 * pc$D, 100 * threshold
        ))
      }
    }
  }

  invisible(list(per_chain = per_chain, threshold = threshold,
                 flagged = flagged))
}
