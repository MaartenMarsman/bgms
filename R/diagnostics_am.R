# ==============================================================================
# Adaptive-Metropolis diagnostics
# ==============================================================================
#
# Post-sampling diagnostics for the componentwise adaptive-Metropolis
# sampler. Mirrors the role of diagnostics_nuts.R for NUTS, but
# necessarily lives in its own file because the underlying sampler
# organisation is different: NUTS does one accept/reject per iteration
# over the whole parameter vector, whereas adaptive-Metropolis does
# one accept/reject per parameter per iteration. The per-parameter move
# summary therefore has a shape (parameter x chain) with no NUTS analogue,
# alongside the per-sweep acceptance trace (chain x iteration) that does
# mirror the NUTS one.
# ==============================================================================


# ------------------------------------------------------------------------------
# summarize_am_diagnostics
# ------------------------------------------------------------------------------
# Combine and summarize adaptive-Metropolis diagnostics across chains.
#
# Two distinct quantities are reported, because the obvious one is not the
# one the tuner targets:
#
#   accept_prob -- the Metropolis acceptance probability the sampler itself
#     records. The C++ sampler averages the acceptance probability over all
#     components updated in a sweep and stores that per iteration as the
#     `am_accept_prob__` trace (metropolis_sampler.h -> ChainResult ->
#     build_output.R). This is the statistic Robbins-Monro adaptation drives
#     towards `target_accept`, so it is the one to compare against it.
#
#   move_rate -- the empirical *move rate*: the fraction of post-warmup
#     iterations on which a parameter's stored value differed from the
#     previous iteration. It is not an acceptance rate whenever a parameter
#     can go un-proposed: under edge selection an excluded pairwise effect is
#     held at exactly 0 and never moves, so its move rate is approximately
#     P(included) * P(accept) and sits well below `target_accept` even when
#     the sampler is tuned correctly. Reported per parameter, which the
#     acceptance trace is not, so it still shows which parameters are moving.
#
# @param out             List of chain outputs. Each element is a named list
#                        that must contain "main_samples" and
#                        "pairwise_samples" matrices (iterations x params),
#                        and optionally the "am_accept_prob__" trace.
# @param names_main      Character vector of main-effect parameter names.
# @param names_pairwise  Character vector of pairwise interaction names.
# @param target_accept   Target acceptance rate the sampler was tuned to
#                        (default: 0.44 for componentwise random-walk MH).
#
# Returns: An invisible named list with:
#   - accept_prob:   Numeric matrix (chains x iterations) of the sampler's
#       own per-sweep mean Metropolis acceptance probability, or NULL when
#       the fit carries no `am_accept_prob__` trace (pre-0.2.0.0 fits, and
#       models that take no Metropolis step). Rows are labelled
#       "chain 1", "chain 2", ...
#   - move_rate:     Numeric matrix (parameter x chain) of per-parameter
#       per-chain empirical move rates. Rows are labelled by parameter
#       (main effects first, then pairwise); columns are labelled
#       "chain 1", "chain 2", ...
#   - target_accept: Numeric scalar; the target acceptance rate. Comparable
#       with accept_prob, not with move_rate.
#   - summary:       List with mean_accept_prob (mean of accept_prob over
#       all iterations and chains, NA when the trace is absent) and
#       mean_move_rate (mean of move_rate over all parameters and chains).
# ------------------------------------------------------------------------------
summarize_am_diagnostics = function(out, names_main, names_pairwise,
                                    target_accept = 0.44) {
  am_chains = Filter(function(chain) {
    all(c("main_samples", "pairwise_samples") %in% names(chain))
  }, out)

  if(length(am_chains) == 0) {
    return(NULL)
  }

  per_chain_move_rate = function(chain) {
    main_acc = apply(
      chain$main_samples, 2,
      function(col) mean(diff(col) != 0)
    )
    pair_acc = apply(
      chain$pairwise_samples, 2,
      function(col) mean(diff(col) != 0)
    )
    c(main_acc, pair_acc)
  }

  move_rate_mat = vapply(
    am_chains, per_chain_move_rate,
    numeric(length(names_main) + length(names_pairwise))
  )

  rownames(move_rate_mat) = c(names_main, names_pairwise)
  colnames(move_rate_mat) = paste0("chain ", seq_along(am_chains))

  # The acceptance trace is required on every chain before it is combined: a
  # mix of chains with and without it has no rectangular form, and a partial
  # summary would be reported under the same name as a complete one.
  # A trace that is present but never written stays at the NaN fill value
  # (ChainResult), which is the same as having no trace at all.
  has_accept_trace = all(vapply(
    am_chains,
    function(chain) {
      trace = chain[["am_accept_prob__"]]
      !is.null(trace) && any(is.finite(trace))
    },
    logical(1)
  ))

  if(has_accept_trace) {
    accept_prob_mat = do.call(rbind, lapply(
      am_chains, function(chain) as.numeric(chain[["am_accept_prob__"]])
    ))
    rownames(accept_prob_mat) = paste0("chain ", seq_along(am_chains))
    mean_accept_prob = mean(accept_prob_mat, na.rm = TRUE)
  } else {
    accept_prob_mat = NULL
    mean_accept_prob = NA_real_
  }

  invisible(list(
    accept_prob = accept_prob_mat,
    move_rate = move_rate_mat,
    target_accept = target_accept,
    summary = list(
      mean_accept_prob = mean_accept_prob,
      mean_move_rate = mean(move_rate_mat, na.rm = TRUE)
    )
  ))
}
