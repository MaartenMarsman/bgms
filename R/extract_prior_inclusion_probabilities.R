# ==============================================================================
# extract_prior_inclusion_probabilities() — parallel to the posterior extractor
# ==============================================================================
#
# Returns a symmetric matrix of prior edge-inclusion probabilities for a bgms
# fit, in the same shape and orientation as
# extract_posterior_inclusion_probabilities(), for use in prior/posterior
# inclusion-odds computations.
#
# Under the joint spike-and-slab prior on a continuous block the graph
# marginal is reweighted by the per-graph normalizer Z(Gamma) — the
# positive-definite-cone mass of the slab under the pattern, further shaped
# by the determinant tilt — so the prior inclusion probability of a
# continuous-continuous edge is NOT the edge-prior marginal at ANY delta,
# including delta = 0: it is the prior edge density of the joint block.
# Discrete-discrete and cross interactions are unconstrained with normalized
# slabs, so their edges keep the edge-prior marginal given the
# hyperparameters. The nodes are exchangeable a priori (within variable
# type), so the marginal is constant within an edge class and the whole
# matrix reduces to per-class values:
#
#   - GGM: one class (all pairs in the joint block).
#   - Mixed MRF: three classes — discrete-discrete and cross pairs carry the
#     edge-prior marginal, continuous-continuous pairs the joint-block
#     density (when at least two variables are continuous).
#   - OMRF: one class, fully factorized.
#
# Per-class values by edge prior:
#   - Bernoulli(theta0): factorized classes theta0; joint-block classes the
#     prior edge density e(theta0) read from the correction table's edens
#     curve.
#   - Beta-Bernoulli(a, b): factorized classes a/(a+b) (the theta marginal
#     is Beta(a, b) under the corrected chain); joint-block classes the
#     Beta(a, b) quadrature of the edens curve.
#   - Stochastic-Block: fully factorized fits (no joint block) have the
#     analytic exchangeable marginal s * mu_within + (1 - s) * mu_between
#     with s = P(same block) =
#     sum_k dpois(k - 1, lambda) (alpha + 1)/(k alpha + 1); fits with a
#     joint block need a prior-only chain (the block reweights the
#     partition), run with the same correction resolution as the fit and
#     pooled within edge class. Chain estimates are cached on fit$cache.
#
# A prior-only chain is also the fallback for Bernoulli/Beta-Bernoulli fits
# whose slab family the table builder does not support (beta-prime).
# ==============================================================================


# ------------------------------------------------------------------
# prior_pip_table (internal)
# ------------------------------------------------------------------
# Get-or-build the correction table for the fit's continuous-block
# cell, whose edens curve is the tilted prior edge density. NULL when
# the table builder does not support the slab family.
# ------------------------------------------------------------------
prior_pip_table = function(prior, num_continuous) {
  interaction_prior = correction_interaction_prior(prior)
  if(is.null(interaction_prior)) {
    return(NULL)
  }
  ggm_correction_table(
    p = num_continuous, delta = prior$delta,
    interaction_prior = interaction_prior,
    precision_scale_prior = correction_scale_prior(prior),
    update_method = "gibbs"
  )
}


# ------------------------------------------------------------------
# tilted_edge_density (internal)
# ------------------------------------------------------------------
# Tilted prior edge density at a fixed inclusion probability theta0
# (linear interpolation on the table grid, constant extension).
# ------------------------------------------------------------------
tilted_edge_density = function(table, theta0) {
  stats::approx(table$theta, table$edens, theta0, rule = 2)$y
}


# ------------------------------------------------------------------
# tilted_bb_edge_marginal (internal)
# ------------------------------------------------------------------
# Beta(a, b) mixture of the tilted edge density: the corrected chain's
# theta marginal is the Beta(a, b) hyperprior, and given theta the
# tilted edge density is the table's edens curve. Trapezoid quadrature
# on the table grid, normalized by the same-grid Beta mass so the
# truncated tails cancel.
# ------------------------------------------------------------------
tilted_bb_edge_marginal = function(table, alpha, beta) {
  theta = table$theta
  w = stats::dbeta(theta, alpha, beta)
  trap = function(y) {
    sum(diff(theta) * (y[-length(y)] + y[-1]) / 2)
  }
  trap(w * table$edens) / trap(w)
}


# ------------------------------------------------------------------
# sbm_factorized_edge_marginal (internal)
# ------------------------------------------------------------------
# Exchangeable pair-inclusion marginal of the MFM-SBM prior when no
# joint continuous block reweights the graph: allocations are iid
# given symmetric Dirichlet(alpha) weights over K = 1 + Poisson(lambda)
# components, so P(same block | K) = (alpha + 1) / (K alpha + 1), and
# the edge probability mixes the within- and between-block Beta means.
# ------------------------------------------------------------------
sbm_factorized_edge_marginal = function(prior) {
  alpha = prior$dirichlet_alpha
  lambda = prior$lambda
  k = seq_len(200L)
  s = sum(stats::dpois(k - 1, lambda) * (alpha + 1) / (k * alpha + 1))
  mu_within = prior$beta_bernoulli_alpha /
    (prior$beta_bernoulli_alpha + prior$beta_bernoulli_beta)
  mu_between = if(is.null(prior$beta_bernoulli_alpha_between)) {
    mu_within
  } else {
    prior$beta_bernoulli_alpha_between /
      (prior$beta_bernoulli_alpha_between + prior$beta_bernoulli_beta_between)
  }
  s * mu_within + (1 - s) * mu_between
}


# ------------------------------------------------------------------
# prior_only_chain_pips (internal)
# ------------------------------------------------------------------
# Run a prior-only chain from the fit's spec (likelihood muted with
# n = 0 data) and return the pooled per-class inclusion proportions:
# a single value for GGM, and c(dd, cc, cross) for the mixed MRF (NA
# for empty classes). The edge-prior correction is resolved exactly as
# at fit time, so the chain targets the same prior the fit used.
# ------------------------------------------------------------------
prior_only_chain_pips = function(spec, iter, warmup) {
  d = spec$data
  p = spec$prior
  s = spec$sampler
  chain_sampler = list(verbose = FALSE, cores = 1L)

  bb_alpha_between = bb_between_or_sentinel(p$beta_bernoulli_alpha_between)
  bb_beta_between = bb_between_or_sentinel(p$beta_bernoulli_beta_between)

  if(identical(spec$model_type, "ggm")) {
    num_vars = d$num_variables
    correction = ggm_edge_prior_correction(
      p, chain_sampler, num_vars
    )
    conjugate_slab = p$interaction_prior_type %in% c("cauchy", "normal") &&
      p$scale_prior_type %in% c("gamma", "exponential")
    results = sample_ggm(
      inputFromR = list(
        n = 0L,
        suf_stat = matrix(0, num_vars, num_vars),
        pairwise_scale = p$pairwise_scale,
        interaction_prior_type = p$interaction_prior_type,
        interaction_alpha = p$interaction_alpha,
        interaction_beta = p$interaction_beta,
        scale_prior_type = p$scale_prior_type,
        scale_shape = p$scale_shape,
        scale_rate = p$scale_rate
      ),
      prior_inclusion_prob = p$inclusion_probability,
      initial_edge_indicators = matrix(1L, num_vars, num_vars),
      no_iter = as.integer(iter),
      no_warmup = as.integer(warmup),
      no_chains = 1L,
      edge_selection = TRUE,
      sampler_type = if(conjugate_slab) "gibbs" else "adaptive-metropolis",
      seed = as.integer(s$seed + 1L),
      no_threads = 1L,
      progress_type = 0L,
      edge_prior = p$edge_prior,
      beta_bernoulli_alpha = p$beta_bernoulli_alpha,
      beta_bernoulli_beta = p$beta_bernoulli_beta,
      beta_bernoulli_alpha_between = bb_alpha_between,
      beta_bernoulli_beta_between = bb_beta_between,
      dirichlet_alpha = p$dirichlet_alpha,
      lambda = p$lambda,
      delta = p$delta,
      edge_prior_correction = correction
    )
    if(length(results) == 0L || isTRUE(results[[1L]]$error)) {
      stop(
        "Prior-only chain failed while estimating prior inclusion ",
        "probabilities."
      )
    }
    # Indicator rows are the full upper triangle including the diagonal in
    # (i <= j) order: pool the off-diagonal rows.
    ind = results[[1L]]$indicator_samples
    is_diag = logical(nrow(ind))
    pos = 0L
    for(i in seq_len(num_vars)) {
      for(j in i:num_vars) {
        pos = pos + 1L
        is_diag[pos] = (i == j)
      }
    }
    return(mean(ind[!is_diag, , drop = FALSE]))
  }

  # Mixed MRF: mute the likelihood with zero-row data of the fit's shape.
  num_disc = d$num_discrete
  num_cont = d$num_continuous
  q_tot = d$num_variables
  correction = ggm_edge_prior_correction(
    p, chain_sampler, q_tot, num_cont
  )
  results = sample_mixed_mrf(
    inputFromR = list(
      discrete_observations = matrix(0L, 0, num_disc),
      continuous_observations = matrix(0, 0, num_cont),
      num_categories = d$num_categories,
      is_ordinal_variable = as.integer(spec$variables$is_ordinal),
      baseline_category = spec$variables$baseline_category,
      interaction_prior_type = p$interaction_prior_type,
      pairwise_scale = p$pairwise_scale,
      interaction_alpha = p$interaction_alpha,
      interaction_beta = p$interaction_beta,
      threshold_prior_type = p$threshold_prior_type,
      main_alpha = p$main_alpha,
      main_beta = p$main_beta,
      threshold_scale = p$threshold_scale,
      means_prior_type = p$means_prior_type,
      means_scale = p$means_scale,
      means_alpha = p$means_alpha,
      means_beta = p$means_beta,
      scale_prior_type = p$scale_prior_type,
      scale_shape = p$scale_shape,
      scale_rate = p$scale_rate
    ),
    prior_inclusion_prob = p$inclusion_probability,
    initial_edge_indicators = matrix(1L, q_tot, q_tot),
    no_iter = as.integer(iter),
    no_warmup = as.integer(warmup),
    no_chains = 1L,
    edge_selection = TRUE,
    seed = as.integer(s$seed + 1L),
    no_threads = 1L,
    progress_type = 0L,
    edge_prior = p$edge_prior,
    beta_bernoulli_alpha = p$beta_bernoulli_alpha,
    beta_bernoulli_beta = p$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = bb_alpha_between,
    beta_bernoulli_beta_between = bb_beta_between,
    dirichlet_alpha = p$dirichlet_alpha,
    lambda = p$lambda,
    sampler_type = "adaptive-metropolis",
    delta = p$delta,
    edge_prior_correction = correction
  )
  if(length(results) == 0L || isTRUE(results[[1L]]$error)) {
    stop(
      "Prior-only chain failed while estimating prior inclusion ",
      "probabilities."
    )
  }
  # Indicator rows are [dd upper-tri | cc upper-tri | cross], all pairwise.
  ind = results[[1L]]$indicator_samples
  n_dd = num_disc * (num_disc - 1L) / 2L
  n_cc = num_cont * (num_cont - 1L) / 2L
  n_cross = num_disc * num_cont
  stopifnot(nrow(ind) == n_dd + n_cc + n_cross)
  pool = function(rows) {
    if(length(rows) == 0L) NA_real_ else mean(ind[rows, , drop = FALSE])
  }
  c(
    dd = pool(seq_len(n_dd)),
    cc = pool(n_dd + seq_len(n_cc)),
    cross = pool(n_dd + n_cc + seq_len(n_cross))
  )
}


# ------------------------------------------------------------------
# prior_pip_class_values (internal)
# ------------------------------------------------------------------
# Per-class prior inclusion probabilities c(dd, cc, cross) for a fit
# spec (GGM and OMRF use the cc slot only). Runs the prior-only chain
# when no closed or table form applies.
# ------------------------------------------------------------------
prior_pip_class_values = function(spec, iter, warmup) {
  p = spec$prior
  model_type = spec$model_type
  num_cont = if(identical(model_type, "mixed_mrf")) {
    spec$data$num_continuous
  } else if(identical(model_type, "ggm")) {
    spec$data$num_variables
  } else {
    0L
  }
  has_joint_block = num_cont >= 2L

  if(identical(p$edge_prior, "Stochastic-Block")) {
    if(!has_joint_block) {
      value = sbm_factorized_edge_marginal(p)
      return(c(dd = value, cc = value, cross = value))
    }
    pips = prior_only_chain_pips(spec, iter, warmup)
    if(identical(model_type, "ggm")) {
      return(c(dd = NA_real_, cc = pips, cross = NA_real_))
    }
    return(pips)
  }

  # Bernoulli / Beta-Bernoulli: classes outside the joint block keep the
  # edge-prior marginal; the joint-block class reads the table's edens
  # curve, or falls back to a prior-only chain when the slab family has
  # no table.
  offdiag = p$inclusion_probability[
    upper.tri(p$inclusion_probability)
  ]
  factorized_value = switch(p$edge_prior,
    Bernoulli = offdiag[1L],
    `Beta-Bernoulli` = p$beta_bernoulli_alpha /
      (p$beta_bernoulli_alpha + p$beta_bernoulli_beta)
  )
  if(identical(p$edge_prior, "Bernoulli") &&
    length(unique(offdiag)) > 1L && has_joint_block) {
    stop(
      "extract_prior_inclusion_probabilities(): per-edge Bernoulli ",
      "probabilities are not supported for models with a continuous ",
      "block; the joint-block marginal is only tabulated for a shared ",
      "probability."
    )
  }
  if(!has_joint_block) {
    return(c(
      dd = factorized_value, cc = factorized_value, cross = factorized_value
    ))
  }

  table = prior_pip_table(p, num_cont)
  if(is.null(table)) {
    pips = prior_only_chain_pips(spec, iter, warmup)
    if(identical(model_type, "ggm")) {
      return(c(dd = NA_real_, cc = pips, cross = NA_real_))
    }
    return(pips)
  }
  cc_value = switch(p$edge_prior,
    Bernoulli = tilted_edge_density(table, factorized_value),
    `Beta-Bernoulli` = tilted_bb_edge_marginal(
      table, p$beta_bernoulli_alpha, p$beta_bernoulli_beta
    )
  )
  c(dd = factorized_value, cc = cc_value, cross = factorized_value)
}


# ------------------------------------------------------------------
# prior_pip_matrix_from_classes (internal)
# ------------------------------------------------------------------
# Assemble the symmetric matrix from per-class values, in the same
# orientation as extract_posterior_inclusion_probabilities (zero
# diagonal, variable names).
# ------------------------------------------------------------------
prior_pip_matrix_from_classes = function(spec, class_values) {
  d = spec$data
  data_columnnames = d$data_columnnames
  if(identical(spec$model_type, "mixed_mrf")) {
    num_disc = d$num_discrete
    num_cont = d$num_continuous
    values = c(
      rep(class_values[["dd"]], num_disc * (num_disc - 1L) / 2L),
      rep(class_values[["cc"]], num_cont * (num_cont - 1L) / 2L),
      rep(class_values[["cross"]], num_disc * num_cont)
    )
    return(fill_mixed_symmetric(
      values, num_disc, num_cont,
      d$discrete_indices, d$continuous_indices,
      list(data_columnnames, data_columnnames)
    ))
  }
  num_vars = d$num_variables
  pip_matrix = matrix(0, num_vars, num_vars)
  pip_matrix[lower.tri(pip_matrix)] = class_values[["cc"]]
  pip_matrix = pip_matrix + t(pip_matrix)
  colnames(pip_matrix) = data_columnnames
  rownames(pip_matrix) = data_columnnames
  pip_matrix
}


# ------------------------------------------------------------------
# extract_prior_inclusion_probabilities() — main extractor
# ------------------------------------------------------------------

#' @title Extract Prior Inclusion Probabilities
#'
#' @description
#' Returns the prior edge-inclusion probabilities of a model fitted with
#' [bgm()] with `edge_selection = TRUE`, as a symmetric matrix in the same
#' shape and orientation as [extract_posterior_inclusion_probabilities()],
#' so the two can be combined element-wise into prior and posterior
#' inclusion odds.
#'
#' @details
#' For models without a continuous block of at least two variables
#' (ordinal MRFs, mixed models with a single continuous variable) the
#' prior inclusion probability is the edge prior's marginal: the fixed
#' probability for `bernoulli_prior()`, `alpha / (alpha + beta)` for
#' `beta_bernoulli_prior()`, and the exchangeable partition mixture of
#' the within- and between-block means for `sbm_prior()`.
#'
#' Under the joint spike-and-slab prior on a continuous block, the graph
#' marginal is reweighted by the per-graph normalizer — the
#' positive-definite-cone mass of the slab under the edge pattern,
#' further shaped by the determinant tilt — so the prior inclusion
#' probability of a continuous-continuous edge differs from the
#' edge-prior marginal at any `delta`, including `delta = 0`. For
#' `bernoulli_prior()` and `beta_bernoulli_prior()` it is read from the
#' same cached normalizing-constant table that corrects the fit's
#' hyperparameter updates (see [bgm()]); a fit whose slab family has no
#' table (beta-prime) falls back to a prior-only chain. For
#' `sbm_prior()` the reweighting also shifts the partition, and the
#' probabilities are estimated by a prior-only chain run with the fit's
#' own prior and correction settings; the estimate is cached on the fit,
#' and `recompute = TRUE` re-runs it.
#'
#' In mixed models only continuous-continuous edges live in the joint
#' block, so the matrix carries up to three distinct values:
#' discrete-discrete, continuous-continuous, and cross edges. Variables
#' are a priori exchangeable within type, which makes the probabilities
#' constant within these classes.
#'
#' @param bgms_object A fitted model object of class `bgms` from [bgm()]
#'   run with `edge_selection = TRUE`.
#' @param iter Integer. Post-warmup iterations for the prior-only chain,
#'   when one is needed. Default `4000`.
#' @param warmup Integer. Warmup iterations for the prior-only chain.
#'   Default `1000`.
#' @param recompute Logical. Re-run the prior-only chain even when a
#'   cached estimate is present. Default `FALSE`.
#'
#' @return A symmetric matrix of prior inclusion probabilities with the
#'   variable names as row and column names and a zero diagonal, matching
#'   [extract_posterior_inclusion_probabilities()].
#'
#' @seealso [extract_posterior_inclusion_probabilities()], [bgm()]
#' @family extractors
#' @export
extract_prior_inclusion_probabilities = function(bgms_object,
                                                 iter = 4000L,
                                                 warmup = 1000L,
                                                 recompute = FALSE) {
  UseMethod("extract_prior_inclusion_probabilities")
}


#' @inheritParams extract_prior_inclusion_probabilities
#' @exportS3Method
#' @noRd
extract_prior_inclusion_probabilities.bgms = function(bgms_object,
                                                      iter = 4000L,
                                                      warmup = 1000L,
                                                      recompute = FALSE) {
  spec = get_fit_spec(bgms_object)
  if(is.null(spec)) {
    stop(
      "extract_prior_inclusion_probabilities(): fit has no embedded ",
      "specification. Re-fit with the current bgms version."
    )
  }
  if(!spec$model_type %in% c("ggm", "omrf", "mixed_mrf")) {
    stop(
      "extract_prior_inclusion_probabilities() supports bgm() fits ",
      "(GGM, ordinal, mixed). Got model_type = '", spec$model_type, "'."
    )
  }
  if(!isTRUE(spec$prior$edge_selection)) {
    stop(
      "To extract prior inclusion probabilities, run bgm() with ",
      "edge_selection = TRUE."
    )
  }

  # Heterogeneous per-edge Bernoulli probabilities pass through untouched
  # when no joint continuous block reweights the graph.
  p = spec$prior
  num_cont = if(identical(spec$model_type, "mixed_mrf")) {
    spec$data$num_continuous
  } else if(identical(spec$model_type, "ggm")) {
    spec$data$num_variables
  } else {
    0L
  }
  if(identical(p$edge_prior, "Bernoulli") && num_cont < 2L) {
    offdiag = p$inclusion_probability[upper.tri(p$inclusion_probability)]
    if(length(unique(offdiag)) > 1L) {
      pip_matrix = p$inclusion_probability
      diag(pip_matrix) = 0
      colnames(pip_matrix) = spec$data$data_columnnames
      rownames(pip_matrix) = spec$data$data_columnnames
      return(pip_matrix)
    }
  }

  cache = get_fit_cache(bgms_object)
  class_values = if(!recompute && !is.null(cache) &&
    !is.null(cache$prior_inclusion_class_values)) {
    cache$prior_inclusion_class_values
  } else {
    values = prior_pip_class_values(spec, iter = iter, warmup = warmup)
    if(!is.null(cache)) {
      cache$prior_inclusion_class_values = values
    }
    values
  }

  prior_pip_matrix_from_classes(spec, class_values)
}
