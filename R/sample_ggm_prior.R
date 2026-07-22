#' @title Sample from the GGM (Partial-Association) Prior
#'
#' @description
#' Draws from the prior of a Gaussian graphical model. The likelihood is
#' omitted (\eqn{n = 0}, \eqn{S = 0}), so the chain targets the prior alone.
#' Three specifications are supported via the \code{spec} argument:
#' \itemize{
#'   \item \code{"conditional"} (default): fix a graph \eqn{\Gamma} and
#'     sample \eqn{K \mid \Gamma} via the same theta-space NUTS sampler
#'     that drives \code{\link{bgm}} for continuous data. The chain
#'     targets \eqn{p(K \mid \Gamma) \propto \mathrm{slab}(K) \cdot
#'     \mathrm{diag}(K) \cdot |K|^{\delta} \cdot \mathbf{1}\{K \in
#'     \mathcal{M}^{+}(\Gamma)\} / Z(\Gamma)}.
#'   \item \code{"joint"}: sample \eqn{(K, \Gamma)} jointly from the
#'     un-normalised joint prior \eqn{p(K, \Gamma) \propto
#'     \mathrm{slab}(K) \cdot \mathrm{diag}(K) \cdot |K|^{\delta} \cdot
#'     \mathbf{1}\{K \in \mathcal{M}^{+}(\Gamma)\} \cdot \pi(\Gamma)}.
#'     Uses the adaptive-Metropolis MH chain from \code{\link{bgm}} with
#'     edge selection on and the likelihood off, so the marginal on
#'     \eqn{\Gamma} is \eqn{\pi(\Gamma) \cdot Z(\Gamma)} (joint
#'     specification, not hierarchical). Useful for simulation-based
#'     calibration of \code{\link{bgm}}'s default sampler.
#'   \item \code{"hierarchical"}: sample \eqn{(K, \Gamma)} from the
#'     hierarchical specification \eqn{p(\Gamma) \, p(K \mid \Gamma)} with
#'     \eqn{p(K \mid \Gamma)} normalized per graph, so the marginal on
#'     \eqn{\Gamma} is exactly the edge prior \eqn{\pi(\Gamma)}. The
#'     per-graph normalizer ratio in each between-edge move is evaluated
#'     by the deterministic local Z-ratio approximation. Requires
#'     \code{normal_prior()} or \code{cauchy_prior()} interactions.
#' }
#'
#' @details
#' The priors are specified on the partial-association scale
#' \eqn{K_{yy} = -K/2}: \code{interaction_prior} acts on
#' \eqn{K_{yy,ij} = -K_{ij}/2}, and \code{precision_scale_prior} acts on
#' \eqn{-K_{yy,ii} = K_{ii}/2}. The same convention is used by
#' \code{\link{bgm}} and by the continuous block of the mixed-MRF model, so
#' a prior argument passed here means the same distribution it would mean
#' there. Output samples are reported as entries of \eqn{K}; convert with
#' \eqn{K_{yy} = -K/2} if you want them on the partial-association scale.
#'
#' When \code{spec = "conditional"} and \code{edge_indicators} is supplied,
#' off-diagonals at excluded positions are constrained to zero throughout
#' the chain. \code{edge_indicators} is ignored when \code{spec = "joint"}
#' (the chain samples \eqn{\Gamma}).
#'
#' When \code{spec = "joint"}, the chain is initialized from an ancestral
#' draw of the edge prior (hyperparameters from their prior, then
#' indicators given the hyperparameters), keyed to \code{seed}. Under a
#' hierarchical edge prior the inclusion parameter and the graph density
#' are coupled, and a full-graph start can pin both near 1 for a large
#' number of sweeps in zero-evidence chains.
#'
#' @param p Integer. Dimension of the precision matrix (\eqn{p \ge 2}).
#' @param n_samples Integer. Number of post-warmup draws to keep.
#' @param n_warmup Integer. NUTS warmup iterations. Default \code{2000}.
#' @param interaction_prior A \code{bgms_parameter_prior} for the
#'   partial-association off-diagonals \eqn{K_{yy,ij} = -K_{ij}/2}. Use
#'   \code{\link{cauchy_prior}()} or \code{\link{normal_prior}()};
#'   \code{\link{beta_prime_prior}()} is not supported here. Default:
#'   \code{cauchy_prior(scale = 2.5)} (i.e. \eqn{K_{ij}} has an implied
#'   \eqn{\textrm{Cauchy}(0, 5)} prior).
#' @param precision_scale_prior A \code{bgms_scale_prior} for
#'   \eqn{K_{ii}/2}. Use \code{\link{gamma_prior}()} or
#'   \code{\link{exponential_prior}()}. Both accept the rate in the raw
#'   frame (\code{rate}) or the standardized frame (\code{eta}; the raw
#'   rate is derived as \code{eta / s} for interaction-prior scale
#'   \code{s}). Default: \code{exponential_prior(eta = 1)}; with the
#'   default \code{cauchy_prior(scale = 2.5)} interaction prior this
#'   resolves to \eqn{K_{ii}/2 \sim \textrm{Exponential}(0.4)}.
#' @param step_size Positive numeric. Initial NUTS step size used to seed
#'   dual-averaging adaptation. Default \code{0.1}. Used only for
#'   \code{spec = "conditional"} (NUTS path); ignored for the
#'   \code{"joint"} MH path.
#' @param max_depth Integer. Maximum NUTS tree depth. Default \code{10}.
#'   Used only for \code{spec = "conditional"}.
#' @param seed Integer. RNG seed for the chain. Default \code{1L}.
#' @param verbose Logical. If \code{TRUE} (default), print a progress bar.
#' @param edge_indicators Optional integer \eqn{p \times p} matrix with
#'   \code{1} = edge included, \code{0} = excluded. Must be symmetric with
#'   \code{1}s on the diagonal. Default: full graph (all edges included).
#'   Used only for \code{spec = "conditional"} (the chain samples
#'   \eqn{K \mid \Gamma}); ignored for \code{spec = "joint"}.
#' @param spec One of \code{"conditional"} (default, sample
#'   \eqn{K \mid \Gamma} at fixed \eqn{\Gamma}), \code{"joint"} (sample
#'   \eqn{(K, \Gamma)} jointly from the un-normalised joint prior), or
#'   \code{"hierarchical"} (sample \eqn{(K, \Gamma)} from the per-graph
#'   normalized specification via the Z-ratio approximation).
#' @param edge_inclusion_prob Probability in \eqn{(0, 1)} for the
#'   Bernoulli edge prior used when \code{spec = "joint"}. Default
#'   \code{0.5}. Ignored when \code{spec = "conditional"}.
#' @param update_method One of \code{"adaptive-metropolis"} (default) or
#'   \code{"gibbs"}. Sampler driving the \code{spec = "joint"} chain; the
#'   Gibbs chain uses the conjugate row and edge updates and needs no
#'   proposal tuning. Ignored when \code{spec = "conditional"} (NUTS).
#' @param edge_prior An edge prior specification object from
#'   \code{\link{bernoulli_prior}()}, \code{\link{beta_bernoulli_prior}()},
#'   or \code{\link{sbm_prior}()}, or \code{NULL} (default) for a Bernoulli
#'   prior with probability \code{edge_inclusion_prob}. Only for
#'   \code{spec = "joint"}.
#' @param apply_correction Logical. For the hierarchical edge priors
#'   (\code{beta_bernoulli_prior()}, \code{sbm_prior()}), apply the
#'   normalizing-constant correction to the hyperparameter updates (default
#'   \code{TRUE}; the correction table is built from the tilted prior
#'   sampler and cached across calls). With \code{FALSE} the plain conjugate
#'   updates are used, whose hyperparameter marginals do not match the
#'   hyperpriors under the determinant tilt.
#' @param zratio_diagnostics Logical (default \code{TRUE}). Only for
#'   \code{spec = "hierarchical"}: run the trust gauge
#'   (\code{\link{summarize_zratio_gauge}}) on the returned chain and attach
#'   the result; detected issues are printed when \code{verbose}. The gauge
#'   redoes a subset of the chain's edge decisions with the exact
#'   calculation and reports two alarms: how often the decision outcome
#'   differs (\code{flip_rate}), and the projected distortion of the mean
#'   inclusion probability from the measured error under the edge prior's
#'   feedback (\code{harm_pred}). Evidence-free sampling is the regime where
#'   the second alarm matters: a small consistent error can shift the graph
#'   marginal without flipping individual decisions.
#' @param delta Non-negative numeric, or \code{NULL} for the dimension-
#'   adaptive default. Determinant-tilt exponent: multiplies the prior
#'   by \eqn{|K|^{\delta}}, softly repelling the chain from the
#'   positive-definite cone boundary. \code{delta = NULL} (default)
#'   auto-resolves to \eqn{0.5 \log(p)}, the simple form of the
#'   dimension-adaptive rule \eqn{\delta(p) = c \log p} with
#'   \eqn{c \in (0.3, 0.6)} (Marsman et al., in preparation). Pass
#'   \code{delta = 0} for the untilted prior or a non-negative numeric to
#'   override.
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{K_offdiag}}{Numeric matrix of size
#'       \code{n_samples} x \code{p * (p - 1) / 2} containing the upper-triangle
#'       off-diagonal entries of \eqn{K} for each draw, in row-major order
#'       (the upper triangle traversed by row)
#'       \eqn{(K_{12}, K_{13}, \ldots, K_{1p}, K_{23}, K_{24}, \ldots, K_{2p}, K_{34}, \ldots)}. Under
#'       \code{spec = "conditional"}, excluded edges are returned as
#'       \code{0}; under \code{spec = "joint"}, off-diagonals at excluded
#'       edges are sampled at \code{0} per the inclusion indicator.}
#'     \item{\code{K_diag}}{Numeric matrix of size
#'       \code{n_samples} x \code{p} containing the diagonal entries
#'       \eqn{K_{11}, \ldots, K_{pp}}.}
#'     \item{\code{offdiag_names}}{Character vector of length
#'       \code{p * (p - 1) / 2} naming the columns of \code{K_offdiag}
#'       (e.g. \code{"K_1_2"}).}
#'     \item{\code{diag_names}}{Character vector of length \code{p} naming
#'       the columns of \code{K_diag}.}
#'     \item{\code{edge_indicators}}{Under \code{spec = "conditional"}, the
#'       \code{p x p} integer matrix of fixed inclusion indicators used
#'       (full graph if not supplied). Under \code{spec = "joint"}, an
#'       \code{n_samples x p(p-1)/2} integer matrix of sampled
#'       \eqn{\Gamma_{ij}} indicators (column order matches
#'       \code{K_offdiag}).}
#'     \item{\code{theta}}{Only with \code{beta_bernoulli_prior()}: numeric
#'       vector of length \code{n_samples} with the sampled inclusion
#'       probability.}
#'     \item{\code{allocations}}{Only with \code{sbm_prior()}: integer
#'       matrix (\code{n_samples x p}) of sampled cluster allocations
#'       (1-based).}
#'     \item{\code{zratio_diagnostics}}{Only with
#'       \code{spec = "hierarchical"} and \code{zratio_diagnostics = TRUE}:
#'       the trust-gauge summary from
#'       \code{\link{summarize_zratio_gauge}}.}
#'   }
#'
#' @seealso \code{\link{cauchy_prior}}, \code{\link{normal_prior}},
#'   \code{\link{gamma_prior}}, \code{\link{exponential_prior}},
#'   \code{\link{bgm}}
#'
#' @examples
#' \donttest{
#' # Default Cauchy(0, 2.5) off-diagonal, Gamma(1, 1) diagonal, p = 4.
#' draws = sample_ggm_prior(
#'   p = 4, n_samples = 200, n_warmup = 200,
#'   verbose = FALSE
#' )
#' dim(draws$K_offdiag) # 200 x 6
#' colnames(draws$K_offdiag) = draws$offdiag_names
#' head(draws$K_offdiag)
#'
#' # Sparser graph: drop the (1, 4) edge.
#' E = matrix(1L, 4, 4)
#' E[1, 4] = E[4, 1] = 0L
#' draws = sample_ggm_prior(
#'   p = 4, n_samples = 200, n_warmup = 200,
#'   edge_indicators = E, verbose = FALSE
#' )
#' colnames(draws$K_offdiag) = draws$offdiag_names
#' all(draws$K_offdiag[, "K_1_4"] == 0) # TRUE
#' }
#' @export
sample_ggm_prior = function(
  p,
  n_samples,
  n_warmup = 2e3,
  interaction_prior = cauchy_prior(scale = 2.5),
  precision_scale_prior = exponential_prior(eta = 1),
  step_size = 0.1,
  max_depth = 10L,
  seed = 1L,
  verbose = TRUE,
  edge_indicators = NULL,
  delta = NULL,
  spec = c("conditional", "joint", "hierarchical"),
  edge_inclusion_prob = 0.5,
  update_method = c("adaptive-metropolis", "gibbs"),
  edge_prior = NULL,
  apply_correction = TRUE,
  zratio_diagnostics = TRUE
) {
  spec = match.arg(spec)
  update_method = match.arg(update_method)
  ep = if(is.null(edge_prior)) {
    NULL
  } else {
    unpack_indicator_prior(edge_prior, num_variables = as.integer(p))
  }
  if(spec == "conditional" && !is.null(ep) &&
    !identical(ep$edge_prior, "Bernoulli")) {
    stop(
      "Hierarchical edge priors require spec = \"joint\" or ",
      "\"hierarchical\" (the conditional spec fixes the graph)."
    )
  }
  if(!is.logical(apply_correction) || length(apply_correction) != 1L ||
    is.na(apply_correction)) {
    stop("'apply_correction' must be TRUE or FALSE.")
  }
  if(!is.logical(zratio_diagnostics) || length(zratio_diagnostics) != 1L ||
    is.na(zratio_diagnostics)) {
    stop("'zratio_diagnostics' must be TRUE or FALSE.")
  }
  validate_integer(p, "p", min_value = 2L)
  validate_integer(n_samples, "n_samples", min_value = 1L)
  validate_integer(n_warmup, "n_warmup", min_value = 0L)
  validate_integer(max_depth, "max_depth", min_value = 1L)
  validate_finite_scalar(step_size, "step_size", positive = TRUE)
  validate_integer(seed, "seed", min_value = 0L)
  if(is.null(delta)) {
    delta = 0.5 * log(p)
  }
  if(!is.numeric(delta) || length(delta) != 1L || is.na(delta) ||
    !is.finite(delta) || delta < 0) {
    stop("'delta' must be a single finite non-negative numeric, or NULL.")
  }
  if(!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.")
  }
  if(!is.numeric(edge_inclusion_prob) || length(edge_inclusion_prob) != 1L ||
    is.na(edge_inclusion_prob) || edge_inclusion_prob <= 0 ||
    edge_inclusion_prob >= 1) {
    stop("'edge_inclusion_prob' must be a single numeric in (0, 1).")
  }

  ip = unpack_interaction_prior(interaction_prior)
  if(identical(ip$interaction_prior_type, "beta-prime")) {
    stop(
      "beta_prime_prior() is not supported for 'interaction_prior' in ",
      "sample_ggm_prior(). Use cauchy_prior() or normal_prior()."
    )
  }
  sp = unpack_scale_prior(precision_scale_prior)
  sp$scale_rate = resolve_scale_rate(
    sp$scale_rate, sp$scale_eta, ip$pairwise_scale
  )

  edge_indicators = validate_ggm_prior_edge_indicators(edge_indicators, p)

  if(spec == "conditional") {
    return(sample_ggm_prior_cpp(
      p                        = as.integer(p),
      n_samples                = as.integer(n_samples),
      n_warmup                 = as.integer(n_warmup),
      pairwise_scale           = ip$pairwise_scale,
      interaction_prior_type   = ip$interaction_prior_type,
      scale_prior_type         = sp$scale_prior_type,
      gamma_shape              = sp$scale_shape,
      gamma_rate               = sp$scale_rate,
      step_size                = step_size,
      max_depth                = as.integer(max_depth),
      seed                     = as.integer(seed),
      verbose                  = verbose,
      edge_indicators_nullable = edge_indicators,
      delta                    = as.numeric(delta)
    ))
  }

  # spec == "joint" / "hierarchical": drive the bgm() MH chain with edge
  # selection on and zero data (n = 0, S = 0). The joint chain targets the
  # un-normalised joint prior; the hierarchical chain adds the per-edge
  # Z-ratio to the between-edge moves so the graph marginal is pi(Gamma).
  inputFromR = list(
    n                      = 0L,
    suf_stat               = matrix(0, p, p),
    pairwise_scale         = ip$pairwise_scale,
    interaction_prior_type = ip$interaction_prior_type,
    scale_prior_type       = sp$scale_prior_type,
    scale_shape            = sp$scale_shape,
    scale_rate             = sp$scale_rate
  )

  if(is.null(ep)) {
    ep = unpack_indicator_prior(
      bernoulli_prior(edge_inclusion_prob),
      num_variables = as.integer(p)
    )
  }
  correction = NULL
  zratio = NULL
  if(spec == "hierarchical") {
    # Hierarchical spec p(K | Gamma) = rho_Gamma(K)/Z(Gamma): the per-edge
    # Z-ratio engine carries the normalizer into the between-edge moves,
    # and the hyperparameter updates are the clean conjugate draws (no
    # C-correction on this path).
    if(!ip$interaction_prior_type %in% c("normal", "cauchy")) {
      stop(sprintf(
        paste0(
          "spec = \"hierarchical\" supports a normal or Cauchy interaction ",
          "(slab) prior. Got %s_prior(). Use interaction_prior = ",
          "normal_prior() or cauchy_prior()."
        ),
        ip$interaction_prior_type
      ))
    }
    zc = zratio_cell_constants(
      delta, ip$pairwise_scale, sp$scale_rate, sp$scale_eta,
      scale_shape = sp$scale_shape,
      slab = ip$interaction_prior_type
    )
    zratio = list(
      addc = zc$addc, tg = zc$tg, ihat = zc$ihat, ghat = zc$ghat,
      wt = zc$wt, psi0 = zc$psi0,
      delta = zc$delta, eta = zc$eta, alpha = zc$alpha, slab = zc$slab,
      gauge_sweeps = if(isTRUE(zratio_diagnostics)) 2L else 0L
    )
    # Deploy the same Option-B surface the posterior chain uses, so the prior
    # chain (SBC reference) carries an identical per-edge correction.
    surf = build_surfaces_allmc(
      zc, max_size = min(p, 44L), cores = zratio_surface_build_cores()
    )
    if(!is.null(surf)) zratio$surface = surf
  } else if(!identical(ep$edge_prior, "Bernoulli") && apply_correction) {
    table = ggm_correction_table(
      p = p, delta = delta,
      interaction_prior = interaction_prior,
      precision_scale_prior = precision_scale_prior,
      update_method = "gibbs",
      verbose = isTRUE(verbose)
    )
    correction = correction_list_from_table(table, ep$edge_prior)
  }

  results = sample_ggm(
    inputFromR = inputFromR,
    prior_inclusion_prob = ep$inclusion_probability,
    initial_edge_indicators = ggm_prior_ancestral_indicators(p, ep, seed),
    no_iter = as.integer(n_samples),
    no_warmup = as.integer(n_warmup),
    no_chains = 1L,
    edge_selection = TRUE,
    sampler_type = update_method,
    seed = as.integer(seed),
    no_threads = 1L,
    progress_type = if(verbose) 2L else 0L,
    edge_prior = ep$edge_prior,
    beta_bernoulli_alpha = ep$beta_bernoulli_alpha,
    beta_bernoulli_beta = ep$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = ep$beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = ep$beta_bernoulli_beta_between,
    dirichlet_alpha = ep$dirichlet_alpha,
    lambda = ep$lambda,
    delta = as.numeric(delta),
    edge_prior_correction = correction,
    zratio_spec = zratio
  )
  if(length(results) == 0L || isTRUE(results[[1L]]$error)) {
    msg = if(length(results) > 0L) results[[1L]]$error_msg else "empty result"
    stop("sample_ggm_prior (joint): chain failed (", msg, ")")
  }

  # Reformat the bgm() output to match the conditional spec's return shape.
  # Both `samples` and `indicator_samples` are emitted as the full upper
  # triangle in (i <= j) order, p(p+1)/2 rows per iteration. Diagonals on
  # the indicator side are always 1 and discarded here.
  upper = results[[1L]]$samples # ((p*(p+1))/2) x n_samples
  inds = results[[1L]]$indicator_samples # ((p*(p+1))/2) x n_samples
  n_edges = as.integer(p * (p - 1) / 2)

  K_offdiag = matrix(0, n_samples, n_edges)
  K_diag = matrix(0, n_samples, p)
  gamma_offdiag = matrix(0L, n_samples, n_edges)
  e = 1L
  off_idx = 1L
  for(i in seq_len(p)) {
    for(j in i:p) {
      if(i == j) {
        K_diag[, i] = upper[e, ]
      } else {
        K_offdiag[, off_idx] = upper[e, ]
        gamma_offdiag[, off_idx] = inds[e, ]
        off_idx = off_idx + 1L
      }
      e = e + 1L
    }
  }

  offdiag_names = character(n_edges)
  idx = 1L
  for(i in seq_len(p - 1L)) {
    for(j in (i + 1L):p) {
      offdiag_names[idx] = paste0("K_", i, "_", j)
      idx = idx + 1L
    }
  }
  diag_names = paste0("K_", seq_len(p), "_", seq_len(p))

  out = list(
    K_offdiag       = K_offdiag,
    K_diag          = K_diag,
    offdiag_names   = offdiag_names,
    diag_names      = diag_names,
    edge_indicators = gamma_offdiag
  )
  if(!is.null(results[[1L]]$inclusion_parameter_samples)) {
    out$theta = as.numeric(results[[1L]]$inclusion_parameter_samples)
  }
  if(!is.null(results[[1L]]$allocation_samples)) {
    out$allocations = t(results[[1L]]$allocation_samples)
  }
  if(spec == "hierarchical" && isTRUE(zratio_diagnostics)) {
    harm_inputs = zratio_harm_inputs(
      list(colMeans(gamma_offdiag)), ep$edge_prior,
      a = ep$beta_bernoulli_alpha, b = ep$beta_bernoulli_beta
    )
    out$zratio_diagnostics = summarize_zratio_gauge(
      results,
      verbose = verbose, harm_inputs = harm_inputs
    )
    if(isTRUE(verbose) && isTRUE(out$zratio_diagnostics$flagged)) {
      cat("See vignette('diagnostics') for guidance.\n")
    }
  }
  out
}

# Ancestral draw of the initial edge-indicator matrix for the joint-spec
# chain: hyperparameters from their prior, then indicators given the
# hyperparameters. Runs in an RNG scope keyed to `seed` and restores the
# caller's RNG state on exit.
ggm_prior_ancestral_indicators = function(p, ep, seed) {
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
  set.seed(seed)

  prob = switch(ep$edge_prior,
    "Bernoulli" = ep$inclusion_probability,
    "Beta-Bernoulli" = matrix(
      rbeta(1, ep$beta_bernoulli_alpha, ep$beta_bernoulli_beta), p, p
    ),
    "Stochastic-Block" = {
      z = ancestral_mfm_sbm_partition(p, ep$lambda, ep$dirichlet_alpha)
      ancestral_sbm_pair_probabilities(
        z,
        ep$beta_bernoulli_alpha, ep$beta_bernoulli_beta,
        ep$beta_bernoulli_alpha_between, ep$beta_bernoulli_beta_between
      )
    }
  )

  g = matrix(0L, p, p)
  upper = upper.tri(g)
  g[upper] = as.integer(runif(sum(upper)) < prob[upper])
  g = g + t(g)
  diag(g) = 1L
  g
}

# Ancestral draw from the MFM-SBM hyperprior: shifted-Poisson component
# count, Dirichlet weights, and allocations.
ancestral_mfm_sbm_partition = function(p, lambda, dirichlet_alpha) {
  num_components = rpois(1, lambda) + 1L
  w = rgamma(num_components, dirichlet_alpha)
  sample.int(num_components, p, replace = TRUE, prob = w)
}

# Per-pair inclusion probabilities implied by an allocation vector, with
# within-block and between-block Beta draws for each block pair.
ancestral_sbm_pair_probabilities = function(
  z, a_within, b_within, a_between, b_between
) {
  labs = sort(unique(z))
  nl = length(labs)
  th_rs = matrix(0, nl, nl)
  for(r in seq_len(nl)) {
    for(s in r:nl) {
      v = if(r == s) {
        rbeta(1, a_within, b_within)
      } else {
        rbeta(1, a_between, b_between)
      }
      th_rs[r, s] = v
      th_rs[s, r] = v
    }
  }
  zi = match(z, labs)
  p = length(z)
  prob = matrix(0.5, p, p)
  for(i in seq_len(p - 1)) {
    for(j in (i + 1):p) {
      prob[i, j] = th_rs[zi[i], zi[j]]
      prob[j, i] = prob[i, j]
    }
  }
  prob
}

validate_integer = function(x, name, min_value = 1L) {
  if(!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("'%s' must be a single finite integer.", name))
  }
  if(x != as.integer(x)) {
    stop(sprintf("'%s' must be an integer (got %s).", name, format(x)))
  }
  if(x < min_value) {
    stop(sprintf("'%s' must be >= %d.", name, as.integer(min_value)))
  }
  invisible(as.integer(x))
}

validate_finite_scalar = function(x, name, positive = FALSE) {
  if(!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("'%s' must be a single finite numeric.", name))
  }
  if(positive && x <= 0) {
    stop(sprintf("'%s' must be positive.", name))
  }
  invisible(x)
}

validate_ggm_prior_edge_indicators = function(edge_indicators, p) {
  if(is.null(edge_indicators)) {
    return(NULL)
  }
  if(!is.matrix(edge_indicators) ||
    nrow(edge_indicators) != p || ncol(edge_indicators) != p) {
    stop("'edge_indicators' must be a p x p matrix.")
  }
  if(any(is.na(edge_indicators))) {
    stop("'edge_indicators' must not contain NA values.")
  }
  vals = as.integer(edge_indicators)
  if(any(!vals %in% c(0L, 1L))) {
    stop("'edge_indicators' must contain only 0 or 1.")
  }
  E = matrix(vals, nrow = p, ncol = p)
  if(!isTRUE(all.equal(E, t(E)))) {
    stop("'edge_indicators' must be symmetric.")
  }
  if(!all(diag(E) == 1L)) {
    stop("'edge_indicators' must have 1s on the diagonal.")
  }
  E
}
