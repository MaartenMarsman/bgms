#' @title Sample from the GGM (Partial-Association) Prior
#'
#' @description
#' Draws from the prior of a Gaussian graphical model. The likelihood is
#' omitted (\eqn{n = 0}, \eqn{S = 0}), so the chain targets the prior alone.
#' Two specifications are supported via the \code{spec} argument:
#' \itemize{
#'   \item \code{"conditional"} (default): fix a graph \eqn{\Gamma} and
#'     sample \eqn{K \mid \Gamma} via the same constrained NUTS sampler
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
#'   \code{s}). Default: \code{gamma_prior(shape = 1, eta = 1)}; with the
#'   default \code{cauchy_prior(scale = 2.5)} interaction prior this
#'   resolves to \eqn{K_{ii}/2 \sim \textrm{Gamma}(1, 0.4)}.
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
#'   \eqn{K \mid \Gamma} at fixed \eqn{\Gamma}) or \code{"joint"} (sample
#'   \eqn{(K, \Gamma)} jointly from the un-normalised joint prior).
#' @param edge_inclusion_prob Probability in \eqn{(0, 1)} for the
#'   Bernoulli edge prior used when \code{spec = "joint"}. Default
#'   \code{0.5}. Ignored when \code{spec = "conditional"}.
#' @param update_method One of \code{"adaptive-metropolis"} (default) or
#'   \code{"gibbs"}. Sampler driving the \code{spec = "joint"} chain; the
#'   Gibbs chain uses the conjugate row and edge updates and needs no
#'   proposal tuning. Ignored when \code{spec = "conditional"} (NUTS).
#' @param delta Non-negative numeric, or \code{NULL} for the dimension-
#'   adaptive default. Determinant-tilt exponent: multiplies the prior
#'   by \eqn{|K|^{\delta}}, softly repelling the chain from the
#'   positive-definite cone boundary. \code{delta = NULL} (default)
#'   auto-resolves to \eqn{0.5 \log(p)}, the simple form of the
#'   dimension-adaptive rule \eqn{\delta(p) = c \log p} with
#'   \eqn{c \in (0.3, 0.6)} discussed in the companion paper on
#'   determinant-tilted spike-and-slab priors (Marsman et al., in
#'   preparation). Pass \code{delta = 0} for the untilted prior (the
#'   companion-paper baseline) or a non-negative numeric to override.
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
  precision_scale_prior = gamma_prior(shape = 1, eta = 1),
  step_size = 0.1,
  max_depth = 10L,
  seed = 1L,
  verbose = TRUE,
  edge_indicators = NULL,
  delta = NULL,
  spec = c("conditional", "joint"),
  edge_inclusion_prob = 0.5,
  update_method = c("adaptive-metropolis", "gibbs")
) {
  spec = match.arg(spec)
  update_method = match.arg(update_method)
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

  # spec == "joint": drive the bgm() MH chain with edge selection on and
  # zero data (n = 0, S = 0). The chain targets the un-normalised joint
  # prior p(K, Gamma) and produces (K, Gamma) draws.
  inputFromR = list(
    n                      = 0L,
    suf_stat               = matrix(0, p, p),
    pairwise_scale         = ip$pairwise_scale,
    interaction_prior_type = ip$interaction_prior_type,
    scale_prior_type       = sp$scale_prior_type,
    scale_shape            = sp$scale_shape,
    scale_rate             = sp$scale_rate
  )

  results = sample_ggm(
    inputFromR              = inputFromR,
    prior_inclusion_prob    = matrix(edge_inclusion_prob, p, p),
    initial_edge_indicators = matrix(1L, p, p),
    no_iter                 = as.integer(n_samples),
    no_warmup               = as.integer(n_warmup),
    no_chains               = 1L,
    edge_selection          = TRUE,
    sampler_type            = update_method,
    seed                    = as.integer(seed),
    no_threads              = 1L,
    progress_type           = if(verbose) 2L else 0L,
    edge_prior              = "Bernoulli",
    delta                   = as.numeric(delta)
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

  list(
    K_offdiag       = K_offdiag,
    K_diag          = K_diag,
    offdiag_names   = offdiag_names,
    diag_names      = diag_names,
    edge_indicators = gamma_offdiag
  )
}


# Internal helpers -------------------------------------------------------------

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
