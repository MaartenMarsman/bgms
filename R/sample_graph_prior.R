# ==============================================================================
# Graph-level prior samplers
# ==============================================================================
#
# sample_graph_prior(): draws (hyperparameters, edge indicators) from the
# graph level of the spike-and-slab prior, under either specification of
# how the precision prior composes with the graph:
#
#   hierarchical: p(hyper) . pi(Gamma | hyper) -- ancestral and exact
#   joint:        p(hyper) . q(Gamma | hyper), with
#                 q(Gamma | hyper) proportional to Z(Gamma) pi(Gamma | hyper)
#                 -- the determinant-tilted law, sampled by the zero-data
#                 (K, Gamma) chain with K discarded from the output
#
# sample_sbm_prior(): ancestral draws of (allocations, pair probabilities)
# from the MFM-SBM edge-prior hyperprior.
#
# Edge indicators are returned in row-major upper-triangle order (i < j),
# matching the K_offdiag column order of sample_ggm_prior().
# ==============================================================================


# Run fn() inside an RNG scope keyed to `seed`, restoring the caller's RNG
# state on exit (same idiom as ggm_prior_ancestral_indicators).
with_graph_prior_seed = function(seed, fn) {
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
  fn()
}


# Row-major upper-triangle pair indices (i < j): a (E x 2) matrix with
# E = p(p-1)/2 rows, ordered (1,2), (1,3), ..., (1,p), (2,3), ...
graph_pair_indices = function(p) {
  ii = rep(seq_len(p - 1L), times = (p - 1L):1L)
  jj = unlist(lapply(seq_len(p - 1L), function(i) (i + 1L):p))
  cbind(ii, jj)
}


# Resolve the optional conditioning arguments to a fixed p x p pair
# probability matrix, or NULL when no conditioning was requested.
graph_prior_conditioning = function(ep, p, theta, allocations, block_probs) {
  if(!is.null(theta) && (!is.null(allocations) || !is.null(block_probs))) {
    stop(
      "Supply either 'theta' or ('allocations', 'block_probs') to condition ",
      "the graph prior, not both."
    )
  }
  if(!is.null(theta)) {
    if(identical(ep$edge_prior, "Stochastic-Block")) {
      stop(
        "Condition the Stochastic-Block prior with 'allocations' and ",
        "'block_probs'; 'theta' conditions the Bernoulli and Beta-Bernoulli ",
        "priors."
      )
    }
    if(!is.numeric(theta) || length(theta) != 1L || is.na(theta) ||
      theta <= 0 || theta >= 1) {
      stop("'theta' must be a single numeric in (0, 1).")
    }
    return(matrix(theta, p, p))
  }
  if(!is.null(allocations) || !is.null(block_probs)) {
    if(is.null(allocations) || is.null(block_probs)) {
      stop(
        "Conditioning the Stochastic-Block prior needs both 'allocations' ",
        "and 'block_probs'."
      )
    }
    if(!identical(ep$edge_prior, "Stochastic-Block")) {
      stop(
        "'allocations'/'block_probs' condition the Stochastic-Block prior; ",
        "use edge_prior = sbm_prior()."
      )
    }
    allocations = as.integer(allocations)
    if(length(allocations) != p || anyNA(allocations) ||
      any(allocations < 1L)) {
      stop("'allocations' must be a length-p vector of 1-based block labels.")
    }
    if(!is.matrix(block_probs) || nrow(block_probs) != ncol(block_probs) ||
      !isSymmetric(unname(block_probs)) ||
      any(block_probs <= 0) || any(block_probs >= 1)) {
      stop(
        "'block_probs' must be a symmetric matrix with entries in (0, 1)."
      )
    }
    if(max(allocations) > nrow(block_probs)) {
      stop("'allocations' labels exceed the size of 'block_probs'.")
    }
    prob = matrix(0.5, p, p)
    for(i in seq_len(p - 1L)) {
      for(j in (i + 1L):p) {
        prob[i, j] = block_probs[allocations[i], allocations[j]]
        prob[j, i] = prob[i, j]
      }
    }
    return(prob)
  }
  NULL
}


#' @title Sample from the Graph Prior
#'
#' @description
#' Draws edge-inclusion indicators, together with any edge-prior
#' hyperparameters, from the graph level of the spike-and-slab prior used by
#' \code{\link{bgm}} for models with continuous variables. The \code{spec}
#' argument selects how the precision prior composes with the graph:
#' \itemize{
#'   \item \code{"hierarchical"} (default): the graph marginal is exactly the
#'     edge prior, \eqn{p(\mathrm{hyper}) \, \pi(\Gamma \mid \mathrm{hyper})}.
#'     Sampling is ancestral and exact: hyperparameters from their prior,
#'     then independent pair flips.
#'   \item \code{"joint"}: the graph marginal is reweighted by the per-graph
#'     normalizer of the determinant-tilted precision prior,
#'     \eqn{q(\Gamma \mid \mathrm{hyper}) \propto Z(\Gamma) \,
#'     \pi(\Gamma \mid \mathrm{hyper})}. Sampling runs the zero-data
#'     \eqn{(K, \Gamma)} chain of \code{\link{sample_ggm_prior}} and discards
#'     \eqn{K}; with a Beta-Bernoulli or Stochastic-Block prior the
#'     hyperparameter updates apply the normalizing-constant correction, so
#'     the first call for a model cell may build the correction table
#'     (cached across fits).
#' }
#'
#' @details
#' The optional conditioning arguments fix the edge-prior hyperparameters
#' instead of sampling them: \code{theta} fixes the inclusion probability of
#' a Bernoulli or Beta-Bernoulli prior, and \code{allocations} plus
#' \code{block_probs} fix the block structure of a Stochastic-Block prior
#' (reducing it to independent pair flips at the given block probabilities).
#'
#' Under \code{spec = "joint"} the tilted graph law depends on the precision
#' prior through its normalizer, so \code{interaction_prior},
#' \code{precision_scale_prior}, and \code{delta} are part of the graph law;
#' they are ignored under \code{spec = "hierarchical"}.
#'
#' @param p Integer. Number of nodes (\eqn{p \ge 2}).
#' @param n_samples Integer. Number of prior draws.
#' @param edge_prior An edge prior specification object:
#'   \code{\link{bernoulli_prior}()}, \code{\link{beta_bernoulli_prior}()},
#'   or \code{\link{sbm_prior}()}. Default \code{bernoulli_prior(0.5)}.
#' @param spec One of \code{"hierarchical"} (default) or \code{"joint"}.
#' @param interaction_prior A \code{\link{cauchy_prior}()} or
#'   \code{\link{normal_prior}()} for the pairwise (slab) part of the
#'   precision prior. Used only when \code{spec = "joint"}.
#' @param precision_scale_prior A \code{\link{gamma_prior}()} or
#'   \code{\link{exponential_prior}()} for the precision diagonal. Used only
#'   when \code{spec = "joint"}.
#' @param delta Non-negative numeric or \code{NULL} (default): determinant
#'   tilt exponent; \code{NULL} resolves to \eqn{0.5 \log(p)}. Used only
#'   when \code{spec = "joint"}.
#' @param theta Optional numeric in (0, 1): fix the inclusion probability of
#'   a Bernoulli or Beta-Bernoulli edge prior instead of sampling it.
#' @param allocations Optional integer vector of length \code{p} with 1-based
#'   block labels: fix the Stochastic-Block allocation instead of sampling it.
#'   Requires \code{block_probs}.
#' @param block_probs Optional symmetric matrix with entries in (0, 1): the
#'   block-pair inclusion probabilities that go with \code{allocations}.
#' @param n_warmup Integer. Warmup iterations of the zero-data chain. Used
#'   only when \code{spec = "joint"}. Default \code{2e3}.
#' @param seed Integer. Seed for the draw; the caller's RNG state is
#'   restored on exit.
#' @param verbose Logical. Print progress of the zero-data chain and of a
#'   correction-table build. Default \code{TRUE}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{edge_indicators}}{Integer matrix
#'       (\code{n_samples x p(p-1)/2}) of edge-inclusion indicators, columns
#'       in row-major upper-triangle order (matching
#'       \code{sample_ggm_prior()}'s \code{K_offdiag}).}
#'     \item{\code{pair_names}}{Character vector labeling the columns as
#'       \code{"i-j"}.}
#'     \item{\code{theta}}{Only with an unconditioned
#'       \code{beta_bernoulli_prior()}: numeric vector of sampled inclusion
#'       probabilities.}
#'     \item{\code{allocations}}{Only with an unconditioned
#'       \code{sbm_prior()}: integer matrix (\code{n_samples x p}) of sampled
#'       block allocations.}
#'     \item{\code{spec}, \code{edge_prior}, \code{p}}{The specification,
#'       edge-prior family, and node count of the draw.}
#'   }
#'
#' @examples
#' # Hierarchical spec: the graph marginal is exactly the edge prior.
#' g = sample_graph_prior(
#'   p = 6, n_samples = 200,
#'   edge_prior = bernoulli_prior(0.3), seed = 11
#' )
#' mean(g$edge_indicators) # about 0.3
#'
#' # Beta-Bernoulli: inclusion probabilities are sampled alongside.
#' g = sample_graph_prior(
#'   p = 6, n_samples = 200,
#'   edge_prior = beta_bernoulli_prior(2, 4), seed = 11
#' )
#' mean(g$theta) # about 1/3
#'
#' \donttest{
#' # Joint spec: the graph law carries the per-graph normalizer Z(Gamma).
#' g = sample_graph_prior(
#'   p = 6, n_samples = 500,
#'   edge_prior = bernoulli_prior(0.3), spec = "joint",
#'   interaction_prior = normal_prior(scale = 0.5),
#'   precision_scale_prior = exponential_prior(rate = 2),
#'   seed = 11, verbose = FALSE
#' )
#' mean(g$edge_indicators) # shifted away from 0.3 by the Z(Gamma) tilt
#' }
#' @seealso \code{\link{sample_ggm_prior}}, \code{\link{sample_sbm_prior}},
#'   \code{\link{bernoulli_prior}}, \code{\link{beta_bernoulli_prior}},
#'   \code{\link{sbm_prior}}, \code{\link{bgm}}
#'
#' @export
sample_graph_prior = function(
  p,
  n_samples,
  edge_prior = bernoulli_prior(0.5),
  spec = c("hierarchical", "joint"),
  interaction_prior = cauchy_prior(scale = 2.5),
  precision_scale_prior = exponential_prior(eta = 1),
  delta = NULL,
  theta = NULL,
  allocations = NULL,
  block_probs = NULL,
  n_warmup = 2e3,
  seed = 1L,
  verbose = TRUE
) {
  spec = match.arg(spec)
  if(!is.numeric(p) || length(p) != 1L || is.na(p) || p < 2 || p != round(p)) {
    stop("'p' must be a single integer >= 2.")
  }
  p = as.integer(p)
  if(!is.numeric(n_samples) || length(n_samples) != 1L || is.na(n_samples) ||
    n_samples < 1 || n_samples != round(n_samples)) {
    stop("'n_samples' must be a single positive integer.")
  }
  n_samples = as.integer(n_samples)

  ep = unpack_indicator_prior(edge_prior, num_variables = p)
  cond_prob = graph_prior_conditioning(ep, p, theta, allocations, block_probs)

  pairs = graph_pair_indices(p)
  pair_names = paste0(pairs[, 1L], "-", pairs[, 2L])

  if(spec == "joint") {
    eff_edge_prior = if(is.null(cond_prob)) {
      edge_prior
    } else {
      bernoulli_prior(inclusion_probability = cond_prob)
    }
    draws = sample_ggm_prior(
      p = p, n_samples = n_samples, n_warmup = as.integer(n_warmup),
      interaction_prior = interaction_prior,
      precision_scale_prior = precision_scale_prior,
      delta = delta,
      spec = "joint",
      edge_prior = eff_edge_prior,
      update_method = "gibbs",
      seed = as.integer(seed),
      verbose = verbose
    )
    out = list(
      edge_indicators = draws$edge_indicators,
      pair_names = pair_names,
      theta = draws$theta,
      allocations = draws$allocations,
      spec = spec,
      edge_prior = ep$edge_prior,
      p = p
    )
    return(out)
  }

  # Hierarchical spec: ancestral. Hyperparameters from their prior, then
  # independent pair flips given the implied pair probabilities.
  res = with_graph_prior_seed(seed, function() {
    n_edges = nrow(pairs)
    theta_out = NULL
    alloc_out = NULL

    if(!is.null(cond_prob) || identical(ep$edge_prior, "Bernoulli")) {
      prob = if(is.null(cond_prob)) ep$inclusion_probability else cond_prob
      pv = prob[pairs]
      ru = matrix(runif(n_samples * n_edges), n_samples, n_edges)
      indicators = 1L * (ru < matrix(pv, n_samples, n_edges, byrow = TRUE))
    } else if(identical(ep$edge_prior, "Beta-Bernoulli")) {
      theta_out = rbeta(
        n_samples, ep$beta_bernoulli_alpha, ep$beta_bernoulli_beta
      )
      ru = matrix(runif(n_samples * n_edges), n_samples, n_edges)
      indicators = 1L * (ru < theta_out)
    } else {
      alloc_out = matrix(0L, n_samples, p)
      indicators = matrix(0L, n_samples, n_edges)
      for(s in seq_len(n_samples)) {
        z = ancestral_mfm_sbm_partition(p, ep$lambda, ep$dirichlet_alpha)
        prob = ancestral_sbm_pair_probabilities(
          z,
          ep$beta_bernoulli_alpha, ep$beta_bernoulli_beta,
          ep$beta_bernoulli_alpha_between, ep$beta_bernoulli_beta_between
        )
        alloc_out[s, ] = as.integer(z)
        indicators[s, ] = as.integer(runif(n_edges) < prob[pairs])
      }
    }
    storage.mode(indicators) = "integer"
    list(indicators = indicators, theta = theta_out, allocations = alloc_out)
  })

  list(
    edge_indicators = res$indicators,
    pair_names = pair_names,
    theta = res$theta,
    allocations = res$allocations,
    spec = spec,
    edge_prior = ep$edge_prior,
    p = p
  )
}


#' @title Sample from the Stochastic-Block Edge-Prior Hyperprior
#'
#' @description
#' Ancestral draws from the MFM-SBM hyperprior of
#' \code{\link{sbm_prior}()}: a partition of the nodes into blocks
#' (shifted-Poisson number of components, Dirichlet-weighted allocation) and
#' Beta-distributed within- and between-block edge-inclusion probabilities.
#' These are the hyperparameters that \code{\link{sample_graph_prior}} and
#' \code{\link{bgm}} integrate over when the edge prior is a Stochastic-Block
#' prior.
#'
#' @param p Integer. Number of nodes (\eqn{p \ge 2}).
#' @param n_samples Integer. Number of prior draws.
#' @param edge_prior An \code{\link{sbm_prior}()} object. Default
#'   \code{sbm_prior()}.
#' @param seed Integer. Seed for the draw; the caller's RNG state is
#'   restored on exit.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{allocations}}{Integer matrix (\code{n_samples x p}) of
#'       block labels.}
#'     \item{\code{pair_probability}}{Numeric matrix
#'       (\code{n_samples x p(p-1)/2}) of implied pair-inclusion
#'       probabilities, columns in row-major upper-triangle order.}
#'     \item{\code{num_blocks}}{Integer vector: number of occupied blocks
#'       per draw.}
#'     \item{\code{pair_names}}{Character vector labeling the pair columns
#'       as \code{"i-j"}.}
#'     \item{\code{p}}{The node count.}
#'   }
#'
#' @examples
#' draws = sample_sbm_prior(p = 8, n_samples = 100, seed = 4)
#' table(draws$num_blocks)
#' range(draws$pair_probability)
#' @seealso \code{\link{sbm_prior}}, \code{\link{sample_graph_prior}},
#'   \code{\link{bgm}}
#'
#' @export
sample_sbm_prior = function(p, n_samples, edge_prior = sbm_prior(),
                            seed = 1L) {
  if(!is.numeric(p) || length(p) != 1L || is.na(p) || p < 2 || p != round(p)) {
    stop("'p' must be a single integer >= 2.")
  }
  p = as.integer(p)
  if(!is.numeric(n_samples) || length(n_samples) != 1L || is.na(n_samples) ||
    n_samples < 1 || n_samples != round(n_samples)) {
    stop("'n_samples' must be a single positive integer.")
  }
  n_samples = as.integer(n_samples)

  ep = unpack_indicator_prior(edge_prior, num_variables = p)
  if(!identical(ep$edge_prior, "Stochastic-Block")) {
    stop("'edge_prior' must be an sbm_prior() object.")
  }

  pairs = graph_pair_indices(p)
  pair_names = paste0(pairs[, 1L], "-", pairs[, 2L])

  res = with_graph_prior_seed(seed, function() {
    allocations = matrix(0L, n_samples, p)
    pair_probability = matrix(0, n_samples, nrow(pairs))
    num_blocks = integer(n_samples)
    for(s in seq_len(n_samples)) {
      z = ancestral_mfm_sbm_partition(p, ep$lambda, ep$dirichlet_alpha)
      prob = ancestral_sbm_pair_probabilities(
        z,
        ep$beta_bernoulli_alpha, ep$beta_bernoulli_beta,
        ep$beta_bernoulli_alpha_between, ep$beta_bernoulli_beta_between
      )
      allocations[s, ] = as.integer(z)
      pair_probability[s, ] = prob[pairs]
      num_blocks[s] = length(unique(z))
    }
    list(
      allocations = allocations,
      pair_probability = pair_probability,
      num_blocks = num_blocks
    )
  })

  list(
    allocations = res$allocations,
    pair_probability = res$pair_probability,
    num_blocks = res$num_blocks,
    pair_names = pair_names,
    p = p
  )
}
