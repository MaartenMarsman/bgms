# Stochastic Block Model (SBM) posterior-summary utilities
#
# Split out of mcmc_summary.R (cleanup S4). These helpers summarize the
# cluster-allocation output of the SBM edge/difference prior: pairwise
# co-appearance diagnostics, a representative clustering (Dahl mean + mode),
# the conditional distribution of the number of blocks, and the top-level
# posterior_summary_SBM wrapper.


# summarize the SBM output -----------------------------------------------------

# Calculate convergence diagnostics on the pairwise cluster co-appearance values
summarize_alloc_pairs = function(allocations, node_names = NULL) {
  n_ch = length(allocations)
  n_iter = nrow(allocations[[1]])
  no_variables = ncol(allocations[[1]])
  for(c in seq_len(n_ch)) {
    stopifnot(nrow(allocations[[c]]) == n_iter, ncol(allocations[[c]]) == no_variables)
  }
  if(!is.null(node_names)) stopifnot(length(node_names) == no_variables)

  # all node pairs
  Pairs = t(combn(seq_len(no_variables), 2))
  nparam = nrow(Pairs)

  # helper to construct a "time-series"
  get_draws_pair = function(i, j) {
    out = matrix(NA, n_iter, n_ch)
    for(c in seq_len(n_ch)) {
      Zc = allocations[[c]]
      out[, c] = as.integer(Zc[, i] == Zc[, j])
    }
    out
  }

  # Pre-build 3D array and batch Rhat via C++
  array3d = array(NA_real_, dim = c(n_iter, n_ch, nparam))
  for(p in seq_len(nparam)) {
    array3d[, , p] = get_draws_pair(Pairs[p, 1], Pairs[p, 2])
  }
  ind_stats = .compute_indicator_ess_cpp(array3d)
  batch_rhat = .compute_rhat_cpp(split_chains(array3d))

  result = cbind(
    ind_stats[, c("mean", "mcse", "sd", "n00", "n01", "n10", "n11", "n_eff_mixt"), drop = FALSE],
    Rhat = batch_rhat
  )
  colnames(result)[4:7] = c("n0->0", "n0->1", "n1->0", "n1->1")
  result[is.na(result[, "n_eff_mixt"]), "Rhat"] = NA_real_
  if(is.null(node_names)) {
    rn = paste0(Pairs[, 1], "-", Pairs[, 2])
    dimn = as.character(seq_len(no_variables))
  } else {
    rn = paste0(node_names[Pairs[, 1]], "-", node_names[Pairs[, 2]])
    dimn = node_names
  }

  sbm_summary = as.data.frame(result, check.names = FALSE)
  rownames(sbm_summary) = rn

  # construct the co-appearance matrix
  co_occur_matrix = matrix(0,
    nrow = no_variables, ncol = no_variables,
    dimnames = list(dimn, dimn)
  )
  diag(co_occur_matrix) = 1
  for(p in seq_len(nparam)) {
    i = Pairs[p, 1]
    j = Pairs[p, 2]
    m = sbm_summary[p, "mean"]
    co_occur_matrix[i, j] = m
    co_occur_matrix[j, i] = m
  }
  list(sbm_summary = sbm_summary, co_occur_matrix = co_occur_matrix)
}

# calculate a representative allocation vector using
# (1) the mean, based on Dahl's method: This part of the code
# was adapted from the R
# code accompanying the paper:
#  Geng, J., Bhattacharya, A., & Pati, D. (2019). Probabilistic Community
#  Detection With Unknown Number of Communities, Journal of the American
#  Statistical Association, 114:526, 893-905, DOI:10.1080/01621459.2018.1458618
# (2) the mode (most frequent co-clustering).
find_representative_clustering = function(cluster_matrix) {
  stopifnot(is.matrix(cluster_matrix))
  n_iter = nrow(cluster_matrix)
  p = ncol(cluster_matrix)

  # Posterior similarity (co-clustering) matrix, streamed one iteration at a
  # time; the per-iteration membership matrices are never stored.
  psm = matrix(0, p, p)
  for(t in seq_len(n_iter)) {
    z = cluster_matrix[t, ]
    psm = psm + (outer(z, z, FUN = "==") * 1L)
  }
  psm = psm / n_iter

  # MEAN representative (Dahl's method): ||M_t - psm||^2 expands to
  # n_pairs_t - 2 * sum(psm over co-clustered pairs) + sum(psm^2), where the
  # co-clustered pairs decompose block-by-block, so M_t is never materialized.
  sum_psm2 = sum(psm^2)
  sqerr = vapply(seq_len(n_iter), function(t) {
    z = cluster_matrix[t, ]
    s = 0
    n_pairs = 0
    for(lab in unique(z)) {
      idx = which(z == lab)
      s = s + sum(psm[idx, idx])
      n_pairs = n_pairs + length(idx)^2
    }
    n_pairs - 2 * s + sum_psm2
  }, numeric(1))
  idx_dahl = which.min(sqerr)
  alloc_dahl = cluster_matrix[idx_dahl, , drop = TRUE]

  # MODE representative: two allocation vectors induce the same partition
  # exactly when their first-occurrence canonical relabelings match, so an
  # O(p) canonical label vector replaces the p x p membership matrix as the
  # hash key.
  keys = vapply(seq_len(n_iter), function(t) {
    z = cluster_matrix[t, ]
    paste(match(z, unique(z)), collapse = ",")
  }, character(1))
  tab = table(keys)
  key_mode = names(tab)[which.max(tab)]
  idx_mode = match(key_mode, keys)
  alloc_mode = cluster_matrix[idx_mode, , drop = TRUE]

  list(
    mean = alloc_dahl,
    mode = alloc_mode
  )
}

# Calculate the conditional probability of the number of blocks given the
# cardinality of a sampled allocation vector based on Equation (3.7) from
# Miller & Harrison (2018). Mixture Models With a Prior on the Number of
# blocks, Journal of the American Statistical Association, 113:521, 340-356,
# DOI:10.1080/01621459.2016.1255636
#
# The prior on the number of components is the shifted Poisson
# K - 1 ~ Poisson(lambda), matching the partition coefficients from
# compute_Vn_mfm_sbm() that drive the sampler.
#' @importFrom stats dpois
compute_p_k_given_t = function(
  t,
  log_Vn,
  dirichlet_alpha,
  num_variables,
  lambda
) {
  # Define the K_values
  K_values = as.numeric(1:num_variables)

  # Initialize vector for probabilities
  p_k_given_t = numeric(length(K_values))

  # Shifted Poisson prior on the number of components (log scale)
  log_poisson_pmf = dpois(K_values - 1, lambda, log = TRUE)

  # Falling factorial K!/(K-t)! and rising factorial
  # prod(alpha*K + 0:(n-1)) = gamma(alpha*K + n)/gamma(alpha*K), both on the
  # log scale so large K/t cannot overflow to Inf.
  valid = K_values >= t
  K = K_values[valid]
  log_falling_factorial = lgamma(K + 1) - lgamma(K - t + 1)
  log_rising_factorial = lgamma(dirichlet_alpha * K + num_variables) -
    lgamma(dirichlet_alpha * K)
  log_p_k = log_falling_factorial - log_rising_factorial +
    log_poisson_pmf[valid] - log_Vn[t]
  p_k_given_t[valid] = exp(log_p_k)

  # Normalize probabilities
  p_k_given_t = p_k_given_t / sum(p_k_given_t)

  return(p_k_given_t)
}

# Wrapper function to compute the posterior summary for the Stochastic Block Model
posterior_summary_SBM = function(
  allocations,
  arguments
) {
  # combine the allocations from the chains
  cluster_allocations = do.call(rbind, allocations)

  dirichlet_alpha = arguments$dirichlet_alpha
  lambda = arguments$lambda
  num_variables = ncol(cluster_allocations)

  # Pre-compute log_Vn for computing the cluster probabilities
  log_Vn = compute_Vn_mfm_sbm(
    num_variables, dirichlet_alpha, num_variables + 10, lambda
  )

  # Compute the number of unique clusters (t) for each iteration, i.e., the
  # cardinality  of the partition z
  clusters = apply(cluster_allocations, 1, function(row) length(unique(row)))

  # Compute the conditional probabilities of the number of clusters once per
  # unique cardinality, then average with the observed frequencies (the
  # per-iteration values only depend on the cardinality t).
  unique_t = sort(unique(clusters))
  p_k_by_t = vapply(unique_t, function(t) {
    compute_p_k_given_t(t, log_Vn, dirichlet_alpha, num_variables, lambda)
  }, numeric(num_variables))
  t_freq = tabulate(match(clusters, unique_t), nbins = length(unique_t))
  p_k_given_t = as.numeric(p_k_by_t %*% (t_freq / length(clusters)))

  # Format the output
  # num_blocks = 1:num_variables
  blocks = cbind(p_k_given_t)
  colnames(blocks) = c("probability")

  # make blocks a data frame
  blocks = as.data.frame(blocks)

  # Compute the mean and mode of the allocations
  allocations = find_representative_clustering(cluster_allocations)

  return(list(
    blocks = blocks,
    allocations_mean = allocations$mean,
    allocations_mode = allocations$mode
  ))
}
