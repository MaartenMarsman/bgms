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
  batch_rhat = .compute_rhat_cpp(array3d)

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

  # Build co-clustering (membership) matrices for all iterations

  Ms = lapply(seq_len(n_iter), function(t) {
    z = cluster_matrix[t, ]
    (outer(z, z, FUN = "==")) * 1L
  })

  # Average (posterior similarity / co-clustering) matrix
  psm = Reduce(`+`, Ms) / n_iter

  # MEAN representative (Dahl's method)
  sqerr = vapply(Ms, function(M) sum((M - psm)^2), numeric(1))
  idx_dahl = which.min(sqerr)
  alloc_dahl = cluster_matrix[idx_dahl, , drop = TRUE]

  #  MODE representative
  hash_mat = function(M) paste(as.integer(t(M)), collapse = ",")
  keys = vapply(Ms, hash_mat, character(1))
  tab = table(keys)
  key_mode = names(tab)[which.max(tab)]
  idx_mode = match(key_mode, keys)
  alloc_mode = cluster_matrix[idx_mode, , drop = TRUE]
  indicator_mode = matrix(
    as.integer(strsplit(key_mode, ",", fixed = TRUE)[[1]]),
    nrow = p, byrow = TRUE
  )
  p_dist = as.numeric(tab) / sum(tab)
  posterior_variance = (1 - sum(p_dist^2)) / (1 - 1 / length(p_dist))

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

  # Normalization constant for t
  log_vn_t = log_Vn[t]

  # Shifted Poisson prior on the number of components
  poisson_pmf = dpois(K_values - 1, lambda)

  # Loop through each value of K
  for(i in seq_along(K_values)) {
    K = K_values[i]
    if(K >= t) {
      # Falling factorial
      falling_factorial = prod(K:(K - t + 1))
      # Rising factorial
      rising_factorial = prod((dirichlet_alpha * K) + 0:(num_variables - 1))
      # Compute log probability
      log_p_k = log(falling_factorial) - log(rising_factorial) +
        log(poisson_pmf[i]) - log_vn_t
      # Convert log probability to probability
      p_k_given_t[i] = exp(log_p_k)
    } else {
      p_k_given_t[i] = 0
    }
  }
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

  # Compute the conditional probabilities of the number of clusters for each
  # row in clusters
  p_k_given_t = matrix(NA, nrow = length(clusters), ncol = num_variables)

  for(i in seq_along(clusters)) {
    p_k_given_t[i, ] = compute_p_k_given_t(
      clusters[i], log_Vn, dirichlet_alpha, num_variables, lambda
    )
  }

  # Average across all iterations
  p_k_given_t = colMeans(p_k_given_t)

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
