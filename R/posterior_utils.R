# Dahl's method to summarize the samples of the cluster_allocations
#  This function was adapted from the R code accompanying the paper:
#  Geng, J., Bhattacharya, A., & Pati, D. (2019). Probabilistic Community
#  Detection With Unknown Number of Communities, Journal of the American
#  Statistical Association, 114:526, 893-905, DOI:10.1080/01621459.2018.1458618
getDahl = function(cluster_allocations) {
  # Dimensions of the input matrix
  niters = nrow(cluster_allocations)  # Number of iterations
  n = ncol(cluster_allocations)       # Number of nodes

  # Compute membership matrices for each iteration
  membershipMatrices = apply(cluster_allocations, 1, function(clusterAssign) {
    outer(clusterAssign, clusterAssign, FUN = "==")
  })

  # Reshape membershipMatrices into a list of matrices
  membershipMatrices = lapply(seq_len(niters), function(i) {
    matrix(membershipMatrices[, i], n, n)
  })

  # Compute the average membership matrix
  membershipAverage = Reduce("+", membershipMatrices) / niters

  # Compute squared error for each iteration
  SqError = sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                   av = membershipAverage)

  # Find the iteration with the minimum squared error
  DahlIndex = which.min(SqError)

  # Extract the cluster assignment corresponding to the best iteration
  DahlAns = cluster_allocations[DahlIndex, , drop = TRUE]

  return(DahlAns)
}


# Calculate the conditional probability of the number of components given the
# cardinality of a sampled allocation vector based on Equation (3.7) from
# Miller & Harrison (2018). Mixture Models With a Prior on the Number of
# Components, Journal of the American Statistical Association, 113:521, 340-356,
# DOI:10.1080/01621459.2016.1255636
#' @importFrom stats dpois
compute_p_k_given_t = function(
    t,
    log_Vn,
    dirichlet_alpha,
    num_variables,
    lambda) {
  # Define the K_values
  K_values = as.numeric(1:num_variables)

  # Initialize vector for probabilities
  p_k_given_t = numeric(length(K_values))

  # Normalization constant for t
  log_vn_t = log_Vn[t]

  # Normalizing factor for the truncated Poisson distribution
  norm_factor = 1 - dpois(0, lambda)
  truncated_poisson_pmf = dpois(K_values, lambda) / norm_factor

  # Loop through each value of K
  for (i in seq_along(K_values)) {
    K = K_values[i]
    if (K >= t) {
      # Falling factorial
      falling_factorial = prod(K:(K - t + 1))
      # Rising factorial
      rising_factorial = prod((dirichlet_alpha * K) + 0:(num_variables - 1))
      # Compute log probability
      log_p_k = log(falling_factorial) - log(rising_factorial) +
        log(truncated_poisson_pmf[i]) - log_vn_t
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

#' Function for summarizing the sampled cluster allocation vectors
#'
#' Th \code{summarySBM} function summarizes the sampled allocation vectors from
#' each iteration of the Gibbs sampler from the output of the  \code{bgm}
#' function ran with \code{edge_prior = "Stochastic-Block"} and
#' \code{save = TRUE}. It also estimates the posterior distribution of the
#' number of clusters.
#'
#' @param bgm_object A fit object created by the bgm function.
#' @param internal_call A logical value indicating whether the function is used
#' within bgms for calculating the posterior probabilities of the number of
#' clusters or by the user. This argument is always set to FALSE.
#' @return Returns a list of two elements: \code{components} and \code{allocations},
#' containing the posterior probabilities for the number of components (clusters)
#' and the estimated cluster allocation of the nodes using Dahl's method.
#' @examples
#' \donttest{
#'   # fit a model with the SBM prior
#'   bgm_object = bgm(
#'     Wenchuan[, c(1:5)],
#'     edge_prior = "Stochastic-Block",
#'     save = TRUE)
#'
#'   summarySBM(bgm_object)
#' }
#' @export
summarySBM = function(
    bgm_object,
    internal_call = FALSE) {

  arguments = extract_arguments(bgm_object)

  if(arguments$edge_prior != "Stochastic-Block")
    stop('The bgm function must be run with edge_prior = "Stochastic-Block".')

  if(arguments$save == FALSE && internal_call == FALSE)
    stop('The bgm function must be run with save = TRUE.')

  cluster_allocations = bgm_object$allocations
  dirichlet_alpha = arguments$dirichlet_alpha
  lambda = arguments$lambda

  # Pre-compute log_Vn for computing the cluster probabilities
  num_variables = ncol(cluster_allocations)
  log_Vn = compute_Vn_mfm_sbm(
    num_variables, dirichlet_alpha, num_variables + 10, lambda)

  # Compute the number of unique clusters (t) for each iteration, i.e., the
  # cardinality  of the partition z
  clusters = apply(cluster_allocations, 1, function(row) length(unique(row)))

  # Compute the conditional probabilities of the number of clusters for each
  # row in clusters
  p_k_given_t = matrix(NA, nrow = length(clusters), ncol = num_variables)

  for (i in 1:length(clusters)) {
    p_k_given_t[i, ] = compute_p_k_given_t(
      clusters[i], log_Vn, dirichlet_alpha, num_variables, lambda)
  }

  # Average across all iterations
  p_k_given_t = colMeans(p_k_given_t)

  # Format the output
  num_components = 1:num_variables
  components = cbind(num_components, p_k_given_t)
  colnames(components) = c("num_components", "probability")

  # Compute the allocations of the nodes based on Dahl's method
  allocations = getDahl(cluster_allocations)

  return(list(components = components,
              allocations = allocations))
}
