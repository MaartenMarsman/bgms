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


# A function that computes the cluster posterior probabilities and allocations

summary_SBM = function(cluster_allocations) {
  # Compute the number of unique clusters for each iteration
  clusters = apply(cluster_allocations, 1, function(row) length(unique(row)))

  # Compute the posterior probabilities of the actual unique clusters
  no_clusters = table(clusters) / length(clusters)

  # Compute the allocations of the nodes based on Dahl's method
  allocations = getDahl(cluster_allocations)

  # Return the results
  return(list(no_clusters = no_clusters,
              allocations = allocations))
}
