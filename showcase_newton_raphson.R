library(bgms)

pseudolikelihood_numerator <- function(thresholds, interactions, suffstats, seen, threshold_counts_without_0, P) {
  result <- 0.0

  for (i in seq_len(P)) {
    for (u in seq_len(seen[i] - 1)) {
      result <- result + threshold_counts_without_0[i, u] * thresholds[i, u]
    }
  }

  result <- result + sum(interactions * suffstats)

  return(result)
}


pseudo_logposterior_full_aij2 <- function(a, i, j, thresholds, interactions, suffstats, seen, threshold_counts_without_0, X, P, N,
                                          prior_cauchy_scale = 2.5) {

  interactions[i, j] <- interactions[j, i] <- a

  pseudolikelihood_numerator(
    thresholds                 = thresholds,
    interactions               = interactions,
    suffstats                  = suffstats,
    seen                       = seen,
    threshold_counts_without_0 = threshold_counts_without_0,
    P                          = P
  ) +
    pseudolikelihood_denominator2(
      thresholds                 = thresholds,
      interactions               = interactions,
      suffstats                  = suffstats,
      seen                       = seen,
      X                          = X,
      P                          = P,
      N                          = N,
      i0                         = i,
      j0                         = j
    ) +
    sum(dcauchy(interactions[lower.tri(interactions)], 0, prior_cauchy_scale, log = TRUE))
}

pseudolikelihood_denominator2 <- function(thresholds, interactions, suffstats, seen, X, P, N, i0, j0) {
  result <- 0.0

  for (v in seq_len(N)) {
    for (i in c(i0, j0)) {
      temp1 <- c(crossprod(interactions[i, ], X[v, ]))

      temp2 <- 1.0
      for (u in seq_len(seen[i] - 1)) {
        temp2 <- temp2 + exp(thresholds[i, u] + u * temp1)
      }

      result <- result - log(temp2)
    }
  }

  return(result)
}

log_pseudolikelihood_full2 <- function(a, i, j, Mu, Sigma, iter, x, suffstats, seen, threshold_counts_without_0) {


  n <- nrow(x)  # Number of observations
  p <- ncol(x)

  MuIter <- Mu[iter, ]
  MuMat <- matrix(0, p, ncol(threshold_counts_without_0))  # Initialize matrix for thresholds
  idx <- 1
  for (ii in 1:p) {
    for (jj in 1:ncol(MuMat)) {
      MuMat[ii, jj] <- MuIter[idx]
      idx <- idx + 1# Fill matrix with threshold values
    }
  }
  SigmaIter <- Sigma[iter, ]
  SigmaMat = matrix(0, p, p)  # Initialize matrix for interactions
  SigmaMat[lower.tri(SigmaMat)] = SigmaIter  # Fill lower triangle with Sigma values
  SigmaMat = SigmaMat + t(SigmaMat)  # Make symmetric

  D = length(a)  # Number of elements in a
  log_pl = numeric(length = D)  # Initialize log pseudolikelihood vector

  # colMax <- unname(matrixStats::colMaxs(x))
  # log_p <- numeric(length = max(colMax) + 1)  # Initialize log probability vector

  for (d in 1:D) {
    log_pl[d] = pseudo_logposterior_full_aij2(a[d], i, j, thresholds = MuMat, interactions = SigmaMat, X = x, N = n, P = p,
                                              seen = seen, suffstats = suffstats, threshold_counts_without_0 = threshold_counts_without_0)
  }

  return(log_pl)  # Return log pseudolikelihood
}


x0 = Wenchuan[1:50, 1:5]  # Select the first 5 columns of Wenchuan dataset
p = ncol(x0)  # Get the number of variables (columns)

samples = bgm(x0, save = TRUE)  # Run the bgm function and save samples
Mu    = samples$main_effect_samples  # Extract threshold estimates
Sigma = samples$pairwise_effect_samples  # Extract interaction estimates

data = bgms:::reformat_data(x = x0,
                            na_action = "listwise",
                            variable_bool = rep(TRUE, p),
                            reference_category = rep(1, p))

x = data$x  # Extract reformatted data
no_categories = data$no_categories  # Get number of categories per variable
no_categories = cumsum(no_categories)  # Cumulative sum for indexing
start = 1 + c(0, no_categories[-length(no_categories)])  # Start indices
stop = no_categories  # Stop indices

K <- max(x)
threshold_counts_wench <- apply(x, 2, \(y) c(table(c(y, 0:K)) - 1))
threshold_counts_without_0_wench <- apply(threshold_counts_wench, 2L, \(y) {
  c(y[y > 0], rep(0, sum(y == 0)))[-1L]
})
threshold_counts_without_0_wench <- t(matrix(threshold_counts_without_0_wench, K, p))

seen_wench <- unname(apply(x, 2, \(y) length(unique(y))))
suffstats_wench <- unname(crossprod(x))


i <- 2; j <- 1
log_pseudolikelihood_full2(c(.2, .5), i, j, Mu, Sigma, iter = 10000, x = x,
                           seen = seen_wench, suffstats = suffstats_wench,
                           threshold_counts_without_0 = threshold_counts_without_0_wench)

optim_res <- optim(Sigma[i, j], function(a) {
  returnVal <- log_pseudolikelihood_full2(a, i, j, Mu, Sigma, iter = 10000, x = x,
                             seen = seen_wench, suffstats = suffstats_wench,
                             threshold_counts_without_0 = threshold_counts_without_0_wench)
  if (any(!is.finite(returnVal))) {
    for (i in seq_along(returnVal)) {
      # if (!is.finite(returnVal[i])) {
      #   print(sprintf("a: %f, returnVal: %f", a[i], returnVal[i]))
      # }
      if (!is.finite(returnVal[i]) && returnVal[i] < 0) {
        returnVal[i] <- -.Machine$double.xmax
      }
    }
  }
  return(returnVal)
}, method = "Brent", lower = -100, upper = 100, control = list(fnscale = -1, trace = 5))


# setup arguments for C++
iter <- 10000
MuMat <- matrix(0, p, ncol(threshold_counts_without_0_wench))  # Initialize matrix for thresholds
idx <- 1
for (ii in 1:p) {
  for (jj in 1:ncol(MuMat)) {
    MuMat[ii, jj] <- Mu[iter, idx]
    idx <- idx + 1# Fill matrix with threshold values
  }
}
SigmaMat = matrix(0, p, p)  # Initialize matrix for interactions
SigmaMat[lower.tri(SigmaMat)] = Sigma[iter, ]  # Fill lower triangle with Sigma values
SigmaMat = SigmaMat + t(SigmaMat)  # Make symmetric


pairwise_effects <- SigmaMat
main_effects     <- MuMat

# const double
initial_value <- Sigma[i, j]
# const arma::mat&
pairwise_effects <- SigmaMat
# const arma::mat&
main_effects <- MuMat
# const arma::imat&
observations <- x
# const arma::ivec&
num_categories <- seen_wench
# const int
num_persons <- nrow(x)
# const int
variable1 <- i
# const int
variable2 <- j
# TODO: these two are unused?
# const double
proposed_state <- 0.0
# const double
current_state <- 0.0
# const arma::mat&
residual_matrix <- matrix(0, nrow(x), p)
# const arma::uvec&
is_ordinal_variable <- rep(1, p)
# const arma::ivec&
reference_category <- data$reference_category
# const double
interaction_scale <- 2.5

newton_raphson_x <- bgms:::optimize_log_pseudoposterior_interaction(
  initial_value       = c(initial_value),
  pairwise_effects    = pairwise_effects,
  main_effects        = main_effects,
  observations        = observations,
  num_categories      = num_categories - 1,
  num_persons         = num_persons,
  variable1           = variable1 - 1,
  variable2           = variable2 - 1,
  proposed_state      = proposed_state,
  current_state       = current_state,
  residual_matrix     = residual_matrix,
  is_ordinal_variable = is_ordinal_variable,
  reference_category  = reference_category,
  interaction_scale   = interaction_scale
)
newton_raphson_fx <- log_pseudolikelihood_full2(newton_raphson, i, j, Mu, Sigma, iter = 10000, x = x,
                                                seen = seen_wench, suffstats = suffstats_wench,
                                                threshold_counts_without_0 = threshold_counts_without_0_wench)
matrix(c(newton_raphson_x, newton_raphson_fx, optim_res$par, optim_res$value),
       nrow = 2, dimnames = list(c("x", "f(x)"), c("Newton-Raphson", "Optim")))
