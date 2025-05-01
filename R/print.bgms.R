#' @name print.bgms
#' @title  Print method for \code{bgms} objects
#'
#' @description Used to prevent bgms output cluttering the console.
#'
#' @param x An object of class \code{bgms}.
#' @param ... Ignored.
#'
#' @importFrom Rdpack reprompt
#' @export
print.bgms <- function(x, ...) {
  arguments = extract_arguments(x)
  if(arguments$edge_selection) {
    if(arguments$edge_prior == "Bernoulli") {
      cat("Bayesian Edge Selection using a Bernoulli prior on edge inclusion\n")
    } else if(arguments$edge_prior == "Beta-Bernoulli") {
      cat("Bayesian Edge Selection using a Beta-Bernoulli prior on edge inclusion\n")
    } else if(arguments$edge_prior == "Stochastic-Block") {
      cat("Bayesian Edge Selection using a Stochastic Block prior on edge inclusion\n")
    }
  } else {
    cat("Bayesian Estimation\n")
  }

  cat(paste0(" Number of variables: ", arguments$num_variables, "\n"))
  if(arguments$na_impute) {
    cat(paste0(" Number of cases: ", arguments$num_cases, " (missings imputed)\n"))
  } else {
    cat(paste0(" Number of cases: ", arguments$num_cases, "\n"))
  }
  if(arguments$save) {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, " (MCMC output saved)\n"))
  } else {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, " (posterior means saved)\n"))
  }
  cat("See the easybgm package for extensive summary and plotting functions\n")
}

#' @name print.bgmCompare
#' @title  Print method for \code{bgms} objects
#'
#' @description Used to prevent bgms output cluttering the console.
#'
#' @param x An object of class \code{bgms}.
#' @param ... Ignored.
#'
#' @export
print.bgmCompare <- function(x, ...) {
  arguments = extract_arguments(x)
  if(arguments$difference_selection) {
    if(arguments$pairwise_difference_prior == "Bernoulli") {
      cat(paste0("Bayesian Variable Selection using a Bernoulli prior on the inclusion of \n",
                 "differences in pairwise interactions\n"))
    } else {
      cat(paste0("Bayesian Variable Selection using a Beta-Bernoulli prior on the inclusion of \n",
                 "differences in pairwise interactions\n"))
    }
    if(arguments$main_difference_model == "Free") {
      cat("Group specific category threshold parameters were estimated")
    } else {
      if(arguments$main_difference_prior == "Bernoulli") {
        cat(paste0("Bayesian Variable Selection using a Bernoulli prior on the inclusion of\n",
                   "differences in the category thresholds\n"))
      } else {
        cat(paste0("Bayesian Variable Selection using a Beta-Bernoulli prior on the inclusion of\n",
                   "differences in the category thresholds\n"))      }
    }
  } else {
    cat("Bayesian Estimation\n")
  }

  cat(paste0(" Number of variables: ", arguments$num_variables, "\n"))
  num_groups = length(unique(arguments$group))
  if(arguments$na_impute) {
    for(group in 1:num_groups) {
      cat(paste0(" Number of cases Group ", group,": ", arguments$num_cases[group], " (missings imputed)\n"))
    }
  } else {
    for(group in 1:num_groups) {
      cat(paste0(" Number of cases Group ", group,": ", arguments$num_cases[group],"\n"))
    }
  }
  if(arguments$save) {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, " (MCMC output saved)\n"))
  } else {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, " (posterior means saved)\n"))
  }
  cat("See the easybgm package for extensive summary and plotting functions \n")
}