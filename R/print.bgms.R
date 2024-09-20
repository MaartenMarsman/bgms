#' @name print.bgms
#' @title  Print method for \code{bgms} objects
#'
#' @description Used to prevent bgms output cluttering the console.
#'
#' @param x An object of class \code{bgms}.
#' @param ... Ignored.
#'
#' @export
print.bgms <- function(x, ...) {
  arguments = extract_arguments(x)
  if(arguments$edge_selection) {
    if(arguments$edge_prior == "Bernoulli") {
      cat("Bayesian Edge Selection using a Bernoulli prior on edge inclusion\n")
    } else {
      cat("Bayesian Edge Selection using a Beta-Bernoulli prior on edge inclusion\n")
    }

  } else {
    cat("Bayesian Estimation\n")
  }

  cat(paste0(" Number of variables: ", arguments$no_variables, "\n"))
  if(arguments$na_impute) {
    cat(paste0(" Number of cases: ", arguments$no_cases, " (missings imputed)\n"))
  } else {
    cat(paste0(" Number of cases: ", arguments$no_cases, "\n"))
  }
  if(arguments$save) {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, " (MCMC output saved)\n"))
  } else {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, " (posterior means saved)\n"))
  }
  cat("See the easybgm package for extensive summary and plotting functions")
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

  cat(paste0(" Number of variables: ", arguments$no_variables, "\n"))
  if(arguments$na_impute) {
    cat(paste0(" Number of cases Group 1: ", arguments$no_cases_gr1, " (missings imputed)\n"))
    cat(paste0(" Number of cases Group 2: ", arguments$no_cases_gr2, " (missings imputed)\n"))
  } else {
    cat(paste0(" Number of cases Group 1: ", arguments$no_cases_gr1, "\n"))
    cat(paste0(" Number of cases Group 2: ", arguments$no_cases_gr2, "\n"))
  }
  if(arguments$save) {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, " (MCMC output saved)\n"))
  } else {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, " (posterior means saved)\n"))
  }
  cat("See the easybgm package for extensive summary and plotting functions")
}