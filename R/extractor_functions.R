extract_arguments <- function(fit_object) {
  if(class(fit_object) != "bgms")
    stop(paste0("Expected an object with class bgms and not one with class ",
      class(fit_object)))

  if(is.null(fit_object$bgm_arguments)) {
    stop(paste0("Extractor functions have been defined for bgms versions 0.1.3 and up but not \n",
                "for older versions. The current fit object predates version 0.1.3."))
  } else {
    return(fit_object$bgm_arguments)
  }
}

extract_edge_indicators <- function(fit_object) {
  arguments = extract_arguments(fit_object)
  if(arguments$save) {
    edge_indicators = fit_object$gamma
    return(edge_indicators)
  } else {
    stop(paste0("To access the sampled edge indicators the bgms package needs to be run using \n",
                "save = TRUE."))
  }
}

extract_posterior_inclusion_probabilities <- function(fit_object) {
  arguments = extract_arguments(fit_object)

  if(arguments$save) {
    edge_means = colMeans(fit_object$gamma)
    no_variables = arguments$no_variables

    posterior_inclusion_probabilities = matrix(0, no_variables, no_variables)
    posterior_inclusion_probabilities[lower.tri(posterior_inclusion_probabilities)] = edge_means
    posterior_inclusion_probabilities = posterior_inclusion_probabilities +
      t(posterior_inclusion_probabilities)

    data_columnnames = arguments$data_columnnames
    colnames(posterior_inclusion_probabilities) = data_columnnames
    rownames(posterior_inclusion_probabilities) = data_columnnames

  } else {

    posterior_inclusion_probabilities = fit_object$gamma

  }
  return(posterior_inclusion_probabilities)
}


extract_edge_priors <- function(fit_object) {
  arguments = extract_arguments(fit_object)

  if(!arguments$edge_selection) {
    stop(paste0("The bgm function did not perform edge selection, so there are no edge priors \n",
                "specified."))
  } else {
    if(arguments$edge_prior == "Bernoulli") {
      edge_prior = list(type = "Bernoulli",
                        prior_inclusion_probability = arguments$inclusion_probability)
    } else {
      edge_prior = list(type = "Beta-Bernoulli",
                        alpha = arguments$beta_bernoulli_alpha,
                        beta = arguments$beta_bernoulli_beta)
    }
  }
  return(edge_prior)
}

extract_pairwise_interactions <- function(fit_object) {
  arguments = extract_arguments(fit_object)

  return(fit_object$interactions)
}

extract_pairwise_thresholds <- function(fit_object) {
  arguments = extract_arguments(fit_object)

  return(fit_object$thresholds)
}