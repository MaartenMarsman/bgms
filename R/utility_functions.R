#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @importFrom methods hasArg

check_model = function(x,
                       variable_type,
                       reference_category,
                       interaction_scale = 2.5,
                       threshold_alpha = 0.5,
                       threshold_beta = 0.5,
                       edge_selection = TRUE,
                       edge_prior = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block"),
                       inclusion_probability = 0.5,
                       beta_bernoulli_alpha = 1,
                       beta_bernoulli_beta = 1,
                       dirichlet_alpha = dirichlet_alpha) {

  #Check variable type input ---------------------------------------------------
  if(length(variable_type) == 1) {
    variable_input = variable_type
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel")),
                        silent = TRUE)
    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  variable_input, "."))
    variable_bool = (variable_type == "ordinal")
    variable_bool = rep(variable_bool, ncol(x))
  } else {
    if(length(variable_type) != ncol(x))
      stop(paste0("The variable type vector variable_type should be either a single character\n",
                  "string or a vector of character strings of length p."))

    variable_input = unique(variable_type)
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel"),
                                  several.ok = TRUE), silent = TRUE)

    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input, collapse = ", "), "."))

    no_types = sapply(variable_input, function(type) {
      tmp = try(match.arg(arg = type,
                          choices = c("ordinal", "blume-capel")),
                silent = TRUE)
      inherits(tmp, what = "try-error")
    })

    if(length(variable_type) != ncol(x))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input[no_types], collapse = ", "), "."))

    variable_bool = (variable_type == "ordinal")
  }

  #Check Blume-Capel variable input --------------------------------------------
  if(any(!variable_bool)) {
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    if(!hasArg("reference_category"))
      stop("The argument reference_category is required for Blume-Capel variables.")

    if(length(reference_category) != ncol(x) && length(reference_category) != 1)
      stop(paste0("The argument reference_category for the Blume-Capel model needs to be a \n",
                  "single integer or a vector of integers of length p."))

    if(length(reference_category) == 1) {
      #Check if the input is integer -------------------------------------------
      integer_check = try(as.integer(reference_category), silent = TRUE)
      if(is.na(integer_check))
        stop(paste0("The reference_category argument for the Blume-Capel model contains either \n",
                    "a missing value or a value that could not be forced into an integer value."))
      integer_check = reference_category - round(reference_category)
      if(integer_check > .Machine$double.eps)
        stop("Reference category needs to an integer value or a vector of integers of length p.")
      reference_category = rep.int(reference_category, times = ncol(x))
    }

    #Check if the input is integer -------------------------------------------
    blume_capel_variables = which(!variable_bool)
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    integer_check = try(as.integer(reference_category[blume_capel_variables]),
                        silent = TRUE)
    if(anyNA(integer_check))
      stop(paste0("The reference_category argument for the Blume-Capel model contains either \n",
                  "missing values or values that could not be forced into an integer value."))

    integer_check = reference_category[blume_capel_variables] -
      round(reference_category[blume_capel_variables])

    if(any(integer_check > .Machine$double.eps)) {
      non_integers = blume_capel_variables[integer_check > .Machine$double.eps]
      if(length(non_integers) > 1) {
        stop(paste0("The entries in reference_category for variables ",
                    paste0(non_integers, collapse = ", "), " need to be integer."))
      } else {
        stop(paste0("The entry in reference_category for variable ",
                    non_integers, " needs to be an integer."))
      }
    }

    variable_lower = apply(x, 2, min, na.rm = TRUE)
    variable_upper = apply(x, 2, max, na.rm = TRUE)

    if(any(reference_category < variable_lower) | any(reference_category > variable_upper)) {
      out_of_range = which(reference_category < variable_lower | reference_category > variable_upper)
      stop(paste0("The Blume-Capel model assumes that the reference category is within the range \n",
                  "of the observed category scores. This was not the case for variable(s) \n",
                  paste0(out_of_range, collapse =", "),
                  "."))
    }

  } else {
    reference_category = rep.int(0, times = ncol(x))
  }

  #Check prior set-up for the interaction parameters ---------------------------
  if(interaction_scale <= 0 || is.na(interaction_scale) || is.infinite(interaction_scale))
    stop("The scale of the Cauchy prior needs to be positive.")

  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0 | !is.finite(threshold_alpha))
    stop("Parameter threshold_alpha needs to be positive.")
  if(threshold_beta <= 0 | !is.finite(threshold_beta))
    stop("Parameter threshold_beta needs to be positive.")

  #Check set-up for the Bayesian edge selection model --------------------------
  edge_selection = as.logical(edge_selection)
  if(is.na(edge_selection))
    stop("The parameter edge_selection needs to be TRUE or FALSE.")
  if(edge_selection == TRUE) {
    #Check prior set-up for the edge indicators --------------------------------
    edge_prior = match.arg(edge_prior)
    if(edge_prior == "Bernoulli") {
      if(length(inclusion_probability) == 1) {
        theta = inclusion_probability[1]
        if(is.na(theta) || is.null(theta))
          stop("There is no value specified for the inclusion probability.")
        if(theta <= 0)
          stop("The inclusion probability needs to be positive.")
        if(theta > 1)
          stop("The inclusion probability cannot exceed the value one.")
        if(theta == 1)
          stop("The inclusion probability cannot equal one.")

        theta = matrix(theta, nrow = ncol(x), ncol = ncol(x))
      } else {
        if(!inherits(inclusion_probability, what = "matrix") &&
           !inherits(inclusion_probability, what = "data.frame"))
          stop("The input for the inclusion probability argument needs to be a single number, matrix, or dataframe.")

        if(inherits(inclusion_probability, what = "data.frame")) {
          theta = data.matrix(inclusion_probability)
        } else {
          theta = inclusion_probability
        }
        if(!isSymmetric(theta))
          stop("The inclusion probability matrix needs to be symmetric.")
        if(ncol(theta) != ncol(x))
          stop("The inclusion probability matrix needs to have as many rows (columns) as there are variables in the data.")

        if(anyNA(theta[lower.tri(theta)]) ||
           any(is.null(theta[lower.tri(theta)])))
          stop("One or more elements of the elements in inclusion probability matrix are not specified.")
        if(any(theta[lower.tri(theta)] <= 0))
          stop(paste0("The inclusion probability matrix contains negative or zero values;\n",
                      "inclusion probabilities need to be positive."))
        if(any(theta[lower.tri(theta)] >= 1))
          stop(paste0("The inclusion probability matrix contains values greater than or equal to one;\n",
                      "inclusion probabilities cannot exceed or equal the value one."))
      }
    }
    if(edge_prior == "Beta-Bernoulli") {
      theta = matrix(0.5, nrow = ncol(x), ncol = ncol(x))
      if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0)
        stop("The scale parameters of the beta distribution need to be positive.")
      if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta))
        stop("The scale parameters of the beta distribution need to be finite.")
      if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
         is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta))
        stop("Values for both scale parameters of the beta distribution need to be specified.")
    }
    if(edge_prior == "Stochastic-Block") {
      theta = matrix(0.5, nrow = ncol(x), ncol = ncol(x))
      if(beta_bernoulli_alpha <= 0 || beta_bernoulli_beta <= 0 || dirichlet_alpha <= 0)
        stop("The scale parameters of the beta and Dirichlet distribution need to be positive.")
      if(!is.finite(beta_bernoulli_alpha) || !is.finite(beta_bernoulli_beta) || !is.finite(dirichlet_alpha))
        stop("The scale parameters of the beta and Dirichlet distribution need to be finite.")
      if(is.na(beta_bernoulli_alpha) || is.na(beta_bernoulli_beta) ||
         is.null(beta_bernoulli_alpha) || is.null(beta_bernoulli_beta) ||
         is.null(dirichlet_alpha) || is.null(dirichlet_alpha))
        stop("Values for both scale parameters of the beta and Dirichlet distribution need to be specified.")
    }
  } else {
    theta = matrix(0.5, nrow = 1, ncol = 1)
    edge_prior = "Not Applicable"
  }

  return(list(variable_bool = variable_bool,
              reference_category = reference_category,
              edge_selection = edge_selection,
              edge_prior = edge_prior,
              theta = theta))
}

check_compare_model = function(x,
                               y,
                               g,
                               difference_selection,
                               variable_type,
                               reference_category,
                               pairwise_difference_scale = 2.5,
                               main_difference_scale = 2.5,
                               pairwise_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                               main_difference_prior = c("Bernoulli", "Beta-Bernoulli"),
                               pairwise_difference_probability = 0.5,
                               main_difference_probability = 0.5,
                               main_beta_bernoulli_alpha = 1,
                               main_beta_bernoulli_beta = 1,
                               pairwise_beta_bernoulli_alpha = 1,
                               pairwise_beta_bernoulli_beta = 1,
                               interaction_scale = 2.5,
                               threshold_alpha = 0.5,
                               threshold_beta = 0.5,
                               main_difference_model = c("Free", "Collapse", "Constrain")) {

  if(!is.null(g)) {
    unique_g = unique(g)
    if(length(unique_g) == 0)
      stop(paste0("The bgmCompare function expects at least two groups, but the input g contains\n",
                  "no group value."))
    if(length(unique_g) == 1)
      stop(paste0("The bgmCompare function expects at least two groups, but the input g contains\n",
                  "only one group value."))
    if(length(unique_g) == length(g))
      stop(paste0("The bgmCompare function expects at least two groups, but the input g contains\n",
                  "only unique group values."))

    group = g
    for(u in unique_g) {
      group[g == u] = which(unique_g == u)
    }
    tab = tabulate(group)

    if(any(tab < 2))
      stop("One or more groups only had one member in the input g.")
  } else {
    group = c(rep.int(1, times = nrow(x)), rep.int(2, times = nrow(y)))
    x = rbind(x, y)
  }

  #Check variable type input ---------------------------------------------------
  if(length(variable_type) == 1) {
    variable_input = variable_type
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel")),
                        silent = TRUE)
    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  variable_input, "."))
    variable_bool = (variable_type == "ordinal")
    variable_bool = rep(variable_bool, ncol(x))
  } else {
    if(length(variable_type) != ncol(x))
      stop(paste0("The variable type vector variable_type should be either a single character\n",
                  "string or a vector of character strings of length p."))

    variable_input = unique(variable_type)
    variable_type = try(match.arg(arg = variable_type,
                                  choices = c("ordinal", "blume-capel"),
                                  several.ok = TRUE), silent = TRUE)

    if(inherits(variable_type, what = "try-error"))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input, collapse = ", "), "."))

    num_types = sapply(variable_input, function(type) {
      tmp = try(match.arg(arg = type,
                          choices = c("ordinal", "blume-capel")),
                silent = TRUE)
      inherits(tmp, what = "try-error")
    })

    if(length(variable_type) != ncol(x))
      stop(paste0("The bgm function supports variables of type ordinal and blume-capel, \n",
                  "but not of type ",
                  paste0(variable_input[num_types], collapse = ", "), "."))

    variable_bool = (variable_type == "ordinal")
  }

  #Check Blume-Capel variable input --------------------------------------------
  if(any(!variable_bool)) {
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    if(!hasArg("reference_category"))
      stop("The argument reference_category is required for Blume-Capel variables.")

    if(length(reference_category) != ncol(x) && length(reference_category) != 1)
      stop(paste0("The argument reference_category for the Blume-Capel model needs to be a \n",
                  "single integer or a vector of integers of length p."))

    if(length(reference_category) == 1) {
      #Check if the input is integer -------------------------------------------
      integer_check = try(as.integer(reference_category), silent = TRUE)
      if(is.na(integer_check))
        stop(paste0("The reference_category argument for the Blume-Capel model contains either \n",
                    "a missing value or a value that could not be forced into an integer value."))
      integer_check = reference_category - round(reference_category)
      if(integer_check > .Machine$double.eps)
        stop("Reference category needs to an integer value or a vector of integers of length p.")
      reference_category = rep.int(reference_category, times = ncol(x))
    }

    #Check if the input is integer -------------------------------------------
    blume_capel_variables = which(!variable_bool)
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)

    integer_check = try(as.integer(reference_category[blume_capel_variables]),
                        silent = TRUE)
    if(anyNA(integer_check))
      stop(paste0("The reference_category argument for the Blume-Capel model contains either \n",
                  "missing values or values that could not be forced into an integer value."))

    integer_check = reference_category[blume_capel_variables] -
      round(reference_category[blume_capel_variables])

    if(any(integer_check > .Machine$double.eps)) {
      non_integers = blume_capel_variables[integer_check > .Machine$double.eps]
      if(length(non_integers) > 1) {
        stop(paste0("The entries in reference_category for variables ",
                    paste0(non_integers, collapse = ", "), " need to be integer."))
      } else {
        stop(paste0("The entry in reference_category for variable ",
                    non_integers, " needs to be an integer."))
      }
    }

    variable_lower = apply(x, 2, min, na.rm = TRUE)
    variable_upper = apply(x, 2, max, na.rm = TRUE)

    if(any(reference_category < variable_lower) | any(reference_category > variable_upper)) {
      out_of_range = which(reference_category < variable_lower | reference_category > variable_upper)
      stop(paste0("The Blume-Capel model assumes that the reference category is within the range \n",
                  "of the observed category scores. This was not the case for variable(s) \n",
                  paste0(out_of_range, collapse =", "),
                  "."))
    }
  } else {
    reference_category = rep.int(0, times = ncol(x))
  }

  #Check prior set-up for the interaction parameters ---------------------------
  if(interaction_scale <= 0 || is.na(interaction_scale) || is.infinite(interaction_scale))
    stop("The scale of the Cauchy prior for the interactions needs to be positive.")

  #Check prior set-up for the interaction differences --------------------------
  if(pairwise_difference_scale <= 0 || is.na(pairwise_difference_scale) || is.infinite(pairwise_difference_scale))
    stop("The scale of the Cauchy prior for the difference in pairwise interactions needs to be positive.")

  #Check prior set-up for the interaction differences --------------------------
  if(main_difference_scale <= 0 || is.na(main_difference_scale) || is.infinite(main_difference_scale))
    stop("The scale of the Cauchy prior for the difference in category thresholds needs to be positive.")

  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0 | !is.finite(threshold_alpha))
    stop("Parameter threshold_alpha needs to be positive.")
  if(threshold_beta <= 0 | !is.finite(threshold_beta))
    stop("Parameter threshold_beta needs to be positive.")

  #Check set-up for the Bayesian difference selection model --------------------
  difference_selection = as.logical(difference_selection)
  if(is.na(difference_selection))
    stop("The parameter difference_selection needs to be TRUE or FALSE.")
  if(difference_selection == TRUE) {
    inclusion_probability_difference = matrix(0,
                                              nrow = ncol(x),
                                              ncol = ncol(x))

    main_difference_prior = match.arg(main_difference_prior)
    if(main_difference_prior == "Bernoulli") {
      if(length(main_difference_probability) == 1) {
        main_difference_inclusion_probability = main_difference_probability[1]
        if(is.na(main_difference_inclusion_probability) || is.null(main_difference_inclusion_probability))
          stop("There is no value specified for the inclusion probability for category threshold differences.")
        if(main_difference_inclusion_probability <= 0)
          stop("The inclusion probability for category threshold differences needs to be positive.")
        if(main_difference_inclusion_probability > 1)
          stop("The inclusion probability for category threshold differences cannot exceed the value one.")
        if(main_difference_inclusion_probability == 1)
          stop("The inclusion probability for category threshold differences cannot equal one.")

        diag(inclusion_probability_difference) = main_difference_inclusion_probability

      } else {
        if(!inherits(main_difference_probability, what = "matrix") &&
           !inherits(main_difference_probability, what = "data.frame"))
          stop(paste0("The input for the inclusion probability argument for category threshold \n",
                      "differences needs to be a single number, matrix, or dataframe."))

        if(length(as.vector(main_difference_probability)) != ncol(x)) {
          stop(paste0("The inclusion probability argument for the category threshold \n",
                      "parameters needs to have either one or as many elements as\n",
                      "there are variables in the data."))
        } else {
          diag(inclusion_probability_difference) = as.vector(main_difference_probability)
        }

        if(anyNA(diag(inclusion_probability_difference)) ||
           any(is.null(diag(inclusion_probability_difference))))
          stop(paste0("One or more of the inclusion probabilities for the category threshold\n",
                      "parameters are not specified."))
        if(any(diag(inclusion_probability_difference) <= 0))
          stop(paste0("One or more of the inclusion probabilities for the category threshold\n",
                      "parameters are non-positive."))
        if(any(diag(inclusion_probability_difference) >= 1))
          stop(paste0("One or more of the inclusion probabilities for the category threshold\n",
                      "parameters are larger than one."))
      }

    } else {
      diag(inclusion_probability_difference) = 0.5
      if(main_beta_bernoulli_alpha <= 0 || main_beta_bernoulli_beta <= 0)
        stop(paste0("The scale parameters of the beta distribution for the differences in category\n",
                    "thresholds need to be positive."))
      if(!is.finite(main_beta_bernoulli_alpha) || !is.finite(main_beta_bernoulli_beta))
        stop(paste0("The scale parameters of the beta distribution for the differences in category\n",
                    "thresholds need to be finite."))
      if(is.na(main_beta_bernoulli_alpha) || is.na(main_beta_bernoulli_beta) ||
         is.null(main_beta_bernoulli_alpha) || is.null(main_beta_bernoulli_beta))
        stop(paste0("The scale parameters of the beta distribution for the differences in category\n",
                    "thresholds need to be specified."))
    }

    pairwise_difference_prior = match.arg(pairwise_difference_prior)
    if(pairwise_difference_prior == "Bernoulli") {
      if(length(pairwise_difference_probability) == 1) {
        pairwise_difference_inclusion_probability = pairwise_difference_probability[1]
        if(is.na(pairwise_difference_inclusion_probability) || is.null(pairwise_difference_inclusion_probability))
          stop(paste0("There is no value specified for the inclusion probability for the differences \n",
                      "in the pairwise interactions."))
        if(pairwise_difference_inclusion_probability <= 0)
          stop(paste0("The inclusion probability for the differences in the pairwise interactions \n",
                      "needs to be positive."))
        if(pairwise_difference_inclusion_probability > 1)
          stop(paste0("The inclusion probability for the differences in the pairwise interactions\n","
                      cannot exceed the value one."))
        if(pairwise_difference_inclusion_probability == 1)
          stop(paste0("The inclusion probability for the differences in the pairwise interactions\n", "
                      cannot equal one."))

        inclusion_probability_difference[lower.tri(inclusion_probability_difference)] =
          pairwise_difference_inclusion_probability
        inclusion_probability_difference[upper.tri(inclusion_probability_difference)] =
          pairwise_difference_inclusion_probability

      } else {
        if(!inherits(pairwise_difference_probability, what = "matrix") &&
           !inherits(pairwise_difference_probability, what = "data.frame"))
          stop(paste0("The input for the inclusion probability argument for differences in the \n",
                      "pairwise interactions needs to be a single number, matrix, or dataframe."))

        if(inherits(pairwise_difference_probability, what = "data.frame")) {
          tmp = data.matrix(pairwise_difference_probability)
          diag(tmp) = diag(inclusion_probability_difference)
          inclusion_probability_difference = tmp
        } else {
          tmp = pairwise_difference_probability
          diag(tmp) = diag(inclusion_probability_difference)
          inclusion_probability_difference = tmp
        }

        if(!isSymmetric(inclusion_probability_difference))
          stop("The inclusion probability matrix needs to be symmetric.")
        if(ncol(inclusion_probability_difference) != ncol(x))
          stop(paste0("The inclusion probability matrix needs to have as many rows (columns) as there\n",
                      " are variables in the data."))

        if(anyNA(inclusion_probability_difference[lower.tri(inclusion_probability_difference)]) ||
           any(is.null(inclusion_probability_difference[lower.tri(inclusion_probability_difference)])))
          stop(paste0("One or more inclusion probabilities for differences in the pairwise \n",
                      "interactions are not specified."))
        if(any(inclusion_probability_difference[lower.tri(inclusion_probability_difference)] <= 0))
          stop(paste0("One or more inclusion probabilities for differences in the pairwise \n",
                      "interactions are negative or equal to zero."))
        if(any(inclusion_probability_difference[lower.tri(inclusion_probability_difference)] >= 1))
          stop(paste0("One or more inclusion probabilities for differences in the pairwise \n",
                      "interactions are larger than one."))
      }
    } else {
      inclusion_probability_difference[lower.tri(inclusion_probability_difference)] =
        0.5
      inclusion_probability_difference[upper.tri(inclusion_probability_difference)] =
        0.5

      if(pairwise_beta_bernoulli_alpha <= 0 || pairwise_beta_bernoulli_beta <= 0)
        stop(paste0("The scale parameters of the beta distribution for the differences in pairwise\n",
                    "interactions need to be positive."))
      if(!is.finite(pairwise_beta_bernoulli_alpha) || !is.finite(pairwise_beta_bernoulli_beta))
        stop(paste0("The scale parameters of the beta distribution for the differences in pairwise\n",
                    "interactions need to be finite."))
      if(is.na(pairwise_beta_bernoulli_alpha) || is.na(pairwise_beta_bernoulli_beta) ||
         is.null(pairwise_beta_bernoulli_alpha) || is.null(pairwise_beta_bernoulli_beta))
        stop(paste0("The scale parameters of the beta distribution for the differences in pairwise\n",
                    "differences need to be specified."))
    }
  } else {
    main_difference_prior = "Not applicable"
    pairwise_difference_prior = "Not applicable"
    inclusion_probability_difference = matrix(0.5, 1, 1)
  }

  main_difference_model = match.arg(main_difference_model)

  u_g = unique(g)
  if(length(u_g) > 2 & main_difference_model == "Constrain")
    stop(paste0("The option Constrain for the main difference model is currently only supported \n",
                "for designs with two groups. Please try again with Free or Collapse."))

  if(main_difference_model == "Constrain")
    warning(paste0("The option Constrain for the main difference model is deprecated and will be\n",
                   "removed in a future version. Please use the Collapse option."))

  return(list(x = x,
              group = group,
              variable_bool = variable_bool,
              reference_category = reference_category,
              main_difference_prior = main_difference_prior,
              pairwise_difference_prior = pairwise_difference_prior,
              inclusion_probability_difference = inclusion_probability_difference,
              main_difference_model = main_difference_model))
}

reformat_data = function(x, na.action, variable_bool, reference_category) {
  if(na.action == "listwise") {
    # Check for missing values ---------------------------------------------------
    missing_values = sapply(1:nrow(x), function(row){anyNA(x[row, ])})
    if(sum(missing_values) == nrow(x))
      stop(paste0("All rows in x contain at least one missing response.\n",
                  "You could try option na.action = impute."))
    if(sum(missing_values) > 1)
      warning(paste0("There were ",
                     sum(missing_values),
                     " rows with missing observations in the input matrix x.\n",
                     "Since na.action = listwise these rows were excluded from the analysis."),
              call. = FALSE)
    if(sum(missing_values) == 1)
      warning(paste0("There was one row with missing observations in the input matrix x.\n",
                     "Since na.action = listwise this row was excluded from \n",
                     "the analysis."),
              call. = FALSE)
    x = x[!missing_values, ]

    if(ncol(x) < 2 || is.null(ncol(x)))
      stop(paste0("After removing missing observations from the input matrix x,\n",
                  "there were less than two columns left in x."))
    if(nrow(x) < 2 || is.null(nrow(x)))
      stop(paste0("After removing missing observations from the input matrix x,\n",
                  "there were less than two rows left in x."))

    missing_index = matrix(NA, nrow = 1, ncol = 1)
    na_impute = FALSE
  } else {
    # Check for missing values -------------------------------------------------
    no_missings = sum(is.na(x))
    no_persons = nrow(x)
    no_variables = ncol(x)
    if(no_missings > 0) {
      missing_index = matrix(0, nrow = no_missings, ncol = 2)
      na_impute = TRUE
      cntr = 0
      for(node in 1:no_variables) {
        mis = which(is.na(x[, node]))
        if(length(mis) > 0) {
          for(i in 1:length(mis)) {
            cntr = cntr + 1
            missing_index[cntr, 1] = mis[i] - 1                                 #c++ index starts at 0
            missing_index[cntr, 2] = node - 1                                   #c++ index starts at 0
            x[mis[i], node] = sample(x[-mis, node],                             #start value for imputation
                                     size = 1)
            #This is non-zero if no zeroes are observed (we then collapse over zero below)
          }
        }
      }
    } else {
      missing_index = matrix(NA, nrow = 1, ncol = 1)
      na_impute = FALSE
    }
  }

  check_fail_zero = FALSE
  no_variables = ncol(x)
  no_categories = vector(length = no_variables)
  for(node in 1:no_variables) {
    unq_vls = sort(unique(x[,  node]))
    mx_vl = max(unq_vls)

    # Check if observed responses are not all unique ---------------------------
    if(mx_vl == nrow(x))
      stop(paste0("Only unique responses observed for variable ",
                  node,
                  ". We expect >= 1 observations per category."))

    # Recode data --------------------------------------------------------------
    if(variable_bool[node]) {#Regular ordinal variable
      # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
      if(length(unq_vls) != mx_vl + 1 || any(unq_vls != 0:mx_vl)) {
        y = x[, node]
        cntr = 0
        for(value in unq_vls) {
          x[y == value, node] = cntr
          cntr = cntr + 1
        }
      }
    } else {#Blume-Capel ordinal variable
      # Check if observations are integer or can be recoded --------------------
      if (any(abs(unq_vls - round(unq_vls)) > .Machine$double.eps)) {
        int_unq_vls = unique(as.integer(unq_vls))
        if(anyNA(int_unq_vls)) {
          stop(paste0(
            "The Blume-Capel model assumes that its observations are coded as integers, but \n",
            "the category scores for node ", node, " were not integer. An attempt to recode \n",
            "them to integer failed. Please inspect the documentation for the base R \n",
            "function as.integer(), which bgm uses for recoding category scores."))
        }

        if(length(int_unq_vls) != length(unq_vls)) {
          stop(paste0("The Blume-Capel model assumes that its observations are coded as integers. The \n",
                      "category scores of the observations for node ", node, " were not integers. An \n",
                      "attempt to recode these observations as integers failed because, after rounding, \n",
                      "a single integer value was used for several observed score categories."))
        }
        x[, node] = as.integer(x[, node])

        if(reference_category[node] < 0 | reference_category[node] > max(x[, node]))
          stop(paste0(
            "The reference category for the Blume-Capel variable ", node, "is outside its \n",
            "range of observations."))
      }

      # Check if observations start at zero and recode otherwise ---------------
      if(min(x[, node]) != 0) {
        reference_category[node] = reference_category[node] - min(x[, node])
        x[, node] = x[, node] - min(x[, node])

        if(check_fail_zero == FALSE) {
          check_fail_zero = TRUE
          failed_zeroes = c(node)
        } else {
          failed_zeroes = c(failed_zeroes, node)
        }
      }

      check_range = length(unique(x[, node]))
      if(check_range < 3)
        stop(paste0("The Blume-Capel is only available for variables with more than one category \n",
                    "observed. There two or less categories observed for variable ",
                    node,
                    "."))
    }

    # Warn that maximum category value is large --------------------------------
    no_categories[node] = max(x[,node])
    if(!variable_bool[node] & no_categories[node] > 10) {
      # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
      warning(paste0("In the (pseudo) likelihood of Blume-Capel variables, the normalization constant \n",
                     "is a sum over all possible values of the ordinal variable. The range of \n",
                     "observed values, possibly after recoding to integers, is assumed to be the \n",
                     "number of possible response categories.  For node ", node,", this range was \n",
                     "equal to ", no_categories[node], "which may cause the analysis to take some \n",
                     "time to run. Note that for the Blume-Capel model, the bgm function does not \n",
                     "collapse the categories that have no observations between zero and the last \n",
                     "category. This may explain the large discrepancy between the first and last \n",
                     "category values."))
    }

    # Check to see if not all responses are in one category --------------------
    if(no_categories[node] == 0)
      stop(paste0("Only one value [",
                  unq_vls,
                  "] was observed for variable ",
                  node,
                  "."))
  }

  if(check_fail_zero == TRUE) {
    if(length(failed_zeroes) == 1) {
      node = failed_zeroes[1]
      warning(paste0("The bgm function assumes that the observed ordinal variables are integers and \n",
                     "that the lowest observed category score is zero. The lowest score for node \n",
                     node, " was recoded to zero for the analysis.\n",
                     "Note that bgm also recoded the corresponding reference category score to ", reference_category[node], "."))
    } else {
      warning(paste0("The bgm function assumes that the observed ordinal variables are integers and \n",
                     "that the lowest observed category score is zero. The lowest score for nodes \n",
                     paste(failed_zeroes, collapse = ","), " were recoded to zero for the analysis.\n",
                     "Note that bgm also recoded the corresponding reference category scores."))
    }
  }

  return(list(x = x,
              no_categories = no_categories,
              reference_category = reference_category,
              missing_index = missing_index,
              na_impute = na_impute))
}

compare_reformat_data = function(x,
                                 group,
                                 na.action,
                                 variable_bool,
                                 reference_category,
                                 main_difference_model) {
  if(na.action == "listwise") {
    # Check for missing values in x --------------------------------------------
    missing_values = sapply(1:nrow(x), function(row){anyNA(x[row, ])})
    if(sum(missing_values) == nrow(x))
      stop(paste0("All rows in x contain at least one missing response.\n",
                  "You could try option na.action = impute."))
    if(sum(missing_values) > 1)
      warning(paste0("There were ",
                     sum(missing_values),
                     " rows with missing observations in the input matrix x.\n",
                     "Since na.action = listwise these rows were excluded from the analysis."),
              call. = FALSE)
    if(sum(missing_values) == 1)
      warning(paste0("There was one row with missing observations in the input matrix x.\n",
                     "Since na.action = listwise this row was excluded from \n",
                     "the analysis."),
              call. = FALSE)

    x = x[!missing_values, ]
    group = group[!missing_values]

    if(nrow(x) < 2 || is.null(nrow(x)))
      stop(paste0("After removing missing observations from the input matrix x,\n",
                  "there were less than two rows left in x."))

    unique_g = unique(group)
    if(length(unique_g) == length(group))
      stop(paste0("After rows with missing observations were excluded, there were no groups, as \n",
                  "there were only unique values in the input g left."))
    if(length(unique_g) == 1)
      stop(paste0("After rows with missing observations were excluded, there were no groups, as \n",
                  "there was only one value in the input g left."))
    g = group
    for(u in unique_g) {
      group[g == u] = which(unique_g == u)
    }
    tab = tabulate(group)

    if(any(tab < 2))
      stop(paste0("After rows with missing observations were excluded, one or more groups, only \n",
                  "had one member in the input g."))

    missing_index = matrix(NA, nrow = 1, ncol = 1)
    na_impute = FALSE
  } else {
    # Check for missing values in x --------------------------------------------
    num_missings = sum(is.na(x))
    num_persons = nrow(x)
    num_variables = ncol(x)
    if(num_missings > 0) {
      missing_index = matrix(0, nrow = num_missings, ncol = 2)
      na_impute = TRUE
      cntr = 0
      for(node in 1:num_variables) {
        mis = which(is.na(x[, node]))
        if(length(mis) > 0) {
          for(i in 1:length(mis)) {
            cntr = cntr + 1
            missing_index[cntr, 1] = mis[i] - 1                             #c++ index starts at 0
            missing_index[cntr, 2] = node - 1                               #c++ index starts at 0
            x[mis[i], node] = sample(x[-mis, node], size = 1)               #start value for imputation
            #This is non-zero if no zeroes are observed (we then collapse over zero below)
          }
        }
      }
    } else {
      missing_index = matrix(NA, nrow = 1, ncol = 1)
      na_impute = FALSE
    }
  }

  check_fail_zero = FALSE
  num_variables = ncol(x)
  num_categories = matrix(0,
                           nrow = num_variables,
                           ncol = max(group))

  for(node in 1:num_variables) {
    unq_vls = sort(unique(x[,  node]))
    mx_vls = length(unq_vls)

    # Check if observed responses are not all unique ---------------------------
    if(mx_vls == nrow(x))
      stop(paste0("Only unique responses observed for variable ",
                  node,
                  " in the matrix x (group 1). We expect >= 1 observations per category."))

    # Recode data --------------------------------------------------------------
    if(variable_bool[node]) {#Regular ordinal variable
      # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
      if(main_difference_model == "Collapse") {

        observed_scores = matrix(NA,
                                 nrow = mx_vls,
                                 ncol = max(group))

        for(value in unq_vls) {
          unique_g = unique(group)
          for(g in unique_g) {
            observed_scores[which(unq_vls == value), g] =
              any(x[group == g, node] == value) * 1
          }
        }

        xx = x[, node]
        cntr = -1
        for(value in unq_vls) {
          #Collapse categories when not observed in one or more groups.
          if(sum(observed_scores[which(unq_vls == value), ]) == max(group)) {
            cntr = cntr + 1 #increment score if category observed in all groups
          }
          x[xx == value, node] = max(0, cntr)
        }

      } else {
        unique_g = unique(group)
        for(g in unique_g) {
          x_g = x[group == g,  node]
          unq_vls_g = sort(unique(x_g))

          cntr = 0
          for(value in unq_vls_g) {
            x_g[x_g == value] = cntr
            cntr = cntr + 1
          }
          x[group == g,  node] = x_g
        }
      }
    } else {#Blume-Capel ordinal variable
      # Check if observations are integer or can be recoded --------------------
      if (any(abs(unq_vls - round(unq_vls)) > .Machine$double.eps)) {
        int_unq_vls = unique(as.integer(unq_vls))
        if(anyNA(int_unq_vls)) {
          stop(paste0(
            "The Blume-Capel model assumes that its observations are coded as integers, but \n",
            "the category scores for node ", node, " were not integer. An attempt to recode \n",
            "them to integer failed. Please inspect the documentation for the base R function \n",
            "as.integer(), which bgmCompare uses for recoding category scores."))
        }

        if(length(int_unq_vls) != length(unq_vls)) {
          stop(paste0("The Blume-Capel model assumes that its observations are coded as integers. The \n",
                      "category scores of the observations for node ", node, " were not integers. An \n",
                      "attempt to recode these observations as integers failed because, after rounding,\n",
                      "a single integer value was used for several observed score categories."))
        }
        x[, node] = as.integer(x[, node])
      }

      mi = min(x[,node])

      ma = max(x[,node])

      if(reference_category[node] < mi | reference_category[node] > ma)
        stop(paste0(
          "The reference category for the Blume-Capel variable ", node, "is outside its \n",
          "range of observations in the matrices x (and y)."))

      # Check if observations start at zero and recode otherwise ---------------
      if(mi != 0) {
        reference_category[node] = reference_category[node] - mi
        x[, node] = x[, node] - mi

        if(check_fail_zero == FALSE) {
          check_fail_zero = TRUE
          failed_zeroes = c(node)
        } else {
          failed_zeroes = c(failed_zeroes, node)
        }
      }

      check_range = length(unique(x[, node]))

      if(check_range < 3)
        stop(paste0("The Blume-Capel is only available for variables with more than two categories \n",
                    "observed. There are two or less categories observed for variable ",
                    node,
                    "."))
    }


    # Warn that maximum category value is large --------------------------------
    if(main_difference_model == "Free") {

      num_categories[node, ] = sapply(1:max(group), function(g) {
        max(x[group == g, node])
      })

      if(!variable_bool[node] & max(num_categories[node, ]) > 10) {
        # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
        warning(paste0("In the (pseudo) likelihood of Blume-Capel variables, the normalization constant \n",
                       "is a sum over all possible values of the ordinal variable. The range of \n",
                       "observed values, possibly after recoding to integers, is assumed to be the \n",
                       "number of possible response categories.  For node ", node,", this range was \n",
                       "equal to ", max(num_categories[node,]), ", which may cause the analysis to take some \n",
                       "time to run. Note that for the Blume-Capel model, the bgm function does not \n",
                       "collapse the categories that have no observations between zero and the last \n",
                       "category. This may explain the large discrepancy between the first and last \n",
                       "category values."))
      }

      if(any(num_categories[node, ] == 0))
        stop(paste0("Only one value was observed for variable ",
                    node,
                    ", in at least one of the groups."))
    } else {
      num_categories[node, ] = max(x[, node])

      if(!variable_bool[node] & max(num_categories[node, ]) > 10) {
        # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
        warning(paste0("In the (pseudo) likelihood of Blume-Capel variables, the normalization constant \n",
                       "is a sum over all possible values of the ordinal variable. The range of \n",
                       "observed values, possibly after recoding to integers, is assumed to be the \n",
                       "number of possible response categories.  For node ", node,", in group 1, this \n",
                       "range was equal to ", max(num_categories[node,]), "which may cause the analysis to take some \n",
                       "time to run. Note that for the Blume-Capel model, the bgm function does not \n",
                       "collapse the categories that have no observations between zero and the last \n",
                       "category. This may explain the large discrepancy between the first and last \n",
                       "category values."))
      }


      # Check to see if not all responses are in one category --------------------
      if(any(num_categories[node, ] == 0))
        stop(paste0("Only one value was observed for variable ",
                    node,
                    ", in at least one of the groups."))
    }
  }

  if(check_fail_zero == TRUE) {
    if(length(failed_zeroes) == 1) {
      node = failed_zeroes[1]
      warning(paste0("The bgm function assumes that the observed ordinal variables are integers and \n",
                     "that the lowest observed category score is zero. The lowest score for node \n",
                     node, " was recoded to zero for the analysis. Note that bgm also recoded the \n",
                     "the corresponding reference category score to ", reference_category[node], "."))
    } else {
      warning(paste0("The bgm function assumes that the observed ordinal variables are integers and \n",
                     "that the lowest observed category score is zero. The lowest score for nodes \n",
                     paste(failed_zeroes, collapse = ","), " were recoded to zero for the analysis. Note that bgm also recoded the \n",
                     "the corresponding reference category scores."))
    }
  }

  return(list(x = x,
              group = group,
              num_categories = num_categories,
              reference_category = reference_category,
              missing_index = missing_index,
              na_impute = na_impute))
}

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


# Helper function for data checks
data_check = function(data, name) {
  if (!inherits(data, c("matrix", "data.frame"))) {
    stop(paste(name, "must be a matrix or data.frame."))
  }
  if (inherits(data, "data.frame")) {
    data = data.matrix(data)
  }
  if (nrow(data) < 2 || ncol(data) < 2) {
    stop(paste(name, "must have at least 2 rows and 2 columns."))
  }
  return(data)
}

check_positive_integer = function(value, name) {
  if (!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value <= 0) {
    stop(sprintf("Parameter `%s` must be a positive integer. Got: %s", name, value))
  }
}

# Helper function for validating non-negative integers
check_non_negative_integer = function(value, name) {
  if (!is.numeric(value) || abs(value - round(value)) > .Machine$double.eps || value < 0) {
    stop(sprintf("Parameter `%s` must be a non-negative integer. Got: %s", name, value))
  }
}

# Helper function for validating logical inputs
check_logical = function(value, name) {
  value = as.logical(value)
  if (is.na(value)) {
    stop(sprintf("Parameter `%s` must be TRUE or FALSE. Got: %s", name, value))
  }
  return(value)
}

# Helper function for computing `num_obs_categories`
compute_num_obs_categories = function(x, num_categories, group = NULL) {
  num_obs_categories = list()
  for (g in unique(group)) {
    num_obs_categories_gr = matrix(0, nrow = max(num_categories[, g]), ncol = ncol(x))
    for (variable in seq_len(ncol(x))) {
      for (category in seq_len(num_categories[variable, g])) {
        num_obs_categories_gr[category, variable] = sum(x[group == g, variable] == category)
      }
    }
    num_obs_categories[[g]] = num_obs_categories_gr
  }
  return(num_obs_categories)

}

# Helper function for computing sufficient statistics for Blume-Capel variables
compute_sufficient_blume_capel = function(x, reference_category, ordinal_variable, group = NULL) {
  if (is.null(group)) {  # Two-group design
    sufficient_stats = matrix(0, nrow = 2, ncol = ncol(x))
    bc_vars = which(!ordinal_variable)
    for (i in bc_vars) {
      sufficient_stats[1, i] = sum(x[, i])
      sufficient_stats[2, i] = sum((x[, i] - reference_category[i]) ^ 2)
    }
    return(sufficient_stats)
  } else {  # Multi-group design
    sufficient_stats = list()
    for (g in unique(group)) {
      sufficient_stats_gr = matrix(0, nrow = 2, ncol = ncol(x))
      bc_vars = which(!ordinal_variable)
      for (i in bc_vars) {
        sufficient_stats_gr[1, i] = sum(x[group == g, i])
        sufficient_stats_gr[2, i] = sum((x[group == g, i] - reference_category[i]) ^ 2)
      }
      sufficient_stats[[g]] = sufficient_stats_gr
    }
    return(sufficient_stats)
  }
}

# Prepare the output
prepare_output_bgmCompare = function(out, x, independent_thresholds,
                                     num_variables, num_categories, group, iter,
                                     data_columnnames, save_options,
                                     difference_selection, na_action,
                                     na_impute, variable_type, burnin,
                                     interaction_scale, threshold_alpha,
                                     threshold_beta, main_difference_model,
                                     pairwise_difference_prior,
                                     main_difference_prior,
                                     inclusion_probability_difference,
                                     pairwise_beta_bernoulli_alpha,
                                     pairwise_beta_bernoulli_beta,
                                     main_beta_bernoulli_alpha,
                                     main_beta_bernoulli_beta,
                                     main_difference_scale,
                                     pairwise_difference_scale,
                                     projection, is_ordinal_variable) {

  save = any(c(save_options$save_main, save_options$save_pairwise, save_options$save_indicator))

  arguments = list(
    num_variables = num_variables, group = group, num_cases = tabulate(group),
    na.action = na.action, na_impute = na_impute, iter = iter, burnin = burnin,
    difference_selection = difference_selection,
    independent_thresholds = independent_thresholds,
    save_main = save_options$save_main,
    save_pairwise = save_options$save_pairwise,
    save_indicator = save_options$save_indicator, save = save,
    variable_type = variable_type, interaction_scale = interaction_scale,
    threshold_alpha = threshold_alpha, threshold_beta = threshold_beta,
    main_difference_model = main_difference_model,
    pairwise_difference_prior = pairwise_difference_prior,
    main_difference_prior = main_difference_prior,
    inclusion_probability_difference = inclusion_probability_difference,
    pairwise_beta_bernoulli_alpha = pairwise_beta_bernoulli_alpha,
    pairwise_beta_bernoulli_beta = pairwise_beta_bernoulli_beta,
    main_beta_bernoulli_alpha = main_beta_bernoulli_alpha,
    main_beta_bernoulli_beta = main_beta_bernoulli_beta,
    main_difference_scale = main_difference_scale,
    pairwise_difference_scale = pairwise_difference_scale,
    version = packageVersion("bgms")
  )

  # Names for pairwise effects
  names_bycol = matrix(rep(data_columnnames, each = num_variables), ncol = num_variables)
  names_byrow = matrix(rep(data_columnnames, each = num_variables), ncol = num_variables, byrow = TRUE)
  names_comb = matrix(paste0(names_byrow, "-", names_bycol), ncol = num_variables)
  names_vec = names_comb[lower.tri(names_comb)]
  names_vec_t = names_comb[lower.tri(names_comb, diag = TRUE)]

  # Prepare output elements
  results = list()

  # Handle output from the new rcpp function (anova)
  num_groups = length(unique(group))
  num_main = nrow(out$posterior_mean_main)

  # Main effects
  tmp = out$posterior_mean_main
  if(independent_thresholds) {
    posterior_mean_main = matrix(NA, nrow = nrow(tmp), ncol = ncol(tmp))
    for(g in 1:num_groups) {
      posterior_mean_main[, g] = tmp[, g]

      #This can probably be done prettier
      if(any(tmp[,g] == 0))
        posterior_mean_main[tmp[,g] == 0, g] = NA
    }
  } else {
    posterior_mean_main = matrix(NA, nrow = nrow(tmp), ncol = ncol(tmp) + 1)
    posterior_mean_main[, 1] = tmp[, 1]
    for (row in 1:nrow(tmp)) {
      posterior_mean_main[row, -1] = projection %*% tmp[row, -1]
    }
  }
  results$posterior_mean_main = posterior_mean_main

  names_variable_categories = vector(length = nrow(tmp))
  counter = 0
  for (v in 1:num_variables) {
    if(is_ordinal_variable[v]) {
      for (c in 1:max(num_categories[v, ])) {
        counter = counter + 1
        names_variable_categories[counter] = paste0(data_columnnames[v], "(", c, ")")
      }
    } else {
      counter = counter + 1
      names_variable_categories[counter] = paste0(data_columnnames[v], "(linear)")
      counter = counter + 1
      names_variable_categories[counter] = paste0(data_columnnames[v], "(quadratic)")
    }

  }
  if(independent_thresholds) {
    dimnames(results$posterior_mean_main) = list(names_variable_categories, paste0("group_", 1:num_groups))
  } else {
    dimnames(results$posterior_mean_main) = list(names_variable_categories, c("overall", paste0("group_", 1:num_groups)))
  }

  # Pairwise effects
  tmp = out$posterior_mean_pairwise
  posterior_mean_pairwise = matrix(0, nrow = nrow(tmp), ncol = ncol(tmp) + 1)
  posterior_mean_pairwise[, 1] = tmp[, 1]
  for (row in 1:nrow(tmp)) {
    posterior_mean_pairwise[row, -1] = projection %*% tmp[row, -1]
  }

  results$posterior_mean_pairwise = posterior_mean_pairwise
  dimnames(results$posterior_mean_pairwise) = list(names_vec, c("overall", paste0("group_", 1:num_groups)))

  if (difference_selection && "posterior_mean_indicator" %in% names(out) ) {
    results$posterior_mean_indicator = out$posterior_mean_indicator
    dimnames(results$posterior_mean_indicator) = list(data_columnnames,
                                                      data_columnnames)
  }


  # Handle main_effect_samples
  if (save_options$save_main && "main_effect_samples" %in% names(out)) {
    main_effect_samples = out$main_effect_samples
    col_names = character()  # Vector to store column names

    if(independent_thresholds) {
      for (gr in 1:num_groups) {
        for (var in 1:num_variables) {
          if(is_ordinal_variable[var]) {
            for (cat in 1:max(num_categories[var, ])) {
              col_names = c(col_names, paste0(data_columnnames[var],"_gr", gr, "(", cat, ")"))
            }
          } else {
            col_names = c(col_names, paste0(data_columnnames[var],"_gr", gr, "(linear)"))
            col_names = c(col_names, paste0(data_columnnames[var],"_gr", gr, "(quadratic)"))
          }
        }
      }

    } else {

      for (gr in 1:num_groups) {
        if(gr == 1) {
          for (var in 1:num_variables) {
            if(is_ordinal_variable[var]) {
              for (cat in 1:max(num_categories[var, ])) {
                col_names = c(col_names, paste0(data_columnnames[var],"_overall", gr, "(", cat, ")"))
              }
            } else {
              col_names = c(col_names, paste0(data_columnnames[var],"_overall", gr, "(linear)"))
              col_names = c(col_names, paste0(data_columnnames[var],"_overall", gr, "(quadratic)"))
            }
          }
        } else {
          for (var in 1:num_variables) {
            if(is_ordinal_variable[var]) {
              for (cat in 1:max(num_categories[var, ])) {
                col_names = c(col_names, paste0(data_columnnames[var],"_contrast_#", gr-1, "(", cat, ")"))
              }
            } else {
              col_names = c(col_names, paste0(data_columnnames[var],"_contrast_#", gr-1, "(linear)"))
              col_names = c(col_names, paste0(data_columnnames[var],"_contrast_#", gr-1, "(quadratic)"))
            }
          }
        }
      }
    }

    dimnames(main_effect_samples) = list(Iter. = 1:nrow(main_effect_samples), col_names)
    results$main_effect_samples = main_effect_samples
  }

  # Handle pairwise_effect_samples
  if (save_options$save_pairwise && "pairwise_effect_samples" %in% names(out)) {
    col_names = character()  # Vector to store column names
    for (v1 in 1:(num_variables - 1)) {
      for (v2 in (v1 + 1):num_variables) {
        col_names = c(col_names, paste0(data_columnnames[v1], "-", data_columnnames[v2]))
        for (gr in 2:num_groups) {
          col_names = c(col_names,
                        paste0("contrast_#",
                               gr-1,
                               "(",
                               data_columnnames[v1],
                               "-",
                               data_columnnames[v2],
                               ")"))
        }
      }
    }

    dimnames(out$pairwise_effect_samples) = list(Iter. = 1:nrow(out$pairwise_effect_samples), col_names)
    results$pairwise_effect_samples = out$pairwise_effect_samples
  }

  # Handle inclusion_indicator_samples
  if (difference_selection && save_options$save_indicator && "inclusion_indicator_samples" %in% names(out)) {
    if(independent_thresholds) {
      #
    } else {
      col_names = character()  # Vector to store column names
      for (v1 in 1:num_variables) {
        for (v2 in v1:num_variables) {
          col_names = c(col_names, paste0(data_columnnames[v1], "-", data_columnnames[v2]))
        }
      }
      dimnames(out$inclusion_indicator_samples) = list(Iter. = 1:nrow(out$inclusion_indicator_samples), col_names)
    }

    results$inclusion_indicator_samples = out$inclusion_indicator_samples
  }

  # Add arguments to the output
  results$arguments = arguments
  class(results) = c("bgmCompare")
  return(results)
}