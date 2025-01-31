#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @importFrom methods hasArg

check_compare_model_dev = function(x,
                                   y,
                                   g,
                                   ttest,
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

    if(ttest == TRUE) {
      variable_lower = apply(rbind(apply(x, 2, min, na.rm = TRUE),
                                   apply(y, 2, min, na.rm = TRUE)), 2, min, na.rm = TRUE)
      variable_upper = apply(rbind(apply(x, 2, max, na.rm = TRUE),
                                   apply(y, 2, max, na.rm = TRUE)), 2, max, na.rm = TRUE)
    } else {
      variable_lower = apply(x, 2, min, na.rm = TRUE)
      variable_upper = apply(x, 2, max, na.rm = TRUE)
    }

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

  return(list(variable_bool = variable_bool,
              reference_category = reference_category,
              main_difference_prior = main_difference_prior,
              pairwise_difference_prior = pairwise_difference_prior,
              inclusion_probability_difference = inclusion_probability_difference,
              main_difference_model = main_difference_model))
}

compare_reformat_data_dev = function(x,
                                     y,
                                     g,
                                     ttest,
                                     na.action,
                                     variable_bool,
                                     reference_category,
                                     main_difference_model) {

  if(ttest == FALSE) {
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
    group = NULL
  }

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
    if(ttest == FALSE)
      group = group[!missing_values]

    if(nrow(x) < 2 || is.null(nrow(x)))
      stop(paste0("After removing missing observations from the input matrix x,\n",
                  "there were less than two rows left in x."))

    if(ttest == TRUE) {
      # Check for missing values in y --------------------------------------------
      missing_values = sapply(1:nrow(y), function(row){anyNA(y[row, ])})
      if(sum(missing_values) == nrow(y))
        stop(paste0("All rows in y contain at least one missing response.\n",
                    "You could try option na.action = impute."))
      if(sum(missing_values) > 1)
        warning(paste0("There were ",
                       sum(missing_values),
                       " rows with missing observations in the input matrix y.\n",
                       "Since na.action = listwise these rows were excluded from the analysis."),
                call. = FALSE)
      if(sum(missing_values) == 1)
        warning(paste0("There was one row with missing observations in the input matrix y.\n",
                       "Since na.action = listwise this row was excluded from \n",
                       "the analysis."),
                call. = FALSE)
      y = y[!missing_values, ]

      if(nrow(y) < 2 || is.null(nrow(y)))
        stop(paste0("After removing missing observations from the input matrix y,\n",
                    "there were less than two rows left in y."))

    } else {
      unique_g = unique(group)
      if(length(unique_g) == length(g))
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
    }

    missing_index_gr1 = matrix(NA, nrow = 1, ncol = 1)
    missing_index_gr2 = matrix(NA, nrow = 1, ncol = 1)
    missing_index = matrix(NA, nrow = 1, ncol = 1)
    na_impute = FALSE
  } else {
    # Check for missing values in x --------------------------------------------
    no_missings_gr1 = sum(is.na(x))
    no_persons_gr1 = nrow(x)
    no_variables = ncol(x)
    if(no_missings_gr1 > 0) {
      missing_index_gr1 = matrix(0, nrow = no_missings_gr1, ncol = 2)
      na_impute = TRUE
      cntr = 0
      for(node in 1:no_variables) {
        mis = which(is.na(x[, node]))
        if(length(mis) > 0) {
          for(i in 1:length(mis)) {
            cntr = cntr + 1
            missing_index_gr1[cntr, 1] = mis[i] - 1                             #c++ index starts at 0
            missing_index_gr1[cntr, 2] = node - 1                               #c++ index starts at 0
            x[mis[i], node] = sample(x[-mis, node],                             #start value for imputation
                                     size = 1)
            #This is non-zero if no zeroes are observed (we then collapse over zero below)
          }
        }
      }
    } else {
      missing_index_gr1 = matrix(NA, nrow = 1, ncol = 1)
      na_impute = FALSE
    }

    if(ttest == TRUE) {
      # Check for missing values in y ------------------------------------------
      no_missings_gr2 = sum(is.na(y))
      no_persons_gr2 = nrow(y)
      no_variables = ncol(y)
      if(no_missings_gr2 > 0) {
        missing_index_gr2 = matrix(0, nrow = no_missings_gr2, ncol = 2)
        na_impute = TRUE
        cntr = 0
        for(node in 1:no_variables) {
          mis = which(is.na(y[, node]))
          if(length(mis) > 0) {
            for(i in 1:length(mis)) {
              cntr = cntr + 1
              missing_index_gr2[cntr, 1] = mis[i] - 1                           #c++ index starts at zero.
              missing_index_gr2[cntr, 2] = node                                 #c++ index starts at zero.
              y[mis[i], node] = sample(y[-mis, node], #start value for imputation
                                       size = 1)
              #This is non-zero if no zeroes are observed (we then collapse over zero below)
            }
          }
        }
      } else {
        missing_index_gr2 = matrix(NA, nrow = 1, ncol = 1)
        if(!na_impute) {
          na_impute = FALSE
        }
      }
      no_missings = NA
      no_persons = NA
      missing_index = matrix(NA, nrow = 1, ncol = 1)
    } else {
      missing_index = missing_index_gr1
      no_missings = no_missings_gr1
      no_persons = no_persons_gr1
      no_missings_gr2 = NA
      no_persons_gr2 = NA
      missing_index_gr2 = matrix(NA, nrow = 1, ncol = 1)
    }
  }

  check_fail_zero = FALSE
  no_variables = ncol(x)
  if(ttest == TRUE) {
    no_categories = matrix(0,
                           nrow = no_variables,
                           ncol = 2)
  } else {
    no_categories = matrix(0,
                           nrow = no_variables,
                           ncol = max(group))
  }

  for(node in 1:no_variables) {
    unq_vls_x = sort(unique(x[,  node]))
    mx_vl_x = length(unq_vls_x)

    # Check if observed responses are not all unique ---------------------------
    if(mx_vl_x == nrow(x))
      stop(paste0("Only unique responses observed for variable ",
                  node,
                  " in the matrix x (group 1). We expect >= 1 observations per category."))

    if(ttest == TRUE) {
      unq_vls_y = sort(unique(y[,  node]))
      mx_vl_y = length(unq_vls_y)

      if(mx_vl_y == nrow(y))
        stop(paste0("Only unique responses observed for variable ",
                    node,
                    " in the matrix y (group 1). We expect >= 1 observations per category."))

      if(main_difference_model != "Free" & variable_bool[node] == TRUE){
        if(sum(is.na(match(unq_vls_x, unq_vls_y))) == mx_vl_x)
          stop(paste0("There is no overlap in response categories between the two groups for variable \n",
                      node,
                      ". As a result, bgmCompare cannot estimate a difference between the groups for any \n",
                      "of its category thresholds. You can run bgmCompare with the option \n",
                      "main_difference_model = ``Free''."))
      }

      unq_vls = sort(unique(c(unq_vls_x, unq_vls_y)))
      mx_vls = max(c(mx_vl_x,mx_vl_y))
    } else {
      unq_vls = sort(unique(unq_vls_x))
      mx_vls = mx_vl_x
    }

    # Recode data --------------------------------------------------------------
    if(variable_bool[node]) {#Regular ordinal variable
      # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
      if(main_difference_model == "Collapse") {

        if(ttest == TRUE) {
          zx = x[, node]
          zy = y[, node]

          cntr = -1
          for(value in unq_vls) {
            #Collapse categories for one group when not observed in the other.
            if(any(zx == value) & any(zy == value))
              cntr = cntr + 1
            x[zx == value, node] = max(0, cntr)
            y[zy == value, node] = max(0, cntr)
          }
        } else {
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
            if(!any(observed_scores[which(unq_vls == value), ] == 0))
              cntr = cntr + 1
            x[xx == value, node] = max(0, cntr)
          }
        }

      } else {
        z = x[, node]
        cntr = 0
        for(value in unq_vls) {
          x[z == value, node] = cntr
          cntr = cntr + 1
        }

        if(ttest == TRUE) {
          z = y[, node]
          cntr = 0
          for(value in unq_vls) {
            y[z == value, node] = cntr
            cntr = cntr + 1
          }
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
        if(ttest == TRUE) {
          y[, node] = as.integer(y[, node])
        }
      }

      mi = min(x[,node])
      if(ttest == TRUE)
        mi = min(c(mi, y[, node]))

      ma = max(x[,node])
      if(ttest == TRUE)
        ma = max(c(ma, y[, node]))

      if(reference_category[node] < mi | reference_category[node] > ma)
        stop(paste0(
          "The reference category for the Blume-Capel variable ", node, "is outside its \n",
          "range of observations in the matrices x (and y)."))

      # Check if observations start at zero and recode otherwise ---------------
      if(mi != 0) {
        reference_category[node] = reference_category[node] - mi
        x[, node] = x[, node] - mi
        if(ttest == TRUE)
          y[, node] = y[, node] - mi

        if(check_fail_zero == FALSE) {
          check_fail_zero = TRUE
          failed_zeroes = c(node)
        } else {
          failed_zeroes = c(failed_zeroes, node)
        }
      }

      if(ttest == TRUE) {
        check_range = length(unique(c(x[, node], y[, node])))
      } else {
        check_range = length(unique(x[, node]))
      }

      if(check_range < 3)
        stop(paste0("The Blume-Capel is only available for variables with more than one category \n",
                    "observed. There are two or less categories observed for variable ",
                    node,
                    "."))
    }


    # Warn that maximum category value is large --------------------------------
    if(main_difference_model == "Free") {
      if(ttest == TRUE) {
        no_categories[node, 1] = max(x[, node])
        no_categories[node, 2] = max(y[, node])
      } else {
        no_categories[node, ] = sapply(1:max(group), function(g) {
          max(x[group == g, node])
        })
      }
      if(!variable_bool[node] & max(no_categories[node, ]) > 10) {
        # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
        warning(paste0("In the (pseudo) likelihood of Blume-Capel variables, the normalization constant \n",
                       "is a sum over all possible values of the ordinal variable. The range of \n",
                       "observed values, possibly after recoding to integers, is assumed to be the \n",
                       "number of possible response categories.  For node ", node,", this range was \n",
                       "equal to ", max(no_categories[node,]), ", which may cause the analysis to take some \n",
                       "time to run. Note that for the Blume-Capel model, the bgm function does not \n",
                       "collapse the categories that have no observations between zero and the last \n",
                       "category. This may explain the large discrepancy between the first and last \n",
                       "category values."))
      }

      if(any(no_categories[node, ] == 0))
        stop(paste0("Only one value was observed for variable ",
                    node,
                    ", in at least one of the groups."))
    } else {
      if(ttest == TRUE) {
        no_categories[node, 1] = max(c(x[, node], y[, node]))
        no_categories[node, 2] = max(c(x[, node], y[, node]))
      } else {
        no_categories[node, ] = max(x[, node])
      }

      if(!variable_bool[node] & max(no_categories[node, ]) > 10) {
        # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
        warning(paste0("In the (pseudo) likelihood of Blume-Capel variables, the normalization constant \n",
                       "is a sum over all possible values of the ordinal variable. The range of \n",
                       "observed values, possibly after recoding to integers, is assumed to be the \n",
                       "number of possible response categories.  For node ", node,", in group 1, this \n",
                       "range was equal to ", max(no_categories[node,]), "which may cause the analysis to take some \n",
                       "time to run. Note that for the Blume-Capel model, the bgm function does not \n",
                       "collapse the categories that have no observations between zero and the last \n",
                       "category. This may explain the large discrepancy between the first and last \n",
                       "category values."))
      }


      # Check to see if not all responses are in one category --------------------
      if(any(no_categories[node, ] == 0))
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
              y = y,
              group = group,
              no_categories = no_categories,
              reference_category = reference_category,
              missing_index = missing_index,
              missing_index_gr1 = missing_index_gr1,
              missing_index_gr2 = missing_index_gr2,
              na_impute = na_impute))
}

# Helper function for computing `n_cat_obs`
compute_n_cat_obs_dev = function(x, no_categories, group = NULL) {
  if (is.null(group)) {  # Two-group design
    n_cat_obs = matrix(0, nrow = max(no_categories) + 1, ncol = ncol(x))
    for (variable in seq_len(ncol(x))) {
      for (category in seq_len(no_categories[variable])) {
        n_cat_obs[category + 1, variable] = sum(x[, variable] == category)
      }
    }
    return(n_cat_obs)
  } else {  # Multi-group design
    n_cat_obs = list()
    for (g in unique(group)) {
      n_cat_obs_gr = matrix(0, nrow = max(no_categories[, g]) + 1, ncol = ncol(x))
      for (variable in seq_len(ncol(x))) {
        for (category in seq_len(no_categories[variable, g])) {
          n_cat_obs_gr[category + 1, variable] = sum(x[group == g, variable] == category)
        }
      }
      n_cat_obs[[g]] = n_cat_obs_gr
    }
    return(n_cat_obs)
  }
}

# Helper function for computing sufficient statistics for Blume-Capel variables
compute_sufficient_blume_capel_dev = function(x, reference_category, ordinal_variable, group = NULL) {
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
prepare_output_bgmCompare_old = function(out, x, t.test, independent_thresholds,
                                         no_variables, no_categories, group, iter,
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
                                         projection,
                                         is_ordinal_variable) {

  save = any(c(save_options$save_main, save_options$save_pairwise, save_options$save_indicator))
  arguments = list(
    no_variables = no_variables, group = group, no_cases = tabulate(group),
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
  names_bycol = matrix(rep(data_columnnames, each = no_variables), ncol = no_variables)
  names_byrow = matrix(rep(data_columnnames, each = no_variables), ncol = no_variables, byrow = TRUE)
  names_comb = matrix(paste0(names_byrow, "-", names_bycol), ncol = no_variables)
  names_vec = names_comb[lower.tri(names_comb)]
  names_vec_t = names_comb[lower.tri(names_comb, diag = TRUE)]

  # Prepare output elements
  results = list()

  # Check output format and process accordingly
  if ("pairwise_difference_indicator" %in% names(out)) {
    no_groups = 2
    projection = matrix(c(-0.5, 0.5), nrow = 2, ncol = 1)

    # Posterior mean pairwise effects
    if(!save) {
      int = out$interactions
      int = int[lower.tri(int)]
      pwd = out$pairwise_difference
      pwd = pwd [lower.tri(pwd)]
    } else {
      int = colMeans(out$interactions)
      pwd = colMeans(out$pairwise_difference)
    }

    posterior_mean_pairwise = matrix(0, nrow = length(int), ncol = 3)
    posterior_mean_pairwise[, 1] = int
    for (row in 1:length(pwd)) {
      posterior_mean_pairwise[row, 2] = projection[1] *  pwd[row]
      posterior_mean_pairwise[row, 3] = projection[2] *  pwd[row]
    }

    results$posterior_mean_pairwise = posterior_mean_pairwise
    colnames(results$posterior_mean_pairwise) = c("overall", "group_1", "group_2")
    rownames(results$posterior_mean_pairwise) = names_vec

    # Posterior mean main effects
    if(independent_thresholds) {
      if(!save) {
        tmpt1 = out$thresholds_gr1
        tmpt2 = out$thresholds_gr2

        n_cats = apply(no_categories, 1, max)
        n_cats[!is_ordinal_variable] = 2
        th1 = vector(length=sum(n_cats))
        th2 = vector(length=sum(n_cats))
        tel2 = tel1 = 0
        for(v in 1:no_variables) {
          if(is_ordinal_variable[v]) {
            for(c in 1:no_categories[v, 1]) {
              tel1 = tel1 + 1
              th1 [tel1] = tmpt1[v, c]
            }
            for(c in 1:no_categories[v, 2]) {
              tel2 = tel2 + 1
              th2 [tel2] = tmpt2[v, c]
            }
          } else {
            for(c in 1:2) {
              tel1 = tel1 + 1
              th1 [tel1] = tmpt1[v, c]
            }
            for(c in 1:2) {
              tel2 = tel2 + 1
              th2 [tel2] = tmpt2[v, c]
            }
          }
          tel2 = tel1 = max(c(tel1, tel2))
        }
      } else {
        th1 = colMeans(out$thresholds_gr1)
        th2 = colMeans(out$thresholds_gr2)
      }

      counter = 0
      main_names = vector(length = length(th1))
      for(v in 1:no_variables) {
        if(is_ordinal_variable[v]) {
          for(c in 1:max(no_categories[v, ])) {
            counter = counter + 1
            main_names[counter] = paste0(data_columnnames[v], "(", c, ")")
          }
        } else {
          counter = counter + 1
          main_names[counter] = paste0(data_columnnames[v], "(linear)")
          counter = counter + 1
          main_names[counter] = paste0(data_columnnames[v], "(quadratic)")
        }
      }

      posterior_mean_main = matrix(0, nrow = length(th1), ncol = 2)
      posterior_mean_main[, 1] = th1
      posterior_mean_main[, 2] = th2

      results$posterior_mean_main = posterior_mean_main
      colnames(results$posterior_mean_main) = c("group_1", "group_2")
      rownames(results$posterior_mean_main) = main_names

    } else {
      if(!save) {
        tmpt = out$thresholds
        tmpm = out$main_difference

        n_cats = apply(no_categories, 1, max)
        n_cats[!is_ordinal_variable] = 2
        th = vector(length=sum(n_cats))
        md = vector(length=sum(n_cats))
        counter = 0
        for(v in 1:no_variables) {
          if(is_ordinal_variable[v]) {
            for(c in 1:max(no_categories[v, ])) {
              counter = counter + 1
              th [counter] = tmpt[v, c]
              md [counter] = tmpm[v, c]
            }
          } else {
            for(c in 1:2) {
              counter = counter + 1
              th [counter] = tmpt[v, c]
              md [counter] = tmpm[v, c]
            }
          }
        }
      } else {
        th = colMeans(out$thresholds)
        md = colMeans(out$main_difference)
      }

      counter = 0
      main_names = vector(length = length(th))
      for(v in 1:no_variables) {
        if(is_ordinal_variable[v]) {
          for(c in 1:max(no_categories[v, ])) {
            counter = counter + 1
            main_names[counter] = paste0(data_columnnames[v], "(", c, ")")
          }
        } else {
          counter = counter + 1
          main_names[counter] = paste0(data_columnnames[v], "(linear)")
          counter = counter + 1
          main_names[counter] = paste0(data_columnnames[v], "(quadratic)")
        }
      }

      posterior_mean_main = matrix(0, nrow = length(th), ncol = 3)
      posterior_mean_main[, 1] = th
      for (row in 1:length(th)) {
        posterior_mean_main[row, 2] = projection[1] *  md[row]
        posterior_mean_main[row, 3] = projection[2] *  md[row]
      }

      results$posterior_mean_main = posterior_mean_main
      colnames(results$posterior_mean_main) = c("overall", "group_1", "group_2")
      rownames(results$posterior_mean_main) = main_names
    }

    if(difference_selection) {
      # Posterior mean indicators
      if(!save) {
        ind = out$pairwise_difference_indicator
        diag(ind) = 1.0
        if(!independent_thresholds)
          diag(ind) = out$main_difference_indicator
      } else {
        tmp = colMeans(out$pairwise_difference_indicator)
        ind = matrix(0, no_variables, no_variables)
        ind[lower.tri(ind)] = tmp
        ind = ind + t(ind)
        if(!independent_thresholds) {
          tmp = colMeans(out$main_difference_indicator)
          diag(ind) = tmp
        } else {
          diag(ind) = 1.0
        }
      }

      results$posterior_mean_indicator = ind
      colnames(results$posterior_mean_indicator) = data_columnnames
      rownames(results$posterior_mean_indicator) = data_columnnames
    }


    if(save) {
      # Raw samples main effects
      if (save) {
        if (independent_thresholds) {
          # Handle thresholds for independent thresholds case
          num_parameters_gr1 = num_parameters_gr2 = 0
          for(v in 1:no_variables) {
            if(is_ordinal_variable[v]) {
              num_parameters_gr1 = num_parameters_gr1 + no_categories[v, 1];
              num_parameters_gr2 = num_parameters_gr2 + no_categories[v, 2];
            } else {
              num_parameters_gr1 = num_parameters_gr1 + 2;
              num_parameters_gr2 = num_parameters_gr2 + 2;
            }
          }

          main_effect_samples = matrix(NA, nrow = nrow(out$thresholds_gr1),
                                       ncol = num_parameters_gr1 + num_parameters_gr2)

          threshold1_counter = 0  # Separate counter for thresholds
          threshold2_counter = 0  # Separate counter for thresholds
          col_counter = 0        # Overall column counter for main_effect_samples
          col_names = character()  # Vector to store column names

          for (var in 1:no_variables) {
            if(is_ordinal_variable[var]) {
              for (cat in 1:max(no_categories[var, ])) {
                # Handle thresholds group 1
                if(cat <= no_categories[var, 1]) {
                  col_counter = col_counter + 1
                  threshold1_counter = threshold1_counter + 1
                  main_effect_samples[, col_counter] = out$thresholds_gr1[, threshold1_counter]
                  col_names = c(col_names, paste0(data_columnnames[var],"_gr1(", cat, ")"))
                }
                # Handle thresholds group 2
                if(cat <= no_categories[var, 2]) {
                  col_counter = col_counter + 1
                  threshold2_counter = threshold2_counter + 1
                  main_effect_samples[, col_counter] = out$thresholds_gr2[, threshold2_counter]
                  col_names = c(col_names, paste0(data_columnnames[var],"_gr2(", cat, ")"))
                }
              }
            } else {
              # Handle thresholds group 1
              col_counter = col_counter + 1
              threshold1_counter = threshold1_counter + 1
              main_effect_samples[, col_counter] = out$thresholds_gr1[, threshold1_counter]
              col_names = c(col_names, paste0(data_columnnames[var],"_gr1(linear)"))

              # Handle thresholds group 2
              col_counter = col_counter + 1
              threshold2_counter = threshold2_counter + 1
              main_effect_samples[, col_counter] = out$thresholds_gr2[, threshold2_counter]
              col_names = c(col_names, paste0(data_columnnames[var],"_gr2(linear)"))

              # Handle thresholds group 1
              col_counter = col_counter + 1
              threshold1_counter = threshold1_counter + 1
              main_effect_samples[, col_counter] = out$thresholds_gr1[, threshold1_counter]
              col_names = c(col_names, paste0(data_columnnames[var],"_gr1(quadratic)"))

              # Handle thresholds group 2
              col_counter = col_counter + 1
              threshold2_counter = threshold2_counter + 1
              main_effect_samples[, col_counter] = out$thresholds_gr2[, threshold2_counter]
              col_names = c(col_names, paste0(data_columnnames[var],"_gr2(quadratic)"))
            }
          }

          colnames(main_effect_samples) = col_names
          rownames(main_effect_samples) = paste0("Iter:", seq_len(nrow(main_effect_samples)))
        } else {
          # Handle thresholds and main differences for non-independent thresholds
          total_categories = sum(apply(no_categories, 1, max))
          main_effect_samples = matrix(NA, nrow = nrow(out$thresholds),
                                       ncol = ncol(out$thresholds) + ncol(out$main_difference)) # Thresholds + main differences

          threshold_counter = 0  # Separate counter for thresholds
          main_diff_counter = 0  # Separate counter for main differences
          col_counter = 0        # Overall column counter for main_effect_samples
          col_names = character()  # Vector to store column names

          for (var in 1:no_variables) {
            if(is_ordinal_variable[var]) {
              for (cat in 1:max(no_categories[var, ])) {
                # Increment overall column counter
                col_counter = col_counter + 1

                # Handle thresholds
                threshold_counter = threshold_counter + 1
                main_effect_samples[, col_counter] = out$thresholds[, threshold_counter]
                col_names = c(col_names, paste0(data_columnnames[var],"(", cat, ")"))

                col_counter = col_counter + 1
                main_diff_counter = main_diff_counter + 1
                main_effect_samples[, col_counter] = out$main_difference[, main_diff_counter]
                col_names = c(col_names, paste0("main_difference_", data_columnnames[var], "(", cat, ")"))
              }
            } else {
              # Increment overall column counter
              col_counter = col_counter + 1

              # Handle thresholds
              threshold_counter = threshold_counter + 1
              main_effect_samples[, col_counter] = out$thresholds[, threshold_counter]
              col_names = c(col_names, paste0(data_columnnames[var],"(linear)"))

              col_counter = col_counter + 1
              main_diff_counter = main_diff_counter + 1
              main_effect_samples[, col_counter] = out$main_difference[, main_diff_counter]
              col_names = c(col_names, paste0("main_difference_", data_columnnames[var], "(linear)"))

              # Increment overall column counter
              col_counter = col_counter + 1

              # Handle thresholds
              threshold_counter = threshold_counter + 1
              main_effect_samples[, col_counter] = out$thresholds[, threshold_counter]
              col_names = c(col_names, paste0(data_columnnames[var],"(quadratic)"))

              col_counter = col_counter + 1
              main_diff_counter = main_diff_counter + 1
              main_effect_samples[, col_counter] = out$main_difference[, main_diff_counter]
              col_names = c(col_names, paste0("main_difference_", data_columnnames[var], "(quadratic)"))

            }

          }

          colnames(main_effect_samples) = col_names
          rownames(main_effect_samples) = paste0("Iter:", seq_len(nrow(main_effect_samples)))
        }

        results$main_effect_samples = main_effect_samples
      }

      # Raw samples pairwise effects
      pairwise_effect_samples = cbind(out$interactions, out$pairwise_difference)
      pairwise_names = unlist(lapply(1:(no_variables - 1), function(v1) {
        sapply((v1 + 1):no_variables, function(v2) {
          c(
            paste0(data_columnnames[v1], "-", data_columnnames[v2]),
            paste0("pairwise_difference(", data_columnnames[v1], "-", data_columnnames[v2], ")")
          )
        })
      }))
      colnames(pairwise_effect_samples) = pairwise_names
      rownames(pairwise_effect_samples) = paste0("Iter:", seq_len(nrow(pairwise_effect_samples)))

      results$pairwise_effect_samples = pairwise_effect_samples

      # Raw samples indicators
      if(difference_selection) {
        if(independent_thresholds) {
          inclusion_indicator_samples = matrix(0,
                                               nrow = nrow(out$pairwise_difference_indicator),
                                               ncol = ncol(out$pairwise_difference_indicator))
        } else {
          inclusion_indicator_samples = matrix(0,
                                               nrow = nrow(out$pairwise_difference_indicator),
                                               ncol = ncol(out$pairwise_difference_indicator) +
                                                 ncol(out$main_difference_indicator))
        }

        for(row in 1:nrow(out$pairwise_difference_indicator)) {
          tmp.p = out$pairwise_difference_indicator[row, ]
          ind = matrix(0, no_variables, no_variables)
          ind[lower.tri(ind)] = tmp.p
          if(!independent_thresholds){
            tmp.t = out$main_difference_indicator[row, ]
            diag(ind) = tmp.t
          }

          if(!independent_thresholds) {
            inclusion_indicator_samples[row, ] = ind[lower.tri(ind, diag = TRUE)]
          } else {
            inclusion_indicator_samples[row, ] = ind[lower.tri(ind)]
          }
        }

        if(independent_thresholds) {
          colnames(inclusion_indicator_samples) = names_vec
        } else {
          colnames(inclusion_indicator_samples) = names_vec_t
        }
        rownames(inclusion_indicator_samples) = paste0("Iter:", 1:nrow(inclusion_indicator_samples))

        results$inclusion_indicator_samples = inclusion_indicator_samples
      }
    }
  }

  # Add arguments to the output
  results$arguments = arguments
  class(results) = c("bgmCompare")
  return(results)
}