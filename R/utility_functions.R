#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @importFrom methods hasArg
#' @importFrom stats qnorm uniroot

check_model = function(x,
                       variable_type,
                       reference_category,
                       interaction_scale = 2.5,
                       threshold_alpha = 0.5,
                       threshold_beta = 0.5,
                       edge_selection = TRUE,
                       edge_prior = c("Bernoulli", "Beta-Bernoulli"),
                       inclusion_probability = 0.5,
                       beta_bernoulli_alpha = 1,
                       beta_bernoulli_beta = 1) {

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
                     "Since na.action = listwise these rows were excluded \n",
                     "from the analysis."),
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
    na.impute = FALSE
  } else {
    # Check for missing values -------------------------------------------------
    no_missings = sum(is.na(x))
    no_persons = nrow(x)
    no_nodes = ncol(x)
    if(no_missings > 0) {
      missing_index = matrix(0, nrow = no_missings, ncol = 2)
      na.impute = TRUE
      cntr = 0
      for(node in 1:no_nodes) {
        mis = which(is.na(x[, node]))
        if(length(mis) > 0) {
          for(i in 1:length(mis)) {
            cntr = cntr + 1
            missing_index[cntr, 1] = mis[i]
            missing_index[cntr, 2] = node
            x[mis[i], node] = stats::median(x[, node], #start value for imputation
                                            na.rm = TRUE)
            #This is non-zero if no zeroes are observed (we then collapse over zero below)
          }
        }
      }
    } else {
      missing_index = matrix(NA, nrow = 1, ncol = 1)
      na.impute = FALSE
    }
  }

  no_nodes = ncol(x)
  no_categories = vector(length = no_nodes)
  for(node in 1:no_nodes) {
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
      if(min(x[, node]) < 0) {
        reference_category[node] = reference_category[node] - min(x[, node])
        x[, node] = x[, node] - min(x[, node])

        warning(paste0("The bgm function assumes that the observed ordinal variables are integers and \n",
                       "that the lowest observed category score is zero. The lowest score for node \n",
                       node, " was recoded to zero for the analysis. Note that bgm also recoded the \n",
                       "the corresponding reference category score to ", reference_category[node], "."))
      }
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

  return(list(x = x,
              no_categories = no_categories,
              reference_category = reference_category,
              missing_index = missing_index,
              na.impute = na.impute))
}