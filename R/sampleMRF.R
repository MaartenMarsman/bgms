#' Sample states of the ordinal MRF
#'
#' @description
#' This function samples states from the ordinal MRF for a one or two-samples
#' design using a Gibbs sampler. The Gibbs sampler is initiated with random
#' values from the response options, after which it proceeds by simulating
#' states for each variable from a logistic model using the other variable
#' states as predictor variables.
#'
#' @details
#' There are two modeling options for the category thresholds. The default
#' option assumes that the category thresholds are free, except that the first
#' threshold is set to zero for identification. The user then only needs to
#' specify the thresholds for the remaining response categories. This option is
#' useful for any type of ordinal variable and gives the user the most freedom
#' in specifying their model.
#'
#' The Blume-Capel option is specifically designed for ordinal variables that
#' have a special type of reference_category category, such as the neutral
#' category in a Likert scale. The Blume-Capel model specifies the following
#' quadratic model for the threshold parameters:
#' \deqn{\mu_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} - \text{r})^2,}{{\mu_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} - \text{r})^2,}}
#' where \eqn{\mu_{\text{c}}}{\mu_{\text{c}}} is the threshold for category c
#' (which now includes zero), \eqn{\alpha}{\alpha} offers a linear trend
#' across categories (increasing threshold values if
#' \eqn{\alpha > 0}{\alpha > 0} and decreasing threshold values if
#' \eqn{\alpha <0}{\alpha <0}), if \eqn{\beta < 0}{\beta < 0}, it offers an
#' increasing penalty for responding in a category further away from the
#' reference_category category r, while \eqn{\beta > 0}{\beta > 0} suggests a
#' preference for responding in the reference_category category.
#'
#' In the two-samples design, the pairwise interactions between the variables \eqn{i}{i}
#' and \eqn{j}{j} are modeled in the first sample as
#' \deqn{\sigma_{\text{ij}} = \theta_{\text{ij}} + \delta_{\text{ij}} / 2,}{\sigma_{\text{ij}} = \theta_{\text{ij}} + \delta_{\text{ij}} / 2,}
#' and in the second samples as
#' \deqn{\sigma_{\text{ij}} = \theta_{\text{ij}} - \delta_{\text{ij}} / 2,}{\sigma_{\text{ij}} = \theta_{\text{ij}} - \delta_{\text{ij}} / 2.}
#' The pairwise interaction parameter \eqn{\theta_{\text{ij}}}{\theta_{\text{ij}}}
#' denotes an overall effect that is considered nuisance, and attention is focused
#' on the pairwise difference parameter \eqn{\delta_{\text{ij}}}{\delta_{\text{ij}}},
#' which reflects the difference in the pairwise interaction between the two groups.
#'
#' In a paired samples design, the pairwise interactions between the samples must
#' be modeled to account for the dependence in the repeated measures. This
#' dependence can also be modeled by assigning a random effect per case to each
#' variable in the model. The matrix of between-sample pairwise interactions can
#' then be viewed as the covariance matrix of the random effects.
#'
#' @param no_states The number of states of the ordinal MRF to be generated.
#'
#' @param no_variables The number of variables in the ordinal MRF.
#'
#' @param no_categories Either a positive integer or a vector of positive
#' integers of length \code{no_variables}. The number of response categories on top
#' of the base category: \code{no_categories = 1} generates binary states.
#'
#' @param interactions A symmetric \code{no_variables} by \code{no_variables} matrix of
#' pairwise interactions. Only its off-diagonal elements are used.
#'
#' @param thresholds A \code{no_variables} by \code{max(no_categories)} matrix of
#' category thresholds. The elements in row \code{i} indicate the thresholds of
#' variable \code{i}. If \code{no_categories} is a vector, only the first
#' \code{no_categories[i]} elements are used in row \code{i}. If the Blume-Capel
#' model is used for the category thresholds for variable \code{i}, then row
#' \code{i} requires two values (details below); the first is
#' \eqn{\alpha}{\alpha}, the linear contribution of the Blume-Capel model and
#' the second is \eqn{\beta}{\beta}, the quadratic contribution.
#'
#' @param variable_type What kind of variables are simulated? Can be a single
#' character string specifying the variable type of all \code{p} variables at
#' once or a vector of character strings of length \code{p} specifying the type
#' for each variable separately. Currently, bgm supports ``ordinal'' and
#' ``blume-capel''. Binary variables are automatically treated as ``ordinal’’.
#' Defaults to \code{variable_type = "ordinal"}.
#'
#' @param reference_category An integer vector of length \code{no_variables} specifying the
#' reference_category category that is used for the Blume-Capel model (details below).
#' Can be any integer value between \code{0} and \code{no_categories} (or
#' \code{no_categories[i]}).
#'
#' @param iter The number of iterations used by the Gibbs sampler.
#' The function provides the last state of the Gibbs sampler as output. By
#' default set to \code{1e3}.
#'
#' @param compare Logical, if \code{TRUE} models the case-specific dependence using
#' a paired-samples design; if \code{FALSE} treats the groups as independent.
#' Default is \code{FALSE}.
#' compare = FALSE,
#'
#' @param no_states_gr1 The number of cases simulated for the first sample.
#'
#' @param no_states_gr2 The number of cases simulated for the second sample. In
#' the paired samples scenario this number is set equal to \code{no_states_gr1}.
#'
#' @param main_difference A \code{no_variables} by \code{max(no_categories)}
#' matrix containing the difference in category thresholds between the two
#' samples.
#'
#' @param pairwise_difference A \code{no_variables} by \code{no_variables}
#' matrix containing the difference in pairwise interactions between the two
#' samples.
#'
#' @param paired Logical, if \code{TRUE} the data come from a paired-samples
#' design; if \code{FALSE} the data come from two independent samples. Default is
#' \code{FALSE}.
#'
#' @param cross_lagged A \code{no_variables} by \code{no_variables} matrix with
#' pairwise interactions between the two samples.
#'
#' @return If \code{paired = FALSE}, a \code{no_states} by \code{no_variables}
#' matrix of simulated states of the ordinal MRF. If \code{paired = TRUE}, a
#' list with a matrix of simulated states for each of the two samples.
#'
#' @examples
#' # Generate responses from a network of five binary and ordinal variables.
#' no_variables = 5
#' no_categories = sample(1:5, size = no_variables, replace = TRUE)
#'
#' Interactions = matrix(0, nrow = no_variables, ncol = no_variables)
#' Interactions[2, 1] = Interactions[4, 1] = Interactions[3, 2] =
#'   Interactions[5, 2] = Interactions[5, 4] = .25
#' Interactions = Interactions + t(Interactions)
#' Thresholds = matrix(0, nrow = no_variables, ncol = max(no_categories))
#'
#' x = mrfSampler(no_states = 1e3,
#'                no_variables = no_variables,
#'                no_categories = no_categories,
#'                interactions = Interactions,
#'                thresholds = Thresholds)
#'
#' # Generate responses from a network of 2 ordinal and 3 Blume-Capel variables.
#' no_variables = 5
#' no_categories = 4
#'
#' Interactions = matrix(0, nrow = no_variables, ncol = no_variables)
#' Interactions[2, 1] = Interactions[4, 1] = Interactions[3, 2] =
#'   Interactions[5, 2] = Interactions[5, 4] = .25
#' Interactions = Interactions + t(Interactions)
#'
#' Thresholds = matrix(NA, no_variables, no_categories)
#' Thresholds[, 1] = -1
#' Thresholds[, 2] = -1
#' Thresholds[3, ] = sort(-abs(rnorm(4)), decreasing = TRUE)
#' Thresholds[5, ] = sort(-abs(rnorm(4)), decreasing = TRUE)
#'
#' x = mrfSampler(no_states = 1e3,
#'                no_variables = no_variables,
#'                no_categories = no_categories,
#'                interactions = Interactions,
#'                thresholds = Thresholds,
#'                variable_type = c("b","b","o","b","o"),
#'                reference_category = 2)
#' @export
mrfSampler = function(no_states,
                      no_variables,
                      no_categories,
                      interactions,
                      thresholds,
                      variable_type = "ordinal",
                      reference_category,
                      iter = 1e3,
                      compare = FALSE,
                      no_states_gr1,
                      no_states_gr2,
                      main_difference,
                      pairwise_difference,
                      paired = FALSE,
                      cross_lagged) {

  # Check if compare and paired are logical variables --------------------------
  compare = as.logical(compare)
  if(is.na(compare))
    stop("The parameter compare needs to be TRUE or FALSE.")
  paired = as.logical(paired)
  if(is.na(paired))
    stop("The parameter paired needs to be TRUE or FALSE.")


  # Check no_states, no_variables, iter ----------------------------------------
  if(compare == FALSE) {
    if(no_states <= 0 ||
       abs(no_states - round(no_states)) > .Machine$double.eps)
      stop("The parameter no_states needs be a positive integer.")
  } else {
    if(no_states_gr1 <= 0 ||
       abs(no_states_gr1 - round(no_states_gr1)) > .Machine$double.eps)
      stop("The parameter no_states_gr1 needs be a positive integer.")

    if(paired == TRUE) {
      no_states_gr2 = no_states_gr1
    } else {
      if(no_states_gr2 <= 0 ||
         abs(no_states_gr2 - round(no_states_gr2)) > .Machine$double.eps)
        stop("The parameter no_states_gr1 needs be a positive integer.")
    }
  }
  if(no_variables <= 0 ||
     abs(no_variables - round(no_variables)) > .Machine$double.eps)
    stop("``no_variables'' needs be a positive integer.")
  if(iter <= 0 ||
     abs(iter - round(iter)) > .Machine$double.eps)
    stop("``iter'' needs be a positive integer.")

  # Check no_categories --------------------------------------------------------
  if(length(no_categories) == 1) {
    if(no_categories <= 0 ||
       abs(no_categories - round(no_categories)) > .Machine$double.eps)
      stop("``no_categories'' needs be a (vector of) positive integer(s).")
    no_categories = rep(no_categories, no_variables)
  } else {
    for(variable in 1:no_variables) {
      if(no_categories[variable] <= 0 ||
         abs(no_categories[variable] - round(no_categories[variable])) >
         .Machine$double.eps)
        stop(paste("For variable", variable, "``no_categories'' was not a positive integer."))
    }
  }

  # Check variable specification -----------------------------------------------
  if(length(variable_type) == 1) {
    variable_type = match.arg(arg = variable_type,
                              choices = c("ordinal", "blume-capel"))

    if(variable_type == "blume-capel" && any(no_categories < 2)) {
      stop(paste0("The Blume-Capel model only works for ordinal variables with more than two \n",
                  "response options. But variables ", which(no_categories < 2), " are binary variables."))
    }
    variable_type = rep(variable_type, no_variables)

  } else {
    if(length(variable_type) != no_variables) {
      stop(paste0("The argument ``variable_type'' should be either a single character string or a \n",
                  "vector of character strings of length ``no_variables''."))
    } else {
      for(variable in 1:no_variables)
        variable_type[variable] = match.arg(arg = variable_type[variable],
                                            choices = c("ordinal", "blume-capel"))
      if(any(variable_type == "blume-capel" & no_categories < 2)) {
        stop(paste0("The Blume-Capel model only works for ordinal variables with more than two \n",
                    "response options. But variables ",
                    which(variable_type == "blume-capel" & no_categories < 2),
                    " are binary variables."))
      }
    }
  }

  # Check the reference_category for Blume-Capel variables ---------------------
  if(any(variable_type == "blume-capel")) {
    if(length(reference_category) == 1) {
      reference_category = rep(reference_category, no_variables)
    }
    if(any(reference_category < 0) || any(abs(reference_category - round(reference_category)) > .Machine$double.eps)) {
      stop(paste0("For variables ",
                  which(reference_category < 0),
                  " ``reference_category'' was either negative or not integer."))
    }
    if(any(reference_category - no_categories > 0)) {
      stop(paste0("For variables ",
                  which(reference_category - no_categories > 0),
                  " the ``reference_category'' category was larger than the maximum category value."))
    }
  }

  # Check interactions ---------------------------------------------------------
  if(!inherits(interactions, what = "matrix"))
    interactions = as.matrix(interactions)
  if(!isSymmetric(interactions))
    stop("The matrix ``interactions'' needs to be symmetric.")
  if(nrow(interactions) != no_variables)
    stop("The matrix ``interactions'' needs to have ``no_variables'' rows and columns.")


  # Check pairwise_difference and cross_lagged ---------------------------------
  if(compare == TRUE) {
    if(!inherits(pairwise_difference, what = "matrix"))
      pairwise_difference = as.matrix(pairwise_difference)
    if(!isSymmetric(pairwise_difference))
      stop("The matrix pairwise_difference needs to be symmetric.")
    if(nrow(pairwise_difference) != no_variables)
      stop("The matrix pairwise_difference needs to have no_variables rows and columns.")

    if(paired == TRUE) {
      if(!inherits(cross_lagged, what = "matrix"))
        cross_lagged = as.matrix(cross_lagged)
      if(!isSymmetric(cross_lagged))
        stop("The matrix cross_lagged needs to be symmetric.")
      if(nrow(cross_lagged) != no_variables)
        stop("The matrix cross_lagged needs to have no_variables rows and columns.")
    }
  }

  # Check the threshold values -------------------------------------------------
  if(!inherits(thresholds, what = "matrix")) {
    if(max(no_categories) == 1) {
      if(length(thresholds) == no_variables) {
        thresholds = matrix(thresholds, ncol = 1)
      } else {
        stop(paste0("The matrix ``thresholds'' has ",
                    length(thresholds),
                    " elements, but requires",
                    no_variables,
                    "."))
      }
    } else {
      stop("The input thresholds needs to be a matrix.")
    }
  }

  if(nrow(thresholds) != no_variables)
    stop("The matrix thresholds needs to be have no_variables rows.")

  if(compare == TRUE) {
    if(inherits(thresholds, what = "matrix")) {
      if(nrow(main_difference) != nrow(thresholds))
        stop('The matrix main_difference needs to have no_variable rows.')
      if(ncol(main_difference) != ncol(thresholds))
        stop('The matrix main_difference needs to have as many columns as the thresholds matrix.')
    } else {
      stop("The input main_difference needs to be a matrix.")
    }
  }

  for(variable in 1:no_variables) {
    if(variable_type[variable] != "blume-capel") {
      if(anyNA(thresholds[variable, 1:no_categories[variable]])) {
        tmp = which(is.na(thresholds[variable, 1:no_categories[variable]]))

        string = paste(tmp, sep = ",")

        stop(paste0("The input thresholds contains NA(s) for variable ",
                    variable,
                    " in category \n",
                    "(categories) ",
                    paste(which(is.na(thresholds[variable, 1:no_categories[variable]])), collapse = ", "),
                    ", where a numeric value is needed."))
        if(compare == TRUE) {
          tmp = which(is.na(main_difference[variable, 1:no_categories[variable]]))

          string = paste(tmp, sep = ",")

          stop(paste0("The input main_difference contains NA(s) for variable ",
                      variable,
                      " in category \n",
                      "(categories) ",
                      paste(which(is.na(main_difference[variable, 1:no_categories[variable]])), collapse = ", "),
                      ", where a numeric value is needed."))
        }
      }
      if(ncol(thresholds) > no_categories[variable]) {
        if(!anyNA(thresholds[variable, (no_categories[variable]+1):ncol(thresholds)])) {
          warning(paste0("The matrix ``thresholds'' contains numeric values for variable ",
                         variable,
                         " for category \n",
                         "(categories, i.e., columns) exceding the maximum of ",
                         no_categories[variable],
                         ". These values will \n",
                         "be ignored."))
        }
        if(compare == TRUE) {
          if(!anyNA(main_difference[variable, (no_categories[variable]+1):ncol(main_difference)])) {
            warning(paste0("The input main_difference contains numeric values for variable ",
                           variable,
                           " for category \n",
                           "(categories, i.e., columns) exceding the maximum of ",
                           no_categories[variable],
                           ". These values will \n",
                           "be ignored."))
          }
        }
      }
    } else {
      if(anyNA(thresholds[variable, 1:2])) {
        stop(paste0("The Blume-Capel model is chosen for the category thresholds of variable ",
                    variable,
                    ". \n",
                    "This model has two parameters that need to be placed in columns 1 and 2, row \n",
                    variable,
                    ", of the ``thresholds'' input matrix. Currently, there are NA(s) in these \n",
                    "entries, where a numeric value is needed."))
      }
      if(ncol(thresholds) > 2) {
        if(!anyNA(thresholds[variable, 3:ncol(thresholds)])) {
          warning(paste0("The Blume-Capel model is chosen for the category thresholds of variable ",
                         variable,
                         ". \n",
                         "This model has two parameters that need to be placed in columns 1 and 2, row \n",
                         variable,
                         ", of the ``thresholds'' input matrix. However, there are numeric values \n",
                         "in higher categories. These values will be ignored."))
        }
      }
      if(paired == TRUE) {
        if(anyNA(main_difference[variable, 1:2])) {
          stop(paste0("The Blume-Capel model is chosen for the category thresholds of variable ",
                      variable,
                      ". \n",
                      "This model has two parameters that need to be placed in columns 1 and 2, row \n",
                      variable,
                      ", of the thresholds input matrix. Since a two-sample model is simulated, the \n",
                      "differences in these thresholds between the two groups are modeled with the \n",
                      "main_difference matrix input, with the differences for the Blume-Capel model \n",
                      "also placed in columns 1 and 2 of the main_difference input matrix. Currently, \n",
                      "there are NA(s) in these entries, where a numeric value is needed."))
        }
        if(ncol(main_difference) > 2) {
          if(!anyNA(main_difference[variable, 3:ncol(thresholds)])) {
            warning(paste0("The Blume-Capel model is chosen for the category thresholds of variable ",
                           variable,
                           ". \n",
                           "This model has two parameters that need to be placed in columns 1 and 2, row \n",
                           variable,
                           ", of the thresholds input matrix. Since a two-sample model is simulated, the \n",
                           "differences in these thresholds between the two groups are modeled with the \n",
                           "main_difference matrix input, with the differences for the Blume-Capel model \n",
                           "also placed in columns 1 and 2 of the main_difference input matrix. However, \n",
                           "there are numeric values in higher categories. These values will be ignored."))
          }
        }
      }
    }
  }

  for(variable in 1:no_variables) {
    if(variable_type[variable] != "blume-capel") {
      for(category in 1:no_categories[variable]) {
        if(!is.finite(thresholds[variable, category]))
          stop(paste("The threshold parameter for variable", variable, "and category",
                     category, "is NA or not finite."))
        if(compare == TRUE) {
          if(!is.finite(main_difference[variable, category]))
            stop(paste("The difference in threshold parameters for variable", variable, "and category",
                       category, "is NA or not finite."))
        }
      }
    } else {
      if(!is.finite(thresholds[variable, 1]))
        stop(paste0(
          "The alpha parameter for the Blume-Capel model for variable ",
          variable,
          " is NA \n",
          " or not finite."))
      if(!is.finite(thresholds[variable, 2]))
        stop(paste0("The beta parameter for the Blume-Capel model for variable",
                    variable,
                    "is NA \n",
                    " or not finite."))
      if(compare == TRUE) {
        if(!is.finite(main_difference[variable, 1]))
          stop(paste0(
            "The between sample difference in the alpha parameter for the Blume-Capel model \n",
            "for variable ",
            variable,
            " is NA \n",
            " or not finite."))
        if(!is.finite(main_difference[variable, 2]))
          stop(paste0(
            "The between sample difference in the beta parameter for the Blume-Capel model \n",
            "for variable",
            variable,
            "is NA \n",
            " or not finite."))
      }
    }
  }

  # The Gibbs sampler ----------------------------------------------------------
  if(compare == FALSE) {
    if(!any(variable_type == "blume-capel")) {
      x <- sample_omrf_gibbs(no_states = no_states,
                             no_variables = no_variables,
                             no_categories = no_categories,
                             interactions = interactions,
                             thresholds = thresholds,
                             iter = iter)
    } else {
      x <- sample_bcomrf_gibbs(no_states = no_states,
                               no_variables = no_variables,
                               no_categories = no_categories,
                               interactions = interactions,
                               thresholds = thresholds,
                               variable_type = variable_type,
                               reference_category = reference_category,
                               iter = iter)
    }
    return(x)
  } else {
    if(paired == FALSE)
      cross_lagged = matrix(0, 1, 1)
    if(!any(variable_type == "blume-capel")) {
      out <- sample_twosample_omrf_gibbs(no_states_gr1 = no_states_gr1,
                                         no_states_gr2 = no_states_gr2,
                                         no_variables = no_variables,
                                         no_categories = no_categories,
                                         interactions = interactions,
                                         thresholds = thresholds,
                                         main_difference = main_difference,
                                         pairwise_difference = pairwise_difference,
                                         cross_lagged = cross_lagged,
                                         paired = paired,
                                         iter = iter)
    } else {
      out <- sample_twosample_bcomrf_gibbs(no_states_gr1 = no_states_gr1,
                                           no_states_gr2 = no_states_gr2,
                                           no_variables = no_variables,
                                           no_categories = no_categories,
                                           interactions = interactions,
                                           thresholds = thresholds,
                                           main_difference = main_difference,
                                           pairwise_difference = pairwise_difference,
                                           cross_lagged = cross_lagged,
                                           paired = paired,
                                           variable_type = variable_type,
                                           reference_category = reference_category,
                                           iter = iter)

    }
    return(list(x = out$x, y = out$y))
  }
}
