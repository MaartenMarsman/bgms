reformat_data = function(x,
                         na_action,
                         variable_bool,
                         reference_category) {
  if(na_action == "listwise") {
    # Check for missing values ---------------------------------------------------
    missing_values = sapply(1:nrow(x), function(row){anyNA(x[row, ])})
    if(sum(missing_values) == nrow(x))
      stop(paste0("All rows in x contain at least one missing response.\n",
                  "You could try option na_action = impute."))
    if(sum(missing_values) > 1)
      warning(paste0("There were ",
                     sum(missing_values),
                     " rows with missing observations in the input matrix x.\n",
                     "Since na_action = listwise these rows were excluded from the analysis."),
              call. = FALSE)
    if(sum(missing_values) == 1)
      warning(paste0("There was one row with missing observations in the input matrix x.\n",
                     "Since na_action = listwise this row was excluded from \n",
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
                                 na_action,
                                 variable_bool,
                                 reference_category,
                                 main_difference_model) {
  if(na_action == "listwise") {
    # Check for missing values in x --------------------------------------------
    missing_values = sapply(1:nrow(x), function(row){anyNA(x[row, ])})
    if(sum(missing_values) == nrow(x))
      stop(paste0("All rows in x contain at least one missing response.\n",
                  "You could try option na_action = impute."))
    if(sum(missing_values) > 1)
      warning(paste0("There were ",
                     sum(missing_values),
                     " rows with missing observations in the input matrix x.\n",
                     "Since na_action = listwise these rows were excluded from the analysis."),
              call. = FALSE)
    if(sum(missing_values) == 1)
      warning(paste0("There was one row with missing observations in the input matrix x.\n",
                     "Since na_action = listwise this row was excluded from \n",
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