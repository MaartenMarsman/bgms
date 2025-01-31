#' @importFrom utils packageVersion

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
    na_action = na_action, na_impute = na_impute, iter = iter, burnin = burnin,
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