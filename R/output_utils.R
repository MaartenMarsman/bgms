#' @importFrom utils packageVersion

prepare_output_bgm = function (
    out, x, num_categories, iter, data_columnnames, is_ordinal_variable,
    save_options, burnin, interaction_scale, threshold_alpha, threshold_beta,
    na_action, na_impute, edge_selection, edge_prior, theta,
    beta_bernoulli_alpha, beta_bernoulli_beta, dirichlet_alpha, lambda,
    variable_type) {

  save = any(c(save_options$save_main, save_options$save_pairwise, save_options$save_indicator))

  arguments = list (
    prepared_data = x, num_variables = ncol(x), num_cases = nrow(x),
    na_impute = na_impute, variable_type = variable_type, iter = iter,
    burnin = burnin, interaction_scale = interaction_scale,
    threshold_alpha = threshold_alpha, threshold_beta = threshold_beta,
    edge_selection = edge_selection, edge_prior = edge_prior,
    inclusion_probability = theta, beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta =  beta_bernoulli_beta,
    dirichlet_alpha = dirichlet_alpha, lambda = lambda, na_action = na_action,
    save = save, version = packageVersion("bgms")
  )

  num_variables = ncol(x)
  results = list()

  # Basic posterior means
  results$posterior_mean_main = out$main
  results$posterior_mean_pairwise = out$pairwise

  # Assign names to posterior mean matrices
  rownames(results$posterior_mean_main) = data_columnnames
  col_names <- character(ncol(results$posterior_mean_main))
  # Default: name by category index
  col_names[] <- paste0("category ", seq_along(col_names))
  # Override first two if only Blume-Capel variables exist
  if (all(!is_ordinal_variable)) {
    results$posterior_mean_main = results$posterior_mean_main[, 1:2]
    col_names <- c("linear", "quadratic")
  } else if (any(!is_ordinal_variable)) {
    # If mixed, prefix first two with dual meaning
    col_names[1:2] <- c("category 1 / linear", "category 2 / quadratic")
  }
  colnames(results$posterior_mean_main) <- col_names


  rownames(results$posterior_mean_pairwise) = data_columnnames
  colnames(results$posterior_mean_pairwise) = data_columnnames

  if(edge_selection) {
    results$posterior_mean_indicator = out$inclusion_indicator
    rownames(results$posterior_mean_indicator) = data_columnnames
    colnames(results$posterior_mean_indicator) = data_columnnames
  }

  # Generate variable × category names
  names_variable_categories = character()
  for (v in seq_len(num_variables)) {
    if (is_ordinal_variable[v]) {
      cats = seq_len(num_categories[v])
      names_variable_categories = c(
        names_variable_categories,
        paste0(data_columnnames[v], "(", cats, ")")
      )
    } else {
      names_variable_categories = c(
        names_variable_categories,
        paste0(data_columnnames[v], "(linear)"),
        paste0(data_columnnames[v], "(quadratic)")
      )
    }
  }

  # Sample storage
  if (save_options$save_main && "main_samples" %in% names(out)) {
    results$main_effect_samples = out$main_samples
    colnames(results$main_effect_samples) = names_variable_categories
  }

  if (save_options$save_pairwise && "pairwise_samples" %in% names(out)) {
    edge_names = character()
    for (i in 1:(num_variables - 1)) {
      for (j in (i + 1):num_variables) {
        edge_names = c(edge_names, paste0(data_columnnames[i], "-", data_columnnames[j]))
      }
    }
    results$pairwise_effect_samples = out$pairwise_samples
    colnames(results$pairwise_effect_samples) = edge_names
  }

  if (edge_selection && save_options$save_indicator && "inclusion_indicator_samples" %in% names(out)) {
    edge_names = character()
    for (i in 1:(num_variables - 1)) {
      for (j in (i + 1):num_variables) {
        edge_names = c(edge_names, paste0(data_columnnames[i], "-", data_columnnames[j]))
      }
    }
    results$inclusion_indicator_samples = out$inclusion_indicator_samples
    colnames(results$inclusion_indicator_samples) = edge_names
  }

  # SBM postprocessing
  if (edge_selection && edge_prior == "Stochastic-Block" && "allocations" %in% names(out)) {
    results$allocations <- out$allocations
    # Requires that summarySBM() is available in namespace
    sbm_summary <- summarySBM(list(
      indicator = results$inclusion_indicator,
      allocations = results$allocations
    ), internal_call = TRUE)
    results$components <- sbm_summary$components
    results$allocations <- sbm_summary$allocations
  }

  results$arguments = arguments
  class(results) = "bgms"
  return(results)
}


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
    if(main_difference_model == "Free"){
      diag(results$posterior_mean_indicator) = NA
    }
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