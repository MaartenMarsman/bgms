# build_spec_*() — per-family bgm_spec builders
#
# Split out of bgm_spec.R (cleanup S4). The user-facing bgm_spec() constructor,
# new_bgm_spec(), validate_bgm_spec(), and print.bgm_spec() remain in bgm_spec.R.
# These builders (+ the shared edge_prior_spec_fields / sampler_sublist helpers)
# turn the validated bgm_spec into the family-specific spec object.


# ==============================================================================
# Internal builders (one per model type)
# ==============================================================================

# The inclusion / edge-prior fields shared by the GGM, OMRF, and mixed-MRF spec
# prior lists, selected by name from the flattened edge prior (ep). Spliced into
# each builder's prior list via c() so the field set stays in sync. bgmCompare
# uses a separate difference-prior block and does not share this.
edge_prior_spec_fields = function(ep) {
  list(
    edge_selection = ep$edge_selection,
    edge_prior = ep$edge_prior,
    inclusion_probability = ep$inclusion_probability,
    beta_bernoulli_alpha = ep$beta_bernoulli_alpha,
    beta_bernoulli_beta = ep$beta_bernoulli_beta,
    beta_bernoulli_alpha_between = ep$beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = ep$beta_bernoulli_beta_between,
    dirichlet_alpha = ep$dirichlet_alpha,
    lambda = ep$lambda
  )
}


# ==============================================================================
# sampler_sublist()  --- extract validated sampler list for new_bgm_spec()
# ==============================================================================
sampler_sublist = function(s) {
  list(
    update_method     = s$update_method,
    target_accept     = s$target_accept,
    iter              = as.integer(s$iter),
    warmup            = as.integer(s$warmup),
    chains            = as.integer(s$chains),
    cores             = as.integer(s$cores),
    nuts_max_depth    = as.integer(s$nuts_max_depth),
    learn_mass_matrix = s$learn_mass_matrix,
    seed              = as.integer(s$seed),
    progress_type     = as.integer(s$progress_type),
    progress_callback = s$progress_callback,
    verbose           = isTRUE(s$verbose)
  )
}


build_spec_ggm = function(x, data_columnnames, num_variables,
                          variable_type, is_ordinal, is_continuous,
                          baseline_category,
                          na_action, sampler,
                          interaction_prior_type, pairwise_scale,
                          interaction_alpha, interaction_beta,
                          scale_prior_type, scale_shape, scale_rate,
                          scale_eta = NA_real_,
                          delta = 0,
                          precision_graph_prior = "joint",
                          zratio_active = FALSE,
                          edge_prior_flat) {
  # Missing data
  md = validate_missing_data(
    x = x, na_action = na_action,
    is_continuous = TRUE
  )
  x = md$x

  # Center continuous data (GGM likelihood assumes zero mean). The column
  # means are kept so prediction can center newdata on the training scale.
  column_means = colMeans(x)
  x = center_continuous_data(x)

  # Standardized-frame scale prior: derive the raw diagonal rate eta / s
  scale_rate = resolve_scale_rate(scale_rate, scale_eta, pairwise_scale)

  ep = edge_prior_flat

  new_bgm_spec(
    model_type = "ggm",
    data = list(
      x                = x,
      data_columnnames = data_columnnames,
      num_variables    = as.integer(ncol(x)),
      num_cases        = as.integer(nrow(x)),
      column_means     = column_means
    ),
    variables = list(
      variable_type     = variable_type,
      is_ordinal        = is_ordinal,
      is_continuous     = TRUE,
      baseline_category = baseline_category
    ),
    missing = list(
      na_action     = na_action,
      na_impute     = md$na_impute,
      missing_index = md$missing_index
    ),
    prior = c(
      list(
        interaction_prior_type = interaction_prior_type,
        pairwise_scale = pairwise_scale,
        interaction_alpha = interaction_alpha,
        interaction_beta = interaction_beta,
        scale_prior_type = scale_prior_type,
        scale_shape = scale_shape,
        scale_rate = scale_rate,
        scale_eta = scale_eta,
        delta = delta,
        precision_graph_prior = precision_graph_prior,
        zratio_active = zratio_active
      ),
      edge_prior_spec_fields(ep)
    ),
    sampler = sampler_sublist(sampler),
    precomputed = list()
  )
}


build_spec_omrf = function(x, data_columnnames, num_variables,
                           variable_type, is_ordinal, is_continuous,
                           baseline_category,
                           na_action, sampler,
                           interaction_prior_type, pairwise_scale,
                           interaction_alpha, interaction_beta,
                           threshold_prior_type, main_alpha, main_beta,
                           threshold_scale,
                           edge_prior_flat) {
  # Baseline category
  bc = validate_baseline_category(
    baseline_category = baseline_category,
    baseline_category_provided = !is.null(baseline_category),
    x = x,
    variable_bool = is_ordinal
  )

  # Missing data + ordinal recoding
  md = validate_missing_data(
    x = x, na_action = na_action,
    is_continuous = FALSE
  )
  x_clean = md$x
  ord = reformat_ordinal_data(
    x = x_clean, is_ordinal = is_ordinal,
    baseline_category = bc
  )
  x_recoded = ord$x
  num_categories = ord$num_categories
  bc_final = ord$baseline_category

  missing_index = md$missing_index

  ep = edge_prior_flat

  num_thresholds = sum(ifelse(is_ordinal, num_categories, 2L))

  new_bgm_spec(
    model_type = "omrf",
    data = list(
      x = x_recoded,
      data_columnnames = data_columnnames,
      num_variables = as.integer(num_variables),
      num_cases = as.integer(nrow(x_recoded)),
      num_categories = as.integer(num_categories),
      # Recode map (sorted original values per ordinal variable) so predict()
      # can recode newdata the same way bgm() recoded the training data.
      category_levels = ord$category_levels,
      # Additive shift to the 0-based scale per Blume-Capel variable.
      blume_capel_shift = ord$blume_capel_shift
    ),
    variables = list(
      variable_type     = variable_type,
      is_ordinal        = is_ordinal,
      is_continuous     = FALSE,
      baseline_category = as.integer(bc_final)
    ),
    missing = list(
      na_action     = na_action,
      na_impute     = md$na_impute,
      missing_index = missing_index
    ),
    prior = c(
      list(
        interaction_prior_type = interaction_prior_type,
        pairwise_scale = pairwise_scale,
        interaction_alpha = interaction_alpha,
        interaction_beta = interaction_beta,
        threshold_prior_type = threshold_prior_type,
        main_alpha = main_alpha,
        main_beta = main_beta,
        threshold_scale = threshold_scale
      ),
      edge_prior_spec_fields(ep)
    ),
    sampler = sampler_sublist(sampler),
    precomputed = list(
      num_thresholds = as.integer(num_thresholds)
    )
  )
}


# ------------------------------------------------------------------
# build_spec_mixed_mrf
# ------------------------------------------------------------------
# Builds a bgm_spec for the mixed MRF model (discrete + continuous).
# Splits the input data matrix into discrete and continuous parts,
# validates and recodes discrete variables (ordinal/BC), and assembles
# the spec with metadata needed by sample_mixed_mrf() and
# build_output_mixed_mrf().
# ------------------------------------------------------------------
build_spec_mixed_mrf = function(x, data_columnnames, num_variables,
                                variable_type, is_ordinal,
                                baseline_category,
                                na_action, sampler,
                                interaction_prior_type, pairwise_scale,
                                interaction_alpha, interaction_beta,
                                threshold_prior_type, main_alpha, main_beta,
                                threshold_scale,
                                means_prior_type, means_scale,
                                means_alpha, means_beta,
                                scale_prior_type, scale_shape, scale_rate,
                                scale_eta = NA_real_,
                                delta = 0,
                                precision_graph_prior = "joint",
                                zratio_active = FALSE,
                                edge_prior_flat) {
  # Standardized-frame scale prior: derive the raw diagonal rate eta / s
  scale_rate = resolve_scale_rate(scale_rate, scale_eta, pairwise_scale)

  # Identify discrete vs continuous columns
  cont_idx = which(variable_type == "continuous")
  disc_idx = which(variable_type != "continuous")
  p = length(disc_idx)
  q = length(cont_idx)

  # Split data
  x_disc = x[, disc_idx, drop = FALSE]
  x_cont = x[, cont_idx, drop = FALSE]

  # Ensure integer matrix for discrete data
  storage.mode(x_disc) = "integer"
  # Ensure numeric matrix for continuous data
  storage.mode(x_cont) = "double"

  # Discrete variable properties (subset to discrete columns)
  is_ordinal_disc = is_ordinal[disc_idx]
  vtype_disc = variable_type[disc_idx]

  # Subset baseline_category to discrete columns when the user supplies a

  # full-length vector (one entry per variable, including continuous ones).
  if(length(baseline_category) == num_variables && num_variables != p) {
    baseline_category = baseline_category[disc_idx]
  }

  # Baseline category for discrete variables
  bc = validate_baseline_category(
    baseline_category = baseline_category,
    baseline_category_provided = !is.null(baseline_category),
    x = x_disc,
    variable_bool = is_ordinal_disc
  )

  # Missing data handling
  na_impute = FALSE
  missing_index_discrete = NULL
  missing_index_continuous = NULL

  if(na_action == "listwise") {
    missing_rows = apply(x_disc, 1, anyNA) | apply(x_cont, 1, anyNA)
    if(all(missing_rows)) {
      stop(paste0(
        "All rows in x contain at least one missing response.\n",
        "You could try option na_action = \"impute\"."
      ))
    }
    n_removed = sum(missing_rows)
    if(n_removed > 0 && isTRUE(getOption("bgms.verbose", TRUE))) {
      n_remaining = nrow(x_disc) - n_removed
      message(
        n_removed, " row", if(n_removed > 1) "s" else "",
        " with missing values excluded (n = ", n_remaining, " remaining).\n",
        "To impute missing values instead, use na_action = \"impute\"."
      )
    }
    x_disc = x_disc[!missing_rows, , drop = FALSE]
    x_cont = x_cont[!missing_rows, , drop = FALSE]
    if(nrow(x_disc) < 2) {
      stop(paste0(
        "After removing missing observations from the input matrix x,\n",
        "there were less than two rows left in x."
      ))
    }
  } else {
    # Impute path: handle discrete and continuous sub-matrices separately
    md_disc = handle_impute(x_disc)
    md_cont = handle_impute(x_cont)
    x_disc = md_disc$x
    x_cont = md_cont$x
    na_impute = md_disc$na_impute || md_cont$na_impute
    if(md_disc$na_impute) missing_index_discrete = md_disc$missing_index
    if(md_cont$na_impute) missing_index_continuous = md_cont$missing_index
  }

  # Ordinal recoding (reformat discrete data)
  ord = reformat_ordinal_data(
    x = x_disc, is_ordinal = is_ordinal_disc,
    baseline_category = bc
  )
  x_disc_recoded = ord$x
  num_categories = ord$num_categories
  bc_final = ord$baseline_category

  ep = edge_prior_flat

  num_thresholds = sum(ifelse(is_ordinal_disc, num_categories, 2L))

  new_bgm_spec(
    model_type = "mixed_mrf",
    data = list(
      x_discrete = x_disc_recoded,
      x_continuous = x_cont,
      data_columnnames = data_columnnames,
      data_columnnames_discrete = data_columnnames[disc_idx],
      data_columnnames_continuous = data_columnnames[cont_idx],
      num_variables = as.integer(num_variables),
      num_discrete = as.integer(p),
      num_continuous = as.integer(q),
      num_cases = as.integer(nrow(x_disc_recoded)),
      num_categories = as.integer(num_categories),
      # Recode map (sorted original values per discrete variable) for original-
      # scale threshold labels; NULL for Blume-Capel.
      category_levels = ord$category_levels,
      # Additive shift to the 0-based scale per Blume-Capel variable.
      blume_capel_shift = ord$blume_capel_shift,
      discrete_indices = disc_idx,
      continuous_indices = cont_idx
    ),
    variables = list(
      variable_type     = variable_type,
      is_ordinal        = is_ordinal_disc,
      is_continuous     = FALSE,
      is_mixed          = TRUE,
      baseline_category = as.integer(bc_final)
    ),
    missing = list(
      na_action = na_action,
      na_impute = na_impute,
      missing_index_discrete = missing_index_discrete,
      missing_index_continuous = missing_index_continuous
    ),
    prior = c(
      list(
        interaction_prior_type = interaction_prior_type,
        pairwise_scale = pairwise_scale,
        interaction_alpha = interaction_alpha,
        interaction_beta = interaction_beta,
        threshold_prior_type = threshold_prior_type,
        main_alpha = main_alpha,
        main_beta = main_beta,
        threshold_scale = threshold_scale,
        means_prior_type = means_prior_type,
        means_scale = means_scale,
        means_alpha = means_alpha,
        means_beta = means_beta,
        scale_prior_type = scale_prior_type,
        scale_shape = scale_shape,
        scale_rate = scale_rate,
        scale_eta = scale_eta,
        delta = delta,
        precision_graph_prior = precision_graph_prior,
        zratio_active = zratio_active
      ),
      edge_prior_spec_fields(ep)
    ),
    sampler = sampler_sublist(sampler),
    precomputed = list(
      num_thresholds = as.integer(num_thresholds)
    )
  )
}


build_spec_compare = function(x, y, group_indicator,
                              data_columnnames, num_variables,
                              variable_type, is_ordinal, is_continuous,
                              baseline_category,
                              na_action, sampler,
                              interaction_prior_type, pairwise_scale,
                              interaction_alpha, interaction_beta,
                              threshold_prior_type, main_alpha, main_beta,
                              threshold_scale,
                              difference_selection, main_difference_selection,
                              difference_prior,
                              difference_scale, difference_probability,
                              difference_prior_type = "normal",
                              beta_bernoulli_alpha, beta_bernoulli_beta,
                              beta_bernoulli_alpha_between = 1,
                              beta_bernoulli_beta_between = 1,
                              dirichlet_alpha = 1,
                              lambda = 1) {
  # --- Combine x/y and create group vector ------------------------------------
  if(!is.null(y)) {
    y = data_check(y, "y")
    if(ncol(x) != ncol(y)) stop("x and y must have the same number of columns.")
  }
  if(is.null(y) && is.null(group_indicator)) {
    stop(paste0(
      "For multi-group designs, the bgmCompare function requires input for\n",
      "either y (group 2 data) or group_indicator (group indicator)."
    ))
  }

  if(!is.null(group_indicator)) {
    group_indicator = as.vector(group_indicator)
    if(anyNA(group_indicator)) {
      stop("group_indicator cannot contain missing values.")
    }
    if(length(group_indicator) != nrow(x)) {
      stop("Length of group_indicator must match number of rows in x.")
    }

    unique_g = unique(group_indicator)
    if(length(unique_g) == 0L) {
      stop("The bgmCompare function expects at least two groups, but the input group_indicator contains no group value.")
    }
    if(length(unique_g) == 1L) {
      stop("The bgmCompare function expects at least two groups, but the input group_indicator contains only one group value.")
    }
    if(length(unique_g) == length(group_indicator)) {
      stop("The input group_indicator contains only unique group values.")
    }

    # Number the groups by first appearance, whatever the indicator's storage
    # type: a character or factor indicator would otherwise coerce the recode
    # back to character and fail tabulate().
    group = match(group_indicator, unique_g)
    tab = tabulate(group)
    if(any(tab < 2L)) {
      stop("One or more groups only had one member in the input group_indicator.")
    }
    # The indicator's own values, in the order the numbering assigned them, so
    # human-facing output can say which group number is which group.
    group_labels = as.character(unique_g)
  } else {
    group = c(rep.int(1L, nrow(x)), rep.int(2L, nrow(y)))
    x = rbind(x, y)
    group_labels = c("x", "y")
  }

  num_variables = ncol(x)

  # --- Baseline category (needs combined x) -----------------------------------
  bc = validate_baseline_category(
    baseline_category = baseline_category,
    baseline_category_provided = !is.null(baseline_category),
    x = x,
    variable_bool = is_ordinal
  )

  # --- Difference prior -------------------------------------------------------
  dp = validate_difference_prior(
    difference_selection = difference_selection,
    difference_prior = difference_prior,
    difference_probability = difference_probability,
    num_variables = num_variables,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = beta_bernoulli_beta_between,
    dirichlet_alpha = dirichlet_alpha,
    lambda = lambda
  )

  # --- Missing data (compare path) --------------------------------------------
  md = validate_missing_data(
    x             = x,
    na_action     = na_action,
    is_continuous = FALSE,
    group         = group
  )
  x = md$x
  na_impute = md$na_impute
  missing_index = md$missing_index
  group = md$group

  # Post-listwise group validation (bgmCompare-specific) -----------------------
  if(na_action == "listwise" && md$n_removed > 0) {
    unique_g = unique(group)
    if(length(unique_g) == length(group)) {
      stop(paste0(
        "After rows with missing observations were excluded, there were no groups, as \n",
        "there were only unique values in the input g left."
      ))
    }
    if(length(unique_g) == 1) {
      stop(paste0(
        "After rows with missing observations were excluded, there were no groups, as \n",
        "there was only one value in the input g left."
      ))
    }
    # Renumbering drops the groups listwise deletion emptied; the labels follow
    # the surviving codes so label i still names group i.
    group_labels = group_labels[unique_g]
    group = match(group, unique_g)
    tab = tabulate(group)
    if(any(tab < 2)) {
      stop(paste0(
        "After rows with missing observations were excluded, one or more groups, only \n",
        "had one member in the input g."
      ))
    }
  }

  # --- Ordinal recoding (compare path) ----------------------------------------
  # Keep the original category values to build the recode map for predict().
  x_original = x
  ord = reformat_ordinal_data(
    x                 = x,
    is_ordinal        = is_ordinal,
    baseline_category = bc
  )
  x = ord$x
  num_categories = ord$num_categories
  bc_final = ord$baseline_category

  # --- Collapse categories across groups (compare-specific) -------------------
  col = collapse_categories_across_groups(
    x                 = x,
    group             = group,
    is_ordinal        = is_ordinal,
    num_categories    = num_categories,
    baseline_category = bc_final
  )
  x_recoded = col$x
  num_categories = col$num_categories
  bc_final = col$baseline_category
  # Per-group observation counts on the final category codes, one matrix per
  # ordinal variable (NULL for Blume-Capel). Carried into the fitted object so
  # a user can see which group-by-category cells are empty.
  category_support = col$category_support
  ordinal_variable = is_ordinal

  # Recode map per ordinal variable: a named vector mapping each original
  # category value to its final (collapsed) 0-based category. reformat + cross-
  # group collapse can merge categories, so this is many-to-one in general --
  # hence a lookup (names = original values) rather than a sorted-value vector.
  category_levels = vector("list", ncol(x_recoded))
  for(vi in seq_len(ncol(x_recoded))) {
    if(ordinal_variable[vi]) {
      pairs = unique(cbind(x_original[, vi], x_recoded[, vi]))
      lookup = pairs[, 2]
      names(lookup) = pairs[, 1]
      category_levels[[vi]] = lookup[order(as.numeric(names(lookup)))]
    }
  }

  # True gaps have to be reported from here, not from
  # collapse_categories_across_groups(): reformat_ordinal_data() runs first and
  # has already closed them, so by the time the cross-group pass sees the data
  # the values are contiguous. Read them off the supplied values instead. A
  # scale that merely starts above 0 (1..7, say) is an offset, not a gap, and
  # is not worth telling anyone about.
  gap_vars = integer(0)
  gap_counts = integer(0)
  for(vi in seq_len(ncol(x_recoded))) {
    if(!ordinal_variable[vi]) next
    observed = as.numeric(names(category_levels[[vi]]))
    missing_values = as.integer(diff(range(observed)) + 1 - length(observed))
    if(missing_values > 0) {
      gap_vars = c(gap_vars, vi)
      gap_counts = c(gap_counts, missing_values)
    }
  }
  if(length(gap_vars) > 0 && isTRUE(getOption("bgms.verbose", TRUE))) {
    labels = data_columnnames[gap_vars]
    message(
      "Some category values were not used by any group. They were dropped ",
      "and the remaining categories renumbered, for ",
      paste0("'", labels, "' (", gap_counts, " dropped)", collapse = ", "),
      ". No observed category was merged."
    )
  }

  num_variables = ncol(x_recoded)
  num_groups = length(unique(group))

  # Compute precomputed structures
  counts_per_category = compute_counts_per_category(
    x_recoded, num_categories, group
  )
  blume_capel_stats = compute_blume_capel_stats(
    x_recoded, bc_final, ordinal_variable, group
  )

  # Center BC variables for pairwise stats
  x_centered = x_recoded
  for(i in which(!ordinal_variable)) {
    x_centered[, i] = x_centered[, i] - bc_final[i]
  }
  pairwise_stats = compute_pairwise_stats(x_centered, group)

  # Index structures
  num_interactions = as.integer(num_variables * (num_variables - 1) / 2)

  main_effect_indices = matrix(NA_integer_, nrow = num_variables, ncol = 2)
  for(variable in seq_len(num_variables)) {
    if(variable > 1) {
      main_effect_indices[variable, 1] = 1L + main_effect_indices[variable - 1, 2]
    } else {
      main_effect_indices[variable, 1] = 0L
    }
    if(ordinal_variable[variable]) {
      main_effect_indices[variable, 2] = main_effect_indices[variable, 1] +
        num_categories[variable] - 1L
    } else {
      main_effect_indices[variable, 2] = main_effect_indices[variable, 1] + 1L
    }
  }

  pairwise_effect_indices = matrix(NA_integer_,
    nrow = num_variables, ncol = num_variables
  )
  tel = 0L
  for(v1 in seq_len(num_variables - 1)) {
    for(v2 in seq(v1 + 1, num_variables)) {
      pairwise_effect_indices[v1, v2] = tel
      pairwise_effect_indices[v2, v1] = tel
      tel = tel + 1L
    }
  }

  # Interaction index matrix (used by C++ to iterate edges in random order)
  interaction_index_matrix = matrix(0L, nrow = num_interactions, ncol = 3)
  counter = 0L
  for(v1 in seq_len(num_variables - 1)) {
    for(v2 in seq(v1 + 1, num_variables)) {
      counter = counter + 1L
      interaction_index_matrix[counter, ] = c(counter, v1 - 1L, v2 - 1L)
    }
  }

  # Group indices and projection
  group_indices = matrix(NA_integer_, nrow = num_groups, ncol = 2)
  observations = x_centered
  sorted_group = sort(group)
  row_permutation = integer(nrow(x_centered))
  for(g in unique(group)) {
    src = which(group == g)
    dst = which(sorted_group == g)
    observations[dst, ] = x_centered[src, ]
    row_permutation[src] = dst
    group_indices[g, 1] = as.integer(min(dst) - 1)
    group_indices[g, 2] = as.integer(max(dst) - 1)
  }

  # missing_index rows index observations, which is in group-sorted order;
  # map the 0-based row column onto that order.
  if(na_impute && nrow(missing_index) > 0) {
    missing_index[, 1] = row_permutation[missing_index[, 1] + 1L] - 1L
  }

  one = matrix(1, nrow = num_groups, ncol = num_groups)
  V = diag(num_groups) - one / num_groups
  projection = eigen(V)$vectors[, -num_groups, drop = FALSE]
  if(num_groups == 2) {
    projection = projection / sqrt(2)
  }

  new_bgm_spec(
    model_type = "compare",
    data = list(
      x = observations,
      data_columnnames = data_columnnames,
      num_variables = as.integer(num_variables),
      num_cases = as.integer(nrow(observations)),
      num_categories = as.integer(num_categories),
      category_levels = category_levels,
      category_support = category_support,
      # Additive shift to the 0-based scale per Blume-Capel variable (the
      # cross-group collapse leaves Blume-Capel columns unchanged).
      blume_capel_shift = ord$blume_capel_shift,
      group = as.integer(group),
      group_labels = group_labels,
      num_groups = as.integer(num_groups),
      group_indices = group_indices,
      projection = projection
    ),
    variables = list(
      variable_type     = variable_type,
      is_ordinal        = ordinal_variable,
      is_continuous     = FALSE,
      baseline_category = as.integer(bc_final)
    ),
    missing = list(
      na_action     = na_action,
      na_impute     = na_impute,
      missing_index = missing_index
    ),
    prior = list(
      interaction_prior_type = interaction_prior_type,
      pairwise_scale = pairwise_scale,
      interaction_alpha = interaction_alpha,
      interaction_beta = interaction_beta,
      threshold_prior_type = threshold_prior_type,
      main_alpha = main_alpha,
      main_beta = main_beta,
      threshold_scale = threshold_scale,
      difference_selection = dp$difference_selection,
      main_difference_selection = main_difference_selection,
      difference_prior = dp$difference_prior,
      difference_scale = difference_scale,
      difference_prior_type = difference_prior_type,
      inclusion_probability_difference = dp$inclusion_probability_difference,
      beta_bernoulli_alpha = dp$beta_bernoulli_alpha,
      beta_bernoulli_beta = dp$beta_bernoulli_beta,
      beta_bernoulli_alpha_between = dp$beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between = dp$beta_bernoulli_beta_between,
      dirichlet_alpha = dp$dirichlet_alpha,
      lambda = dp$lambda
    ),
    sampler = sampler_sublist(sampler),
    precomputed = list(
      counts_per_category      = counts_per_category,
      blume_capel_stats        = blume_capel_stats,
      pairwise_stats           = pairwise_stats,
      main_effect_indices      = main_effect_indices,
      pairwise_effect_indices  = pairwise_effect_indices,
      interaction_index_matrix = interaction_index_matrix
    )
  )
}
