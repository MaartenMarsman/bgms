xi_delta_matching = function(xi, delta, n) {
  n * log(n / xi) / (n / xi - 1) - delta ^ 2
}

set_slab = function(x, no_categories, thresholds, interactions) {
  no_persons = nrow(x)
  no_nodes = ncol(x)
  no_thresholds = sum(no_categories)
  no_interactions = no_nodes * (no_nodes - 1) / 2
  no_parameters = no_thresholds + no_interactions
  
  hessian = matrix(data = NA, 
                   nrow = no_parameters,
                   ncol = no_parameters)
  slab_var = matrix(0, 
                    nrow = no_nodes, 
                    ncol = no_nodes)
  
  # Determine asymptotic covariance matrix (inverse hessian) ------------------
  if(!hasArg("thresholds") || !hasArg("interactions")) {
    fit = try(mple(x = x, no_categories = no_categories), 
              silent = TRUE)
    if(class(fit) != "try-error") {
      thresholds = fit$thresholds
      interactions = fit$interactions
    } else {
      stop("Pseudolikelihood optimization failed. Please check the data. If the 
           data checks out, please try different starting values.")
    }
  }
  
  # Compute Hessian  -----------------------------------------------------------
  hessian[1:no_thresholds, 1:no_thresholds] = 
    hessian_thresholds_pl(interactions = interactions, 
                          thresholds = thresholds, 
                          observations = x, 
                          no_categories = no_categories)
  hessian[-(1:no_thresholds), -(1:no_thresholds)] = 
    hessian_interactions(interactions = interactions, 
                         thresholds = thresholds, 
                         observations = x, 
                         no_categories = no_categories, 
                         interaction_var = matrix(Inf, 
                                                  nrow = no_nodes, 
                                                  ncol = no_nodes))
  hessian[-(1:no_thresholds), 1:no_thresholds] = 
    hessian_crossparameters(interactions = interactions, 
                            thresholds = thresholds, 
                            observations = x, 
                            no_categories = no_categories)
  
  hessian[1:no_thresholds, -(1:no_thresholds)] = 
    t(hessian[-(1:no_thresholds), 1:no_thresholds])
  
  asymptotic_covariance = -solve(hessian)[-(1:no_thresholds), 
                                          -(1:no_thresholds)]
  
  # Specify spike and slab matrices -------------------------------------------
  cntr = 0
  for(node in 1:(no_nodes - 1)) {
    for(node_2 in (node + 1):no_nodes) {
      cntr = cntr + 1
      slab_var[node, node_2] = slab_var[node_2, node] = 
        no_persons * asymptotic_covariance[cntr, cntr]
    }
  }
  return(slab_var)
}
