// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the threshold parameters
// ----------------------------------------------------------------------------|
NumericMatrix metropolis_thresholds(NumericMatrix interactions, 
                                    NumericMatrix thresholds,
                                    IntegerMatrix observations,
                                    IntegerVector no_categories,
                                    IntegerMatrix n_cat_obs, 
                                    double threshold_alpha,
                                    double threshold_beta) {
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  
  NumericVector g(no_persons);  
  NumericVector q(no_persons);  
  
  double c = 1.0;
  double log_prob = 0.0;  
  double rest_score;
  double a;
  double b;
  double tmp;
  double proposed_state;
  double current_state;
  double U;
  
  for(int node = 0; node < no_nodes; node++) {
    for(int category = 0; category < no_categories[node]; category++) {
      current_state = thresholds(node, category);
      //The pseudolikelihood has the form exp(m) / (g + q * exp(m)) 
      //The proposal has the form [c * exp(m)] ^ a / (1 + c * exp(m)) ^ (a + b)
      //First, we compute g, q, and c
      c = (threshold_alpha + threshold_beta) / (1 + std::exp(current_state));
      for(int person = 0; person < no_persons; person++) {
        g[person] = 1.0;
        q[person] = 1.0;
        rest_score = 0.0;
        for(int s = 0; s <  no_nodes; s++) {
          rest_score += observations(person, s) *                              
            interactions(s, node);
        }
        for(int cat = 0; cat < no_categories[node]; cat++) {
          if(cat != category) {
            g[person] += std::exp(thresholds(node, cat) + 
              (cat + 1) * rest_score);  
          }
        } 
        q[person] = std::exp((category + 1) * rest_score);    
        c +=  q[person] / (g[person] + q[person] * std::exp(current_state));
      }
      c = c / ((no_persons + threshold_alpha + threshold_beta) - 
        std::exp(current_state) * c);
      
      //Proposal is generalized beta-prime. 
      a = n_cat_obs(category + 1, node) + threshold_alpha;
      b = no_persons + threshold_beta - n_cat_obs(category + 1, node);
      tmp = R::rbeta(a, b);
      proposed_state = std::log(tmp / (1  - tmp) / c);
      
      //Compute log_acceptance probability for Metropolis.
      //First, we use g and q above to compute the ratio of pseudolikelihoods
      log_prob = 0;
      for(int person = 0; person < no_persons; person++) {
        log_prob += std::log(g[person] + q[person] * std::exp(current_state));
        log_prob -= std::log(g[person] + q[person] * std::exp(proposed_state));
      } 
      //Second, we add the ratio of prior probabilities
      log_prob -= (threshold_alpha + threshold_beta) * 
        std::log(1 + std::exp(proposed_state));
      log_prob += (threshold_alpha + threshold_beta) * 
        std::log(1 + std::exp(current_state)); 
      //Third, we add the ratio of proposals
      log_prob -= (a + b) * std::log(1 + c * std::exp(current_state));
      log_prob += (a + b) * std::log(1 + c * std::exp(proposed_state));
      
      U = std::log(R::runif(0, 1));      
      if(U < log_prob) {
        thresholds(node, category) = proposed_state;
      }
    }
  }
  return thresholds;
}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for an interaction
// ----------------------------------------------------------------------------|
double log_pseudolikelihood_ratio(NumericMatrix interactions,
                                  NumericMatrix thresholds,
                                  IntegerMatrix observations,
                                  IntegerVector no_categories,
                                  int node1,
                                  int node2,
                                  double proposed_state,
                                  double current_state) {
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  double rest_score = 0.0;
  double bound =  0.0;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop = 0.0;
  double denominator_curr = 0.0;
  double exponent = 0.0;
  int score = 0;
  int obs_score1 = 0;
  int obs_score2 = 0;
  
  for(int person = 0; person < no_persons; person++) {
    obs_score1 = observations(person, node1);                                   
    obs_score2 = observations(person, node2);                                   
    
    pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 *                         
      (proposed_state - current_state);                                         

        //Node 1 log pseudolikelihood ratio
    rest_score = 0.0;
    for(int node = 0; node < no_nodes; node++) {
      rest_score += observations(person, node) * interactions(node, node1);  
    }
    rest_score -= obs_score2 * interactions(node2, node1);
    
    if(rest_score > 0) {
      bound = no_categories[node1] * rest_score;
    } else {
      bound = 0.0;
    }
    
    denominator_prop = std::exp(-bound);
    denominator_curr = std::exp(-bound);
    for(int category = 0; category < no_categories[node1]; category++) {
      score = category + 1;
      exponent = thresholds(node1, category) + 
        score * rest_score - 
        bound;
      denominator_prop += 
        std::exp(exponent + score * obs_score2 * proposed_state);
      denominator_curr += 
        std::exp(exponent + score * obs_score2 * current_state);
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);
    
    //Node 2 log pseudolikelihood ratio
    rest_score = 0.0;
    for(int node = 0; node < no_nodes; node++) {
      rest_score += observations(person, node) * interactions(node, node2);  
    }
    rest_score -= obs_score1 * interactions(node1, node2);
    
    if(rest_score > 0) {
      bound = no_categories[node2] * rest_score;
    } else {
      bound = 0.0;
    }
    
    denominator_prop = std::exp(-bound);
    denominator_curr = std::exp(-bound);
    for(int category = 0; category < no_categories[node2]; category++) {
      score = category + 1;
      exponent = thresholds(node2, category) + 
        score * rest_score - 
        bound;
      denominator_prop += 
        std::exp(exponent + score * obs_score1 * proposed_state);
      denominator_curr += 
        std::exp(exponent + score * obs_score1 * current_state);
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);
  }
  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of the active interaction 
//  parameters (using a cauchy prior)
// ----------------------------------------------------------------------------|
NumericMatrix metropolis_interactions_cauchy(NumericMatrix interactions, 
                                             NumericMatrix thresholds,
                                             IntegerMatrix gamma,
                                             IntegerMatrix observations,
                                             IntegerVector no_categories,
                                             NumericMatrix proposal_sd,
                                             double cauchy_scale) {
  int no_nodes = observations.ncol();
  
  double proposed_state;
  double current_state;
  double log_prob;
  double U;
  
  for(int node1 = 0; node1 <  no_nodes - 1; node1++) {
    for(int node2 = node1 + 1; node2 <  no_nodes; node2++)
      if(gamma(node1, node2) == 1) {
        current_state = interactions(node1, node2);
        proposed_state = R::rnorm(current_state,
                                  proposal_sd(node1, node2));
        
        log_prob = log_pseudolikelihood_ratio(interactions,
                                              thresholds,
                                              observations,
                                              no_categories,
                                              node1,
                                              node2,
                                              proposed_state,
                                              current_state);
        log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
        log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale, true);

        U = R::runif(0, 1);
        if(std::log(U) < log_prob) {
          interactions(node1, node2) = proposed_state;
          interactions(node2, node1) = proposed_state;
        }
      }
  }
  return interactions;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of the active interaction 
//  parameters (using a unit information prior)
// ----------------------------------------------------------------------------|
NumericMatrix metropolis_interactions_unitinfo(NumericMatrix interactions, 
                                               NumericMatrix thresholds,
                                               IntegerMatrix gamma,
                                               IntegerMatrix observations,
                                               IntegerVector no_categories,
                                               NumericMatrix proposal_sd,
                                               NumericMatrix unit_info) {
  int no_nodes = observations.ncol();
  
  double proposed_state;
  double current_state;
  double log_prob;
  double U;
  
  for(int node1 = 0; node1 <  no_nodes - 1; node1++) {
    for(int node2 = node1 + 1; node2 <  no_nodes; node2++)
      if(gamma(node1, node2) == 1) {
        current_state = interactions(node1, node2);
        proposed_state = R::rnorm(current_state,
                                  proposal_sd(node1, node2));
        
        log_prob = log_pseudolikelihood_ratio(interactions,
                                              thresholds,
                                              observations,
                                              no_categories,
                                              node1,
                                              node2,
                                              proposed_state,
                                              current_state);
        log_prob += R::dnorm(proposed_state, 
                             0.0, 
                             unit_info(node1, node2), 
                             true);
        log_prob -= R::dnorm(current_state, 
                             0.0, 
                             unit_info(node1, node2), 
                             true);
        
        U = R::runif(0, 1);
        if(std::log(U) < log_prob) {
          interactions(node1, node2) = proposed_state;
          interactions(node2, node1) = proposed_state;
        }
      }
  }
  return interactions;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of an edge + interaction 
//  pair (using a cauchy prior)
// ----------------------------------------------------------------------------|
List metropolis_edge_interaction_pair_cauchy(NumericMatrix interactions, 
                                             NumericMatrix thresholds,
                                             IntegerMatrix gamma,
                                             IntegerMatrix observations,
                                             IntegerVector no_categories,
                                             NumericMatrix proposal_sd,
                                             double cauchy_scale,
                                             IntegerMatrix index,
                                             int no_moves) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;
  
//  int no_moves = index.nrow();
  int node1;
  int node2;
  
  for(int cntr = 0; cntr < no_moves; cntr ++) {
    node1 = index(cntr, 1) - 1;
    node2 = index(cntr, 2) - 1;
    if(gamma(node1, node2) == 0) {
      current_state = interactions(node1, node2);
      proposed_state = R::rnorm(current_state,
                                proposal_sd(node1, node2));
      
      log_prob = log_pseudolikelihood_ratio(interactions,
                                            thresholds,
                                            observations,
                                            no_categories,
                                            node1,
                                            node2,
                                            proposed_state,
                                            current_state);
      log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_sd(node1, node2),
                           true);
      
      U = R::runif(0, 1);    
      if(std::log(U) < log_prob) {
        gamma(node1, node2) = 1;
        gamma(node2, node1) = 1;
        interactions(node1, node2) = proposed_state;
        interactions(node2, node1) = proposed_state;
      }
    } else {
      current_state = interactions(node1, node2);
      proposed_state = 0.0;
      
      log_prob = log_pseudolikelihood_ratio(interactions,
                                            thresholds,
                                            observations,
                                            no_categories,
                                            node1,
                                            node2,
                                            proposed_state,
                                            current_state);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_sd(node1, node2),
                           true);
      log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale,  true);
      
      U = R::runif(0, 1);    
      if(std::log(U) < log_prob) {
        gamma(node1, node2) = 0;
        gamma(node2, node1) = 0;
        interactions(node1, node2) = proposed_state;
        interactions(node2, node1) = proposed_state;
      }
    }
  }
  return List::create(Named("interactions") = interactions,
                      Named("gamma") = gamma);
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of an edge + interaction 
//  pair (using a unit information prior)
// ----------------------------------------------------------------------------|
List metropolis_edge_interaction_pair_unitinfo(NumericMatrix interactions, 
                                               NumericMatrix thresholds,
                                               IntegerMatrix gamma,
                                               IntegerMatrix observations,
                                               IntegerVector no_categories,
                                               NumericMatrix proposal_sd,
                                               NumericMatrix unit_info,
                                               IntegerMatrix index,
                                               int no_moves) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;
  
//  int no_moves = index.nrow();
  int node1;
  int node2;
  
  for(int cntr = 0; cntr < no_moves; cntr ++) {
    node1 = index(cntr, 1) - 1;
    node2 = index(cntr, 2) - 1;
    if(gamma(node1, node2) == 0) {
      current_state = interactions(node1, node2);
      proposed_state = R::rnorm(current_state,
                                proposal_sd(node1, node2));
      
      log_prob = log_pseudolikelihood_ratio(interactions,
                                            thresholds,
                                            observations,
                                            no_categories,
                                            node1,
                                            node2,
                                            proposed_state,
                                            current_state);
      log_prob += R::dnorm(proposed_state, 
                           0.0, 
                           unit_info(node1, node2),
                           true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_sd(node1, node2),
                           true);
      
      U = R::runif(0, 1);    
      if(std::log(U) < log_prob) {
        gamma(node1, node2) = 1;
        gamma(node2, node1) = 1;
        interactions(node1, node2) = proposed_state;
        interactions(node2, node1) = proposed_state;
      }
    } else {
      current_state = interactions(node1, node2);
      proposed_state = 0.0;
      
      log_prob = log_pseudolikelihood_ratio(interactions,
                                            thresholds,
                                            observations,
                                            no_categories,
                                            node1,
                                            node2,
                                            proposed_state,
                                            current_state);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_sd(node1, node2),
                           true);
      log_prob -= R::dnorm(current_state, 
                           0.0, 
                           unit_info(node1, node2),
                           true);
      
      U = R::runif(0, 1);    
      if(std::log(U) < log_prob) {
        gamma(node1, node2) = 0;
        gamma(node2, node1) = 0;
        interactions(node1, node2) = proposed_state;
        interactions(node2, node1) = proposed_state;
      }
    }
  }
  return List::create(Named("interactions") = interactions,
                      Named("gamma") = gamma);
}

// ----------------------------------------------------------------------------|
// The Gibbs sampler; returning raw samples
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List gibbs_samples(IntegerMatrix observations,
                   IntegerVector no_categories,
                   String interaction_prior,
                   double cauchy_scale,
                   NumericMatrix unit_info,
                   NumericMatrix proposal_sd,
                   IntegerMatrix Index,
                   int no_iterations,
                   IntegerMatrix n_cat_obs, 
                   double threshold_alpha,
                   double threshold_beta,
                   String moms_method,
                   bool display_progress = false){
  
  int cntr;
  int no_nodes = observations.ncol();
  int no_interactions = Index.nrow();
  int no_thresholds = sum(no_categories);
  int max_no_categories = max(no_categories);
  
  int no_moves = no_interactions;
  if(moms_method == "AddDelete") 
    no_moves = 1;
  
  IntegerVector v = seq(0, no_interactions - 1);
  IntegerVector order(no_interactions);
  IntegerMatrix index(no_interactions, 3);
  
  IntegerMatrix gamma(no_nodes, no_nodes);
  NumericMatrix interactions(no_nodes, no_nodes);
  NumericMatrix thresholds(no_nodes, max_no_categories);  
  
  NumericMatrix samples_gamma(no_iterations, no_interactions);
  NumericMatrix samples_int(no_iterations, no_interactions);
  NumericMatrix samples_tre(no_iterations, no_thresholds);
  
  //Progress bar
  Progress p(no_iterations, display_progress);
  //p.increment();
  
  for(int iteration = 0; iteration < no_iterations; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("samples.gamma") = samples_gamma,
                          Named("samples.interactions") = samples_int,
                          Named("samples.thresholds") = samples_tre);
    }
    p.increment();
    
    //Update interactions and model (between model move)
    if(moms_method == "AddDelete") {
      order = sample(v, 
                     1, 
                     false, 
                     R_NilValue);
      
      for(int cntr = 0; cntr < 1; cntr++) {
        index(cntr, 0) = Index(order[cntr], 0);
        index(cntr, 1) = Index(order[cntr], 1);
        index(cntr, 2) = Index(order[cntr], 2);
      }  
    } else {
      order = sample(v, 
                     no_interactions, 
                     false, 
                     R_NilValue);
      
      for(int cntr = 0; cntr < no_interactions; cntr++) {
        index(cntr, 0) = Index(order[cntr], 0);
        index(cntr, 1) = Index(order[cntr], 1);
        index(cntr, 2) = Index(order[cntr], 2);
      }  
    }
    
    if(interaction_prior == "Cauchy") {
      List out = metropolis_edge_interaction_pair_cauchy(interactions,
                                                         thresholds,
                                                         gamma,
                                                         observations,
                                                         no_categories,
                                                         proposal_sd,
                                                         cauchy_scale,
                                                         index, 
                                                         no_moves);
      IntegerMatrix gamma = out["gamma"];
      NumericMatrix interactions = out["interactions"];
    }
    if(interaction_prior ==  "UnitInfo") {
      List out = metropolis_edge_interaction_pair_unitinfo(interactions,
                                                           thresholds,
                                                           gamma,
                                                           observations,
                                                           no_categories,
                                                           proposal_sd,
                                                           unit_info,
                                                           index, 
                                                           no_moves);
      IntegerMatrix gamma = out["gamma"];
      NumericMatrix interactions = out["interactions"];
    }
    
    //Update interactions (within model move)
    if(moms_method != "Gibbs") {
      if(interaction_prior == "Cauchy") {
        interactions = metropolis_interactions_cauchy(interactions,
                                                      thresholds,
                                                      gamma,
                                                      observations,
                                                      no_categories,
                                                      proposal_sd,
                                                      cauchy_scale);
      }
      if(interaction_prior == "UnitInfo") {
        interactions = metropolis_interactions_unitinfo(interactions,
                                                        thresholds,
                                                        gamma,
                                                        observations,
                                                        no_categories,
                                                        proposal_sd,
                                                        unit_info);
      }  
    }
    
    
    //Update thresholds
    thresholds = metropolis_thresholds(interactions,
                                       thresholds,
                                       observations,
                                       no_categories,
                                       n_cat_obs,
                                       threshold_alpha,
                                       threshold_beta);
    
    //Update estimators
    cntr = 0;
    for(int node1 = 0; node1 < no_nodes - 1; node1++) {
      for(int node2 = node1 + 1; node2 < no_nodes;node2++) {
        samples_gamma(iteration, cntr) = gamma(node1, node2);
        samples_int(iteration, cntr) = interactions(node1, node2);
        cntr++;
      }
    }
    cntr = 0;
    for(int node = 0; node < no_nodes; node++) {
      for(int category = 0; category < no_categories[node]; category++) {
        samples_tre(iteration, cntr) = thresholds(node, category);
        cntr++;
      }
    }
  }
  
  return List::create(Named("samples.gamma") = samples_gamma,
                      Named("samples.interactions") = samples_int,
                      Named("samples.thresholds") = samples_tre);
}

// ----------------------------------------------------------------------------|
// The Gibbs sampler; returning eap estimates
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List gibbs_eap(IntegerMatrix observations,
               IntegerVector no_categories,
               String interaction_prior,
               double cauchy_scale,
               NumericMatrix unit_info,
               NumericMatrix proposal_sd,
               IntegerMatrix Index,
               int no_iterations,
               IntegerMatrix n_cat_obs, 
               double threshold_alpha,
               double threshold_beta,
               String moms_method,
               bool display_progress = false) {
  int no_nodes = observations.ncol();
  int no_interactions = Index.nrow();
  int max_no_categories = max(no_categories);
  
  int no_moves = no_interactions;
  if(moms_method == "AddDelete") 
    no_moves = 1;
  
  IntegerVector v = seq(0, no_interactions - 1);
  IntegerVector order(no_interactions);
  IntegerMatrix index(no_interactions, 3);
  
  IntegerMatrix gamma(no_nodes, no_nodes);
  NumericMatrix interactions(no_nodes, no_nodes);
  NumericMatrix thresholds(no_nodes, max_no_categories);
  
  NumericMatrix edge_estimates(no_nodes, no_nodes);
  NumericMatrix category_estimates(no_nodes, max(no_categories));
  
  //Progress bar
  Progress p(no_iterations, display_progress);
  
  for(int iteration = 0; iteration < no_iterations; iteration++) {
    p.increment();
    
    if (Progress::check_abort()){
      return List::create(Named("eap.edges") = edge_estimates,
                          Named("eap.thresholds") = category_estimates);
    }
    
    //Update interactions and model (between model move)
    if(moms_method == "AddDelete") {
      order = sample(v, 
                     1, 
                     false, 
                     R_NilValue);
      
      for(int cntr = 0; cntr < 1; cntr++) {
        index(cntr, 0) = Index(order[cntr], 0);
        index(cntr, 1) = Index(order[cntr], 1);
        index(cntr, 2) = Index(order[cntr], 2);
      }  
    } else {
      order = sample(v, 
                     no_interactions, 
                     false, 
                     R_NilValue);
      
      for(int cntr = 0; cntr < no_interactions; cntr++) {
        index(cntr, 0) = Index(order[cntr], 0);
        index(cntr, 1) = Index(order[cntr], 1);
        index(cntr, 2) = Index(order[cntr], 2);
      }  
    }
    
    if(interaction_prior == "Cauchy") {
      List out = metropolis_edge_interaction_pair_cauchy(interactions,
                                                         thresholds,
                                                         gamma,
                                                         observations,
                                                         no_categories,
                                                         proposal_sd,
                                                         cauchy_scale,
                                                         index,
                                                         no_moves);
      IntegerMatrix gamma = out["gamma"];
      NumericMatrix interactions = out["interactions"];
    }
    if(interaction_prior ==  "UnitInfo") {
      List out = metropolis_edge_interaction_pair_unitinfo(interactions,
                                                           thresholds,
                                                           gamma,
                                                           observations,
                                                           no_categories,
                                                           proposal_sd,
                                                           unit_info,
                                                           index,
                                                           no_moves);
      IntegerMatrix gamma = out["gamma"];
      NumericMatrix interactions = out["interactions"];
    }
    
    //Update interactions (within model move)
    if(moms_method != "Gibbs") {
      if(interaction_prior == "Cauchy") {
        interactions = metropolis_interactions_cauchy(interactions,
                                                      thresholds,
                                                      gamma,
                                                      observations,
                                                      no_categories,
                                                      proposal_sd,
                                                      cauchy_scale);
      }
      if(interaction_prior == "UnitInfo") {
        interactions = metropolis_interactions_unitinfo(interactions,
                                                        thresholds,
                                                        gamma,
                                                        observations,
                                                        no_categories,
                                                        proposal_sd,
                                                        unit_info);
      }      
    }
    
    //Update thresholds
    thresholds = metropolis_thresholds(interactions,
                                       thresholds,
                                       observations,
                                       no_categories,
                                       n_cat_obs,
                                       threshold_alpha,
                                       threshold_beta);
    
    //Update estimators
    for(int node1 = 0; node1 < no_nodes - 1; node1++) {
      for(int node2 = node1 + 1; node2 < no_nodes; node2++) {
        //Upper diagonal contains the eap estimate of inclusion probabilities
        edge_estimates(node1, node2) *= iteration;
        if(gamma(node1, node2) > 0) {
          edge_estimates(node1, node2) += 1.0;
        }
        edge_estimates(node1, node2) /= (iteration + 1);
        
        //Lower diagonal contains the eap estimate of interaction parameters
        edge_estimates(node2, node1) *= iteration;
        edge_estimates(node2, node1) += interactions(node1, node2);
        edge_estimates(node2, node1) /= (iteration + 1);
      }
      for(int category = 0; category < no_categories[node1];category++) {
        category_estimates(node1, category) *= iteration;
        category_estimates(node1, category) += thresholds(node1, category);
        category_estimates(node1, category) /= iteration + 1;
      }
    }
    for(int category = 0; category < no_categories[no_nodes - 1];category++) {
      category_estimates(no_nodes - 1, category) *= iteration;
      category_estimates(no_nodes - 1, category) += 
        thresholds(no_nodes - 1, category);
      category_estimates(no_nodes - 1, category) /= iteration + 1;
    }
  }
  
  return List::create(Named("eap.edges") = edge_estimates,
                      Named("eap.thresholds") = category_estimates);
}
