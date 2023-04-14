// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "gibbs_graphical_model.h"
#include "gibbs_edge_model.h"
#include "gibbs_graphical_and_edge_model.h"
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// Gibbs step for graphical model parameters
// ----------------------------------------------------------------------------|
List gibbs_step_graphical_model(IntegerMatrix observations,
                                IntegerVector no_categories,
                                String interaction_prior,
                                double cauchy_scale,
                                NumericMatrix unit_info,
                                NumericMatrix proposal_sd,
                                IntegerMatrix index,
                                IntegerMatrix n_cat_obs,
                                double threshold_alpha,
                                double threshold_beta,
                                int no_persons,
                                int no_nodes,
                                int no_interactions,
                                int no_thresholds,
                                int max_no_categories,
                                IntegerMatrix gamma,
                                NumericMatrix interactions,
                                NumericMatrix thresholds,
                                NumericMatrix rest_matrix,
                                NumericMatrix inclusion) {

  if(interaction_prior == "Cauchy") {
    List out = metropolis_edge_interaction_pair_cauchy(interactions,
                                                       thresholds,
                                                       gamma,
                                                       observations,
                                                       no_categories,
                                                       proposal_sd,
                                                       cauchy_scale,
                                                       index,
                                                       no_interactions,
                                                       no_persons,
                                                       rest_matrix,
                                                       inclusion);
    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix rest_matrix = out["rest_matrix"];
  } else if(interaction_prior ==  "UnitInfo") {
    List out = metropolis_edge_interaction_pair_unitinfo(interactions,
                                                         thresholds,
                                                         gamma,
                                                         observations,
                                                         no_categories,
                                                         proposal_sd,
                                                         unit_info,
                                                         index,
                                                         no_interactions,
                                                         no_persons,
                                                         rest_matrix,
                                                         inclusion);
    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix rest_matrix = out["rest_matrix"];
  }

  //Update interactions (within model move)
  if(interaction_prior == "Cauchy") {
    List out = metropolis_interactions_cauchy(interactions,
                                              thresholds,
                                              gamma,
                                              observations,
                                              no_categories,
                                              proposal_sd,
                                              cauchy_scale,
                                              no_persons,
                                              no_nodes,
                                              rest_matrix);
    NumericMatrix interactions = out["interactions"];
    NumericMatrix rest_matrix = out["rest_matrix"];
  }
  if(interaction_prior == "UnitInfo") {
    List out = metropolis_interactions_unitinfo(interactions,
                                                thresholds,
                                                gamma,
                                                observations,
                                                no_categories,
                                                proposal_sd,
                                                unit_info,
                                                no_persons,
                                                no_nodes,
                                                rest_matrix);
    NumericMatrix interactions = out["interactions"];
    NumericMatrix rest_matrix = out["rest_matrix"];
  }

  //Update thresholds
  thresholds = metropolis_thresholds(interactions,
                                     thresholds,
                                     observations,
                                     no_categories,
                                     n_cat_obs,
                                     no_persons,
                                     no_nodes,
                                     threshold_alpha,
                                     threshold_beta,
                                     rest_matrix);

  return List::create(Named("gamma") = gamma,
                      Named("interactions") = interactions,
                      Named("thresholds") = thresholds,
                      Named("rest_matrix") = rest_matrix);
}

// ----------------------------------------------------------------------------|
// The Gibbs sampler without mfm-SBM
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List gibbs_sampler_no_sbm(IntegerMatrix observations,
                          IntegerMatrix gamma,
                          NumericMatrix interactions,
                          NumericMatrix thresholds,
                          IntegerVector no_categories,
                          String interaction_prior,
                          double cauchy_scale,
                          NumericMatrix unit_info,
                          String edge_prior,
                          NumericMatrix proposal_sd,
                          IntegerMatrix Index,
                          int iter,
                          int burnin,
                          IntegerMatrix n_cat_obs,
                          double threshold_alpha,
                          double threshold_beta,
                          double beta_alpha,
                          double beta_beta,
                          double theta = 0.5,
                          bool save = false,
                          bool display_progress = false){
  int cntr;
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  int no_interactions = Index.nrow();
  int no_thresholds = sum(no_categories);
  int max_no_categories = max(no_categories);

  // Used to randomly update edge-association pairs
  IntegerVector v = seq(0, no_interactions - 1);
  IntegerVector order(no_interactions);
  IntegerMatrix index(no_interactions, 3);

  //Size of output depends on ``save''
  int nrow = no_nodes;
  int ncol_edges = no_nodes;
  int ncol_thresholds = max_no_categories;

  if(save == true) {
    nrow = iter;
    ncol_edges= no_interactions;
    ncol_thresholds = no_thresholds;
  }

  NumericMatrix out_gamma(nrow, ncol_edges);
  NumericMatrix out_interactions(nrow, ncol_edges);
  NumericMatrix out_thresholds(nrow, ncol_thresholds);

  //Precompute the matrix of rest scores
  NumericMatrix rest_matrix(no_persons, no_nodes);
  for(int node1 = 0; node1 < no_nodes; node1++) {
    for(int person = 0; person < no_persons; person++) {
      for(int node2 = 0; node2 < no_nodes; node2++) {
        rest_matrix(person, node1) +=
          observations(person, node2) * interactions(node2, node1);
      }
    }
  }

  //Variable declaration edge priors
  NumericMatrix inclusion(no_nodes, no_nodes);

  for(int i = 0; i < no_nodes - 1; i++) {
    for(int j = i + 1; j < no_nodes; j++) {
      inclusion(i, j) = theta;
      inclusion(j, i) = theta;
    }
  }

  //Progress bar
  Progress p(iter + burnin, display_progress);

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------
  for(int iteration = 0; iteration < burnin; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("gamma") = out_gamma,
                          Named("interactions") = out_interactions,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_interactions,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    List out = gibbs_step_graphical_model(observations,
                                          no_categories,
                                          interaction_prior,
                                          cauchy_scale,
                                          unit_info,
                                          proposal_sd,
                                          index,
                                          n_cat_obs,
                                          threshold_alpha,
                                          threshold_beta,
                                          no_persons,
                                          no_nodes,
                                          no_interactions,
                                          no_thresholds,
                                          max_no_categories,
                                          gamma,
                                          interactions,
                                          thresholds,
                                          rest_matrix,
                                          inclusion);

    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix rest_matrix = out["rest_matrix"];

    if(edge_prior == "Beta-Binomial") {
      int sumG = 0;
      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          sumG += gamma(i, j);
        }
      }
      double theta = R::rbeta(beta_alpha + sumG,
                              beta_beta + no_nodes * no_nodes / 2 - sumG);

      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          inclusion(i, j) = theta;
          inclusion(j, i) = theta;
        }
      }
    }
  }

  //The post burn-in iterations ------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("gamma") = out_gamma,
                          Named("interactions") = out_interactions,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_interactions,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }


    List out = gibbs_step_graphical_model(observations,
                                          no_categories,
                                          interaction_prior,
                                          cauchy_scale,
                                          unit_info,
                                          proposal_sd,
                                          index,
                                          n_cat_obs,
                                          threshold_alpha,
                                          threshold_beta,
                                          no_persons,
                                          no_nodes,
                                          no_interactions,
                                          no_thresholds,
                                          max_no_categories,
                                          gamma,
                                          interactions,
                                          thresholds,
                                          rest_matrix,
                                          inclusion);

    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix rest_matrix = out["rest_matrix"];

    if(edge_prior == "Beta-Binomial") {
      int sumG = 0;
      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          sumG += gamma(i, j);
        }
      }
      double theta = R::rbeta(beta_alpha + sumG,
                              beta_beta + no_nodes * no_nodes / 2 - sumG);

      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          inclusion(i, j) = theta;
          inclusion(j, i) = theta;
        }
      }
    }

    //Output -------------------------------------------------------------------
    if(save == TRUE) {
      //Save raw samples -------------------------------------------------------
      cntr = 0;
      for(int node1 = 0; node1 < no_nodes - 1; node1++) {
        for(int node2 = node1 + 1; node2 < no_nodes;node2++) {
          out_gamma(iteration, cntr) = gamma(node1, node2);
          out_interactions(iteration, cntr) = interactions(node1, node2);
          cntr++;
        }
      }
      cntr = 0;
      for(int node = 0; node < no_nodes; node++) {
        for(int category = 0; category < no_categories[node]; category++) {
          out_thresholds(iteration, cntr) = thresholds(node, category);
          cntr++;
        }
      }
    } else {
      //Compute running averages -----------------------------------------------
      for(int node1 = 0; node1 < no_nodes - 1; node1++) {
        for(int node2 = node1 + 1; node2 < no_nodes; node2++) {
          out_gamma(node1, node2) *= iteration;
          out_gamma(node1, node2) += gamma(node1, node2);
          out_gamma(node1, node2) /= iteration + 1;
          out_gamma(node2, node1) = out_gamma(node1, node2);

          out_interactions(node1, node2) *= iteration;
          out_interactions(node1, node2) += interactions(node1, node2);
          out_interactions(node1, node2) /= iteration + 1;
          out_interactions(node2, node1) = out_interactions(node1, node2);
        }

        for(int category = 0; category < no_categories[node1]; category++) {
          out_thresholds(node1, category) *= iteration;
          out_thresholds(node1, category) += thresholds(node1, category);
          out_thresholds(node1, category) /= iteration + 1;
        }
      }
      for(int category = 0; category < no_categories[no_nodes - 1]; category++) {
        out_thresholds(no_nodes - 1, category) *= iteration;
        out_thresholds(no_nodes - 1, category) +=
          thresholds(no_nodes - 1, category);
        out_thresholds(no_nodes - 1, category) /= iteration + 1;
      }
    }
  }

  return List::create(Named("gamma") = out_gamma,
                      Named("interactions") = out_interactions,
                      Named("thresholds") = out_thresholds);
}

// ----------------------------------------------------------------------------|
// The Gibbs sampler with mfm-SBM
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List gibbs_sampler_sbm(IntegerMatrix observations,
                       IntegerMatrix gamma,
                       NumericMatrix interactions,
                       NumericMatrix thresholds,
                       IntegerVector no_categories,
                       String interaction_prior,
                       double cauchy_scale,
                       NumericMatrix unit_info,
                       NumericMatrix proposal_sd,
                       IntegerMatrix Index,
                       int iter,
                       int burnin,
                       IntegerMatrix n_cat_obs,
                       double threshold_alpha,
                       double threshold_beta,
                       double dirichlet_gamma,
                       double beta_alpha,
                       double beta_beta,
                       bool save = false,
                       bool display_progress = false){
  int cntr;
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  int no_interactions = Index.nrow();
  int no_thresholds = sum(no_categories);
  int max_no_categories = max(no_categories);

  // Used to randomly update edge-association pairs
  IntegerVector v = seq(0, no_interactions - 1);
  IntegerVector order(no_interactions);
  IntegerMatrix index(no_interactions, 3);

  //Size of output depends on ``save''
  int nrow = no_nodes;
  int ncol_edges = no_nodes;
  int ncol_thresholds = max_no_categories;

  if(save == true) {
    nrow = iter;
    ncol_edges= no_interactions;
    ncol_thresholds = no_thresholds;
  }

  NumericMatrix out_gamma(nrow, ncol_edges);
  NumericMatrix out_interactions(nrow, ncol_edges);
  NumericMatrix out_thresholds(nrow, ncol_thresholds);

  //Precompute the matrix of rest scores
  NumericMatrix rest_matrix(no_persons, no_nodes);
  for(int node1 = 0; node1 < no_nodes; node1++) {
    for(int person = 0; person < no_persons; person++) {
      for(int node2 = 0; node2 < no_nodes; node2++) {
        rest_matrix(person, node1) +=
          observations(person, node2) * interactions(node2, node1);
      }
    }
  }

  //Variable declaration edge prior
  NumericMatrix inclusion(no_nodes, no_nodes);
  IntegerVector cluster_allocation(no_nodes);
  cluster_allocation[0] = 0;
  cluster_allocation[1] = 1;
  for(int i = 2; i < no_nodes; i++) {
    double U = R::unif_rand();
    if(U > 0.5){
      cluster_allocation[i] = 1;
    } else {
      cluster_allocation[i] = 0;
    }
  }
  NumericMatrix cluster_prob = block_probs_mfm_sbm(cluster_allocation,
                                                   gamma,
                                                   no_nodes,
                                                   beta_alpha,
                                                   beta_beta);

  for(int i = 0; i < no_nodes - 1; i++) {
    for(int j = i + 1; j < no_nodes; j++) {
      inclusion(i, j) = cluster_prob(cluster_allocation[i], cluster_allocation[j]);
      inclusion(j, i) = cluster_prob(cluster_allocation[i], cluster_allocation[j]);
    }
  }

  NumericVector log_Vn = compute_Vn_mfm_sbm(no_nodes,
                                            dirichlet_gamma,
                                            no_nodes + 10);
  NumericMatrix out_allocations(iter, no_nodes);


  //Progress bar
  Progress p(iter + burnin, display_progress);

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------
  for(int iteration = 0; iteration < burnin; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("gamma") = out_gamma,
                          Named("interactions") = out_interactions,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_interactions,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    List out = gibbs_step_graphical_model(observations,
                                          no_categories,
                                          interaction_prior,
                                          cauchy_scale,
                                          unit_info,
                                          proposal_sd,
                                          index,
                                          n_cat_obs,
                                          threshold_alpha,
                                          threshold_beta,
                                          no_persons,
                                          no_nodes,
                                          no_interactions,
                                          no_thresholds,
                                          max_no_categories,
                                          gamma,
                                          interactions,
                                          thresholds,
                                          rest_matrix,
                                          inclusion);

    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix rest_matrix = out["rest_matrix"];


    cluster_allocation = block_allocations_mfm_sbm(cluster_allocation,
                                                   no_nodes,
                                                   log_Vn,
                                                   cluster_prob,
                                                   gamma,
                                                   dirichlet_gamma,
                                                   beta_alpha,
                                                   beta_beta);


    cluster_prob = block_probs_mfm_sbm(cluster_allocation,
                                       gamma,
                                       no_nodes,
                                       beta_alpha,
                                       beta_beta);

    for(int i = 0; i < no_nodes - 1; i++) {
      for(int j = i + 1; j < no_nodes; j++) {
        inclusion(i, j) = cluster_prob(cluster_allocation[i], cluster_allocation[j]);
        inclusion(j, i) = cluster_prob(cluster_allocation[i], cluster_allocation[j]);
      }
    }
  }
  //The post burn-in iterations ------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("gamma") = out_gamma,
                          Named("interactions") = out_interactions,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_interactions,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }


    List out = gibbs_step_graphical_model(observations,
                                          no_categories,
                                          interaction_prior,
                                          cauchy_scale,
                                          unit_info,
                                          proposal_sd,
                                          index,
                                          n_cat_obs,
                                          threshold_alpha,
                                          threshold_beta,
                                          no_persons,
                                          no_nodes,
                                          no_interactions,
                                          no_thresholds,
                                          max_no_categories,
                                          gamma,
                                          interactions,
                                          thresholds,
                                          rest_matrix,
                                          inclusion);

    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix rest_matrix = out["rest_matrix"];

    cluster_allocation = block_allocations_mfm_sbm(cluster_allocation,
                                                   no_nodes,
                                                   log_Vn,
                                                   cluster_prob,
                                                   gamma,
                                                   dirichlet_gamma,
                                                   beta_alpha,
                                                   beta_beta);


    cluster_prob = block_probs_mfm_sbm(cluster_allocation,
                                       gamma,
                                       no_nodes,
                                       beta_alpha,
                                       beta_beta);

    for(int i = 0; i < no_nodes - 1; i++) {
      for(int j = i + 1; j < no_nodes; j++) {
        inclusion(i, j) = cluster_prob(cluster_allocation[i], cluster_allocation[j]);
        inclusion(j, i) = cluster_prob(cluster_allocation[i], cluster_allocation[j]);
      }
    }

    //Output -------------------------------------------------------------------
    if(save == TRUE) {
      //Save raw samples -------------------------------------------------------
      cntr = 0;
      for(int node1 = 0; node1 < no_nodes - 1; node1++) {
        for(int node2 = node1 + 1; node2 < no_nodes;node2++) {
          out_gamma(iteration, cntr) = gamma(node1, node2);
          out_interactions(iteration, cntr) = interactions(node1, node2);
          cntr++;
        }
      }
      cntr = 0;
      for(int node = 0; node < no_nodes; node++) {
        for(int category = 0; category < no_categories[node]; category++) {
          out_thresholds(iteration, cntr) = thresholds(node, category);
          cntr++;
        }
      }
    } else {
      //Compute running averages -----------------------------------------------
      for(int node1 = 0; node1 < no_nodes - 1; node1++) {
        for(int node2 = node1 + 1; node2 < no_nodes; node2++) {
          out_gamma(node1, node2) *= iteration;
          out_gamma(node1, node2) += gamma(node1, node2);
          out_gamma(node1, node2) /= iteration + 1;
          out_gamma(node2, node1) = out_gamma(node1, node2);

          out_interactions(node1, node2) *= iteration;
          out_interactions(node1, node2) += interactions(node1, node2);
          out_interactions(node1, node2) /= iteration + 1;
          out_interactions(node2, node1) = out_interactions(node1, node2);
        }

        for(int category = 0; category < no_categories[node1]; category++) {
          out_thresholds(node1, category) *= iteration;
          out_thresholds(node1, category) += thresholds(node1, category);
          out_thresholds(node1, category) /= iteration + 1;
        }
      }
      for(int category = 0; category < no_categories[no_nodes - 1]; category++) {
        out_thresholds(no_nodes - 1, category) *= iteration;
        out_thresholds(no_nodes - 1, category) +=
          thresholds(no_nodes - 1, category);
        out_thresholds(no_nodes - 1, category) /= iteration + 1;
      }
    }
    for(int i = 0; i < no_nodes; i++) {
      out_allocations(iteration, i) = cluster_allocation[i] + 1;
    }
  }

  return List::create(Named("gamma") = out_gamma,
                      Named("interactions") = out_interactions,
                      Named("thresholds") = out_thresholds,
                      Named("allocations") = out_allocations);
}