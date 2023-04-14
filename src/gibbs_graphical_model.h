#include <Rcpp.h>
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the threshold parameters
// ----------------------------------------------------------------------------|
Rcpp::NumericMatrix metropolis_thresholds(Rcpp::NumericMatrix interactions,
                                          Rcpp::NumericMatrix thresholds,
                                          Rcpp::IntegerMatrix observations,
                                          Rcpp::IntegerVector no_categories,
                                          Rcpp::IntegerMatrix n_cat_obs,
                                          int no_persons,
                                          int no_nodes,
                                          double threshold_alpha,
                                          double threshold_beta,
                                          Rcpp::NumericMatrix rest_matrix);

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for an interaction
// ----------------------------------------------------------------------------|
double log_pseudolikelihood_ratio(Rcpp::NumericMatrix interactions,
                                  Rcpp::NumericMatrix thresholds,
                                  Rcpp::IntegerMatrix observations,
                                  Rcpp::IntegerVector no_categories,
                                  int no_persons,
                                  int node1,
                                  int node2,
                                  double proposed_state,
                                  double current_state,
                                  Rcpp::NumericMatrix rest_matrix);

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of the active interaction
//  parameters (using a cauchy prior)
// ----------------------------------------------------------------------------|
Rcpp::List metropolis_interactions_cauchy(Rcpp::NumericMatrix interactions,
                                          Rcpp::NumericMatrix thresholds,
                                          Rcpp::IntegerMatrix gamma,
                                          Rcpp::IntegerMatrix observations,
                                          Rcpp::IntegerVector no_categories,
                                          Rcpp::NumericMatrix proposal_sd,
                                          double cauchy_scale,
                                          int no_persons,
                                          int no_nodes,
                                          Rcpp::NumericMatrix rest_matrix);

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of the active interaction
//  parameters (using a unit information prior)
// ----------------------------------------------------------------------------|
Rcpp::List metropolis_interactions_unitinfo(Rcpp::NumericMatrix interactions,
                                            Rcpp::NumericMatrix thresholds,
                                            Rcpp::IntegerMatrix gamma,
                                            Rcpp::IntegerMatrix observations,
                                            Rcpp::IntegerVector no_categories,
                                            Rcpp::NumericMatrix proposal_sd,
                                            Rcpp::NumericMatrix unit_info,
                                            int no_persons,
                                            int no_nodes,
                                            Rcpp::NumericMatrix rest_matrix);