#include <Rcpp.h>
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of an edge + interaction
//  pair (using a cauchy prior)
// ----------------------------------------------------------------------------|
Rcpp::List metropolis_edge_interaction_pair_cauchy(Rcpp::NumericMatrix interactions,
                                                   Rcpp::NumericMatrix thresholds,
                                                   Rcpp::IntegerMatrix gamma,
                                                   Rcpp::IntegerMatrix observations,
                                                   Rcpp::IntegerVector no_categories,
                                                   Rcpp::NumericMatrix proposal_sd,
                                                   double cauchy_scale,
                                                   Rcpp::IntegerMatrix index,
                                                   int no_interactions,
                                                   int no_persons,
                                                   Rcpp::NumericMatrix rest_matrix,
                                                   Rcpp::NumericMatrix inclusion);

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of an edge + interaction
//  pair (using a unit information prior)
// ----------------------------------------------------------------------------|
Rcpp::List metropolis_edge_interaction_pair_unitinfo(Rcpp::NumericMatrix interactions,
                                                     Rcpp::NumericMatrix thresholds,
                                                     Rcpp::IntegerMatrix gamma,
                                                     Rcpp::IntegerMatrix observations,
                                                     Rcpp::IntegerVector no_categories,
                                                     Rcpp::NumericMatrix proposal_sd,
                                                     Rcpp::NumericMatrix unit_info,
                                                     Rcpp::IntegerMatrix index,
                                                     int no_interactions,
                                                     int no_persons,
                                                     Rcpp::NumericMatrix rest_matrix,
                                                     Rcpp::NumericMatrix inclusion);