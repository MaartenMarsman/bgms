#include <Rcpp.h>
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the threshold parameters
//   for a regular binary or ordinal variable
// ----------------------------------------------------------------------------|
void metropolis_thresholds_regular(Rcpp::NumericMatrix thresholds,
                                   Rcpp::IntegerMatrix observations,
                                   Rcpp::IntegerVector no_categories,
                                   Rcpp::IntegerMatrix n_cat_obs,
                                   int no_persons,
                                   int variable,
                                   double threshold_alpha,
                                   double threshold_beta,
                                   Rcpp::NumericMatrix rest_matrix);

// ----------------------------------------------------------------------------|
// Adaptive Metropolis algorithm to sample from the full-conditional of the
//   threshold parameters for a Blume-Capel ordinal variable
// ----------------------------------------------------------------------------|
void metropolis_thresholds_blumecapel(Rcpp::NumericMatrix thresholds,
                                      Rcpp::IntegerMatrix observations,
                                      Rcpp::IntegerVector no_categories,
                                      Rcpp::IntegerMatrix sufficient_blume_capel,
                                      int no_persons,
                                      int variable,
                                      Rcpp::IntegerVector reference_category,
                                      double threshold_alpha,
                                      double threshold_beta,
                                      Rcpp::NumericMatrix rest_matrix,
                                      Rcpp::NumericMatrix proposal_sd_blumecapel,
                                      double phi,
                                      double target_ar,
                                      int t,
                                      double epsilon_lo,
                                      double epsilon_hi);

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for an interaction
// ----------------------------------------------------------------------------|
double log_pseudolikelihood_ratio(Rcpp::NumericMatrix interactions,
                                  Rcpp::NumericMatrix thresholds,
                                  Rcpp::IntegerMatrix observations,
                                  Rcpp::IntegerVector no_categories,
                                  int no_persons,
                                  int variable1,
                                  int variable2,
                                  double proposed_state,
                                  double current_state,
                                  Rcpp::NumericMatrix rest_matrix,
                                  Rcpp::LogicalVector variable_bool,
                                  Rcpp::IntegerVector reference_category);
