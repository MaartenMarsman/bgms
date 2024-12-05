#include <Rcpp.h>
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// Compute partition coefficient for the MFM - SBM
// ----------------------------------------------------------------------------|
Rcpp::NumericVector compute_Vn_mfm_sbm(int no_variables,
                                       double dirichlet_alpha,
                                       int t_max);

// ----------------------------------------------------------------------------|
// Sample the block allocations for the MFM - SBM
// ----------------------------------------------------------------------------|
Rcpp::IntegerVector block_allocations_mfm_sbm(Rcpp::IntegerVector cluster_assign,
                                              int no_variables,
                                              Rcpp::NumericVector log_Vn,
                                              Rcpp::NumericMatrix block_probs,
                                              Rcpp::IntegerMatrix indicator,
                                              int dirichlet_alpha,
                                              double beta_bernoulli_alpha,
                                              double beta_bernoulli_beta);

// ----------------------------------------------------------------------------|
// Sample the number of clusters (K) based on Equation (3.7) from Miller & Harrison
// ----------------------------------------------------------------------------|
int sample_K_mfm_sbm(IntegerVector cluster_assign,
                     double dirichlet_alpha,
                     NumericVector log_Vn,
                     int max_K = 12);

// ----------------------------------------------------------------------------|
// Sample the block parameters for the MFM - SBM
// ----------------------------------------------------------------------------|
Rcpp::NumericMatrix block_probs_mfm_sbm(Rcpp::IntegerVector cluster_assign,
                                        Rcpp::IntegerMatrix indicator,
                                        int no_variables,
                                        double beta_bernoulli_alpha,
                                        double beta_bernoulli_beta);