#include <Rcpp.h>
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// Compute partition coefficient for the MFM - SBM
// ----------------------------------------------------------------------------|
Rcpp::NumericVector compute_Vn_mfm_sbm(int no_nodes,
                                 double dirichlet_gamma,
                                 int t_max);

// ----------------------------------------------------------------------------|
// Sample the block allocations for the MFM - SBM
// ----------------------------------------------------------------------------|
Rcpp::IntegerVector block_allocations_mfm_sbm(Rcpp::IntegerVector cluster_assign,
                                              int no_nodes,
                                              Rcpp::NumericVector log_Vn,
                                              Rcpp::NumericMatrix block_probs,
                                              Rcpp::IntegerMatrix gamma,
                                              int dirichlet_gamma,
                                              double beta_alpha,
                                              double beta_beta);

// ----------------------------------------------------------------------------|
// Sample the block parameters for the MFM - SBM
// ----------------------------------------------------------------------------|
Rcpp::NumericMatrix block_probs_mfm_sbm(Rcpp::IntegerVector cluster_assign,
                                        Rcpp::IntegerMatrix gamma,
                                        int no_nodes,
                                        double beta_alpha,
                                        double beta_beta);
