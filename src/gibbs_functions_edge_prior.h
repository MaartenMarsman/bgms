#include <RcppArmadillo.h>

// ----------------------------------------------------------------------------|
// Compute partition coefficient for the MFM - SBM
// ----------------------------------------------------------------------------|
arma::vec compute_Vn_mfm_sbm(arma::uword no_variables,
                                        double dirichlet_alpha,
                                        arma::uword t_max,
                                        double lambda);

// ----------------------------------------------------------------------------|
// Sample the block allocations for the MFM - SBM
// ----------------------------------------------------------------------------|
arma::uvec block_allocations_mfm_sbm(arma::uvec cluster_assign,
                                                arma::uword no_variables,
                                                arma::vec log_Vn,
                                                arma::mat block_probs,
                                                arma::umat indicator,
                                                arma::uword dirichlet_alpha,
                                                double beta_bernoulli_alpha,
                                                double beta_bernoulli_beta);

// ----------------------------------------------------------------------------|
// Sample the block parameters for the MFM - SBM
// ----------------------------------------------------------------------------|
arma::mat block_probs_mfm_sbm(arma::uvec cluster_assign,
                                        arma::umat indicator,
                                        arma::uword no_variables,
                                        double beta_bernoulli_alpha,
                                        double beta_bernoulli_beta);