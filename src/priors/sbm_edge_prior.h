#pragma once

/**
 * @file sbm_edge_prior.h
 * @brief Gibbs samplers for the Mixture of Finite Mixtures
 *        Stochastic Block Model (MFM-SBM) edge prior.
 *
 * Implements the block-allocation and block-probability updates
 * described in Geng, Bhattacharya & Pati (2019, JASA 114:526).
 * Called by StochasticBlockEdgePrior::update() in edge_prior.h.
 */

#include <RcppArmadillo.h>
#include "edge_prior_correction.h"
struct SafeRNG;


/**
 * Sample block allocations for the MFM-SBM.
 *
 * Reassigns each variable to a cluster via a collapsed Gibbs step,
 * integrating over block-level inclusion probabilities.
 *
 * @param cluster_assign   Current cluster assignment (length p).
 * @param no_variables      Number of variables p.
 * @param log_Vn            Log partition coefficients from compute_Vn_mfm_sbm().
 * @param block_probs       Current block-level inclusion probability matrix.
 * @param indicator         Edge indicator matrix (p x p, upper-triangular).
 * @param dirichlet_alpha   Dirichlet concentration parameter for cluster sizes.
 * @param beta_bernoulli_alpha   Beta-Bernoulli alpha for within-block edges.
 * @param beta_bernoulli_beta    Beta-Bernoulli beta for within-block edges.
 * @param beta_bernoulli_alpha_between  Beta-Bernoulli alpha for between-block edges.
 * @param beta_bernoulli_beta_between   Beta-Bernoulli beta for between-block edges.
 * @param rng               Random number generator.
 * @return Updated cluster assignment vector (length p).
 */
arma::uvec block_allocations_mfm_sbm(arma::uvec cluster_assign,
                                                arma::uword no_variables,
                                                arma::vec log_Vn,
                                                arma::mat block_probs,
                                                arma::umat indicator,
                                                arma::uword dirichlet_alpha,
                                                double beta_bernoulli_alpha,
                                                double beta_bernoulli_beta,
                                                double beta_bernoulli_alpha_between,
                                                double beta_bernoulli_beta_between,
                                                SafeRNG& rng);

/**
 * Sample block-level inclusion probabilities for the MFM-SBM.
 *
 * Draws within-block probabilities from Beta(alpha + included, beta + excluded)
 * and between-block probabilities from separate Beta hyperparameters.
 *
 * @param cluster_assign   Current cluster assignment (length p).
 * @param indicator         Edge indicator matrix (p x p, upper-triangular).
 * @param no_variables      Number of variables p.
 * @param beta_bernoulli_alpha   Beta-Bernoulli alpha for within-block edges.
 * @param beta_bernoulli_beta    Beta-Bernoulli beta for within-block edges.
 * @param beta_bernoulli_alpha_between  Beta-Bernoulli alpha for between-block edges.
 * @param beta_bernoulli_beta_between   Beta-Bernoulli beta for between-block edges.
 * @return Block-level inclusion probability matrix (K x K).
 */
arma::mat block_probs_mfm_sbm(arma::uvec cluster_assign,
                                        arma::umat indicator,
                                        arma::uword no_variables,
                                        double beta_bernoulli_alpha,
                                        double beta_bernoulli_beta,
                                        double beta_bernoulli_alpha_between,
                                        double beta_bernoulli_beta_between,
                                        SafeRNG& rng);
/**
 * Corrected block-allocation sweep under the determinant-tilted GGM prior.
 *
 * The collapsed Gibbs weights gain the normalizing-constant ratio via mini
 * thermodynamic integration along the block-morph, and the new-cluster
 * weight uses the corrected (homogeneous whole-pair) collapsed marginal.
 * Grows and shrinks block_probs in place so the subsequent corrected
 * block-probability draw centers on the carried values.
 *
 * @param cluster_assign   Current cluster assignment (length p).
 * @param no_variables     Number of variables p.
 * @param log_Vn           Log partition coefficients.
 * @param block_probs      Block probability matrix (modified in place).
 * @param indicator        Edge indicator matrix (p x p).
 * @param dirichlet_alpha  Dirichlet concentration parameter.
 * @param beta_bernoulli_alpha        Within-block Beta alpha.
 * @param beta_bernoulli_beta         Within-block Beta beta.
 * @param beta_bernoulli_alpha_between Between-block Beta alpha.
 * @param beta_bernoulli_beta_between  Between-block Beta beta.
 * @param correction       Normalizing-constant correction curves.
 * @param rng              Random number generator.
 * @return Updated cluster assignment vector (length p).
 */
arma::uvec block_allocations_mfm_sbm_corrected(arma::uvec cluster_assign,
                                               arma::uword no_variables,
                                               const arma::vec& log_Vn,
                                               arma::mat& block_probs,
                                               const arma::umat& indicator,
                                               arma::uword dirichlet_alpha,
                                               double beta_bernoulli_alpha,
                                               double beta_bernoulli_beta,
                                               double beta_bernoulli_alpha_between,
                                               double beta_bernoulli_beta_between,
                                               const SBMCorrection& correction,
                                               SafeRNG& rng);

/**
 * Corrected block-probability draw under the determinant-tilted GGM prior.
 *
 * Each pair's conjugate density is multiplied by the per-pair correction
 * with self-consistent per-edge slopes, and sampled by inverse CDF on a
 * fine grid centered at the pair's current value.
 *
 * @param cluster_assign       Current cluster assignment (length p).
 * @param block_probs_current  Current block probability matrix (window centers).
 * @param indicator            Edge indicator matrix (p x p).
 * @param no_variables         Number of variables p.
 * @param beta_bernoulli_alpha        Within-block Beta alpha.
 * @param beta_bernoulli_beta         Within-block Beta beta.
 * @param beta_bernoulli_alpha_between Between-block Beta alpha.
 * @param beta_bernoulli_beta_between  Between-block Beta beta.
 * @param correction           Normalizing-constant correction curves.
 * @param rng                  Random number generator.
 * @return Block-level inclusion probability matrix (K x K).
 */
arma::mat block_probs_mfm_sbm_corrected(const arma::uvec& cluster_assign,
                                        const arma::mat& block_probs_current,
                                        const arma::umat& indicator,
                                        arma::uword no_variables,
                                        double beta_bernoulli_alpha,
                                        double beta_bernoulli_beta,
                                        double beta_bernoulli_alpha_between,
                                        double beta_bernoulli_beta_between,
                                        const SBMCorrection& correction,
                                        SafeRNG& rng);

/**
 * Self-consistent per-edge tilt slopes given labels and block probabilities.
 */
arma::mat compute_ce_sbm(const arma::uvec& cluster_assign,
                         const arma::mat& block_probs,
                         arma::uword no_variables,
                         const SBMCorrection& correction);

/**
 * Baseline expected degree densities under (labels, block probabilities).
 */
arma::vec degrees_ld_sbm(const arma::uvec& cluster_assign,
                         const arma::mat& block_probs,
                         arma::uword no_variables,
                         const SBMCorrection& correction);

/**
 * Mini-TI normalizing-constant difference for one node's label move.
 */
double miniti_node_sbm(arma::uword node,
                       const arma::uvec& cluster_assign,
                       const arma::mat& block_probs,
                       arma::uword no_variables,
                       arma::uword cur,
                       arma::uword cand,
                       const arma::vec& deg_base,
                       const SBMCorrection& correction);

/**
 * Corrected collapsed log-marginal for assigning a node to a new cluster.
 */
double corrected_log_marginal_mfm_sbm(const arma::uvec& cluster_assign,
                                      const arma::umat& indicator,
                                      arma::uword node,
                                      arma::uword no_variables,
                                      double beta_bernoulli_alpha_between,
                                      double beta_bernoulli_beta_between,
                                      const SBMCorrection& correction);

/**
 * Mini-TI normalizing-constant contribution of removing one node's edges,
 * anchoring the new-cluster weight to the existing-cluster baseline.
 */
double miniti_removal_sbm(arma::uword node,
                          const arma::uvec& cluster_assign,
                          const arma::mat& block_probs,
                          arma::uword no_variables,
                          arma::uword cur,
                          const arma::vec& deg_base,
                          const SBMCorrection& correction);
