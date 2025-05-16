#include <RcppArmadillo.h>

void update_indicator_interaction_pair_with_metropolis_and_marginal_moms (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::mat& proposal_sd,
    const double interaction_scale,
    const arma::imat& index,
    const int num_interactions,
    const int num_persons,
    arma::mat& residual_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category
);