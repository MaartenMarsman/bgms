#pragma once

#include <memory>
#include <RcppArmadillo.h>
#include "rng/rng_utils.h"
#include "utils/common_helpers.h"
#include "edge_prior_correction.h"
#include "sbm_edge_prior.h"
#include "sbm_edge_prior_interface.h"


/**
 * Abstract base class for edge inclusion priors.
 *
 * The edge prior updates the inclusion probability matrix based on the
 * current edge indicators. This is independent of the model type (GGM, OMRF,
 * etc.), so it is implemented as a separate class hierarchy.
 *
 * The MCMC runner calls update() after each edge indicator update, passing
 * the current edge indicators and inclusion probability matrix. The edge
 * prior modifies inclusion_probability in place.
 */
class BaseEdgePrior {
public:
    virtual ~BaseEdgePrior() = default;

    /**
     * Resample the prior's latent state given the current edge indicators
     * and write the resulting per-edge inclusion probabilities into
     * inclusion_probability (symmetric, modified in place).
     *
     * @param edge_indicators       Current edge inclusion matrix (p x p)
     * @param inclusion_probability Inclusion probability matrix (p x p), updated in place
     * @param num_variables         Number of variables p
     * @param num_pairwise          Number of variable pairs p(p-1)/2
     * @param rng                   Random number generator
     */
    virtual void update(
        const arma::imat& edge_indicators,
        arma::mat& inclusion_probability,
        int num_variables,
        int num_pairwise,
        SafeRNG& rng
    ) = 0;

    /** Deep copy for parallel chains. */
    virtual std::unique_ptr<BaseEdgePrior> clone() const = 0;

    /** Whether the prior maintains cluster allocations (SBM). */
    virtual bool has_allocations() const { return false; }

    /** Cluster allocations, 1-based; empty when has_allocations() is false. */
    virtual arma::ivec get_allocations() const { return arma::ivec(); }

    /** Whether the prior carries a sampled inclusion parameter (BB theta). */
    virtual bool has_inclusion_parameter() const { return false; }

    /** Current inclusion parameter; NA_REAL when has_inclusion_parameter() is false. */
    virtual double get_inclusion_parameter() const {
        return NA_REAL;
    }
};


/**
 * Bernoulli edge prior (fixed inclusion probabilities, no update needed).
 */
class BernoulliEdgePrior : public BaseEdgePrior {
public:
    /** No-op: the inclusion probabilities are fixed. */
    void update(
        const arma::imat& /*edge_indicators*/,
        arma::mat& /*inclusion_probability*/,
        int /*num_variables*/,
        int /*num_pairwise*/,
        SafeRNG& /*rng*/
    ) override {
        // No-op: inclusion probabilities are fixed
    }

    std::unique_ptr<BaseEdgePrior> clone() const override {
        return std::make_unique<BernoulliEdgePrior>(*this);
    }
};


/**
 * Beta-Bernoulli edge prior.
 *
 * Draws a shared inclusion probability from Beta(alpha + #included,
 * beta + #excluded) and assigns it to all edges. With a normalizing-constant
 * correction attached (GGM path), the draw targets the corrected conditional
 * that carries the 1/C(theta) factor instead.
 */
class BetaBernoulliEdgePrior : public BaseEdgePrior {
public:
    BetaBernoulliEdgePrior(double alpha = 1.0, double beta = 1.0)
        : alpha_(alpha), beta_(beta),
          current_prob_(alpha / (alpha + beta)) {}

    /** Attach the whole-graph logC(theta) correction table (GGM path). */
    void set_correction(const EdgePriorCorrection& correction) {
        correction_ = correction;
    }

    void update(
        const arma::imat& edge_indicators,
        arma::mat& inclusion_probability,
        int num_variables,
        int num_pairwise,
        SafeRNG& rng
    ) override {
        int num_edges_included = 0;
        for (int i = 0; i < num_variables - 1; i++) {
            for (int j = i + 1; j < num_variables; j++) {
                num_edges_included += edge_indicators(i, j);
            }
        }

        double a_post = alpha_ + num_edges_included;
        double b_post = beta_ + num_pairwise - num_edges_included;
        double prob = correction_.active()
            ? correction_.draw_theta(rng, a_post, b_post, current_prob_)
            : rbeta(rng, a_post, b_post);
        current_prob_ = prob;

        for (int i = 0; i < num_variables - 1; i++) {
            for (int j = i + 1; j < num_variables; j++) {
                inclusion_probability(i, j) = prob;
                inclusion_probability(j, i) = prob;
            }
        }
    }

    std::unique_ptr<BaseEdgePrior> clone() const override {
        return std::make_unique<BetaBernoulliEdgePrior>(*this);
    }

    /** Always true: the shared inclusion probability theta is sampled. */
    bool has_inclusion_parameter() const override { return true; }
    double get_inclusion_parameter() const override { return current_prob_; }

private:
    double alpha_;
    double beta_;
    double current_prob_;
    EdgePriorCorrection correction_;
};


/**
 * Stochastic Block Model (MFM-SBM) edge prior.
 *
 * Maintains cluster allocations and block-level inclusion probabilities.
 * Each edge's inclusion probability depends on its endpoints' cluster
 * assignments.
 */
class StochasticBlockEdgePrior : public BaseEdgePrior {
public:
    StochasticBlockEdgePrior(
        double beta_bernoulli_alpha,
        double beta_bernoulli_beta,
        double beta_bernoulli_alpha_between,
        double beta_bernoulli_beta_between,
        double dirichlet_alpha,
        double lambda
    ) : beta_bernoulli_alpha_(beta_bernoulli_alpha),
        beta_bernoulli_beta_(beta_bernoulli_beta),
        beta_bernoulli_alpha_between_(beta_bernoulli_alpha_between),
        beta_bernoulli_beta_between_(beta_bernoulli_beta_between),
        dirichlet_alpha_(dirichlet_alpha),
        lambda_(lambda),
        initialized_(false)
    {}

    /** Attach the per-pair correction curves and slope (mixed-MRF/GGM path). */
    void set_correction(const SBMCorrection& correction) {
        correction_ = correction;
    }

    /**
     * Initialize SBM state from the current edge indicators. Called
     * automatically on first update().
     */
    void initialize(
        const arma::imat& edge_indicators,
        arma::mat& inclusion_probability,
        int num_variables,
        SafeRNG& rng
    ) {
        cluster_allocations_.set_size(num_variables);
        cluster_allocations_[0] = 0;
        cluster_allocations_[1] = 1;
        for (int i = 2; i < num_variables; i++) {
            cluster_allocations_[i] = (runif(rng) > 0.5) ? 1 : 0;
        }

        cluster_prob_ = block_probs_mfm_sbm(
            cluster_allocations_,
            arma::conv_to<arma::umat>::from(edge_indicators),
            num_variables,
            beta_bernoulli_alpha_, beta_bernoulli_beta_,
            beta_bernoulli_alpha_between_, beta_bernoulli_beta_between_,
            rng
        );

        for (int i = 0; i < num_variables - 1; i++) {
            for (int j = i + 1; j < num_variables; j++) {
                inclusion_probability(i, j) = cluster_prob_(cluster_allocations_[i], cluster_allocations_[j]);
                inclusion_probability(j, i) = inclusion_probability(i, j);
            }
        }

        log_Vn_ = compute_Vn_mfm_sbm(
            num_variables, dirichlet_alpha_, num_variables + 10, lambda_);

        initialized_ = true;
    }

    void update(
        const arma::imat& edge_indicators,
        arma::mat& inclusion_probability,
        int num_variables,
        int /*num_pairwise*/,
        SafeRNG& rng
    ) override {
        if (!initialized_) {
            initialize(edge_indicators, inclusion_probability, num_variables, rng);
        }

        if (correction_.active()) {
            arma::umat indicator =
                arma::conv_to<arma::umat>::from(edge_indicators);
            cluster_allocations_ = block_allocations_mfm_sbm_corrected(
                cluster_allocations_, num_variables, log_Vn_, cluster_prob_,
                indicator, dirichlet_alpha_,
                beta_bernoulli_alpha_, beta_bernoulli_beta_,
                beta_bernoulli_alpha_between_, beta_bernoulli_beta_between_,
                correction_, rng
            );
            cluster_prob_ = block_probs_mfm_sbm_corrected(
                cluster_allocations_, cluster_prob_, indicator, num_variables,
                beta_bernoulli_alpha_, beta_bernoulli_beta_,
                beta_bernoulli_alpha_between_, beta_bernoulli_beta_between_,
                correction_, rng
            );
        } else {
            cluster_allocations_ = block_allocations_mfm_sbm(
                cluster_allocations_, num_variables, log_Vn_, cluster_prob_,
                arma::conv_to<arma::umat>::from(edge_indicators), dirichlet_alpha_,
                beta_bernoulli_alpha_, beta_bernoulli_beta_,
                beta_bernoulli_alpha_between_, beta_bernoulli_beta_between_, rng
            );

            cluster_prob_ = block_probs_mfm_sbm(
                cluster_allocations_,
                arma::conv_to<arma::umat>::from(edge_indicators), num_variables,
                beta_bernoulli_alpha_, beta_bernoulli_beta_,
                beta_bernoulli_alpha_between_, beta_bernoulli_beta_between_, rng
            );
        }

        for (int i = 0; i < num_variables - 1; i++) {
            for (int j = i + 1; j < num_variables; j++) {
                inclusion_probability(i, j) = cluster_prob_(cluster_allocations_[i], cluster_allocations_[j]);
                inclusion_probability(j, i) = inclusion_probability(i, j);
            }
        }
    }

    std::unique_ptr<BaseEdgePrior> clone() const override {
        return std::make_unique<StochasticBlockEdgePrior>(*this);
    }

    bool has_allocations() const override { return initialized_; }

    arma::ivec get_allocations() const override {
        return arma::conv_to<arma::ivec>::from(cluster_allocations_) + 1; // 1-based
    }

private:
    double beta_bernoulli_alpha_;
    double beta_bernoulli_beta_;
    double beta_bernoulli_alpha_between_;
    double beta_bernoulli_beta_between_;
    double dirichlet_alpha_;
    double lambda_;

    bool initialized_;
    arma::uvec cluster_allocations_;
    arma::mat cluster_prob_;
    arma::vec log_Vn_;
    SBMCorrection correction_;
};


/**
 * Attach a normalizing-constant correction list (assembled by
 * R/correction_tables.R) to a hierarchical edge prior. The beta-bernoulli
 * prior reads the whole-graph logC(theta) curve; the stochastic block prior
 * reads the slope and per-pair curves plus the optional continuous-block
 * node mask used on the mixed-MRF path.
 */
inline void attach_edge_prior_correction(
    BaseEdgePrior* edge_prior_obj,
    const Rcpp::Nullable<Rcpp::List>& edge_prior_correction,
    const char* caller
) {
    if (edge_prior_correction.isNull()) return;
    Rcpp::List correction(edge_prior_correction.get());
    if (auto* bb = dynamic_cast<BetaBernoulliEdgePrior*>(edge_prior_obj)) {
        bb->set_correction(EdgePriorCorrection(
            Rcpp::as<arma::vec>(correction["theta"]),
            Rcpp::as<arma::vec>(correction["logC"])
        ));
    } else if (auto* sbm = dynamic_cast<StochasticBlockEdgePrior*>(edge_prior_obj)) {
        arma::uvec is_continuous;
        if (correction.containsElementNamed("is_continuous")) {
            is_continuous = Rcpp::as<arma::uvec>(correction["is_continuous"]);
        }
        sbm->set_correction(SBMCorrection(
            Rcpp::as<arma::vec>(correction["fprime_density"]),
            Rcpp::as<arma::vec>(correction["fprime"]),
            Rcpp::as<arma::vec>(correction["quad_theta"]),
            Rcpp::as<arma::vec>(correction["quad_f"]),
            is_continuous
        ));
    } else {
        Rcpp::stop("%s: a correction table was supplied for an edge prior "
                   "that does not support it.", caller);
    }
}


/**
 * Factory: create an edge prior from an EdgePrior enum and hyperparameters.
 */
inline std::unique_ptr<BaseEdgePrior> create_edge_prior(
    EdgePrior type,
    double beta_bernoulli_alpha = 1.0,
    double beta_bernoulli_beta = 1.0,
    double beta_bernoulli_alpha_between = 1.0,
    double beta_bernoulli_beta_between = 1.0,
    double dirichlet_alpha = 1.0,
    double lambda = 1.0
) {
    switch (type) {
    case Beta_Bernoulli:
        return std::make_unique<BetaBernoulliEdgePrior>(
            beta_bernoulli_alpha, beta_bernoulli_beta);
    case Stochastic_Block:
        return std::make_unique<StochasticBlockEdgePrior>(
            beta_bernoulli_alpha, beta_bernoulli_beta,
            beta_bernoulli_alpha_between, beta_bernoulli_beta_between,
            dirichlet_alpha, lambda);
    case Bernoulli:
    case Not_Applicable:
    default:
        return std::make_unique<BernoulliEdgePrior>();
    }
}
