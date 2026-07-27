#pragma once

#include <memory>
#include <string>
#include <cmath>
#include <Rmath.h>
#include "math/explog_macros.h"


/**
 * BaseParameterPrior — abstract base for real-valued parameter priors.
 *
 * Follows the same polymorphic pattern as BaseEdgePrior. Each subclass stores
 * its own hyperparameters and provides logp/grad evaluated at a point x.
 *
 * Used for interaction parameters, threshold parameters, and continuous means.
 */
class BaseParameterPrior {
public:
    virtual ~BaseParameterPrior() = default;

    /** Log-density log p(x) up to an additive constant. */
    virtual double logp(double x) const = 0;

    /** Gradient d/dx log p(x). */
    virtual double grad(double x) const = 0;

    /** Deep copy for parallel chains. */
    virtual std::unique_ptr<BaseParameterPrior> clone() const = 0;

    /**
     * Current slab scale, for priors that carry one (Normal, Cauchy).
     *
     * @return the scale; throws for families with no scale to report.
     */
    virtual double scale() const {
        Rf_error("This parameter prior family has no scale.");
        return 0.0; // unreachable
    }

    /**
     * Set the slab scale in place, for priors that carry one.
     *
     * The random interaction-scale hyperprior mutates the owned prior object
     * each iteration so that all downstream logp/grad and edge birth-death
     * evaluations see the current scale without any call-site change. Families
     * with no scale throw, since randomizing them is unsupported.
     *
     * @param scale  New scale; must be positive.
     */
    virtual void set_scale(double scale) {
        (void) scale;
        Rf_error("This parameter prior family has no settable scale.");
    }
};


/**
 * CauchyPrior — Cauchy(0, scale).
 */
class CauchyPrior final : public BaseParameterPrior {
public:
    explicit CauchyPrior(double scale) : scale_(scale) {}

    double logp(double x) const override {
        return R::dcauchy(x, 0.0, scale_, true);
    }

    double grad(double x) const override {
        double s2 = scale_ * scale_;
        return -2.0 * x / (s2 + x * x);
    }

    std::unique_ptr<BaseParameterPrior> clone() const override {
        return std::make_unique<CauchyPrior>(*this);
    }

    /** @return the slab scale. */
    double scale() const override { return scale_; }

    /** Set the slab scale in place. @param scale New scale; must be positive. */
    void set_scale(double scale) override { scale_ = scale; }

private:
    double scale_;
};


/**
 * NormalPrior — Normal(0, scale).
 */
class NormalPrior final : public BaseParameterPrior {
public:
    explicit NormalPrior(double scale) : scale_(scale) {}

    double logp(double x) const override {
        return R::dnorm(x, 0.0, scale_, true);
    }

    double grad(double x) const override {
        double s2 = scale_ * scale_;
        return -x / s2;
    }

    std::unique_ptr<BaseParameterPrior> clone() const override {
        return std::make_unique<NormalPrior>(*this);
    }

    /** @return the slab standard deviation. */
    double scale() const override { return scale_; }

    /** Set the slab standard deviation in place. @param scale New scale; must be positive. */
    void set_scale(double scale) override { scale_ = scale; }

private:
    double scale_;
};


/**
 * BetaPrimePrior — logit-Beta(alpha, beta) prior.
 *
 * If sigma(x) ~ Beta(alpha, beta), then x = logit(Y) where Y ~ Beta(a, b).
 * log p(x) = alpha * x - (alpha + beta) * log(1 + exp(x)) + const
 */
class BetaPrimePrior final : public BaseParameterPrior {
public:
    BetaPrimePrior(double alpha, double beta)
        : alpha_(alpha), beta_(beta) {}

    double logp(double x) const override {
        return x * alpha_ - MY_LOG1P(MY_EXP(x)) * (alpha_ + beta_);
    }

    double grad(double x) const override {
        // alpha - (alpha + beta) * sigmoid(x)
        double p = 1.0 / (1.0 + MY_EXP(-x));
        return alpha_ - (alpha_ + beta_) * p;
    }

    std::unique_ptr<BaseParameterPrior> clone() const override {
        return std::make_unique<BetaPrimePrior>(*this);
    }

private:
    double alpha_;
    double beta_;
};


/**
 * GammaScalePrior — Gamma(shape, rate) prior for positive parameters.
 *
 * Used for precision matrix diagonal elements.
 * log p(x) = (shape - 1) * log(x) - rate * x + const
 */
class GammaScalePrior final : public BaseParameterPrior {
public:
    GammaScalePrior(double shape, double rate)
        : shape_(shape), rate_(rate) {}

    double logp(double x) const override {
        return R::dgamma(x, shape_, 1.0 / rate_, true);
    }

    double grad(double x) const override {
        // d/dx log Gamma(x; shape, rate) = (shape - 1) / x - rate
        return (shape_ - 1.0) / x - rate_;
    }

    std::unique_ptr<BaseParameterPrior> clone() const override {
        return std::make_unique<GammaScalePrior>(*this);
    }

    /** @return the Gamma shape parameter. */
    double shape() const { return shape_; }
    /** @return the Gamma rate parameter. */
    double rate() const { return rate_; }

private:
    double shape_;
    double rate_;
};


// =============================================================================
// Factory functions
// =============================================================================

/**
 * Create a parameter prior from a type string and hyperparameters.
 *
 * @param type    One of "cauchy", "normal", "beta-prime"
 * @param scale   Scale for Cauchy/Normal (ignored for beta-prime)
 * @param alpha   Alpha for beta-prime (ignored for Cauchy/Normal)
 * @param beta    Beta for beta-prime (ignored for Cauchy/Normal)
 */
inline std::unique_ptr<BaseParameterPrior> create_parameter_prior(
    const std::string& type,
    double scale = 1.0,
    double alpha = 0.5,
    double beta = 0.5
) {
    if (type == "cauchy") {
        return std::make_unique<CauchyPrior>(scale);
    } else if (type == "normal") {
        return std::make_unique<NormalPrior>(scale);
    } else if (type == "beta-prime") {
        return std::make_unique<BetaPrimePrior>(alpha, beta);
    }
    Rf_error("Unknown parameter prior type: '%s'", type.c_str());
    return nullptr; // unreachable
}


/**
 * Create a scale prior from a type string and hyperparameters.
 *
 * @param type   One of "gamma", "exponential"
 * @param shape  Shape for Gamma (ignored for exponential, set to 1)
 * @param rate   Rate for Gamma/Exponential
 */
inline std::unique_ptr<BaseParameterPrior> create_scale_prior(
    const std::string& type,
    double shape = 1.0,
    double rate = 1.0
) {
    if (type == "gamma" || type == "exponential") {
        double s = (type == "exponential") ? 1.0 : shape;
        return std::make_unique<GammaScalePrior>(s, rate);
    }
    Rf_error("Unknown scale prior type: '%s'", type.c_str());
    return nullptr; // unreachable
}
