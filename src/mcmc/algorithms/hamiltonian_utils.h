#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include <utility>
#include "mcmc/execution/step_result.h"
#include "mcmc/algorithms/leapfrog.h"
struct SafeRNG;


/**
 * Kinetic energy for Hamiltonian Monte Carlo
 *
 * Computes 0.5 * r^T * M^{-1} * r where M^{-1} is a diagonal mass matrix.
 *
 * @param r              Momentum vector
 * @param inv_mass_diag  Diagonal of the inverse mass matrix
 * @return Scalar kinetic energy
 */
double kinetic_energy(const arma::vec& r, const arma::vec& inv_mass_diag);


/**
 * Heuristic initial step size for HMC/NUTS (identity mass)
 *
 * Iteratively doubles or halves a candidate step size until a single leapfrog
 * step yields an acceptance probability near the target. Delegates to the
 * mass-matrix overload with inv_mass_diag = ones.
 *
 * @param theta             Initial parameter vector
 * @param joint             Joint log-posterior + gradient function
 * @param rng               Random number generator
 * @param target_acceptance Target acceptance probability
 * @param init_step         Starting step size
 * @param max_attempts      Maximum doubling/halving iterations
 * @return Step size yielding acceptance probability near target
 */
double heuristic_initial_step_size(
    const arma::vec& theta,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    SafeRNG& rng,
    double target_acceptance = 0.8,
    double init_step = 1.0,
    int max_attempts = 20
);


/**
 * Heuristic initial step size for HMC/NUTS (with mass matrix)
 *
 * Same algorithm as the identity-mass overload, but samples momentum from
 * N(0, M) and evaluates kinetic energy with the supplied diagonal M^{-1}.
 *
 * @param theta             Initial parameter vector
 * @param joint             Joint log-posterior + gradient function
 * @param inv_mass_diag     Diagonal of the inverse mass matrix
 * @param rng               Random number generator
 * @param target_acceptance Target acceptance probability
 * @param init_step         Starting step size
 * @param max_attempts      Maximum doubling/halving iterations
 * @return Step size yielding acceptance probability near target
 */
double heuristic_initial_step_size(
    const arma::vec& theta,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng,
    double target_acceptance = 0.8,
    double init_step = 1.0,
    int max_attempts = 20
);
