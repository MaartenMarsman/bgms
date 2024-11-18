#include <RcppArmadillo.h>
#include <cmath>
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

inline bool stoppingCriterion(const arma::vec& batch_means,
                              double delta,
                              int batch_num,
                              double cumulative_batch_sum) {
  if (batch_num > 1) {  // Ensure at least 2 batches to calculate variance
    double variance = arma::var(batch_means.head(batch_num));
    double standard_error = std::sqrt(variance / batch_num);
    double running_mean = cumulative_batch_sum / batch_num;

    // Stop if the standard error is less than delta * |running mean|
    return (standard_error < delta * std::abs(running_mean));
  }
  return false;
}


// [[Rcpp::export]]
double sample_moment(const arma::mat& theta,
                     int index_1,
                     int index_2,
                     double atheta,
                     int batch_size,
                     double delta,
                     int max_num_batches) {

  int no_variables = theta.n_cols;
  arma::vec x(no_variables, arma::fill::zeros);
  arma::vec chain(batch_size);
  arma::vec batch_means(max_num_batches);
  double cumulative_batch_sum = 0.0;

  int i = index_1 - 1;
  int j = index_2 - 1;
  double diff = atheta - theta(i, j);
  arma::vec exp_diff(2);
  exp_diff[0] = 1.0;
  exp_diff[1] = std::exp(diff);
  //if(i != j)
  //  exp_diff[1] = std::exp(2 * diff);

  // Random starting values for the chain
  for (int variable = 0; variable < no_variables; variable++) {
    double u = R::unif_rand();
    x[variable] = (u > 0.5) ? 1.0 : 0.0;
  }

  // MCMC Sampling Loop
  for (int batch = 0; batch < max_num_batches; ++batch) {
    Rcpp::checkUserInterrupt();

    // Run a batch of MCMC samples
    for (int iter = 0; iter < batch_size; ++iter) {
      for (int variable = 0; variable < no_variables; variable++) {
        //Set x[variable] to 1 for efficient computation of dot product.
        x[variable] = 1;

        // Compute dot product with Armadillo
        double exponent = arma::dot(x, theta.col(variable));

        double z = 1.0 + std::exp(exponent);
        double u = z * R::unif_rand();

        if (u < 1.0) {
          //x[variable] was set to 1, needs to be 0 if u < 1.0.
          x[variable] = 0;
        }
      }

      // Calculate the chain element based on the updated x vector
      chain[iter] = exp_diff[x[i] * x[j]];  // Lookup precomputed values
    }

    // Calculate the mean of this batch and add it to cumulative sum
    double current_batch_mean = arma::accu(chain) / batch_size;
    batch_means[batch] = current_batch_mean;
    cumulative_batch_sum += current_batch_mean;

    // Check stopping criterion
    if (stoppingCriterion(batch_means, delta, batch + 1, cumulative_batch_sum)) {
      return cumulative_batch_sum / (batch + 1);  // Mean of batch means so far
    }
  }

  // Return the mean of all batch means if stopping criterion not met
  return cumulative_batch_sum / max_num_batches;
}


// [[Rcpp::export]]
double raoblackwellized_sample_moment(const arma::mat& theta,
                                      int index_1,
                                      int index_2,
                                      double atheta,
                                      int batch_size,
                                      double delta,
                                      int max_num_batches) {

  int no_variables = theta.n_cols;
  arma::vec x(no_variables, arma::fill::zeros);
  arma::vec chain(batch_size);
  arma::vec batch_means(max_num_batches);
  double cumulative_batch_sum = 0.0;

  double exponent, z, u, e;

  int i = index_1 - 1;
  int j = index_2 - 1;
  int x_old_i, x_old_j;
  double diff = atheta - theta(i, j);
  double exp_diff = std::exp(diff);

  // Random starting values for the chain
  for (int variable = 0; variable < no_variables; variable++) {
    double u = R::unif_rand();
    x[variable] = (u > 0.5) ? 1.0 : 0.0;
  }

  // MCMC Sampling Loop
  for (int batch = 0; batch < max_num_batches; ++batch) {
    Rcpp::checkUserInterrupt();

    // Run a batch of MCMC samples
    for (int iter = 0; iter < batch_size; ++iter) {
      for (int variable = 0; variable < no_variables; variable++) {
        //Set x[variable] to 1 for efficient computation of dot product.
        x[variable] = 1;

        // Compute dot product with Armadillo
        exponent = arma::dot(x, theta.col(variable));

        z = 1.0 + std::exp(exponent);
        u = z * R::unif_rand();

        if (u < 1.0) {
          //x[variable] was set to 1, needs to be 0 if u < 1.0.
          x[variable] = 0;
        }
      }

      // Calculate the chain element based on the updated x vector
      if(i == j) {
        //Rao-Blackwellization.
        x_old_i = x[i];
        x[i] = 1;

        // Compute dot product with Armadillo
        exponent = arma::dot(x, theta.col(i));
        e = std::exp(exponent);
        z = 1.0 + e;
        double p0 = 1.0 / z;
        double p1 = e / z;
        chain[iter] = p0 + exp_diff * p1;

        x[i] = x_old_i;
      } else {
        x_old_i = x[i];
        x_old_j = x[j];

        //x_i = 0, x_j = 0
        double p00 = 1.0;
        z = 1.0;

        //x_i = 0, x_j = 1
        x[i] = 0;
        x[j] = 1;
        double e01 = arma::dot(x, theta.col(j)) ;
        double p01 = std::exp(e01);
        z += p01;

        //x_i = 1, x_j = 0
        x[i] = 1;
        x[j] = 0;
        double e10 = arma::dot(x, theta.col(i)) ;
        double p10 = std::exp(e10);
        z += p10;

        //x_i = 1, x_j = 1
        x[i] = 1;
        x[j] = 1;
        exponent = e10 + e01 + theta(i,j);
        double p11 = std::exp(exponent);
        z += p11;

        //Rao-Blackwellization.
        e = p00 / z;
        e += p01 / z;
        e += p10 / z;
        e += p11 / z * exp_diff;

        chain[iter] = e;

        x[i] = x_old_i;
        x[j] = x_old_j;
      }
    }

    // Calculate the mean of this batch and add it to cumulative sum
    double current_batch_mean = arma::accu(chain) / batch_size;
    batch_means[batch] = current_batch_mean;
    cumulative_batch_sum += current_batch_mean;

    // Check stopping criterion
    if (stoppingCriterion(batch_means, delta, batch + 1, cumulative_batch_sum)) {
      return cumulative_batch_sum / (batch + 1);  // Mean of batch means so far
    }
  }

  // Return the mean of all batch means if stopping criterion not met
  return cumulative_batch_sum / max_num_batches;
}


// [[Rcpp::export]]
double taylor_raoblackwellized_sample_moment(const arma::mat& theta,
                                             int index_1,
                                             int index_2,
                                             int batch_size,
                                             double delta,
                                             int max_num_batches) {

  int no_variables = theta.n_cols;
  arma::vec x(no_variables, arma::fill::zeros);
  arma::vec chain(batch_size);
  arma::vec batch_means(max_num_batches);
  double cumulative_batch_sum = 0.0;

  double exponent, z, u, e;

  int i = index_1 - 1;
  int j = index_2 - 1;
  int x_old_i, x_old_j;

  // Random starting values for the chain
  for (int variable = 0; variable < no_variables; variable++) {
    double u = R::unif_rand();
    x[variable] = (u > 0.5) ? 1.0 : 0.0;
  }

  // MCMC Sampling Loop
  for (int batch = 0; batch < max_num_batches; ++batch) {
    Rcpp::checkUserInterrupt();

    // Run a batch of MCMC samples
    for (int iter = 0; iter < batch_size; ++iter) {
      for (int variable = 0; variable < no_variables; variable++) {
        //Set x[variable] to 1 for efficient computation of dot product.
        x[variable] = 1;

        // Compute dot product with Armadillo
        exponent = arma::dot(x, theta.col(variable));

        z = 1.0 + std::exp(exponent);
        u = z * R::unif_rand();

        if (u < 1.0) {
          //x[variable] was set to 1, needs to be 0 if u < 1.0.
          x[variable] = 0;
        }
      }

      // Calculate the chain element based on the updated x vector
      if(i == j) {
        //Rao-Blackwellization.
        x_old_i = x[i];
        x[i] = 1;

        // Compute dot product with Armadillo
        exponent = arma::dot(x, theta.col(i));
        e = std::exp(exponent);
        chain[iter] = e / (1.0 + e);

        x[i] = x_old_i;
      } else {
        x_old_i = x[i];
        x_old_j = x[j];

        //x_i = 0, x_j = 0
        z = 1.0;

        //x_i = 0, x_j = 1
        x[i] = 0;
        x[j] = 1;
        double e01 = arma::dot(x, theta.col(j)) ;
        double p01 = std::exp(e01);
        z += p01;

        //x_i = 1, x_j = 0
        x[i] = 1;
        x[j] = 0;
        double e10 = arma::dot(x, theta.col(i)) ;
        double p10 = std::exp(e10);
        z += p10;

        //x_i = 1, x_j = 1
        x[i] = 1;
        x[j] = 1;
        exponent = e10 + e01 + theta(i,j);
        double p11 = std::exp(exponent);
        z += p11;

        //Rao-Blackwellization.
        e = p11 / z;

        chain[iter] = e;

        x[i] = x_old_i;
        x[j] = x_old_j;
      }
    }

    // Calculate the mean of this batch and add it to cumulative sum
    double current_batch_mean = arma::accu(chain) / batch_size;
    batch_means[batch] = current_batch_mean;
    cumulative_batch_sum += current_batch_mean;

    // Check stopping criterion
    if (stoppingCriterion(batch_means, delta, batch + 1, cumulative_batch_sum)) {
      return cumulative_batch_sum / (batch + 1);  // Mean of batch means so far
    }
  }

  // Return the mean of all batch means if stopping criterion not met
  return cumulative_batch_sum / max_num_batches;
}

// [[Rcpp::export]]
Rcpp::List mcmc_mh(arma::mat x,
                   int n_iter,
                   int n_burnin,
                   double target_moms,
                   double target_mh,
                   double epsilon) {
  int n = x.n_rows;
  int p = x.n_cols;

  arma::mat s = x.t() * x;
  arma::mat gamma(p, p, arma::fill::ones);
  arma::mat theta(p, p, arma::fill::zeros);
  arma::mat sd_theta(p, p, arma::fill::ones);
  arma::mat eap_theta = theta;
  arma::mat eap_gamma = gamma;

  double prop_rej_moms = 0;
  double prop_rej_mh = 0;
  double delta_moms = 0.01;
  double delta_mh = 0.01;

  double a = 0.5, b = 0.5;
  int tel = 0;

  Progress q(n_iter + n_burnin, true);

  //burnin
  for (int iter = 1; iter <= n_burnin; ++iter) {
    q.increment();

    prop_rej_moms = 0;
    prop_rej_mh = 0;
    tel = 0;

    for (int i = 0; i < p; ++i) {
      for (int j = i; j < p; ++j) {
        // Between Model Move
        if (i < j) {
          double agamma = 1 - gamma(i, j);
          double atheta = (agamma == 1) ? R::rnorm(theta(i, j), sd_theta(i, j)) : 0.0;

          arma::mat theta_prime = theta;
          theta_prime(i, j) = atheta;
          theta_prime(j, i) = atheta;

          double E1 = raoblackwellized_sample_moment(theta_prime, i + 1, j + 1, theta(i, j), 50, delta_moms, 1e5);
          double E2 = raoblackwellized_sample_moment(theta, i + 1, j + 1, atheta, 50, delta_moms, 1e5);

          if (std::abs(1 - E1 * E2) < epsilon) {
            double lnE = 0.5 * std::log(E1) - 0.5 * std::log(E2);
            double ln_prob = s(i, j) * (atheta - theta(i, j)) + n * lnE;

            if (agamma == 1) {
              ln_prob += R::dcauchy(atheta, 0, 2.5, true) - R::dnorm(atheta, theta(i, j), sd_theta(i, j), true);
            } else {
              ln_prob -= R::dcauchy(theta(i, j), 0, 2.5, true) - R::dnorm(theta(i, j), atheta, sd_theta(i, j), true);
            }

            if (std::log(R::runif(0, 1)) < ln_prob) {
              theta(i, j) = atheta;
              theta(j, i) = atheta;
              gamma(i, j) = agamma;
              gamma(j, i) = agamma;
            }
          } else {
            prop_rej_moms += 1;
          }
        }

        // Within Model Move
        if (gamma(i, j) == 1) {
          tel += 1;
          double atheta = R::rnorm(theta(i, j), sd_theta(i, j));

          arma::mat theta_prime = theta;
          theta_prime(i, j) = atheta;
          theta_prime(j, i) = atheta;

          double E1 = raoblackwellized_sample_moment(theta_prime, i + 1, j + 1, theta(i, j), 50, delta_mh, 1e5);
          double E2 = raoblackwellized_sample_moment(theta, i + 1, j + 1, atheta, 50, delta_mh, 1e5);

          if (std::abs(1 - E1 * E2) < epsilon) {
            double lnE = 0.5 * std::log(E1) - 0.5 * std::log(E2);
            double ln_prob = s(i, j) * (atheta - theta(i, j)) + n * lnE;

            if (i != j) {
              ln_prob += R::dcauchy(atheta, 0, 2.5, true) - R::dcauchy(theta(i, j), 0, 2.5, true);
            } else {
              ln_prob += a * atheta - (a + b) * std::log(1 + std::exp(atheta));
              ln_prob -= a * theta(i, j) - (a + b) * std::log(1 + std::exp(theta(i, j)));
            }

            if (std::log(R::runif(0, 1)) < ln_prob) {
              theta(i, j) = atheta;
              theta(j, i) = atheta;
            }

            double prob_accept = (ln_prob > 0) ? 1.0 : std::exp(ln_prob);
            sd_theta(i, j) += (prob_accept - 0.234) * std::exp(-std::log(iter) * 0.75);
            sd_theta(i, j) = std::clamp(sd_theta(i, j), 1.0 / n, 20.0);
            sd_theta(j, i) = sd_theta(i, j);
          } else {
            prop_rej_mh += 1;
          }
        }
      }
    }

    prop_rej_moms /= (p * (p - 1) / 2);
    prop_rej_mh /= tel;

    // Adapt hyperparameters
    delta_moms = std::exp(std::log(delta_moms) + std::exp(-std::log(iter) * 0.51) * (target_moms - prop_rej_moms));
    delta_moms = std::clamp(delta_moms, 0.001, 0.5);

    delta_mh = std::exp(std::log(delta_mh) + std::exp(-std::log(iter) * 0.51) * (target_mh - prop_rej_mh));
    delta_mh = std::clamp(delta_mh, 0.001, 0.5);
  }

  //The Gibbs sampler
  for (int iter = 1; iter <= n_iter; ++iter) {
    q.increment();

    prop_rej_moms = 0;
    prop_rej_mh = 0;
    tel = 0;

    for (int i = 0; i < p; ++i) {
      for (int j = i; j < p; ++j) {
        // Between Model Move
        if (i < j) {
          double agamma = 1 - gamma(i, j);
          double atheta = (agamma == 1) ? R::rnorm(theta(i, j), sd_theta(i, j)) : 0.0;

          arma::mat theta_prime = theta;
          theta_prime(i, j) = atheta;
          theta_prime(j, i) = atheta;

          double E1 = raoblackwellized_sample_moment(theta_prime, i + 1, j + 1, theta(i, j), 50, delta_moms, 1e5);
          double E2 = raoblackwellized_sample_moment(theta, i + 1, j + 1, atheta, 50, delta_moms, 1e5);

          if (std::abs(1 - E1 * E2) < epsilon) {
            double lnE = 0.5 * std::log(E1) - 0.5 * std::log(E2);
            double ln_prob = s(i, j) * (atheta - theta(i, j)) + n * lnE;

            if (agamma == 1) {
              ln_prob += R::dcauchy(atheta, 0, 2.5, true) - R::dnorm(atheta, theta(i, j), sd_theta(i, j), true);
            } else {
              ln_prob -= R::dcauchy(theta(i, j), 0, 2.5, true) - R::dnorm(theta(i, j), atheta, sd_theta(i, j), true);
            }

            if (std::log(R::runif(0, 1)) < ln_prob) {
              theta(i, j) = atheta;
              theta(j, i) = atheta;
              gamma(i, j) = agamma;
              gamma(j, i) = agamma;
            }
          } else {
            prop_rej_moms += 1;
          }
        }

        // Within Model Move
        if (gamma(i, j) == 1) {
          tel += 1;
          double atheta = R::rnorm(theta(i, j), sd_theta(i, j));

          arma::mat theta_prime = theta;
          theta_prime(i, j) = atheta;
          theta_prime(j, i) = atheta;

          double E1 = raoblackwellized_sample_moment(theta_prime, i + 1, j + 1, theta(i, j), 50, delta_mh, 1e5);
          double E2 = raoblackwellized_sample_moment(theta, i + 1, j + 1, atheta, 50, delta_mh, 1e5);

          if (std::abs(1 - E1 * E2) < epsilon) {
            double lnE = 0.5 * std::log(E1) - 0.5 * std::log(E2);
            double ln_prob = s(i, j) * (atheta - theta(i, j)) + n * lnE;

            if (i != j) {
              ln_prob += R::dcauchy(atheta, 0, 2.5, true) - R::dcauchy(theta(i, j), 0, 2.5, true);
            } else {
              ln_prob += a * atheta - (a + b) * std::log(1 + std::exp(atheta));
              ln_prob -= a * theta(i, j) - (a + b) * std::log(1 + std::exp(theta(i, j)));
            }

            if (std::log(R::runif(0, 1)) < ln_prob) {
              theta(i, j) = atheta;
              theta(j, i) = atheta;
            }

            double prob_accept = (ln_prob > 0) ? 1.0 : std::exp(ln_prob);
            sd_theta(i, j) += (prob_accept - 0.234) * std::exp(-std::log(iter) * 0.75);
            sd_theta(i, j) = std::clamp(sd_theta(i, j), 1.0 / n, 20.0);
            sd_theta(j, i) = sd_theta(i, j);
          } else {
            prop_rej_mh += 1;
          }
        }
      }
    }

    // Update EAP estimates
    eap_theta = (iter - 1) * eap_theta / iter + theta / iter;
    eap_gamma = (iter - 1) * eap_gamma / iter + gamma / iter;

    prop_rej_moms /= (p * (p - 1) / 2);
    prop_rej_mh /= tel;

    // Adapt hyperparameters
    delta_moms = std::exp(std::log(delta_moms) + std::exp(-std::log(iter) * 0.51) * (0.05 - prop_rej_moms));
    delta_moms = std::clamp(delta_moms, 0.001, 0.5);

    delta_mh = std::exp(std::log(delta_mh) + std::exp(-std::log(iter) * 0.51) * (0.05 - prop_rej_mh));
    delta_mh = std::clamp(delta_mh, 0.001, 0.5);
  }

  return Rcpp::List::create(
    Rcpp::Named("eap_theta") = eap_theta,
    Rcpp::Named("eap_gamma") = eap_gamma
  );
}
