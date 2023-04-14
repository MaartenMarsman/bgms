// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the threshold parameters
// ----------------------------------------------------------------------------|
NumericMatrix metropolis_thresholds(NumericMatrix interactions,
                                    NumericMatrix thresholds,
                                    IntegerMatrix observations,
                                    IntegerVector no_categories,
                                    IntegerMatrix n_cat_obs,
                                    int no_persons,
                                    int no_nodes,
                                    double threshold_alpha,
                                    double threshold_beta,
                                    NumericMatrix rest_matrix) {
  NumericVector g(no_persons);
  NumericVector q(no_persons);

  double log_prob, rest_score;
  double a, b, c;
  double tmp;
  double current_state, proposed_state;
  double U;
  double exp_current, exp_proposed;

  for(int node = 0; node < no_nodes; node++) {
    for(int category = 0; category < no_categories[node]; category++) {
      current_state = thresholds(node, category);
      exp_current = std::exp(current_state);
      c = (threshold_alpha + threshold_beta) / (1 + exp_current);
      for(int person = 0; person < no_persons; person++) {
        g[person] = 1.0;
        q[person] = 1.0;
        rest_score = rest_matrix(person, node);
        for(int cat = 0; cat < no_categories[node]; cat++) {
          if(cat != category) {
            g[person] += std::exp(thresholds(node, cat) +
              (cat + 1) * rest_score);
          }
        }
        q[person] = std::exp((category + 1) * rest_score);
        c +=  q[person] / (g[person] + q[person] * exp_current);
      }
      c = c / ((no_persons + threshold_alpha + threshold_beta) -
        exp_current * c);

      //Proposal is generalized beta-prime.
      a = n_cat_obs(category + 1, node) + threshold_alpha;
      b = no_persons + threshold_beta - n_cat_obs(category + 1, node);
      tmp = R::rbeta(a, b);
      proposed_state = std::log(tmp / (1  - tmp) / c);
      exp_proposed = exp(proposed_state);

      //Compute log_acceptance probability for Metropolis.
      //First, we use g and q above to compute the ratio of pseudolikelihoods
      log_prob = 0;
      for(int person = 0; person < no_persons; person++) {
        log_prob += std::log(g[person] + q[person] * exp_current);
        log_prob -= std::log(g[person] + q[person] * exp_proposed);
      }
      //Second, we add the ratio of prior probabilities
      log_prob -= (threshold_alpha + threshold_beta) *
        std::log(1 + exp_proposed);
      log_prob += (threshold_alpha + threshold_beta) *
        std::log(1 + exp_current);
      //Third, we add the ratio of proposals
      log_prob -= (a + b) * std::log(1 + c * exp_current);
      log_prob += (a + b) * std::log(1 + c * exp_proposed);

      U = std::log(R::unif_rand());
      if(U < log_prob) {
        thresholds(node, category) = proposed_state;
      }
    }
  }
  return thresholds;
}


// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for an interaction
// ----------------------------------------------------------------------------|
double log_pseudolikelihood_ratio(NumericMatrix interactions,
                                  NumericMatrix thresholds,
                                  IntegerMatrix observations,
                                  IntegerVector no_categories,
                                  int no_persons,
                                  int node1,
                                  int node2,
                                  double proposed_state,
                                  double current_state,
                                  NumericMatrix rest_matrix) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  double delta_state = proposed_state - current_state;

  for(int person = 0; person < no_persons; person++) {
    obs_score1 = observations(person, node1);
    obs_score2 = observations(person, node2);

    pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * delta_state;

    //Node 1 log pseudolikelihood ratio
    rest_score = rest_matrix(person, node1) -
      obs_score2 * interactions(node2, node1);

    if(rest_score > 0) {
      bound = no_categories[node1] * rest_score;
    } else {
      bound = 0.0;
    }

    denominator_prop = std::exp(-bound);
    denominator_curr = std::exp(-bound);
    for(int category = 0; category < no_categories[node1]; category++) {
      score = category + 1;
      exponent = thresholds(node1, category) +
        score * rest_score -
        bound;
      denominator_prop +=
        std::exp(exponent + score * obs_score2 * proposed_state);
      denominator_curr +=
        std::exp(exponent + score * obs_score2 * current_state);
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);

    //Node 2 log pseudolikelihood ratio
    rest_score = rest_matrix(person, node2) -
      obs_score1 * interactions(node1, node2);

    if(rest_score > 0) {
      bound = no_categories[node2] * rest_score;
    } else {
      bound = 0.0;
    }

    denominator_prop = std::exp(-bound);
    denominator_curr = std::exp(-bound);
    for(int category = 0; category < no_categories[node2]; category++) {
      score = category + 1;
      exponent = thresholds(node2, category) +
        score * rest_score -
        bound;
      denominator_prop +=
        std::exp(exponent + score * obs_score1 * proposed_state);
      denominator_curr +=
        std::exp(exponent + score * obs_score1 * current_state);
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);
  }
  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of the active interaction
//  parameters (using a cauchy prior)
// ----------------------------------------------------------------------------|
List metropolis_interactions_cauchy(NumericMatrix interactions,
                                    NumericMatrix thresholds,
                                    IntegerMatrix gamma,
                                    IntegerMatrix observations,
                                    IntegerVector no_categories,
                                    NumericMatrix proposal_sd,
                                    double cauchy_scale,
                                    int no_persons,
                                    int no_nodes,
                                    NumericMatrix rest_matrix) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  for(int node1 = 0; node1 <  no_nodes - 1; node1++) {
    for(int node2 = node1 + 1; node2 <  no_nodes; node2++)
      if(gamma(node1, node2) == 1) {
        current_state = interactions(node1, node2);
        proposed_state = R::rnorm(current_state, proposal_sd(node1, node2));

        log_prob = log_pseudolikelihood_ratio(interactions,
                                              thresholds,
                                              observations,
                                              no_categories,
                                              no_persons,
                                              node1,
                                              node2,
                                              proposed_state,
                                              current_state,
                                              rest_matrix);
        log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
        log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale, true);

        U = R::unif_rand();
        if(std::log(U) < log_prob) {
          double state_diff = proposed_state - current_state;
          interactions(node1, node2) = proposed_state;
          interactions(node2, node1) = proposed_state;

          //Update the matrix of rest scores
          for(int person = 0; person < no_persons; person++) {
            rest_matrix(person, node1) += observations(person, node2) *
              state_diff;
            rest_matrix(person, node2) += observations(person, node1) *
              state_diff;
          }
        }
      }
  }
  return List::create(Named("interactions") = interactions,
                      Named("rest_matrix") = rest_matrix);
}


// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of the active interaction
//  parameters (using a unit information prior)
// ----------------------------------------------------------------------------|
List metropolis_interactions_unitinfo(NumericMatrix interactions,
                                      NumericMatrix thresholds,
                                      IntegerMatrix gamma,
                                      IntegerMatrix observations,
                                      IntegerVector no_categories,
                                      NumericMatrix proposal_sd,
                                      NumericMatrix unit_info,
                                      int no_persons,
                                      int no_nodes,
                                      NumericMatrix rest_matrix) {
  //NumericMatrix theta // no_nodes x no_nodes
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  for(int node1 = 0; node1 <  no_nodes - 1; node1++) {
    for(int node2 = node1 + 1; node2 <  no_nodes; node2++)
      if(gamma(node1, node2) == 1) {
        current_state = interactions(node1, node2);
        proposed_state = R::rnorm(current_state, proposal_sd(node1, node2));

        log_prob = log_pseudolikelihood_ratio(interactions,
                                              thresholds,
                                              observations,
                                              no_categories,
                                              no_persons,
                                              node1,
                                              node2,
                                              proposed_state,
                                              current_state,
                                              rest_matrix);

        log_prob += R::dnorm(proposed_state,
                             0.0,
                             unit_info(node1, node2),
                             true);
        log_prob -= R::dnorm(current_state,
                             0.0,
                             unit_info(node1, node2),
                             true);

        U = R::unif_rand();
        if(std::log(U) < log_prob) {
          double state_diff = proposed_state - current_state;
          interactions(node1, node2) = proposed_state;
          interactions(node2, node1) = proposed_state;

          //Update the matrix of rest scores
          for(int person = 0; person < no_persons; person++) {
            rest_matrix(person, node1) += observations(person, node2) *
              state_diff;
            rest_matrix(person, node2) += observations(person, node1) *
              state_diff;
          }
        }
      }
  }
  return List::create(Named("interactions") = interactions,
                      Named("rest_matrix") = rest_matrix);
}
