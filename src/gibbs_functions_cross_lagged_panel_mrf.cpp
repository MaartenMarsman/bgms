// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the threshold parameters
// in a ordinal cross-lagged panel mrf
// ----------------------------------------------------------------------------|
NumericMatrix metropolis_thresholds_cross_lagged_panel_mrf(NumericMatrix crsec_interactions,
                                                           NumericMatrix crlag_interactions,
                                                           NumericMatrix thresholds,
                                                           IntegerMatrix observations,
                                                           IntegerVector no_categories,
                                                           IntegerVector start,
                                                           IntegerMatrix n_cat_obs,
                                                           int no_persons,
                                                           int no_nodes,
                                                           int no_timepoints,
                                                           double threshold_alpha,
                                                           double threshold_beta) {
  NumericVector g(no_persons);
  NumericVector q(no_persons);

  double log_prob, rest_score;
  double a, b, c;
  double tmp;
  double current_state, proposed_state;
  double U;
  double exp_current, exp_proposed;

  for(int t = 1; t <= no_timepoints; t++) {
    for(int node = 0; node < no_nodes; node++) {
      for(int category = 0; category < no_categories[node]; category++) {
        current_state = thresholds(node + start[t - 1], category);
        exp_current = std::exp(current_state);
        c = (threshold_alpha + threshold_beta) / (1 + exp_current);
        for(int person = 0; person < no_persons; person++) {
          g[person] = 1.0;
          q[person] = 1.0;
          rest_score = 0.0;
          for(int node2 = 0; node2 < no_nodes; node2++) {
            rest_score += observations(person, node2 + start[t]) *
              crsec_interactions(node, node2);
            rest_score += observations(person, node2 + start[t - 1]) *
              crlag_interactions(node, node2);
          }
          for(int cat = 0; cat < no_categories[node]; cat++) {
            if(cat != category) {
              g[person] += std::exp(thresholds(node + start[t - 1], cat) +
                (cat + 1) * rest_score);
            }
          }
          q[person] = std::exp((category + 1) * rest_score);
          c +=  q[person] / (g[person] + q[person] * exp_current);
        }
        c = c / ((no_persons + threshold_alpha + threshold_beta) -
          exp_current * c);

        //Proposal is generalized beta-prime.
        a = n_cat_obs(category + 1, node + start[t - 1]) +
          threshold_alpha;
        b = no_persons +
          threshold_beta -
          n_cat_obs(category + 1, node + start[t - 1]);
        tmp = R::rbeta(a, b);
        proposed_state = std::log(tmp / (1 - tmp) / c);
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
        if(U < log_prob)
          thresholds(node + start[t - 1], category) = proposed_state;
      }
    }
  }

  return thresholds;
}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for a
//  cross-sectional interaction
// ----------------------------------------------------------------------------|
double log_pseudolikelihood_ratio_cross_sectional(NumericMatrix crsec_interactions,
                                                  NumericMatrix crlag_interactions,
                                                  NumericMatrix thresholds,
                                                  IntegerMatrix observations,
                                                  IntegerVector no_categories,
                                                  IntegerVector start,
                                                  int no_persons,
                                                  int no_nodes,
                                                  int no_timepoints,
                                                  int node1,
                                                  int node2,
                                                  double proposed_state,
                                                  double current_state) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  double diff_state = proposed_state - current_state;

  for(int t = 1; t <= no_timepoints; t++) {
    for(int person = 0; person < no_persons; person++) {
      obs_score1 = observations(person, node1 + start[t]);
      obs_score2 = observations(person, node2 + start[t]);

      pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * diff_state;

      //Node 1 log pseudolikelihood ratio
      rest_score = 0.0;
      for(int node = 0; node < no_nodes; node++) {
        rest_score += crsec_interactions(node1, node) *
          observations(person, node + start[t]);
        rest_score += crlag_interactions(node1, node) *
          observations(person, node + start[t - 1]);
      }
      rest_score -= obs_score2 *
        crsec_interactions(node1, node2);

      if(rest_score > 0) {
        bound = no_categories[node1] * rest_score;
      } else {
        bound = 0.0;
      }

      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories[node1]; category++) {
        score = category + 1;
        exponent = thresholds(node1 + start[t - 1], category) +
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
      rest_score = 0.0;
      for(int node = 0; node < no_nodes; node++) {
        rest_score += crsec_interactions(node2, node) *
          observations(person, node + start[t]);
        rest_score += crlag_interactions(node2, node) *
          observations(person, node + start[t - 1]);
      }
      rest_score -= obs_score1 *
        crsec_interactions(node1, node2);

      if(rest_score > 0) {
        bound = no_categories[node2] * rest_score;
      } else {
        bound = 0.0;
      }

      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories[node2]; category++) {
        score = category + 1;
        exponent = thresholds(node2 + start[t - 1], category) +
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
  }

  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for a
//  cross-lagged interaction
// ----------------------------------------------------------------------------|
double log_pseudolikelihood_ratio_cross_lagged(NumericMatrix crsec_interactions,
                                               NumericMatrix crlag_interactions,
                                               NumericMatrix thresholds,
                                               IntegerMatrix observations,
                                               IntegerVector no_categories,
                                               IntegerVector start,
                                               int no_persons,
                                               int no_nodes,
                                               int no_timepoints,
                                               int node1,
                                               int node2,
                                               double proposed_state,
                                               double current_state) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;
  double diff_state = proposed_state - current_state;

  for(int t = 1; t <= no_timepoints; t++) {
    for(int person = 0; person < no_persons; person++) {
      obs_score1 = observations(person, node1 + start[t]);
      obs_score2 = observations(person, node2 + start[t - 1]);

      pseudolikelihood_ratio += obs_score1 * obs_score2 * diff_state;

      //Node 1 log pseudolikelihood ratio
      rest_score = 0.0;
      for(int node = 0; node < no_nodes; node++) {
        rest_score += crsec_interactions(node1, node) *
          observations(person, node + start[t]);
        rest_score += crlag_interactions(node1, node) *
          observations(person, node + start[t - 1]);
      }
      rest_score -= obs_score2 *
        crlag_interactions(node1, node2);

      if(rest_score > 0) {
        bound = no_categories[node1] * rest_score;
      } else {
        bound = 0.0;
      }

      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories[node1]; category++) {
        score = category + 1;
        exponent = thresholds(node1 + start[t - 1], category) +
          score * rest_score -
          bound;
        denominator_prop +=
          std::exp(exponent + score * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + score * obs_score2 * current_state);
      }
      pseudolikelihood_ratio -= std::log(denominator_prop);
      pseudolikelihood_ratio += std::log(denominator_curr);
    }
  }

  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of the active cross-sectional
//  interaction parameters (using a cauchy prior)
// ----------------------------------------------------------------------------|
List metropolis_cross_sectional_interactions(NumericMatrix crsec_interactions,
                                             NumericMatrix crlag_interactions,
                                             NumericMatrix thresholds,
                                             IntegerMatrix gamma,
                                             IntegerMatrix observations,
                                             IntegerVector no_categories,
                                             IntegerVector start,
                                             NumericMatrix crsec_proposal_sd,
                                             double cauchy_scale,
                                             int no_persons,
                                             int no_nodes,
                                             int no_timepoints,
                                             double phi,
                                             double target_ar,
                                             int tt,
                                             double epsilon_lo,
                                             double epsilon_hi) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  for(int node1 = 0; node1 <  no_nodes - 1; node1++) {
    for(int node2 = node1 + 1; node2 <  no_nodes; node2++) {
      if(gamma(node1, node2) == 1) {
        current_state = crsec_interactions(node1, node2);
        proposed_state = R::rnorm(current_state,
                                  crsec_proposal_sd(node1, node2));

        log_prob = log_pseudolikelihood_ratio_cross_sectional(crsec_interactions,
                                                              crlag_interactions,
                                                              thresholds,
                                                              observations,
                                                              no_categories,
                                                              start,
                                                              no_persons,
                                                              no_nodes,
                                                              no_timepoints,
                                                              node1,
                                                              node2,
                                                              proposed_state,
                                                              current_state);

        log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
        log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale, true);

        U = R::unif_rand();
        if(std::log(U) < log_prob) {
          crsec_interactions(node1, node2) = proposed_state;
          crsec_interactions(node2, node1) = proposed_state;
        }

        if(log_prob > 0) {
          log_prob = 1;
        } else {
          log_prob = std::exp(log_prob);
        }
        crsec_proposal_sd(node1, node2) = crsec_proposal_sd(node1, node2) +
          (log_prob - target_ar) * std::exp(-log(tt) * phi);
        if(crsec_proposal_sd(node1, node2) < epsilon_lo) {
          crsec_proposal_sd(node1, node2) = epsilon_lo;
        } else if (crsec_proposal_sd(node1, node2) > epsilon_hi) {
          crsec_proposal_sd(node1, node2) = epsilon_hi;
        }
        crsec_proposal_sd(node2, node1) = crsec_proposal_sd(node1, node2);
      }
    }
  }

  return List::create(Named("crsec_interactions") = crsec_interactions,
                      Named("crsec_proposal_sd") = crsec_proposal_sd);
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of the active cross-lagged
//  interaction parameters (using a cauchy prior)
// ----------------------------------------------------------------------------|
List metropolis_cross_lagged_interactions(NumericMatrix crsec_interactions,
                                          NumericMatrix crlag_interactions,
                                          NumericMatrix thresholds,
                                          IntegerMatrix delta,
                                          IntegerMatrix observations,
                                          IntegerVector no_categories,
                                          IntegerVector start,
                                          NumericMatrix crlag_proposal_sd,
                                          double cauchy_scale,
                                          int no_persons,
                                          int no_nodes,
                                          int no_timepoints,
                                          double phi,
                                          double target_ar,
                                          int tt,
                                          double epsilon_lo,
                                          double epsilon_hi) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  for(int node1 = 0; node1 <  no_nodes; node1++) {
    for(int node2 = 0; node2 <  no_nodes; node2++) {
      if(delta(node1, node2) == 1) {
        current_state = crlag_interactions(node1, node2);
        proposed_state = R::rnorm(current_state,
                                  crlag_proposal_sd(node1, node2));

        log_prob = log_pseudolikelihood_ratio_cross_lagged(crsec_interactions,
                                                           crlag_interactions,
                                                           thresholds,
                                                           observations,
                                                           no_categories,
                                                           start,
                                                           no_persons,
                                                           no_nodes,
                                                           no_timepoints,
                                                           node1,
                                                           node2,
                                                           proposed_state,
                                                           current_state);

        log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
        log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale, true);

        U = R::unif_rand();
        if(std::log(U) < log_prob)
          crlag_interactions(node1, node2) = proposed_state;

        if(log_prob > 0) {
          log_prob = 1;
        } else {
          log_prob = std::exp(log_prob);
        }
        crlag_proposal_sd(node1, node2) = crlag_proposal_sd(node1, node2) +
          (log_prob - target_ar) * std::exp(-log(tt) * phi);
        if(crlag_proposal_sd(node1, node2) < epsilon_lo) {
          crlag_proposal_sd(node1, node2) = epsilon_lo;
        } else if (crlag_proposal_sd(node1, node2) > epsilon_hi) {
          crlag_proposal_sd(node1, node2) = epsilon_hi;
        }
      }
    }
  }

  return List::create(Named("crlag_interactions") = crlag_interactions,
                      Named("crlag_proposal_sd") = crlag_proposal_sd);
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of a cross-sectional edge +
//  interaction pair (using a cauchy prior)
// ----------------------------------------------------------------------------|
List metropolis_cross_sectional_edge_interaction_pair(NumericMatrix crsec_interactions,
                                                      NumericMatrix crlag_interactions,
                                                      NumericMatrix thresholds,
                                                      IntegerMatrix gamma,
                                                      IntegerMatrix observations,
                                                      IntegerVector no_categories,
                                                      IntegerVector start,
                                                      NumericMatrix crsec_proposal_sd,
                                                      double cauchy_scale,
                                                      IntegerMatrix crsec_index,
                                                      int no_crsec_interactions,
                                                      int no_persons,
                                                      int no_nodes,
                                                      int no_timepoints,
                                                      NumericMatrix crsec_theta) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  int node1;
  int node2;

  for(int cntr = 0; cntr < no_crsec_interactions; cntr ++) {
    node1 = crsec_index(cntr, 1) - 1;
    node2 = crsec_index(cntr, 2) - 1;

    current_state = crsec_interactions(node1, node2);

    if(gamma(node1, node2) == 0) {
      proposed_state = R::rnorm(current_state, crsec_proposal_sd(node1, node2));
    } else {
      proposed_state = 0.0;
    }

    log_prob = log_pseudolikelihood_ratio_cross_sectional(crsec_interactions,
                                                          crlag_interactions,
                                                          thresholds,
                                                          observations,
                                                          no_categories,
                                                          start,
                                                          no_persons,
                                                          no_nodes,
                                                          no_timepoints,
                                                          node1,
                                                          node2,
                                                          proposed_state,
                                                          current_state);

    if(gamma(node1, node2) == 0) {
      log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           crsec_proposal_sd(node1, node2),
                           true);

      log_prob += log(crsec_theta(node1, node2) /
        (1 - crsec_theta(node1, node2)));
    } else {
      log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale,  true);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           crsec_proposal_sd(node1, node2),
                           true);

      log_prob -= log(crsec_theta(node1, node2) /
        (1 - crsec_theta(node1, node2)));
    }

    U = R::unif_rand();
    if(std::log(U) < log_prob) {
      gamma(node1, node2) = 1 - gamma(node1, node2);
      gamma(node2, node1) = 1 - gamma(node2, node1);

      crsec_interactions(node1, node2) = proposed_state;
      crsec_interactions(node2, node1) = proposed_state;
    }
  }
  return List::create(Named("crsec_interactions") = crsec_interactions,
                      Named("gamma") = gamma);
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of a cross-lagged edge +
//  interaction pair (using a cauchy prior)
// ----------------------------------------------------------------------------|
List metropolis_cross_lagged_edge_interaction_pair(NumericMatrix crsec_interactions,
                                                   NumericMatrix crlag_interactions,
                                                   NumericMatrix thresholds,
                                                   IntegerMatrix delta,
                                                   IntegerMatrix observations,
                                                   IntegerVector no_categories,
                                                   IntegerVector start,
                                                   NumericMatrix crlag_proposal_sd,
                                                   double cauchy_scale,
                                                   IntegerMatrix crlag_index,
                                                   int no_crlag_interactions,
                                                   int no_persons,
                                                   int no_nodes,
                                                   int no_timepoints,
                                                   NumericMatrix crlag_theta) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  int node1;
  int node2;

  for(int cntr = 0; cntr < no_crlag_interactions; cntr ++) {
    node1 = crlag_index(cntr, 1) - 1;
    node2 = crlag_index(cntr, 2) - 1;

    current_state = crlag_interactions(node1, node2);

    if(delta(node1, node2) == 0) {
      proposed_state = R::rnorm(current_state, crlag_proposal_sd(node1, node2));
    } else {
      proposed_state = 0.0;
    }

    log_prob = log_pseudolikelihood_ratio_cross_lagged(crsec_interactions,
                                                       crlag_interactions,
                                                       thresholds,
                                                       observations,
                                                       no_categories,
                                                       start,
                                                       no_persons,
                                                       no_nodes,
                                                       no_timepoints,
                                                       node1,
                                                       node2,
                                                       proposed_state,
                                                       current_state);

    if(delta(node1, node2) == 0) {
      log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           crlag_proposal_sd(node1, node2),
                           true);

      log_prob += log(crlag_theta(node1, node2) /
        (1 - crlag_theta(node1, node2)));
    } else {
      log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale,  true);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           crlag_proposal_sd(node1, node2),
                           true);

      log_prob -= log(crlag_theta(node1, node2) /
        (1 - crlag_theta(node1, node2)));
    }

    U = R::unif_rand();
    if(std::log(U) < log_prob) {
      delta(node1, node2) = 1 - delta(node1, node2);
      crlag_interactions(node1, node2) = proposed_state;
    }
  }
  return List::create(Named("crlag_interactions") = crlag_interactions,
                      Named("delta") = delta);
}

// ----------------------------------------------------------------------------|
// Gibbs step for graphical model parameters
// ----------------------------------------------------------------------------|
List gibbs_step_cross_lagged_mrf(IntegerMatrix observations,
                                 IntegerVector no_categories,
                                 IntegerVector start,
                                 double cauchy_scale,
                                 NumericMatrix crsec_proposal_sd,
                                 NumericMatrix crlag_proposal_sd,
                                 IntegerMatrix crsec_index,
                                 IntegerMatrix crlag_index,
                                 IntegerMatrix n_cat_obs,
                                 double threshold_alpha,
                                 double threshold_beta,
                                 int no_persons,
                                 int no_nodes,
                                 int no_timepoints,
                                 int no_crsec_interactions,
                                 int no_crlag_interactions,
                                 int no_thresholds,
                                 int max_no_categories,
                                 IntegerMatrix gamma,
                                 IntegerMatrix delta,
                                 NumericMatrix crsec_interactions,
                                 NumericMatrix crlag_interactions,
                                 NumericMatrix thresholds,
                                 NumericMatrix crsec_theta,
                                 NumericMatrix crlag_theta,
                                 double phi,
                                 double target_ar,
                                 int tt,
                                 double epsilon_lo,
                                 double epsilon_hi) {

  //Block I: Thresholds
  thresholds = metropolis_thresholds_cross_lagged_panel_mrf(crsec_interactions,
                                                            crlag_interactions,
                                                            thresholds,
                                                            observations,
                                                            no_categories,
                                                            start,
                                                            n_cat_obs,
                                                            no_persons,
                                                            no_nodes,
                                                            no_timepoints,
                                                            threshold_alpha,
                                                            threshold_beta);

  // //Block II: Cross-Sectional
  List out = metropolis_cross_sectional_edge_interaction_pair(crsec_interactions,
                                                              crlag_interactions,
                                                              thresholds,
                                                              gamma,
                                                              observations,
                                                              no_categories,
                                                              start,
                                                              crsec_proposal_sd,
                                                              cauchy_scale,
                                                              crsec_index,
                                                              no_crsec_interactions,
                                                              no_persons,
                                                              no_nodes,
                                                              no_timepoints,
                                                              crsec_theta);
  IntegerMatrix G = out["gamma"];
  NumericMatrix THETA = out["crsec_interactions"];
  gamma = G;
  crsec_interactions = THETA;

  // //Block III: Cross-Lagged
  out = metropolis_cross_lagged_edge_interaction_pair(crsec_interactions,
                                                      crlag_interactions,
                                                      thresholds,
                                                      delta,
                                                      observations,
                                                      no_categories,
                                                      start,
                                                      crlag_proposal_sd,
                                                      cauchy_scale,
                                                      crlag_index,
                                                      no_crlag_interactions,
                                                      no_persons,
                                                      no_nodes,
                                                      no_timepoints,
                                                      crlag_theta);
  IntegerMatrix D = out["delta"];
  NumericMatrix PHI = out["crlag_interactions"];
  delta = D;
  crlag_interactions = PHI;

  //Update cross-sectional interactions (within model move)
  out = metropolis_cross_sectional_interactions(crsec_interactions,
                                                crlag_interactions,
                                                thresholds,
                                                gamma,
                                                observations,
                                                no_categories,
                                                start,
                                                crsec_proposal_sd,
                                                cauchy_scale,
                                                no_persons,
                                                no_nodes,
                                                no_timepoints,
                                                phi,
                                                target_ar,
                                                tt,
                                                epsilon_lo,
                                                epsilon_hi);
  NumericMatrix THETAs = out["crsec_interactions"];
  crsec_interactions = THETAs;
  NumericMatrix NU = out["crsec_proposal_sd"];
  crsec_proposal_sd = NU;

  //Update cross-lagged interactions (within model move)
  out = metropolis_cross_lagged_interactions(crsec_interactions,
                                             crlag_interactions,
                                             thresholds,
                                             delta,
                                             observations,
                                             no_categories,
                                             start,
                                             crlag_proposal_sd,
                                             cauchy_scale,
                                             no_persons,
                                             no_nodes,
                                             no_timepoints,
                                             phi,
                                             target_ar,
                                             tt,
                                             epsilon_lo,
                                             epsilon_hi);
  NumericMatrix PHIs = out["crlag_interactions"];
  crlag_interactions = PHIs;
  NumericMatrix NUl = out["crlag_proposal_sd"];
  crlag_proposal_sd = NUl;

  return List::create(Named("gamma") = gamma,
                      Named("delta") = delta,
                      Named("crsec_interactions") = crsec_interactions,
                      Named("crlag_interactions") = crlag_interactions,
                      Named("thresholds") = thresholds,
                      Named("crsec_proposal_sd") = crsec_proposal_sd,
                      Named("crlag_proposal_sd") = crlag_proposal_sd);
}

// ----------------------------------------------------------------------------|
// The Gibbs sampler for the Cross-Lagged Network Model
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List gibbs_sampler_cross_lagged_mrf(IntegerMatrix observations,
                                    int no_persons,
                                    int no_nodes,
                                    int no_timepoints,
                                    IntegerMatrix gamma,
                                    IntegerMatrix delta,
                                    NumericMatrix crsec_interactions,
                                    NumericMatrix crlag_interactions,
                                    NumericMatrix thresholds,
                                    IntegerVector no_categories,
                                    IntegerVector start,
                                    double cauchy_scale,
                                    NumericMatrix crsec_proposal_sd,
                                    NumericMatrix crlag_proposal_sd,
                                    String crsec_edge_prior,
                                    String crlag_edge_prior,
                                    NumericMatrix crsec_theta,
                                    NumericMatrix crlag_theta,
                                    double crsec_beta_bernoulli_alpha,
                                    double crsec_beta_bernoulli_beta,
                                    double crlag_beta_bernoulli_alpha,
                                    double crlag_beta_bernoulli_beta,
                                    IntegerMatrix crsec_Index,
                                    IntegerMatrix crlag_Index,
                                    int iter,
                                    int burnin,
                                    IntegerMatrix n_cat_obs,
                                    double threshold_alpha,
                                    double threshold_beta,
                                    bool save = false,
                                    bool display_progress = false) {
  int cntr;
  int no_crsec_interactions = no_nodes * (no_nodes - 1) / 2;
  int no_crlag_interactions = no_nodes * no_nodes;
  int no_thresholds = sum(no_categories) * no_timepoints;
  int max_no_categories = max(no_categories);

  IntegerVector crsec_v = seq(0, no_crsec_interactions - 1);
  IntegerVector crlag_v = seq(0, no_crlag_interactions - 1);
  IntegerVector crsec_order(no_crsec_interactions);
  IntegerVector crlag_order(no_crlag_interactions);
  IntegerMatrix crsec_index(no_crsec_interactions, 3);
  IntegerMatrix crlag_index(no_crlag_interactions, 3);

  //Parameters of adaptive proposals -------------------------------------------
  double phi = .75;
  double target_ar = 0.234;
  double epsilon_lo = 1 / no_persons;
  double epsilon_hi = 2.0;

  //The resizing based on ``save'' could probably be prettier ------------------
  int nrow = no_nodes;
  int nrow_thresholds = no_nodes * no_timepoints;
  int ncol_crsec_edges = no_nodes;
  int ncol_crlag_edges = no_nodes;
  int ncol_thresholds = max_no_categories;

  if(save == true) {
    nrow = iter;
    ncol_crsec_edges = no_crsec_interactions;
    ncol_crlag_edges = no_crlag_interactions;
    nrow_thresholds = iter;
    ncol_thresholds = no_thresholds;
  }

  NumericMatrix out_gamma(nrow, ncol_crsec_edges);
  NumericMatrix out_delta(nrow, ncol_crlag_edges);
  NumericMatrix out_crsec_interactions(nrow, ncol_crsec_edges);
  NumericMatrix out_crlag_interactions(nrow, ncol_crlag_edges);
  NumericMatrix out_thresholds(nrow_thresholds, ncol_thresholds);

  //Progress bar
  Progress p(iter + burnin, display_progress);

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------
  for(int iteration = 0; iteration < burnin; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("gamma") = out_gamma,
                          Named("delta") = out_delta,
                          Named("crsec_interactions") = out_crsec_interactions,
                          Named("crlag_interactions") = out_crlag_interactions,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update cross-sectional interactions and model (between model move) -------
    //Use a random order to update the edge - interaction pairs ----------------
    crsec_order = sample(crsec_v,
                         no_crsec_interactions,
                         false,
                         R_NilValue);

    for(int cntr = 0; cntr < no_crsec_interactions; cntr++) {
      crsec_index(cntr, 0) = crsec_Index(crsec_order[cntr], 0);
      crsec_index(cntr, 1) = crsec_Index(crsec_order[cntr], 1);
      crsec_index(cntr, 2) = crsec_Index(crsec_order[cntr], 2);
    }

    //Update cross-lagged interactions and model (between model move) ----------
    //Use a random order to update the edge - interaction pairs ----------------
    crlag_order = sample(crlag_v,
                         no_crlag_interactions,
                         false,
                         R_NilValue);

    for(int cntr = 0; cntr < no_crlag_interactions; cntr++) {
      crlag_index(cntr, 0) = crlag_Index(crlag_order[cntr], 0);
      crlag_index(cntr, 1) = crlag_Index(crlag_order[cntr], 1);
      crlag_index(cntr, 2) = crlag_Index(crlag_order[cntr], 2);
    }

    List out = gibbs_step_cross_lagged_mrf(observations,
                                           no_categories,
                                           start,
                                           cauchy_scale,
                                           crsec_proposal_sd,
                                           crlag_proposal_sd,
                                           crsec_index,
                                           crlag_index,
                                           n_cat_obs,
                                           threshold_alpha,
                                           threshold_beta,
                                           no_persons,
                                           no_nodes,
                                           no_timepoints,
                                           no_crsec_interactions,
                                           no_crlag_interactions,
                                           no_thresholds,
                                           max_no_categories,
                                           gamma,
                                           delta,
                                           crsec_interactions,
                                           crlag_interactions,
                                           thresholds,
                                           crsec_theta,
                                           crlag_theta,
                                           phi,
                                           target_ar,
                                           iteration + 1,
                                           epsilon_lo,
                                           epsilon_hi);

    IntegerMatrix gamma = out["gamma"];
    IntegerMatrix delta = out["delta"];
    NumericMatrix crsec_interactions = out["crsec_interactions"];
    NumericMatrix crlag_interactions = out["crlag_interactions"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix crsec_proposal_sd = out["crsec_proposal_sd"];
    NumericMatrix crlag_proposal_sd = out["crlag_proposal_sd"];

    if(crsec_edge_prior == "Beta-Bernoulli") {
      int sumG = 0;
      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          sumG += gamma(i, j);
        }
      }
      double probability = R::rbeta(crsec_beta_bernoulli_alpha + sumG,
                                    crsec_beta_bernoulli_beta + no_crsec_interactions - sumG);

      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          crsec_theta(i, j) = probability;
          crsec_theta(j, i) = probability;
        }
      }
    }

    if(crlag_edge_prior == "Beta-Bernoulli") {
      int sumG = 0;
      for(int i = 0; i < no_nodes; i++) {
        for(int j = 0; j < no_nodes; j++) {
          sumG += delta(i, j);
        }
      }
      double probability = R::rbeta(crlag_beta_bernoulli_alpha + sumG,
                                    crlag_beta_bernoulli_beta + no_crlag_interactions - sumG);

      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          crlag_theta(i, j) = probability;
          crlag_theta(j, i) = probability;
        }
      }
    }
  }

  //The post burn-in iterations ------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("gamma") = out_gamma,
                          Named("delta") = out_delta,
                          Named("crsec_interactions") = out_crsec_interactions,
                          Named("crlag_interactions") = out_crlag_interactions,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update cross-sectional interactions and model (between model move) -------
    //Use a random order to update the edge - interaction pairs ----------------
    crsec_order = sample(crsec_v,
                         no_crsec_interactions,
                         false,
                         R_NilValue);

    for(int cntr = 0; cntr < no_crsec_interactions; cntr++) {
      crsec_index(cntr, 0) = crsec_Index(crsec_order[cntr], 0);
      crsec_index(cntr, 1) = crsec_Index(crsec_order[cntr], 1);
      crsec_index(cntr, 2) = crsec_Index(crsec_order[cntr], 2);
    }

    //Update cross-lagged interactions and model (between model move) ----------
    //Use a random order to update the edge - interaction pairs ----------------
    crlag_order = sample(crlag_v,
                         no_crlag_interactions,
                         false,
                         R_NilValue);

    for(int cntr = 0; cntr < no_crlag_interactions; cntr++) {
      crlag_index(cntr, 0) = crlag_Index(crlag_order[cntr], 0);
      crlag_index(cntr, 1) = crlag_Index(crlag_order[cntr], 1);
      crlag_index(cntr, 2) = crlag_Index(crlag_order[cntr], 2);
    }

    List out = gibbs_step_cross_lagged_mrf(observations,
                                           no_categories,
                                           start,
                                           cauchy_scale,
                                           crsec_proposal_sd,
                                           crlag_proposal_sd,
                                           crsec_index,
                                           crlag_index,
                                           n_cat_obs,
                                           threshold_alpha,
                                           threshold_beta,
                                           no_persons,
                                           no_nodes,
                                           no_timepoints,
                                           no_crsec_interactions,
                                           no_crlag_interactions,
                                           no_thresholds,
                                           max_no_categories,
                                           gamma,
                                           delta,
                                           crsec_interactions,
                                           crlag_interactions,
                                           thresholds,
                                           crsec_theta,
                                           crlag_theta,
                                           phi,
                                           target_ar,
                                           iteration + 1,
                                           epsilon_lo,
                                           epsilon_hi);


    IntegerMatrix gamma = out["gamma"];
    IntegerMatrix delta = out["delta"];
    NumericMatrix crsec_interactions = out["crsec_interactions"];
    NumericMatrix crlag_interactions = out["crlag_interactions"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix crsec_proposal_sd = out["crsec_proposal_sd"];
    NumericMatrix crlag_proposal_sd = out["crlag_proposal_sd"];

    if(crsec_edge_prior == "Beta-Bernoulli") {
      int sumG = 0;
      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          sumG += gamma(i, j);
        }
      }
      double probability = R::rbeta(crsec_beta_bernoulli_alpha + sumG,
                                    crsec_beta_bernoulli_beta + no_crsec_interactions - sumG);

      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          crsec_theta(i, j) = probability;
          crsec_theta(j, i) = probability;
        }
      }
    }

    if(crlag_edge_prior == "Beta-Bernoulli") {
      int sumG = 0;
      for(int i = 0; i < no_nodes; i++) {
        for(int j = 0; j < no_nodes; j++) {
          sumG += delta(i, j);
        }
      }
      double probability = R::rbeta(crlag_beta_bernoulli_alpha + sumG,
                                    crlag_beta_bernoulli_beta + no_crlag_interactions - sumG);

      for(int i = 0; i < no_nodes - 1; i++) {
        for(int j = i + 1; j < no_nodes; j++) {
          crlag_theta(i, j) = probability;
          crlag_theta(j, i) = probability;
        }
      }
    }

    //Output -------------------------------------------------------------------
    if(save == TRUE) {
      //Save raw samples -------------------------------------------------------
      cntr = 0;
      for(int node1 = 0; node1 < no_nodes - 1; node1++) {
        for(int node2 = node1 + 1; node2 < no_nodes;node2++) {
          out_gamma(iteration, cntr) = gamma(node1, node2);
          out_crsec_interactions(iteration, cntr) = crsec_interactions(node1, node2);
          cntr++;
        }
      }
      cntr = 0;
      for(int node1 = 0; node1 < no_nodes; node1++) {
        for(int node2 = 0; node2 < no_nodes;node2++) {
          out_delta(iteration, cntr) = delta(node1, node2);
          out_crlag_interactions(iteration, cntr) = crlag_interactions(node1, node2);
          cntr++;
        }
      }
      cntr = 0;
      for(int t = 0; t < no_timepoints; t++) {
        for(int node = 0; node < no_nodes; node++) {
          for(int category = 0; category < no_categories[node]; category++) {
            out_thresholds(iteration, cntr) = thresholds(node + start[t], category);
            cntr++;
          }
        }
      }

    } else {
      //Compute running averages -----------------------------------------------
      for(int node1 = 0; node1 < no_nodes - 1; node1++) {
        for(int node2 = node1 + 1; node2 < no_nodes; node2++) {
          out_gamma(node1, node2) *= iteration;
          out_gamma(node1, node2) += gamma(node1, node2);
          out_gamma(node1, node2) /= iteration + 1;
          out_gamma(node2, node1) = out_gamma(node1, node2);

          out_crsec_interactions(node1, node2) *= iteration;
          out_crsec_interactions(node1, node2) += crsec_interactions(node1, node2);
          out_crsec_interactions(node1, node2) /= iteration + 1;
          out_crsec_interactions(node2, node1) = out_crsec_interactions(node1, node2);
        }
      }

      for(int node1 = 0; node1 < no_nodes; node1++) {
        for(int node2 = 0; node2 < no_nodes; node2++) {
          out_delta(node1, node2) *= iteration;
          out_delta(node1, node2) += delta(node1, node2);
          out_delta(node1, node2) /= iteration + 1;

          out_crlag_interactions(node1, node2) *= iteration;
          out_crlag_interactions(node1, node2) += crlag_interactions(node1, node2);
          out_crlag_interactions(node1, node2) /= iteration + 1;
        }
      }

      for(int t = 0; t < no_timepoints; t++) {
        for(int node = 0; node < no_nodes; node++) {
          for(int category = 0; category < no_categories[node]; category++) {
            out_thresholds(node + start[t], category) *= iteration;
            out_thresholds(node + start[t], category) += thresholds(node + start[t], category);
            out_thresholds(node + start[t], category) /= iteration + 1;
          }
        }
      }
    }
  }

  return List::create(Named("gamma") = out_gamma,
                      Named("delta") = out_delta,
                      Named("crsec_interactions") = out_crsec_interactions,
                      Named("crlag_interactions") = out_crlag_interactions,
                      Named("thresholds") = out_thresholds);
}