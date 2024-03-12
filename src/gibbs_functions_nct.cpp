// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "gibbs_functions.h"
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// Impute missing data from full-conditional
// ----------------------------------------------------------------------------|
List nct_impute_missing_data(NumericMatrix interactions,
                             NumericMatrix thresholds,
                             NumericMatrix interactions_difference,
                             IntegerVector group_indicator,
                             IntegerMatrix observations,
                             IntegerMatrix n_cat_obs,
                             IntegerMatrix sufficient_blume_capel,
                             IntegerVector no_categories,
                             NumericMatrix rest_matrix,
                             IntegerMatrix missing_index,
                             LogicalVector variable_bool,
                             IntegerVector reference_category) {

  int no_variables = observations.ncol();
  int no_missings = missing_index.nrow();
  int max_no_categories = 0;
  for(int variable = 0; variable < no_variables; variable++) {
    if(no_categories[variable] > max_no_categories) {
      max_no_categories = no_categories[variable];
    }
  }
  NumericVector probabilities(max_no_categories + 1);
  double exponent, rest_score, cumsum, u, interaction, effect;
  int score, person, variable, new_observation, old_observation;

  for(int missing = 0; missing < no_missings; missing++) {
    //Which observation to impute? ---------------------------------------------
    person = missing_index(missing, 0) - 1; //R to C++ indexing
    variable = missing_index(missing, 1) - 1; //R to C++ indexing
    effect = group_indicator[person] - 0.5;

    //Generate new observation -------------------------------------------------
    rest_score = rest_matrix(person, variable);

    //Two distinct (ordinal) variable types ------------------------------------
    if(variable_bool[variable] == true) {

      //Regular binary or ordinal MRF variable ---------------------------------
      cumsum = 1.0;
      probabilities[0] = 1.0;
      for(int category = 0; category < no_categories[variable]; category++) {
        exponent = thresholds(variable, category);
        exponent += (category + 1) * rest_score;
        cumsum += std::exp(exponent);
        probabilities[category + 1] = cumsum;
      }

    } else {

      //Blume-Capel ordinal MRF variable ---------------------------------------
      exponent = thresholds(variable, 1) *
        reference_category[variable] *
        reference_category[variable];
      cumsum = std::exp(exponent);
      probabilities[0] = cumsum;
      for(int category = 0; category < no_categories[variable]; category++) {
        exponent = thresholds(variable, 0) * (category + 1);
        exponent += thresholds(variable, 1) *
          (category + 1 - reference_category[variable]) *
          (category + 1 - reference_category[variable]);
        exponent += (category + 1) * rest_score;
        cumsum += std::exp(exponent);
        probabilities[category + 1] = cumsum;
      }
    }

    u = cumsum * R::unif_rand();
    score = 0;
    while (u > probabilities[score]) {
      score++;
    }

    //Update observations
    new_observation = score;
    old_observation = observations(person, variable);
    if(old_observation != new_observation) {
      observations(person, variable) = new_observation;
      if(variable_bool[variable] == true) {
        //Regular binary or ordinal MRF variable -------------------------------
        n_cat_obs(old_observation, variable)--;
        n_cat_obs(new_observation, variable)++;
      } else {
        //Regular binary or ordinal MRF variable -------------------------------
        sufficient_blume_capel(0, variable) -= old_observation;
        sufficient_blume_capel(0, variable) += new_observation;
        sufficient_blume_capel(1, variable) -=
          (old_observation - reference_category[variable]) *
          (old_observation - reference_category[variable]);
        sufficient_blume_capel(1, variable) +=
          (new_observation - reference_category[variable]) *
          (new_observation - reference_category[variable]);
      }

      for(int vertex = 0; vertex < no_variables; vertex++) {
        interaction = interactions(vertex, variable) +
          effect * interactions_difference(vertex, variable);
        //interactions(i, i) = 0
        rest_matrix(person, vertex) -= old_observation *
          interaction;
        rest_matrix(person, vertex) += new_observation *
          interaction;
      }
    }
  }

  return List::create(Named("observations") = observations,
                      Named("n_cat_obs") = n_cat_obs,
                      Named("sufficient_blume_capel") =
                        sufficient_blume_capel,
                        Named("rest_matrix") = rest_matrix);
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the active interaction
//  parameters for Bayesian edge selection
// ----------------------------------------------------------------------------|
void nct_metropolis_interactions(NumericMatrix interactions,
                                 NumericMatrix thresholds,
                                 IntegerMatrix observations,
                                 IntegerVector no_categories,
                                 NumericMatrix proposal_sd,
                                 double interaction_scale,
                                 int no_persons,
                                 int no_variables,
                                 NumericMatrix rest_matrix,
                                 double phi,
                                 double target_ar,
                                 int t,
                                 double epsilon_lo,
                                 double epsilon_hi,
                                 LogicalVector variable_bool,
                                 IntegerVector reference_category) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  for(int variable1 = 0; variable1 <  no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  no_variables; variable2++) {
      current_state = interactions(variable1, variable2);
      proposed_state = R::rnorm(current_state, proposal_sd(variable1, variable2));

      log_prob = log_pseudolikelihood_ratio(interactions,
                                            thresholds,
                                            observations,
                                            no_categories,
                                            no_persons,
                                            variable1,
                                            variable2,
                                            proposed_state,
                                            current_state,
                                            rest_matrix,
                                            variable_bool,
                                            reference_category);
      log_prob += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_prob -= R::dcauchy(current_state, 0.0, interaction_scale, true);

      U = R::unif_rand();
      if(std::log(U) < log_prob) {
        double state_diff = proposed_state - current_state;
        interactions(variable1, variable2) = proposed_state;
        interactions(variable2, variable1) = proposed_state;

        //Update the matrix of rest scores
        for(int person = 0; person < no_persons; person++) {
          rest_matrix(person, variable1) += observations(person, variable2) *
            state_diff;
          rest_matrix(person, variable2) += observations(person, variable1) *
            state_diff;
        }
      }

      if(log_prob > 0) {
        log_prob = 1;
      } else {
        log_prob = std::exp(log_prob);
      }
      proposal_sd(variable1, variable2) = proposal_sd(variable1, variable2) +
        (log_prob - target_ar) * std::exp(-log(t) * phi);
      if(proposal_sd(variable1, variable2) < epsilon_lo) {
        proposal_sd(variable1, variable2) = epsilon_lo;
      } else if (proposal_sd(variable1, variable2) > epsilon_hi) {
        proposal_sd(variable1, variable2) = epsilon_hi;
      }
    }
  }
}


// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for anm,. interaction
// ----------------------------------------------------------------------------|
double nct_log_pseudolikelihood_ratio_int_diff(NumericMatrix thresholds,
                                               NumericMatrix interactions_difference,
                                               IntegerVector group_indicator,
                                               IntegerMatrix observations,
                                               IntegerVector no_categories,
                                               int no_persons,
                                               int variable1,
                                               int variable2,
                                               double proposed_state,
                                               double current_state,
                                               NumericMatrix rest_matrix,
                                               LogicalVector variable_bool,
                                               IntegerVector reference_category) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;
  double delta_state = proposed_state - current_state;

  for(int person = 0; person < no_persons; person++) {
    double effect = group_indicator[person] - 0.5;

    obs_score1 = observations(person, variable1);
    obs_score2 = observations(person, variable2);

    pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * delta_state * effect;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix(person, variable1) -
      effect * obs_score2 * interactions_difference(variable2, variable1);
    if(rest_score > 0) {
      bound = no_categories[variable1] * rest_score;
    } else {
      bound = 0.0;
    }

    if(variable_bool[variable1] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories[variable1]; category++) {
        score = category + 1;
        exponent = thresholds(variable1, category) +
          score * rest_score -
          bound;
        denominator_prop +=
          std::exp(exponent + score * obs_score2 * proposed_state * effect);
        denominator_curr +=
          std::exp(exponent + score * obs_score2 * current_state * effect);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable1] + 1; category++) {
        exponent = thresholds(variable1, 0) * category;
        exponent += thresholds(variable1, 1) *
          (category - reference_category[variable1]) *
          (category - reference_category[variable1]);
        exponent+= category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * obs_score2 * proposed_state * effect);
        denominator_curr +=
          std::exp(exponent + category * obs_score2 * current_state * effect);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);

    //variable 2 log pseudolikelihood ratio
    rest_score = rest_matrix(person, variable2) -
      effect * obs_score1 * interactions_difference(variable1, variable2);

    if(rest_score > 0) {
      bound = no_categories[variable2] * rest_score;
    } else {
      bound = 0.0;
    }

    if(variable_bool[variable2] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories[variable2]; category++) {
        score = category + 1;
        exponent = thresholds(variable2, category) +
          score * rest_score -
          bound;
        denominator_prop +=
          std::exp(exponent + score * obs_score1 * proposed_state * effect);
        denominator_curr +=
          std::exp(exponent + score * obs_score1 * current_state * effect);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable2] + 1; category++) {
        exponent = thresholds(variable2, 0) * category;
        exponent += thresholds(variable2, 1) *
          (category - reference_category[variable2]) *
          (category - reference_category[variable2]);
        exponent +=  category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * obs_score1 * proposed_state * effect);
        denominator_curr +=
          std::exp(exponent + category * obs_score1 * current_state * effect);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);
  }
  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the active interaction
//  parameters for Bayesian edge selection
// ----------------------------------------------------------------------------|
void nct_metropolis_int_diff(NumericMatrix thresholds,
                             NumericMatrix interactions_difference,
                             IntegerVector group_indicator,
                             IntegerMatrix gamma,
                             IntegerMatrix observations,
                             IntegerVector no_categories,
                             NumericMatrix proposal_sd_int_diff,
                             double difference_scale,
                             int no_persons,
                             int no_variables,
                             NumericMatrix rest_matrix,
                             double phi,
                             double target_ar,
                             int t,
                             double epsilon_lo,
                             double epsilon_hi,
                             LogicalVector variable_bool,
                             IntegerVector reference_category) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  for(int variable1 = 0; variable1 <  no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  no_variables; variable2++) {
      if(gamma(variable1, variable2) == 1) {
        current_state = interactions_difference(variable1, variable2);
        proposed_state = R::rnorm(current_state, proposal_sd_int_diff(variable1, variable2));

        log_prob = nct_log_pseudolikelihood_ratio_int_diff(thresholds,
                                                           interactions_difference,
                                                           group_indicator,
                                                           observations,
                                                           no_categories,
                                                           no_persons,
                                                           variable1,
                                                           variable2,
                                                           proposed_state,
                                                           current_state,
                                                           rest_matrix,
                                                           variable_bool,
                                                           reference_category);

        log_prob += R::dcauchy(proposed_state, 0.0, difference_scale, true);
        log_prob -= R::dcauchy(current_state, 0.0, difference_scale, true);

        U = R::unif_rand();
        if(std::log(U) < log_prob) {
          double state_diff = proposed_state - current_state;
          interactions_difference(variable1, variable2) = proposed_state;
          interactions_difference(variable2, variable1) = proposed_state;

          //Update the matrix of rest scores
          for(int person = 0; person < no_persons; person++) {
            double effect = group_indicator[person] - 0.5;
            rest_matrix(person, variable1) += observations(person, variable2) *
              state_diff * effect;
            rest_matrix(person, variable2) += observations(person, variable1) *
              state_diff * effect;
          }
        }

        if(log_prob > 0) {
          log_prob = 1;
        } else {
          log_prob = std::exp(log_prob);
        }
        proposal_sd_int_diff(variable1, variable2) = proposal_sd_int_diff(variable1, variable2) +
          (log_prob - target_ar) * std::exp(-log(t) * phi);
        if(proposal_sd_int_diff(variable1, variable2) < epsilon_lo) {
          proposal_sd_int_diff(variable1, variable2) = epsilon_lo;
        } else if (proposal_sd_int_diff(variable1, variable2) > epsilon_hi) {
          proposal_sd_int_diff(variable1, variable2) = epsilon_hi;
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of an edge + interaction
//  pair for Bayesian edge selection
// ----------------------------------------------------------------------------|
void nct_metropolis_effect_int_diff_pair(NumericMatrix thresholds,
                                         NumericMatrix interactions_difference,
                                         IntegerVector group_indicator,
                                         IntegerMatrix gamma,
                                         IntegerMatrix observations,
                                         IntegerVector no_categories,
                                         NumericMatrix proposal_sd_int_diff,
                                         double difference_scale,
                                         IntegerMatrix index,
                                         int no_interactions,
                                         int no_persons,
                                         NumericMatrix rest_matrix,
                                         NumericMatrix theta,
                                         LogicalVector variable_bool,
                                         IntegerVector reference_category) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  int variable1;
  int variable2;

  for(int cntr = 0; cntr < no_interactions; cntr ++) {
    variable1 = index(cntr, 1) - 1;
    variable2 = index(cntr, 2) - 1;

    current_state = interactions_difference(variable1, variable2);

    if(gamma(variable1, variable2) == 0) {
      proposed_state = R::rnorm(current_state, proposal_sd_int_diff(variable1, variable2));
    } else {
      proposed_state = 0.0;
    }

    log_prob = nct_log_pseudolikelihood_ratio_int_diff(thresholds,
                                                       interactions_difference,
                                                       group_indicator,
                                                       observations,
                                                       no_categories,
                                                       no_persons,
                                                       variable1,
                                                       variable2,
                                                       proposed_state,
                                                       current_state,
                                                       rest_matrix,
                                                       variable_bool,
                                                       reference_category);

    if(gamma(variable1, variable2) == 0) {
      log_prob += R::dcauchy(proposed_state, 0.0, difference_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_sd_int_diff(variable1, variable2),
                           true);

      log_prob += log(theta(variable1, variable2) / (1 - theta(variable1, variable2)));
    } else {
      log_prob -= R::dcauchy(current_state, 0.0, difference_scale, true);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_sd_int_diff(variable1, variable2),
                           true);

      log_prob -= log(theta(variable1, variable2) / (1 - theta(variable1, variable2)));
    }

    U = R::unif_rand();
    if(std::log(U) < log_prob) {
      gamma(variable1, variable2) = 1 - gamma(variable1, variable2);
      gamma(variable2, variable1) = 1 - gamma(variable2, variable1);

      interactions_difference(variable1, variable2) = proposed_state;
      interactions_difference(variable2, variable1) = proposed_state;

      double state_diff = proposed_state - current_state;

      //Update the matrix of rest scores ---------------------------------------
      for(int person = 0; person < no_persons; person++) {
        double effect = group_indicator[person] - 0.5;
        rest_matrix(person, variable1) += observations(person, variable2) *
          state_diff * effect;
        rest_matrix(person, variable2) += observations(person, variable1) *
          state_diff * effect;
      }
    }
  }
}

// ----------------------------------------------------------------------------|
// A Gibbs step for graphical model parameters for Bayesian edge selection
// ----------------------------------------------------------------------------|
List nct_gibbs_step_gm(IntegerMatrix observations,
                       IntegerVector no_categories,
                       double interaction_scale,
                       double difference_scale,
                       NumericMatrix proposal_sd,
                       NumericMatrix proposal_sd_blumecapel,
                       NumericMatrix proposal_sd_int_diff,
                       IntegerMatrix index,
                       IntegerMatrix n_cat_obs,
                       IntegerMatrix sufficient_blume_capel,
                       double threshold_alpha,
                       double threshold_beta,
                       int no_persons,
                       int no_variables,
                       int no_interactions,
                       int no_thresholds,
                       int max_no_categories,
                       IntegerMatrix gamma,
                       NumericMatrix interactions,
                       NumericMatrix thresholds,
                       NumericMatrix interactions_difference,
                       IntegerVector group_indicator,
                       NumericMatrix rest_matrix,
                       NumericMatrix theta,
                       double phi,
                       double target_ar,
                       int t,
                       double epsilon_lo,
                       double epsilon_hi,
                       LogicalVector variable_bool,
                       IntegerVector reference_category,
                       bool edge_selection) {

  //Between model move (update effect indicators and difference parameters)
  nct_metropolis_effect_int_diff_pair(thresholds,
                                      interactions_difference,
                                      group_indicator,
                                      gamma,
                                      observations,
                                      no_categories,
                                      proposal_sd_int_diff,
                                      difference_scale,
                                      index,
                                      no_interactions,
                                      no_persons,
                                      rest_matrix,
                                      theta,
                                      variable_bool,
                                      reference_category);

  //Within model move (update interaction parameters)
  nct_metropolis_interactions(interactions,
                          thresholds,
                          observations,
                          no_categories,
                          proposal_sd,
                          interaction_scale,
                          no_persons,
                          no_variables,
                          rest_matrix,
                          phi,
                          target_ar,
                          t,
                          epsilon_lo,
                          epsilon_hi,
                          variable_bool,
                          reference_category);

  //Within model move (update interaction parameters)
  nct_metropolis_int_diff(thresholds,
                          interactions_difference,
                          group_indicator,
                          gamma,
                          observations,
                          no_categories,
                          proposal_sd_int_diff,
                          difference_scale,
                          no_persons,
                          no_variables,
                          rest_matrix,
                          phi,
                          target_ar,
                          t,
                          epsilon_lo,
                          epsilon_hi,
                          variable_bool,
                          reference_category);

  //Update threshold parameters
  for(int variable = 0; variable < no_variables; variable++) {
    if(variable_bool[variable] == true) {
      metropolis_thresholds_regular(thresholds,
                                    observations,
                                    no_categories,
                                    n_cat_obs,
                                    no_persons,
                                    variable,
                                    threshold_alpha,
                                    threshold_beta,
                                    rest_matrix);
    } else {
      metropolis_thresholds_blumecapel(thresholds,
                                       observations,
                                       no_categories,
                                       sufficient_blume_capel,
                                       no_persons,
                                       variable,
                                       reference_category,
                                       threshold_alpha,
                                       threshold_beta,
                                       rest_matrix,
                                       proposal_sd_blumecapel,
                                       phi,
                                       target_ar,
                                       t,
                                       epsilon_lo,
                                       epsilon_hi);
    }
  }

  return List::create(Named("gamma") = gamma,
                      Named("interactions") = interactions,
                      Named("thresholds") = thresholds,
                      Named("interactions_difference") = interactions_difference,
                      Named("rest_matrix") = rest_matrix,
                      Named("proposal_sd") = proposal_sd);
}

// ----------------------------------------------------------------------------|
// The Gibbs sampler for Bayesian edge selection
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List nct_gibbs_sampler(IntegerMatrix observations,
                       IntegerMatrix gamma,
                       NumericMatrix interactions,
                       NumericMatrix thresholds,
                       NumericMatrix interactions_difference,
                       IntegerVector group_indicator,
                       IntegerVector no_categories,
                       double interaction_scale,
                       double difference_scale,
                       NumericMatrix proposal_sd,
                       NumericMatrix proposal_sd_blumecapel,
                       NumericMatrix proposal_sd_int_diff,
                       String difference_prior,
                       NumericMatrix theta,
                       double beta_bernoulli_alpha,
                       double beta_bernoulli_beta,
                       IntegerMatrix Index,
                       int iter,
                       int burnin,
                       IntegerMatrix n_cat_obs,
                       IntegerMatrix sufficient_blume_capel,
                       double threshold_alpha,
                       double threshold_beta,
                       bool na_impute,
                       IntegerMatrix missing_index,
                       LogicalVector variable_bool,
                       IntegerVector reference_category,
                       bool save = false,
                       bool display_progress = false,
                       bool edge_selection = true) {
  int cntr;
  int no_variables = observations.ncol();
  int no_persons = observations.nrow();
  int no_interactions = Index.nrow();
  int no_thresholds = sum(no_categories);
  int max_no_categories = max(no_categories);

  IntegerVector v = seq(0, no_interactions - 1);
  IntegerVector order(no_interactions);
  IntegerMatrix index(no_interactions, 3);

  //Parameters of adaptive proposals -------------------------------------------
  double phi = .75;
  double target_ar = 0.234;
  double epsilon_lo = 1 / no_persons;
  double epsilon_hi = 2.0;

  //The resizing based on ``save'' could probably be prettier ------------------
  int nrow = no_variables;
  int ncol_edges = no_variables;
  int ncol_thresholds = max_no_categories;

  if(save == true) {
    nrow = iter;
    ncol_edges= no_interactions;
    ncol_thresholds = no_thresholds;
  }

  NumericMatrix out_interactions(nrow, ncol_edges);
  NumericMatrix out_interactions_difference(nrow, ncol_edges);
  NumericMatrix out_gamma(nrow, ncol_edges);
  NumericMatrix out_thresholds(nrow, ncol_thresholds);

  NumericMatrix rest_matrix(no_persons, no_variables);
  for(int variable1 = 0; variable1 < no_variables; variable1++) {
    for(int person = 0; person < no_persons; person++) {
      double effect = group_indicator[person] - 0.5;
      for(int variable2 = 0; variable2 < no_variables; variable2++) {
        rest_matrix(person, variable1) +=
          observations(person, variable2) * (interactions(variable2, variable1) +
            effect * interactions_difference(variable2, variable1));
      }
    }
  }

  //Progress bar
  Progress p(iter + burnin, display_progress);

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------
  for(int iteration = 0; iteration < burnin; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("gamma") = out_gamma,
                          Named("interactions") = out_interactions,
                          Named("interactions_difference") = out_interactions_difference,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_interactions,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    if(na_impute == true) {
      List out = nct_impute_missing_data(interactions,
                                     thresholds,
                                     interactions_difference,
                                     group_indicator,
                                     observations,
                                     n_cat_obs,
                                     sufficient_blume_capel,
                                     no_categories,
                                     rest_matrix,
                                     missing_index,
                                     variable_bool,
                                     reference_category);

      IntegerMatrix observations = out["observations"];
      IntegerMatrix n_cat_obs = out["n_cat_obs"];
      IntegerMatrix sufficient_blume_capel = out["sufficient_blume_capel"];
      NumericMatrix rest_matrix = out["rest_matrix"];
    }

    List out = nct_gibbs_step_gm(observations,
                                 no_categories,
                                 interaction_scale,
                                 difference_scale,
                                 proposal_sd,
                                 proposal_sd_blumecapel,
                                 proposal_sd_int_diff,
                                 index,
                                 n_cat_obs,
                                 sufficient_blume_capel,
                                 threshold_alpha,
                                 threshold_beta,
                                 no_persons,
                                 no_variables,
                                 no_interactions,
                                 no_thresholds,
                                 max_no_categories,
                                 gamma,
                                 interactions,
                                 thresholds,
                                 interactions_difference,
                                 group_indicator,
                                 rest_matrix,
                                 theta,
                                 phi,
                                 target_ar,
                                 iteration + 1,
                                 epsilon_lo,
                                 epsilon_hi,
                                 variable_bool,
                                 reference_category,
                                 edge_selection);

    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix interactions_difference = out["interactions_difference"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix rest_matrix = out["rest_matrix"];
    NumericMatrix proposal_sd = out["proposal_sd"];

    if(edge_selection == true) {
      if(difference_prior == "Beta-Bernoulli") {
        int sumG = 0;
        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            sumG += gamma(i, j);
          }
        }
        double probability = R::rbeta(beta_bernoulli_alpha + sumG,
                                      beta_bernoulli_beta + no_interactions - sumG);

        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            theta(i, j) = probability;
            theta(j, i) = probability;
          }
        }
      }
    }
  }

  //The post burn-in iterations ------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    if (Progress::check_abort()) {
      return List::create(Named("gamma") = out_gamma,
                          Named("interactions") = out_interactions,
                          Named("interactions_difference") = out_interactions_difference,
                          Named("thresholds") = out_thresholds);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_interactions,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_interactions; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    if(na_impute == true) {
      List out = nct_impute_missing_data(interactions,
                                         thresholds,
                                         interactions_difference,
                                         group_indicator,
                                         observations,
                                         n_cat_obs,
                                         sufficient_blume_capel,
                                         no_categories,
                                         rest_matrix,
                                         missing_index,
                                         variable_bool,
                                         reference_category);

      IntegerMatrix observations = out["observations"];
      IntegerMatrix n_cat_obs = out["n_cat_obs"];
      IntegerMatrix sufficient_blume_capel = out["sufficient_blume_capel"];
      NumericMatrix rest_matrix = out["rest_matrix"];
    }

    List out = nct_gibbs_step_gm(observations,
                                 no_categories,
                                 interaction_scale,
                                 difference_scale,
                                 proposal_sd,
                                 proposal_sd_blumecapel,
                                 proposal_sd_int_diff,
                                 index,
                                 n_cat_obs,
                                 sufficient_blume_capel,
                                 threshold_alpha,
                                 threshold_beta,
                                 no_persons,
                                 no_variables,
                                 no_interactions,
                                 no_thresholds,
                                 max_no_categories,
                                 gamma,
                                 interactions,
                                 thresholds,
                                 interactions_difference,
                                 group_indicator,
                                 rest_matrix,
                                 theta,
                                 phi,
                                 target_ar,
                                 iteration + 1,
                                 epsilon_lo,
                                 epsilon_hi,
                                 variable_bool,
                                 reference_category,
                                 edge_selection);

    IntegerMatrix gamma = out["gamma"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix interactions_difference = out["interactions_difference"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix rest_matrix = out["rest_matrix"];
    NumericMatrix proposal_sd = out["proposal_sd"];

    if(edge_selection == true) {
      if(difference_prior == "Beta-Bernoulli") {
        int sumG = 0;
        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            sumG += gamma(i, j);
          }
        }
        double probability = R::rbeta(beta_bernoulli_alpha + sumG,
                                      beta_bernoulli_beta + no_interactions - sumG);

        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            theta(i, j) = probability;
            theta(j, i) = probability;
          }
        }
      }
    }


    //Output -------------------------------------------------------------------
    if(save == true) {
      //Save raw samples -------------------------------------------------------
      cntr = 0;
      for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
        for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
          out_gamma(iteration, cntr) = gamma(variable1, variable2);
          out_interactions(iteration, cntr) = interactions(variable1, variable2);
          out_interactions_difference(iteration, cntr) =
            interactions_difference(variable1, variable2);
          cntr++;
        }
      }
      cntr = 0;
      for(int variable = 0; variable < no_variables; variable++) {
        if(variable_bool[variable] == true) {
          for(int category = 0; category < no_categories[variable]; category++) {
            out_thresholds(iteration, cntr) = thresholds(variable, category);
            cntr++;
          }
        } else {
          out_thresholds(iteration, cntr) = thresholds(variable, 0);
          cntr++;
          out_thresholds(iteration, cntr) = thresholds(variable, 1);
          cntr++;
        }
      }
    } else {
      //Compute running averages -----------------------------------------------
      for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
        for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
          if(edge_selection == true) {
            out_gamma(variable1, variable2) *= iteration;
            out_gamma(variable1, variable2) += gamma(variable1, variable2);
            out_gamma(variable1, variable2) /= iteration + 1;
            out_gamma(variable2, variable1) = out_gamma(variable1, variable2);
          }

          out_interactions(variable1, variable2) *= iteration;
          out_interactions(variable1, variable2) += interactions(variable1, variable2);
          out_interactions(variable1, variable2) /= iteration + 1;
          out_interactions(variable2, variable1) = out_interactions(variable1, variable2);

          out_interactions_difference(variable1, variable2) *= iteration;
          out_interactions_difference(variable1, variable2) += interactions_difference(variable1, variable2);
          out_interactions_difference(variable1, variable2) /= iteration + 1;
          out_interactions_difference(variable2, variable1) = out_interactions_difference(variable1, variable2);
        }

        if(variable_bool[variable1] == true) {
          for(int category = 0; category < no_categories[variable1]; category++) {
            out_thresholds(variable1, category) *= iteration;
            out_thresholds(variable1, category) += thresholds(variable1, category);
            out_thresholds(variable1, category) /= iteration + 1;
          }
        } else {
          out_thresholds(variable1, 0) *= iteration;
          out_thresholds(variable1, 0) += thresholds(variable1, 0);
          out_thresholds(variable1, 0) /= iteration + 1;
          out_thresholds(variable1, 1) *= iteration;
          out_thresholds(variable1, 1) += thresholds(variable1, 1);
          out_thresholds(variable1, 1) /= iteration + 1;
        }
      }
      if(variable_bool[no_variables - 1] == true) {
        for(int category = 0; category < no_categories[no_variables - 1]; category++) {
          out_thresholds(no_variables - 1, category) *= iteration;
          out_thresholds(no_variables - 1, category) += thresholds(no_variables - 1, category);
          out_thresholds(no_variables - 1, category) /= iteration + 1;
        }
      } else {
        out_thresholds(no_variables - 1, 0) *= iteration;
        out_thresholds(no_variables - 1, 0) += thresholds(no_variables - 1, 0);
        out_thresholds(no_variables - 1, 0) /= iteration + 1;
        out_thresholds(no_variables - 1, 1) *= iteration;
        out_thresholds(no_variables - 1, 1) += thresholds(no_variables - 1, 1);
        out_thresholds(no_variables - 1, 1) /= iteration + 1;
      }
    }
  }

  return List::create(Named("gamma") = out_gamma,
                      Named("interactions") = out_interactions,
                      Named("interactions_difference") = out_interactions_difference,
                      Named("thresholds") = out_thresholds);
}