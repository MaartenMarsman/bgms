// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "gibbs_functions.h"
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// Impute missing data for two-samples designs
// ----------------------------------------------------------------------------|
List compare_impute_missing_data(NumericMatrix interactions,
                                 NumericMatrix thresholds,
                                 NumericMatrix cross_lagged,
                                 NumericMatrix pairwise_difference,
                                 NumericMatrix main_difference,
                                 IntegerMatrix observations_gr1,
                                 IntegerMatrix observations_gr2,
                                 IntegerMatrix n_cat_obs_gr1,
                                 IntegerMatrix n_cat_obs_gr2,
                                 IntegerMatrix sufficient_blume_capel_gr1,
                                 IntegerMatrix sufficient_blume_capel_gr2,
                                 IntegerVector no_categories,
                                 NumericMatrix rest_matrix_gr1,
                                 NumericMatrix rest_matrix_gr2,
                                 IntegerMatrix missing_index_gr1,
                                 IntegerMatrix missing_index_gr2,
                                 LogicalVector variable_bool,
                                 IntegerVector reference_category,
                                 bool paired) {

  int no_variables = observations_gr1.ncol();
  int no_missings_gr1 = missing_index_gr1.nrow();
  int no_missings_gr2 = missing_index_gr2.nrow();
  int max_no_categories = max(no_categories);
  NumericVector probabilities(max_no_categories + 1);
  double exponent, rest_score, cumsum, u, interaction;
  int score, person, variable, new_observation, old_observation;

  if(no_missings_gr1 > 1) {
    for(int missing = 0; missing < no_missings_gr1; missing++) {
      //Which observation to impute? ---------------------------------------------
      person = missing_index_gr1(missing, 0) - 1; //R to C++ indexing
      variable = missing_index_gr1(missing, 1) - 1; //R to C++ indexing

      //Generate new observation -------------------------------------------------
      rest_score = rest_matrix_gr1(person, variable);

      //Two distinct (ordinal) variable types ------------------------------------
      if(variable_bool[variable] == true) {

        //Regular binary or ordinal MRF variable ---------------------------------
        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < no_categories[variable]; category++) {
          exponent = thresholds(variable, category) -
            .5 * main_difference(variable, category);
          exponent += (category + 1) * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }

      } else {

        //Blume-Capel ordinal MRF variable ---------------------------------------
        exponent = (thresholds(variable, 1) - .5 * main_difference(variable, 1)) *
          reference_category[variable] *
          reference_category[variable];
        cumsum = std::exp(exponent);
        probabilities[0] = cumsum;
        for(int category = 0; category < no_categories[variable]; category++) {
          exponent = (thresholds(variable, 0) - .5 * main_difference(variable, 0))  *
            (category + 1);
          exponent += (thresholds(variable, 1) - .5 * main_difference(variable, 1))  *
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
      old_observation = observations_gr1(person, variable);
      if(old_observation != new_observation) {
        observations_gr1(person, variable) = new_observation;
        if(variable_bool[variable] == true) {
          //Regular binary or ordinal MRF variable -------------------------------
          n_cat_obs_gr1(old_observation, variable)--;
          n_cat_obs_gr1(new_observation, variable)++;
        } else {
          //Regular binary or ordinal MRF variable -------------------------------
          sufficient_blume_capel_gr1(0, variable) -= old_observation;
          sufficient_blume_capel_gr1(0, variable) += new_observation;
          sufficient_blume_capel_gr1(1, variable) -=
            (old_observation - reference_category[variable]) *
            (old_observation - reference_category[variable]);
          sufficient_blume_capel_gr1(1, variable) +=
            (new_observation - reference_category[variable]) *
            (new_observation - reference_category[variable]);
        }

        for(int vertex = 0; vertex < no_variables; vertex++) {
          interaction = interactions(vertex, variable) -
            .5 * pairwise_difference(vertex, variable);
          //interactions(i, i) = 0
          rest_matrix_gr1(person, vertex) -= old_observation *
          interaction;
          rest_matrix_gr1(person, vertex) += new_observation *
            interaction;
          if(paired == true) {
            rest_matrix_gr2(person, vertex) -=
              .25 * old_observation *
              cross_lagged(vertex, variable);
            rest_matrix_gr2(person, vertex) +=
              .25 * new_observation *
              cross_lagged(vertex, variable);
          }
        }
      }
    }
  }

  if(no_missings_gr2 > 1) {
    for(int missing = 0; missing < no_missings_gr2; missing++) {
      //Which observation to impute? ---------------------------------------------
      person = missing_index_gr2(missing, 0) - 1; //R to C++ indexing
      variable = missing_index_gr2(missing, 1) - 1; //R to C++ indexing

      //Generate new observation -------------------------------------------------
      rest_score = rest_matrix_gr2(person, variable);

      //Two distinct (ordinal) variable types ------------------------------------
      if(variable_bool[variable] == true) {

        //Regular binary or ordinal MRF variable ---------------------------------
        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < no_categories[variable]; category++) {
          exponent = thresholds(variable, category) -
            .5 * main_difference(variable, category);
          exponent += (category + 1) * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }

      } else {

        //Blume-Capel ordinal MRF variable ---------------------------------------
        exponent = (thresholds(variable, 1) - .5 * main_difference(variable, 1)) *
          reference_category[variable] *
          reference_category[variable];
        cumsum = std::exp(exponent);
        probabilities[0] = cumsum;
        for(int category = 0; category < no_categories[variable]; category++) {
          exponent = (thresholds(variable, 0) - .5 * main_difference(variable, 0))  *
            (category + 1);
          exponent += (thresholds(variable, 1) - .5 * main_difference(variable, 1))  *
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
      old_observation = observations_gr2(person, variable);
      if(old_observation != new_observation) {
        observations_gr2(person, variable) = new_observation;
        if(variable_bool[variable] == true) {
          //Regular binary or ordinal MRF variable -------------------------------
          n_cat_obs_gr2(old_observation, variable)--;
          n_cat_obs_gr2(new_observation, variable)++;
        } else {
          //Regular binary or ordinal MRF variable -------------------------------
          sufficient_blume_capel_gr2(0, variable) -= old_observation;
          sufficient_blume_capel_gr2(0, variable) += new_observation;
          sufficient_blume_capel_gr2(1, variable) -=
            (old_observation - reference_category[variable]) *
            (old_observation - reference_category[variable]);
          sufficient_blume_capel_gr2(1, variable) +=
            (new_observation - reference_category[variable]) *
            (new_observation - reference_category[variable]);
        }

        for(int vertex = 0; vertex < no_variables; vertex++) {
          interaction = interactions(vertex, variable) -
            .5 * pairwise_difference(vertex, variable);
          //interactions(i, i) = 0
          rest_matrix_gr2(person, vertex) -= old_observation *
          interaction;
          rest_matrix_gr2(person, vertex) += new_observation *
            interaction;
          if(paired == true) {
            rest_matrix_gr1(person, vertex) -=
              .25 *
              old_observation *
              cross_lagged(vertex, variable);
            rest_matrix_gr1(person, vertex) +=
              .25 *
              new_observation *
              cross_lagged(vertex, variable);
          }
        }
      }
    }
  }

  return List::create(Named("observations_gr1") = observations_gr1,
                      Named("observations_gr2") = observations_gr2,
                      Named("n_cat_obs_gr1") = n_cat_obs_gr1,
                      Named("n_cat_obs_gr2") = n_cat_obs_gr2,
                      Named("sufficient_blume_capel_gr1") = sufficient_blume_capel_gr1,
                      Named("sufficient_blume_capel_gr2") = sufficient_blume_capel_gr2,
                      Named("rest_matrix_gr1") = rest_matrix_gr1,
                      Named("rest_matrix_gr2") = rest_matrix_gr2);
}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for an interaction
//  for the two independent samples design
// ----------------------------------------------------------------------------|
double compare_log_pseudolikelihood_ratio_interaction(NumericMatrix thresholds,
                                                      NumericMatrix main_difference,
                                                      IntegerMatrix observations_gr1,
                                                      IntegerMatrix observations_gr2,
                                                      IntegerVector no_categories,
                                                      int no_persons_gr1,
                                                      int no_persons_gr2,
                                                      int variable1,
                                                      int variable2,
                                                      double proposed_state,
                                                      double current_state,
                                                      NumericMatrix rest_matrix_gr1,
                                                      NumericMatrix rest_matrix_gr2,
                                                      LogicalVector variable_bool,
                                                      IntegerVector reference_category) {

  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  double delta_state = proposed_state - current_state;

  for(int person = 0; person < no_persons_gr1; person++) {
    obs_score1 = observations_gr1(person, variable1);
    obs_score2 = observations_gr1(person, variable2);

    pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * delta_state;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix_gr1(person, variable1) -
      obs_score2 * current_state;

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
        exponent = thresholds(variable1, category) -
          .5 * main_difference(variable1, category) +
        score * rest_score -
        bound;
        denominator_prop +=
          std::exp(exponent + score * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + score * obs_score2 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable1] + 1; category++) {
        exponent = (thresholds(variable1, 0) -
          .5 * main_difference(variable1, 0)) *
        category;
        exponent += (thresholds(variable1, 1) -
          .5 * main_difference(variable1, 1)) *
        (category - reference_category[variable1]) *
        (category - reference_category[variable1]);
        exponent+= category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + category * obs_score2 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);

    //variable 2 log pseudolikelihood ratio
    rest_score = rest_matrix_gr1(person, variable2) -
      obs_score1 * current_state;

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
        exponent = (thresholds(variable2, category) -
          .5 * main_difference(variable2, category)) +
        score * rest_score -
        bound;
        denominator_prop +=
          std::exp(exponent + score * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent + score * obs_score1 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable2] + 1; category++) {
        exponent = (thresholds(variable2, 0) -
          .5 * main_difference(variable2, 0)) *
        category;
        exponent += (thresholds(variable2, 1) -
          .5 * main_difference(variable2, 1)) *
        (category - reference_category[variable2]) *
        (category - reference_category[variable2]);
        exponent+=  category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent + category * obs_score1 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);
  }

  for(int person = 0; person < no_persons_gr2; person++) {
    obs_score1 = observations_gr2(person, variable1);
    obs_score2 = observations_gr2(person, variable2);

    pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * delta_state;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix_gr2(person, variable1) -
      obs_score2 * current_state;

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
          .5 * main_difference(variable1, category) +
        score * rest_score -
        bound;
        denominator_prop +=
          std::exp(exponent + score * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + score * obs_score2 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable1] + 1; category++) {
        exponent = (thresholds(variable1, 0) +
          .5 * main_difference(variable1, 0)) *
        category;
        exponent += (thresholds(variable1, 1) +
          .5 * main_difference(variable1, 1)) *
        (category - reference_category[variable1]) *
        (category - reference_category[variable1]);
        exponent+= category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + category * obs_score2 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);

    //variable 2 log pseudolikelihood ratio
    rest_score = rest_matrix_gr2(person, variable2) -
      obs_score1 * current_state;

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
        exponent = (thresholds(variable2, category) +
          .5 * main_difference(variable2, category)) +
        score * rest_score -
        bound;
        denominator_prop +=
          std::exp(exponent + score * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent + score * obs_score1 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable2] + 1; category++) {
        exponent = (thresholds(variable2, 0) +
          .5 * main_difference(variable2, 0)) *
        category;
        exponent += (thresholds(variable2, 1) +
          .5 * main_difference(variable2, 1)) *
        (category - reference_category[variable2]) *
        (category - reference_category[variable2]);
        exponent+=  category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent + category * obs_score1 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);
  }

  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of an overall pairwise
//  interaction parameter (nuisance)
// ----------------------------------------------------------------------------|
void compare_metropolis_interaction(NumericMatrix interactions,
                                    NumericMatrix thresholds,
                                    NumericMatrix main_difference,
                                    IntegerMatrix observations_gr1,
                                    IntegerMatrix observations_gr2,
                                    IntegerVector no_categories,
                                    NumericMatrix proposal_sd_interaction,
                                    double interaction_scale,
                                    int no_persons_gr1,
                                    int no_persons_gr2,
                                    int no_variables,
                                    NumericMatrix rest_matrix_gr1,
                                    NumericMatrix rest_matrix_gr2,
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
      proposed_state = R::rnorm(current_state, proposal_sd_interaction(variable1, variable2));

      log_prob = compare_log_pseudolikelihood_ratio_interaction(thresholds,
                                                                main_difference,
                                                                observations_gr1,
                                                                observations_gr2,
                                                                no_categories,
                                                                no_persons_gr1,
                                                                no_persons_gr2,
                                                                variable1,
                                                                variable2,
                                                                proposed_state,
                                                                current_state,
                                                                rest_matrix_gr1,
                                                                rest_matrix_gr2,
                                                                variable_bool,
                                                                reference_category);

      log_prob += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_prob -= R::dcauchy(current_state, 0.0, interaction_scale, true);

      U = R::unif_rand();
      if(std::log(U) < log_prob) {
        double state_difference = proposed_state - current_state;
        interactions(variable1, variable2) = proposed_state;
        interactions(variable2, variable1) = proposed_state;

        //Update the matrix of rest scores
        for(int person = 0; person < no_persons_gr1; person++) {
          rest_matrix_gr1(person, variable1) += observations_gr1(person, variable2) *
            state_difference;
          rest_matrix_gr1(person, variable2) += observations_gr1(person, variable1) *
            state_difference;
        }
        //Update the matrix of rest scores
        for(int person = 0; person < no_persons_gr2; person++) {
          rest_matrix_gr2(person, variable1) += observations_gr2(person, variable2) *
            state_difference;
          rest_matrix_gr2(person, variable2) += observations_gr2(person, variable1) *
            state_difference;
        }
      }

      if(log_prob > 0) {
        log_prob = 1.0;
      } else {
        log_prob = std::exp(log_prob);
      }

      double update_proposal_sd =
        proposal_sd_interaction(variable1, variable2) +
        (log_prob - target_ar) * std::exp(-log(t) * phi);

      if(std::isnan(update_proposal_sd) == true) {
        update_proposal_sd = 1.0;
      }

      if(update_proposal_sd < epsilon_lo) {
        update_proposal_sd = epsilon_lo;
      } else if (update_proposal_sd > epsilon_hi) {
        update_proposal_sd = epsilon_hi;
      }

      proposal_sd_interaction(variable1, variable2) = update_proposal_sd;
    }
  }
}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for the difference
//  in a pairwise interaction for a two independent samples design
// ----------------------------------------------------------------------------|
double compare_log_pseudolikelihood_ratio_pairwise_difference(NumericMatrix thresholds,
                                                              NumericMatrix main_difference,
                                                              IntegerMatrix observations_gr1,
                                                              IntegerMatrix observations_gr2,
                                                              IntegerVector no_categories,
                                                              int no_persons_gr1,
                                                              int no_persons_gr2,
                                                              int variable1,
                                                              int variable2,
                                                              double proposed_state,
                                                              double current_state,
                                                              NumericMatrix rest_matrix_gr1,
                                                              NumericMatrix rest_matrix_gr2,
                                                              LogicalVector variable_bool,
                                                              IntegerVector reference_category) {

  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  double delta_state = proposed_state - current_state;

  for(int person = 0; person < no_persons_gr1; person++) {
    obs_score1 = observations_gr1(person, variable1);
    obs_score2 = observations_gr1(person, variable2);

    pseudolikelihood_ratio -= obs_score1 * obs_score2 * delta_state;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix_gr1(person, variable1) +
      .5 * obs_score2 * current_state;

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
        exponent = thresholds(variable1, category) -
          .5 * main_difference(variable1, category) +
        score * rest_score -
        bound;
        denominator_prop +=
          std::exp(exponent - score * .5 * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent - score * .5 * obs_score2 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable1] + 1; category++) {
        exponent = (thresholds(variable1, 0) -
          .5 * main_difference(variable1, 0)) *
        category;
        exponent += (thresholds(variable1, 1) -
          .5 * main_difference(variable1, 1)) *
        (category - reference_category[variable1]) *
        (category - reference_category[variable1]);
        exponent+= category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent - category * .5 * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent - category * .5 * obs_score2 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);

    //variable 2 log pseudolikelihood ratio
    rest_score = rest_matrix_gr1(person, variable2) +
      .5 * obs_score1 * current_state;

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
        exponent = (thresholds(variable2, category) -
          .5 * main_difference(variable2, category)) +
        score * rest_score -
        bound;
        denominator_prop +=
          std::exp(exponent - score * .5 * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent - score * .5 * obs_score1 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable2] + 1; category++) {
        exponent = (thresholds(variable2, 0) -
          .5 * main_difference(variable2, 0)) *
        category;
        exponent += (thresholds(variable2, 1) -
          .5 * main_difference(variable2, 1)) *
        (category - reference_category[variable2]) *
        (category - reference_category[variable2]);
        exponent +=  category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent - category * .5 * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent - category * .5 * obs_score1 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);
  }


  for(int person = 0; person < no_persons_gr2; person++) {
    obs_score1 = observations_gr2(person, variable1);
    obs_score2 = observations_gr2(person, variable2);

    pseudolikelihood_ratio += obs_score1 * obs_score2 * delta_state;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix_gr2(person, variable1) -
      .5 * obs_score2 * current_state;

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
          .5 * main_difference(variable1, category) +
        score * rest_score -
        bound;
        denominator_prop +=
          std::exp(exponent + score * .5 * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + score * .5 * obs_score2 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable1] + 1; category++) {
        exponent = (thresholds(variable1, 0) +
          .5 * main_difference(variable1, 0)) *
        category;
        exponent += (thresholds(variable1, 1) +
          .5 * main_difference(variable1, 1)) *
        (category - reference_category[variable1]) *
        (category - reference_category[variable1]);
        exponent+= category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * .5 * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + category * .5 * obs_score2 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);

    //variable 2 log pseudolikelihood ratio
    rest_score = rest_matrix_gr2(person, variable2) -
      .5 * obs_score1 * current_state;

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
        exponent = (thresholds(variable2, category) +
          .5 * main_difference(variable2, category)) +
        score * rest_score -
        bound;
        denominator_prop +=
          std::exp(exponent + score * .5 * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent + score * .5 * obs_score1 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable2] + 1; category++) {
        exponent = (thresholds(variable2, 0) +
          .5 * main_difference(variable2, 0)) *
        category;
        exponent += (thresholds(variable2, 1) +
          .5 * main_difference(variable2, 1)) *
        (category - reference_category[variable2]) *
        (category - reference_category[variable2]);
        exponent+=  category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * .5 * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent + category * .5 * obs_score1 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);
  }

  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the active differences
//  in the pairwise interaction parameters
// ----------------------------------------------------------------------------|
void compare_metropolis_pairwise_difference(NumericMatrix thresholds,
                                            NumericMatrix pairwise_difference,
                                            NumericMatrix main_difference,
                                            IntegerMatrix observations_gr1,
                                            IntegerMatrix observations_gr2,
                                            IntegerVector no_categories,
                                            IntegerMatrix indicator,
                                            NumericMatrix proposal_sd_pairwise_difference,
                                            double pairwise_difference_scale,
                                            int no_persons_gr1,
                                            int no_persons_gr2,
                                            int no_variables,
                                            NumericMatrix rest_matrix_gr1,
                                            NumericMatrix rest_matrix_gr2,
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
      if(indicator(variable1, variable2) == 1) {
        current_state = pairwise_difference(variable1, variable2);
        proposed_state = R::rnorm(current_state, proposal_sd_pairwise_difference(variable1, variable2));

        log_prob = compare_log_pseudolikelihood_ratio_pairwise_difference(thresholds,
                                                                          main_difference,
                                                                          observations_gr1,
                                                                          observations_gr2,
                                                                          no_categories,
                                                                          no_persons_gr1,
                                                                          no_persons_gr2,
                                                                          variable1,
                                                                          variable2,
                                                                          proposed_state,
                                                                          current_state,
                                                                          rest_matrix_gr1,
                                                                          rest_matrix_gr2,
                                                                          variable_bool,
                                                                          reference_category);


        log_prob += R::dcauchy(proposed_state, 0.0, pairwise_difference_scale, true);
        log_prob -= R::dcauchy(current_state, 0.0, pairwise_difference_scale, true);

        U = R::unif_rand();
        if(std::log(U) < log_prob) {
          double state_difference = .5 * (proposed_state - current_state);
          pairwise_difference(variable1, variable2) = proposed_state;
          pairwise_difference(variable2, variable1) = proposed_state;

          //Update the matrices of rest scores
          for(int person = 0; person < no_persons_gr1; person++) {
            rest_matrix_gr1(person, variable1) -= observations_gr1(person, variable2) *
              state_difference;
            rest_matrix_gr1(person, variable2) -= observations_gr1(person, variable1) *
              state_difference;
          }
          for(int person = 0; person < no_persons_gr2; person++) {
            rest_matrix_gr2(person, variable1) += observations_gr2(person, variable2) *
              state_difference;
            rest_matrix_gr2(person, variable2) += observations_gr2(person, variable1) *
              state_difference;
          }
        }

        // Update the Robbins-Monro parameters
        if(log_prob > 0) {
          log_prob = 1;
        } else {
          log_prob = std::exp(log_prob);
        }

        double update_proposal_sd =
          proposal_sd_pairwise_difference(variable1, variable2) +
          (log_prob - target_ar) * std::exp(-log(t) * phi);

        if(std::isnan(update_proposal_sd) == true) {
          update_proposal_sd = 1.0;
        }

        if(update_proposal_sd < epsilon_lo) {
          update_proposal_sd = epsilon_lo;
        } else if (update_proposal_sd > epsilon_hi) {
          update_proposal_sd = epsilon_hi;
        }

        proposal_sd_pairwise_difference(variable1, variable2) = update_proposal_sd;
      }
    }
  }
}

// ----------------------------------------------------------------------------|
// Between model MH algorithm for the pairwise differences
// ----------------------------------------------------------------------------|
void compare_metropolis_pairwise_difference_between_model(NumericMatrix thresholds,
                                                          NumericMatrix pairwise_difference,
                                                          NumericMatrix main_difference,
                                                          IntegerMatrix observations_gr1,
                                                          IntegerMatrix observations_gr2,
                                                          IntegerVector no_categories,
                                                          IntegerVector indicator,
                                                          NumericMatrix proposal_sd_pairwise_difference,
                                                          double pairwise_difference_scale,
                                                          NumericMatrix inclusion_probability_difference,
                                                          int no_persons_gr1,
                                                          int no_persons_gr2,
                                                          int no_interactions,
                                                          IntegerMatrix index,
                                                          NumericMatrix rest_matrix_gr1,
                                                          NumericMatrix rest_matrix_gr2,
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

  int variable1;
  int variable2;

  for(int cntr = 0; cntr < no_interactions; cntr ++) {
    variable1 = index(cntr, 1) - 1;
    variable2 = index(cntr, 2) - 1;

    current_state = pairwise_difference(variable1, variable2);

    if(indicator(variable1, variable2) == 0) {
      proposed_state = R::rnorm(current_state, proposal_sd_pairwise_difference(variable1, variable2));
    } else {
      proposed_state = 0.0;
    }

    log_prob = compare_log_pseudolikelihood_ratio_pairwise_difference(thresholds,
                                                                      main_difference,
                                                                      observations_gr1,
                                                                      observations_gr2,
                                                                      no_categories,
                                                                      no_persons_gr1,
                                                                      no_persons_gr2,
                                                                      variable1,
                                                                      variable2,
                                                                      proposed_state,
                                                                      current_state,
                                                                      rest_matrix_gr1,
                                                                      rest_matrix_gr2,
                                                                      variable_bool,
                                                                      reference_category);

    if(indicator(variable1, variable2) == 0) {
      log_prob += R::dcauchy(proposed_state, 0.0, pairwise_difference_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_sd_pairwise_difference(variable1, variable2),
                           true);

      log_prob += log(inclusion_probability_difference(variable1, variable2));
      log_prob -= log(1 - inclusion_probability_difference(variable1, variable2));

    } else {
      log_prob -= R::dcauchy(current_state, 0.0, pairwise_difference_scale, true);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_sd_pairwise_difference(variable1, variable2),
                           true);

      log_prob -= log(inclusion_probability_difference(variable1, variable2));
      log_prob += log(1 - inclusion_probability_difference(variable1, variable2));
    }

    U = R::unif_rand();
    if(std::log(U) < log_prob) {
      indicator(variable1, variable2) = 1 - indicator(variable1, variable2);
      indicator(variable2, variable1) = 1 - indicator(variable2, variable1);

      pairwise_difference(variable1, variable2) = proposed_state;
      pairwise_difference(variable2, variable1) = proposed_state;

      double state_difference = .5 * (proposed_state - current_state);

      //Update the matrices of rest scores
      for(int person = 0; person < no_persons_gr1; person++) {
        rest_matrix_gr1(person, variable1) -= observations_gr1(person, variable2) *
          state_difference;
        rest_matrix_gr1(person, variable2) -= observations_gr1(person, variable1) *
          state_difference;
      }
      for(int person = 0; person < no_persons_gr2; person++) {
        rest_matrix_gr2(person, variable1) += observations_gr2(person, variable2) *
          state_difference;
        rest_matrix_gr2(person, variable2) += observations_gr2(person, variable1) *
          state_difference;
      }
    }
  }
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the overall category
//  threshold parameters (nuisance) -- regular ordinal variable
// ----------------------------------------------------------------------------|
void compare_metropolis_threshold_regular(NumericMatrix thresholds,
                                          NumericMatrix main_difference,
                                          IntegerVector no_categories,
                                          IntegerMatrix n_cat_obs_gr1,
                                          IntegerMatrix n_cat_obs_gr2,
                                          int no_persons_gr1,
                                          int no_persons_gr2,
                                          int variable,
                                          double threshold_alpha,
                                          double threshold_beta,
                                          NumericMatrix rest_matrix_gr1,
                                          NumericMatrix rest_matrix_gr2) {

  NumericVector g(no_persons_gr1 + no_persons_gr2);
  NumericVector q(no_persons_gr1 + no_persons_gr2);

  double log_prob, rest_score;
  double a, b, c;
  double tmp;
  double current_state, proposed_state;
  double U;
  double exp_current, exp_proposed;

  for(int category = 0; category < no_categories[variable]; category++) {
    current_state = thresholds(variable, category);
    exp_current = std::exp(current_state);
    c = (threshold_alpha + threshold_beta) / (1 + exp_current);

    for(int person = 0; person < no_persons_gr1; person++) {
      g[person] = 1.0;
      q[person] = 1.0;
      rest_score = rest_matrix_gr1(person, variable);
      for(int cat = 0; cat < no_categories[variable]; cat++) {
        if(cat != category) {
          g[person] +=
            std::exp(thresholds(variable, cat) -
            .5 * main_difference(variable, cat) +
            (cat + 1) * rest_score);
        }
      }
      q[person] = std::exp((category + 1) * rest_score);
      c +=  q[person] / (g[person] + q[person] * exp_current);
    }
    for(int person = 0; person < no_persons_gr2; person++) {
      g[person + no_persons_gr1] = 1.0;
      q[person + no_persons_gr1] = 1.0;
      rest_score = rest_matrix_gr2(person, variable);
      for(int cat = 0; cat < no_categories[variable]; cat++) {
        if(cat != category) {
          g[person + no_persons_gr1] +=
            std::exp(thresholds(variable, cat) +
            .5 * main_difference(variable, cat) +
            (cat + 1) * rest_score);
        }
      }
      q[person + no_persons_gr1] = std::exp((category + 1) * rest_score);
      c +=  q[person + no_persons_gr1] /
        (g[person + no_persons_gr1] + q[person + no_persons_gr1] * exp_current);
    }
    c = c / ((no_persons_gr1 + no_persons_gr2 + threshold_alpha + threshold_beta) -
      exp_current * c);

    //Proposal is generalized beta-prime.
    a = n_cat_obs_gr1(category + 1, variable) +
      n_cat_obs_gr2(category + 1, variable) +
      threshold_alpha;
    b = no_persons_gr1 +
      no_persons_gr2 +
      threshold_beta -
      n_cat_obs_gr1(category + 1, variable) -
      n_cat_obs_gr2(category + 1, variable);
    tmp = R::rbeta(a, b);
    proposed_state = std::log(tmp / (1  - tmp) / c);
    exp_proposed = exp(proposed_state);

    //Compute log_acceptance probability for Metropolis.
    //First, we use g and q above to compute the ratio of pseudolikelihoods
    log_prob = 0;
    for(int person = 0; person < no_persons_gr2 + no_persons_gr1; person++) {
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
      thresholds(variable, category) = proposed_state;
    }
  }

}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for a category
//  threshold difference for the two independent samples design
// ----------------------------------------------------------------------------|
double compare_log_pseudolikelihood_ratio_main_difference(NumericMatrix thresholds,
                                                          NumericMatrix main_difference,
                                                          IntegerMatrix n_cat_obs_gr1,
                                                          IntegerMatrix n_cat_obs_gr2,
                                                          IntegerVector no_categories,
                                                          int no_persons_gr1,
                                                          int no_persons_gr2,
                                                          int variable,
                                                          int category,
                                                          double proposed_state,
                                                          double current_state,
                                                          NumericMatrix rest_matrix_gr1,
                                                          NumericMatrix rest_matrix_gr2) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_proposed, denominator_current, exponent;
  int score;
  double delta_state = proposed_state - current_state;

  NumericVector proposed_thresholds(no_categories[variable]);
  NumericVector current_thresholds(no_categories[variable]);

  //Pre-compute vector of threshold differences for Group 1:
  for(int cat = 0; cat < no_categories[variable]; cat++) {
    proposed_thresholds[cat] =
      thresholds(variable, cat) -
      .5 * main_difference(variable, cat);
      current_thresholds[cat] = proposed_thresholds[cat];
  }
  proposed_thresholds[category] =
    thresholds(variable, category) -
    .5 * proposed_state;

  //Compute the pseudolikelihood ratio
  pseudolikelihood_ratio -= delta_state * .5 * n_cat_obs_gr1(category + 1, variable);
  pseudolikelihood_ratio += delta_state * .5 * n_cat_obs_gr2(category + 1, variable);

  //Compute the contribution for group 1
  for(int person = 0; person < no_persons_gr1; person++) {
    rest_score = rest_matrix_gr1(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score;
    } else {
      bound = 0.0;
    }

    denominator_proposed = std::exp(-bound);
    denominator_current = std::exp(-bound);
    for(int cat = 0; cat < no_categories[variable]; cat++) {
      score = cat + 1;
      exponent = score * rest_score - bound;
      denominator_proposed += std::exp(exponent + proposed_thresholds[cat]);
      denominator_current += std::exp(exponent + current_thresholds[cat]);
    }

    pseudolikelihood_ratio -= std::log(denominator_proposed);
    pseudolikelihood_ratio += std::log(denominator_current);
  }

  //Pre-compute vector of threshold differences for Group 2:
  for(int cat = 0; cat < no_categories[variable]; cat++) {
    proposed_thresholds[cat] =
      thresholds(variable, cat) +
      .5 * main_difference(variable, cat);
      current_thresholds[cat] = proposed_thresholds[cat];
  }
  proposed_thresholds[category] =
    thresholds(variable, category) +
    .5 * proposed_state;

  //Compute the contribution for group 2
  for(int person = 0; person < no_persons_gr2; person++) {
    rest_score = rest_matrix_gr2(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score;
    } else {
      bound = 0.0;
    }

    denominator_proposed = std::exp(-bound);
    denominator_current = std::exp(-bound);
    for(int cat = 0; cat < no_categories[variable]; cat++) {
      score = cat + 1;
      exponent = score * rest_score - bound;
      denominator_proposed += std::exp(exponent + proposed_thresholds[cat]);
      denominator_current += std::exp(exponent + current_thresholds[cat]);
    }

    pseudolikelihood_ratio -= std::log(denominator_proposed);
    pseudolikelihood_ratio += std::log(denominator_current);
  }
  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the category threshold
//  difference parameter -- regular ordinal variable
// ----------------------------------------------------------------------------|
void compare_metropolis_main_difference_regular(NumericMatrix thresholds,
                                                NumericMatrix main_difference,
                                                IntegerMatrix n_cat_obs_gr1,
                                                IntegerMatrix n_cat_obs_gr2,
                                                IntegerVector no_categories,
                                                IntegerMatrix indicator,
                                                NumericMatrix proposal_sd_main_difference,
                                                double main_difference_scale,
                                                int no_persons_gr1,
                                                int no_persons_gr2,
                                                int variable,
                                                NumericMatrix rest_matrix_gr1,
                                                NumericMatrix rest_matrix_gr2,
                                                double phi,
                                                double target_ar,
                                                int t,
                                                double epsilon_lo,
                                                double epsilon_hi) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  if(indicator(variable, variable) == 1) {
    for(int category = 0; category < no_categories[variable]; category++) {
      current_state = main_difference(variable, category);
      proposed_state = R::rnorm(current_state,
                                proposal_sd_main_difference(variable, category));

      log_prob = compare_log_pseudolikelihood_ratio_main_difference(thresholds,
                                                                      main_difference,
                                                                      n_cat_obs_gr1,
                                                                      n_cat_obs_gr2,
                                                                      no_categories,
                                                                      no_persons_gr1,
                                                                      no_persons_gr2,
                                                                      variable,
                                                                      category,
                                                                      proposed_state,
                                                                      current_state,
                                                                      rest_matrix_gr1,
                                                                      rest_matrix_gr2);
      log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
      log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);

      U = R::unif_rand();
      if(std::log(U) < log_prob) {
        main_difference(variable, category) = proposed_state;
      }

     // Update the Robbins-Monro parameters
     if(log_prob > 0) {
       log_prob = 1;
     } else {
       log_prob = std::exp(log_prob);
     }

     double update_proposal_sd =
       proposal_sd_main_difference(variable, category) +
       (log_prob - target_ar) * std::exp(-log(t) * phi);

     if(std::isnan(update_proposal_sd) == true) {
       update_proposal_sd = 1.0;
     }

     if(update_proposal_sd < epsilon_lo) {
       update_proposal_sd = epsilon_lo;
     } else if (update_proposal_sd > epsilon_hi) {
       update_proposal_sd = epsilon_hi;
     }

     proposal_sd_main_difference(variable, category) = update_proposal_sd;
    }
  }
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the overall category
//  threshold parameters (nuisance) -- Blume-Capel ordinal variable
// ----------------------------------------------------------------------------|
void compare_metropolis_threshold_blumecapel(NumericMatrix thresholds,
                                             NumericMatrix main_difference,
                                             IntegerVector no_categories,
                                             IntegerMatrix sufficient_blume_capel_gr1,
                                             IntegerMatrix sufficient_blume_capel_gr2,
                                             int no_persons_gr1,
                                             int no_persons_gr2,
                                             int variable,
                                             IntegerVector reference_category,
                                             double threshold_alpha,
                                             double threshold_beta,
                                             NumericMatrix rest_matrix_gr1,
                                             NumericMatrix rest_matrix_gr2,
                                             NumericMatrix proposal_sd_blumecapel,
                                             double phi,
                                             double target_ar,
                                             int t,
                                             double epsilon_lo,
                                             double epsilon_hi) {

  double log_prob, U;
  double current_state, proposed_state;
  double numerator, denominator;
  double lbound, bound, exponent, rest_score;
  NumericVector constant_numerator (no_categories[variable] + 1);
  NumericVector constant_denominator (no_categories[variable] + 1);

  //--------------------------------------------------------------------------
  // Adaptive Metropolis for the linear Blume-Capel parameter
  //--------------------------------------------------------------------------
  current_state = thresholds(variable, 0);
  proposed_state = R::rnorm(current_state,
                            proposal_sd_blumecapel(variable, 0));

  //--------------------------------------------------------------------------
  //Compute the log acceptance probability -----------------------------------
  //--------------------------------------------------------------------------

  //Compute the prior ratio --------------------------------------------------
  log_prob = threshold_alpha * proposed_state;
  log_prob -= threshold_alpha * current_state;
  log_prob -= (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(proposed_state));
  log_prob += (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(current_state));

  //Precompute bounds for group 1 for numeric stability and efficiency -------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    exponent = (thresholds(variable, 1) - .5 * main_difference(variable, 1)) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_numerator[category] =
      (current_state - .5 * main_difference(variable, 0)) *
      category + exponent;
    constant_denominator[category] =
      (proposed_state - .5 * main_difference(variable, 0)) *
      category + exponent;
  }

  double tmp_n = max(constant_numerator);
  double tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log likelihood ratio group 1----------------------------------
  log_prob += sufficient_blume_capel_gr1(0, variable) * proposed_state;
  log_prob -= sufficient_blume_capel_gr1(0, variable) * current_state;

  for(int person = 0; person < no_persons_gr1; person++) {
    rest_score = rest_matrix_gr1(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  //Precompute bounds for group 2 for numeric stability and efficiency -------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    exponent = (thresholds(variable, 1) + .5 * main_difference(variable, 1)) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_numerator[category] =
      (current_state + .5 * main_difference(variable, 0)) *
      category + exponent;
    constant_denominator[category] =
      (proposed_state + .5 * main_difference(variable, 0)) *
      category + exponent;
  }

  tmp_n = max(constant_numerator);
  tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log likelihood ratio group 2----------------------------------
  log_prob += sufficient_blume_capel_gr2(0, variable) * proposed_state;
  log_prob -= sufficient_blume_capel_gr2(0, variable) * current_state;

  for(int person = 0; person < no_persons_gr2; person++) {
    rest_score = rest_matrix_gr2(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    thresholds(variable, 0) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ----------------------------
  if(log_prob > 0) {
    log_prob = 1;
  } else {
    log_prob = std::exp(log_prob);
  }

  double update_proposal_sd = proposal_sd_blumecapel(variable, 0) +
    (log_prob - target_ar) * std::exp(-log(t) * phi);

  if(std::isnan(update_proposal_sd) == true) {
    update_proposal_sd = 1.0;
  }

  if(update_proposal_sd < epsilon_lo) {
    update_proposal_sd = epsilon_lo;
  } else if (update_proposal_sd > epsilon_hi) {
    update_proposal_sd = epsilon_hi;
  }

  proposal_sd_blumecapel(variable, 0) = update_proposal_sd;


  //--------------------------------------------------------------------------
  // Adaptive Metropolis for the quadratic Blume-Capel parameter
  //--------------------------------------------------------------------------
  current_state = thresholds(variable, 1);
  proposed_state = R::rnorm(current_state,
                            proposal_sd_blumecapel(variable, 1));

  //--------------------------------------------------------------------------
  // Compute the log acceptance probability ----------------------------------
  //--------------------------------------------------------------------------

  //Compute the prior ratio --------------------------------------------------
  log_prob = threshold_alpha * proposed_state;
  log_prob -= threshold_alpha * current_state;
  log_prob -= (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(proposed_state));
  log_prob += (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(current_state));


  //Precompute bounds for group 1 for numeric stability and efficiency -------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    exponent =  thresholds(variable, 0) * category;
    exponent -= .5 * main_difference(variable, 0) * category;
    constant_numerator[category] = exponent +
      (current_state - .5 * main_difference(variable, 1)) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_denominator[category] = exponent +
      (proposed_state - .5 * main_difference(variable, 1)) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
  }

  tmp_n = max(constant_numerator);
  tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log likelihood ratio group 1----------------------------------
  log_prob += sufficient_blume_capel_gr1(1, variable) * proposed_state;
  log_prob -= sufficient_blume_capel_gr1(1, variable) * current_state;

  for(int person = 0; person < no_persons_gr1; person++) {
    rest_score = rest_matrix_gr1(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }


  //Precompute bounds for group 2 for numeric stability and efficiency -------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    exponent =  thresholds(variable, 0) * category;
    exponent += .5 * main_difference(variable, 0) * category;
    constant_numerator[category] = exponent +
      (current_state + .5 * main_difference(variable, 1)) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_denominator[category] = exponent +
      (proposed_state + .5 * main_difference(variable, 1)) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
  }

  tmp_n = max(constant_numerator);
  tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log likelihood ratio group 2----------------------------------
  log_prob += sufficient_blume_capel_gr2(1, variable) * proposed_state;
  log_prob -= sufficient_blume_capel_gr2(1, variable) * current_state;

  for(int person = 0; person < no_persons_gr2; person++) {
    rest_score = rest_matrix_gr2(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    thresholds(variable, 1) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ----------------------------
  if(log_prob > 0) {
    log_prob = 1;
  } else {
    log_prob = std::exp(log_prob);
  }

  update_proposal_sd = proposal_sd_blumecapel(variable, 1) +
    (log_prob - target_ar) * std::exp(-log(t) * phi);

  if(std::isnan(update_proposal_sd) == true) {
    update_proposal_sd = 1.0;
  }

  if(update_proposal_sd < epsilon_lo) {
    update_proposal_sd = epsilon_lo;
  } else if (update_proposal_sd > epsilon_hi) {
    update_proposal_sd = epsilon_hi;
  }

  proposal_sd_blumecapel(variable, 1) = update_proposal_sd;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the category threshold
//  difference parameters-- Blume-Capel ordinal variable
// ----------------------------------------------------------------------------|
void compare_metropolis_main_difference_blumecapel(NumericMatrix thresholds,
                                                   NumericMatrix main_difference,
                                                   IntegerVector no_categories,
                                                   IntegerMatrix sufficient_blume_capel_gr1,
                                                   IntegerMatrix sufficient_blume_capel_gr2,
                                                   int no_persons_gr1,
                                                   int no_persons_gr2,
                                                   int variable,
                                                   IntegerVector reference_category,
                                                   double main_difference_scale,
                                                   NumericMatrix rest_matrix_gr1,
                                                   NumericMatrix rest_matrix_gr2,
                                                   NumericMatrix proposal_sd_main_difference,
                                                   double phi,
                                                   double target_ar,
                                                   int t,
                                                   double epsilon_lo,
                                                   double epsilon_hi) {

  double log_prob, U;
  double current_state, proposed_state;
  double numerator, denominator;
  double lbound, bound, exponent, rest_score;
  NumericVector constant_numerator (no_categories[variable] + 1);
  NumericVector constant_denominator (no_categories[variable] + 1);

  //--------------------------------------------------------------------------
  // Adaptive Metropolis for the difference in the linear Blume-Capel parameter
  //--------------------------------------------------------------------------
  current_state = main_difference(variable, 0);
  proposed_state = R::rnorm(current_state,
                            proposal_sd_main_difference(variable, 0));

  //--------------------------------------------------------------------------
  // Compute the log acceptance probability -----------------------------------
  //--------------------------------------------------------------------------

  //Compute the prior ratio --------------------------------------------------
  log_prob = R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
  log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);

  //Precompute bounds for group 1 for numeric stability and efficiency -------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    exponent = (thresholds(variable, 1) - .5 * main_difference(variable, 1)) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_numerator[category] =
      (thresholds(variable, 0) - .5 * current_state) *
      category + exponent;
    constant_denominator[category] =
      (thresholds(variable, 0) - .5 * proposed_state) *
      category + exponent;
  }

  double tmp_n = max(constant_numerator);
  double tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log likelihood ratio group 1----------------------------------
  log_prob -= .5 * sufficient_blume_capel_gr1(0, variable) * proposed_state;
  log_prob += .5 * sufficient_blume_capel_gr1(0, variable) * current_state;

  for(int person = 0; person < no_persons_gr1; person++) {
    rest_score = rest_matrix_gr1(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  //Precompute bounds for group 2 for numeric stability and efficiency -------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    exponent = (thresholds(variable, 1) + .5 * main_difference(variable, 1)) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_numerator[category] =
      (thresholds(variable, 0) + .5 * current_state) *
      category + exponent;
    constant_denominator[category] =
      (thresholds(variable, 0) + .5 * proposed_state) *
      category + exponent;
  }

  tmp_n = max(constant_numerator);
  tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log likelihood ratio group 2 ---------------------------------
  log_prob += .5 * sufficient_blume_capel_gr2(0, variable) * proposed_state;
  log_prob -= .5 * sufficient_blume_capel_gr2(0, variable) * current_state;

  for(int person = 0; person < no_persons_gr2; person++) {
    rest_score = rest_matrix_gr2(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    main_difference(variable, 0) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ----------------------------
  if(log_prob > 0) {
    log_prob = 1;
  } else {
    log_prob = std::exp(log_prob);
  }

  double update_proposal_sd = proposal_sd_main_difference(variable, 0) +
    (log_prob - target_ar) * std::exp(-log(t) * phi);

  if(std::isnan(update_proposal_sd) == true) {
    update_proposal_sd = 1.0;
  }

  if(update_proposal_sd < epsilon_lo) {
    update_proposal_sd = epsilon_lo;
  } else if (update_proposal_sd > epsilon_hi) {
    update_proposal_sd = epsilon_hi;
  }

  proposal_sd_main_difference(variable, 0) = update_proposal_sd;

  //--------------------------------------------------------------------------
  // Adaptive Metropolis for the difference in quadratic Blume-Capel parameter
  //--------------------------------------------------------------------------
  current_state = main_difference(variable, 1);
  proposed_state = R::rnorm(current_state,
                            proposal_sd_main_difference(variable, 1));

  //--------------------------------------------------------------------------
  // Compute the log acceptance probability ----------------------------------
  //--------------------------------------------------------------------------

  // Compute the prior ratio -------------------------------------------------
  log_prob = R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
  log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);

  //Precompute bounds for group 1 for numeric stability and efficiency -------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    exponent =  thresholds(variable, 0) * category;
    exponent -= .5 * main_difference(variable, 0) * category;
    constant_numerator[category] = exponent +
      (thresholds(variable, 1) - .5 * current_state) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_denominator[category] = exponent +
      (thresholds(variable, 1) - .5 * current_state) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
  }

  tmp_n = max(constant_numerator);
  tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log likelihood ratio group 1----------------------------------
  log_prob -= .5 * sufficient_blume_capel_gr1(1, variable) * proposed_state;
  log_prob += .5 * sufficient_blume_capel_gr1(1, variable) * current_state;

  for(int person = 0; person < no_persons_gr1; person++) {
    rest_score = rest_matrix_gr1(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  //Precompute bounds for group 2 for numeric stability and efficiency -------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    exponent =  thresholds(variable, 0) * category;
    exponent += .5 * main_difference(variable, 0) * category;
    constant_numerator[category] = exponent +
      (thresholds(variable, 1) + .5 * current_state) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_denominator[category] = exponent +
      (thresholds(variable, 1) + .5 * current_state) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
  }

  tmp_n = max(constant_numerator);
  tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log likelihood ratio group 2 ---------------------------------
  log_prob += .5 * sufficient_blume_capel_gr2(1, variable) * proposed_state;
  log_prob -= .5 * sufficient_blume_capel_gr2(1, variable) * current_state;

  for(int person = 0; person < no_persons_gr2; person++) {
    rest_score = rest_matrix_gr2(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    main_difference(variable, 1) = proposed_state;
  }
  //Robbins-Monro update of the proposal variance ----------------------------
  if(log_prob > 0) {
    log_prob = 1;
  } else {
    log_prob = std::exp(log_prob);
  }

  update_proposal_sd = proposal_sd_main_difference(variable, 1) +
    (log_prob - target_ar) * std::exp(-log(t) * phi);

  if(std::isnan(update_proposal_sd) == true) {
    update_proposal_sd = 1.0;
  }

  if(update_proposal_sd < epsilon_lo) {
    update_proposal_sd = epsilon_lo;
  } else if (update_proposal_sd > epsilon_hi) {
    update_proposal_sd = epsilon_hi;
  }

  proposal_sd_main_difference(variable, 1) = update_proposal_sd;
}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for a vector of
// category threshold differences for the two independent samples design
// ----------------------------------------------------------------------------|
double compare_log_pseudolikelihood_ratio_main_differences(NumericMatrix thresholds,
                                                           IntegerMatrix n_cat_obs_gr1,
                                                           IntegerMatrix n_cat_obs_gr2,
                                                           IntegerVector no_categories,
                                                           int no_persons_gr1,
                                                           int no_persons_gr2,
                                                           int variable,
                                                           NumericVector proposed_states,
                                                           NumericVector current_states,
                                                           NumericMatrix rest_matrix_gr1,
                                                           NumericMatrix rest_matrix_gr2) {

  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_proposed, denominator_current, exponent;
  int score;

  //Compute the pseudolikelihood ratio
  for(int category = 0; category < no_categories[variable]; category++) {
    pseudolikelihood_ratio -= .5 * n_cat_obs_gr1(category + 1, variable) *
      (proposed_states[category] - current_states[category]);
    pseudolikelihood_ratio += .5 * n_cat_obs_gr2(category + 1, variable) *
      (proposed_states[category] - current_states[category]);
  }

  //Compute the contribution for group 1
  for(int person = 0; person < no_persons_gr1; person++) {
    rest_score = rest_matrix_gr1(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score;
    } else {
      bound = 0.0;
    }

    denominator_proposed = std::exp(-bound);
    denominator_current = std::exp(-bound);
    for(int category = 0; category < no_categories[variable]; category++) {
      score = category + 1;
      exponent = score * rest_score - bound;
      denominator_proposed += std::exp(exponent +
        thresholds(variable, category) -
        .5 * proposed_states[category]);
        denominator_current += std::exp(exponent +
        thresholds(variable, category) -
        .5 * current_states[category]);
    }

    pseudolikelihood_ratio -= std::log(denominator_proposed);
    pseudolikelihood_ratio += std::log(denominator_current);
  }

  //Compute the contribution for group 2
  for(int person = 0; person < no_persons_gr2; person++) {
    rest_score = rest_matrix_gr2(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score;
    } else {
      bound = 0.0;
    }

    denominator_proposed = std::exp(-bound);
    denominator_current = std::exp(-bound);
    for(int category = 0; category < no_categories[variable]; category++) {
      score = category + 1;
      exponent = score * rest_score - bound;
      denominator_proposed += std::exp(exponent +
        thresholds(variable, category) +
        .5 * proposed_states[category]);
        denominator_current += std::exp(exponent +
        thresholds(variable, category) +
        .5 * current_states[category]);
    }

    pseudolikelihood_ratio -= std::log(denominator_proposed);
    pseudolikelihood_ratio += std::log(denominator_current);
  }

  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// Between model MH algorithm for the threshold differences
//    -- regular ordinal variable
// ----------------------------------------------------------------------------|
void compare_metropolis_main_difference_regular_between_model(NumericMatrix thresholds,
                                                              NumericMatrix main_difference,
                                                              IntegerMatrix n_cat_obs_gr1,
                                                              IntegerMatrix n_cat_obs_gr2,
                                                              IntegerVector no_categories,
                                                              IntegerMatrix indicator,
                                                              NumericMatrix proposal_sd_main_difference,
                                                              double main_difference_scale,
                                                              int no_persons_gr1,
                                                              int no_persons_gr2,
                                                              int variable,
                                                              NumericMatrix rest_matrix_gr1,
                                                              NumericMatrix rest_matrix_gr2,
                                                              NumericMatrix inclusion_probability_difference) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;
  int max_no_categories = max(no_categories);
  NumericVector proposed_states(max_no_categories);
  NumericVector current_states(max_no_categories);

  log_prob = 0.0;

  for(int category = 0; category < no_categories[variable]; category++) {
    current_state = main_difference(variable, category);
    current_states[category] = current_state;

    if(indicator(variable, variable) == 0) {
      proposed_state = R::rnorm(current_state, proposal_sd_main_difference(variable, category));
      proposed_states[category] = proposed_state;
    } else {
      proposed_state = 0.0;
      proposed_states[category] = proposed_state;
    }
    if(indicator(variable, variable) == 0) {
      log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_sd_main_difference(variable, category),
                           true);
      log_prob += std::log(inclusion_probability_difference(variable, variable));
      log_prob -= std::log(1 - inclusion_probability_difference(variable, variable));
    } else {
      log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_sd_main_difference(variable, category),
                           true);
      log_prob -= log(inclusion_probability_difference(variable, variable));
      log_prob += log(1 - inclusion_probability_difference(variable, variable));
    }
  }

  if(indicator(variable, variable) == 0) {
    log_prob += std::log(inclusion_probability_difference(variable, variable));
    log_prob -= std::log(1 - inclusion_probability_difference(variable, variable));
  } else {
    log_prob -= log(inclusion_probability_difference(variable, variable));
    log_prob += log(1 - inclusion_probability_difference(variable, variable));
  }

  log_prob += compare_log_pseudolikelihood_ratio_main_differences(thresholds,
                                                                  n_cat_obs_gr1,
                                                                  n_cat_obs_gr2,
                                                                  no_categories,
                                                                  no_persons_gr1,
                                                                  no_persons_gr2,
                                                                  variable,
                                                                  proposed_states,
                                                                  current_states,
                                                                  rest_matrix_gr1,
                                                                  rest_matrix_gr2);

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    indicator(variable, variable) = 1 - indicator(variable, variable);
    for(int category = 0; category < no_categories[variable]; category++) {
      main_difference(variable, category) = proposed_states[category];
    }
  }
}

// ----------------------------------------------------------------------------|
// Between model MH algorithm for the threshold differences
//    -- blume capel ordinal variable
// ----------------------------------------------------------------------------|
void compare_metropolis_main_difference_blumecapel_between_model(NumericMatrix thresholds,
                                                                 NumericMatrix main_difference,
                                                                 IntegerMatrix sufficient_blume_capel_gr1,
                                                                 IntegerMatrix sufficient_blume_capel_gr2,
                                                                 IntegerVector no_categories,
                                                                 IntegerMatrix indicator,
                                                                 NumericMatrix proposal_sd_main_difference,
                                                                 double main_difference_scale,
                                                                 int no_persons_gr1,
                                                                 int no_persons_gr2,
                                                                 int variable,
                                                                 NumericMatrix rest_matrix_gr1,
                                                                 NumericMatrix rest_matrix_gr2,
                                                                 NumericMatrix inclusion_probability_difference,
                                                                 IntegerVector reference_category) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;
  double lbound, bound, rest_score;
  double numerator, denominator, exponent;
  NumericVector proposed_states(2);
  NumericVector current_states(2);
  NumericVector constant_numerator (no_categories[variable] + 1);
  NumericVector constant_denominator (no_categories[variable] + 1);

  log_prob = 0.0;

  for(int parameter = 0; parameter < 2; parameter++) {
    //------------------------------------------------------------------------
    // Adaptive Metropolis for the difference in Blume-Capel parameters
    //------------------------------------------------------------------------
    current_state = main_difference(variable, parameter);
    current_states[parameter] = current_state;

    if(indicator(variable, variable) == 0) {
      proposed_state = R::rnorm(current_state,
                                proposal_sd_main_difference(variable, parameter));
      proposed_states[parameter] = proposed_state;
    } else {
      proposed_state = 0.0;
      proposed_states[parameter] = proposed_state;
    }

    //------------------------------------------------------------------------
    // Compute the log acceptance probability --------------------------------
    //------------------------------------------------------------------------

    //Compute the prior and proposal ratios ----------------------------------
    if(indicator(variable, variable) == 0) {
      log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_sd_main_difference(variable, parameter),
                           true);
    } else {
      log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_sd_main_difference(variable, parameter),
                           true);
    }
  }

  //Compute the prior inclusion ratios ---------------------------------------
  if(indicator(variable, variable) == 0) {
    log_prob += std::log(inclusion_probability_difference(variable, variable));
    log_prob -= std::log(1 - inclusion_probability_difference(variable, variable));
  } else {
    log_prob -= log(inclusion_probability_difference(variable, variable));
    log_prob += log(1 - inclusion_probability_difference(variable, variable));
  }

  //Precompute bounds for group 1 for numeric stability and efficiency -------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    constant_numerator[category] =
      (thresholds(variable, 0) - .5 * current_states[0]) *
      category +
      (thresholds(variable, 1) - .5 * current_states[1]) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_denominator[category] =
      (thresholds(variable, 0) - .5 * proposed_states[0]) *
      category +
      (thresholds(variable, 1) - .5 * proposed_states[1]) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
  }

  double tmp_n = max(constant_numerator);
  double tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log likelihood ratio group 1----------------------------------
  log_prob -= .5 * sufficient_blume_capel_gr1(0, variable) * proposed_states[0];
  log_prob += .5 * sufficient_blume_capel_gr1(0, variable) * current_states[0];
  log_prob -= .5 * sufficient_blume_capel_gr1(1, variable) * proposed_states[1];
  log_prob += .5 * sufficient_blume_capel_gr1(1, variable) * current_states[1];

  for(int person = 0; person < no_persons_gr1; person++) {
    rest_score = rest_matrix_gr1(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  //Precompute bounds for group 2 for numeric stability and efficiency -------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    constant_numerator[category] =
      (thresholds(variable, 0) + .5 * current_states[0]) *
      category +
      (thresholds(variable, 1) + .5 * current_states[1]) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_denominator[category] =
      (thresholds(variable, 0) + .5 * proposed_states[0]) *
      category +
      (thresholds(variable, 1) + .5 * proposed_states[1]) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
  }

  tmp_n = max(constant_numerator);
  tmp_d = max(constant_denominator);
  if(tmp_n > 0) {
    if(tmp_n > tmp_d) {
      lbound = tmp_n;
    } else {
      lbound = tmp_d;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log likelihood ratio group 2----------------------------------
  log_prob += .5 * sufficient_blume_capel_gr2(0, variable) * proposed_states[0];
  log_prob -= .5 * sufficient_blume_capel_gr2(0, variable) * current_states[0];
  log_prob += .5 * sufficient_blume_capel_gr2(1, variable) * proposed_states[1];
  log_prob -= .5 * sufficient_blume_capel_gr2(1, variable) * current_states[1];

  for(int person = 0; person < no_persons_gr2; person++) {
    rest_score = rest_matrix_gr2(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories[variable]; category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    indicator(variable, variable) = 1 - indicator(variable, variable);
    main_difference(variable, 0) = proposed_states[0];
    main_difference(variable, 1) = proposed_states[1];
  }
}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for a cross-lagged
//  interaction for the paired samples design
// ----------------------------------------------------------------------------|
double compare_log_pseudolikelihood_ratio_cross_lagged(NumericMatrix thresholds,
                                                       NumericMatrix main_difference,
                                                       IntegerMatrix observations_gr1,
                                                       IntegerMatrix observations_gr2,
                                                       IntegerVector no_categories,
                                                       int no_persons_gr1,
                                                       int no_persons_gr2,
                                                       int variable1,
                                                       int variable2,
                                                       double proposed_state,
                                                       double current_state,
                                                       NumericMatrix rest_matrix_gr1,
                                                       NumericMatrix rest_matrix_gr2,
                                                       LogicalVector variable_bool,
                                                       IntegerVector reference_category) {

  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  double delta_state = proposed_state - current_state;

  // The cross_lagged effect [i, j] occurs in four full-conditionals:
  //  X[i] | Y
  //  Y[i] | X
  //  X[j] | Y
  //  Y[j] | X
  // But if i = j, then it occurs in:
  //  X[i] | Y
  //  Y[i] | X
  // Also, all cross-lagged effects occur twice, so need to be rescaled.

  for(int person = 0; person < no_persons_gr1; person++) {
    //X[i] | Y
    obs_score1 = observations_gr1(person, variable1);
    obs_score2 = observations_gr2(person, variable2);

    pseudolikelihood_ratio += obs_score1 * obs_score2 * delta_state * .25;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix_gr1(person, variable1) -
      .25 * obs_score2 * current_state;

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
        exponent = thresholds(variable1, category) -
          .5 * main_difference(variable1, category) +
        score * rest_score -
        bound;
        denominator_prop +=
          std::exp(exponent + .25 * score * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + .25 * score * obs_score2 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable1] + 1; category++) {
        exponent = (thresholds(variable1, 0) -
          .5 * main_difference(variable1, 0)) *
        category;
        exponent += (thresholds(variable1, 1) -
          .5 * main_difference(variable1, 1)) *
        (category - reference_category[variable1]) *
        (category - reference_category[variable1]);
        exponent+= category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + .25 * category * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + .25 * category * obs_score2 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);

    if(variable1 != variable2) {
      //X[j] | Y
      obs_score1 = observations_gr1(person, variable2);
      obs_score2 = observations_gr2(person, variable1);

      pseudolikelihood_ratio += obs_score1 * obs_score2 * delta_state * .25;

      rest_score = rest_matrix_gr1(person, variable2) -
        .25 * obs_score2 * current_state;

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
          exponent = (thresholds(variable2, category) -
            .5 * main_difference(variable2, category)) +
          score * rest_score -
          bound;
          denominator_prop +=
            std::exp(exponent + .25 * score * obs_score2 * proposed_state);
          denominator_curr +=
            std::exp(exponent + .25 * score * obs_score2 * current_state);
        }
      } else {
        //Blume-Capel ordinal MRF variable ---------------------------------------
        denominator_prop = 0.0;
        denominator_curr = 0.0;
        for(int category = 0; category < no_categories[variable2] + 1; category++) {
          exponent = (thresholds(variable2, 0) -
            .5 * main_difference(variable2, 0)) *
          category;
          exponent += (thresholds(variable2, 1) -
            .5 * main_difference(variable2, 1)) *
          (category - reference_category[variable2]) *
          (category - reference_category[variable2]);
          exponent+=  category * rest_score - bound;
          denominator_prop +=
            std::exp(exponent + .25 * category * obs_score2 * proposed_state);
          denominator_curr +=
            std::exp(exponent + .25 * category * obs_score2 * current_state);
        }
      }
      pseudolikelihood_ratio -= std::log(denominator_prop);
      pseudolikelihood_ratio += std::log(denominator_curr);
    }
  }

  for(int person = 0; person < no_persons_gr2; person++) {
    //Y[i] | X
    obs_score1 = observations_gr2(person, variable1);
    obs_score2 = observations_gr1(person, variable2);

    pseudolikelihood_ratio += obs_score1 * obs_score2 * delta_state * .25;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix_gr2(person, variable1) -
      .25 * obs_score2 * current_state;

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
        exponent = thresholds(variable1, category) -
          .5 * main_difference(variable1, category) +
        score * rest_score -
        bound;
        denominator_prop +=
          std::exp(exponent + .25 * score * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + .25 * score * obs_score2 * current_state);
      }
    } else {
      //Blume-Capel ordinal MRF variable ---------------------------------------
      denominator_prop = 0.0;
      denominator_curr = 0.0;
      for(int category = 0; category < no_categories[variable1] + 1; category++) {
        exponent = (thresholds(variable1, 0) -
          .5 * main_difference(variable1, 0)) *
        category;
        exponent += (thresholds(variable1, 1) -
          .5 * main_difference(variable1, 1)) *
        (category - reference_category[variable1]) *
        (category - reference_category[variable1]);
        exponent+= category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + .25 * category * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + .25 * category * obs_score2 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);

    if(variable1 != variable2) {
      //Y[j] | X
      obs_score1 = observations_gr2(person, variable2);
      obs_score2 = observations_gr1(person, variable1);

      pseudolikelihood_ratio += obs_score1 * obs_score2 * delta_state * .25;

      rest_score = rest_matrix_gr2(person, variable2) -
        .25 * obs_score2 * current_state;

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
          exponent = (thresholds(variable2, category) -
            .5 * main_difference(variable2, category)) +
          score * rest_score -
          bound;
          denominator_prop +=
            std::exp(exponent + .25 * score * obs_score2 * proposed_state);
          denominator_curr +=
            std::exp(exponent + .25 * score * obs_score2 * current_state);
        }
      } else {
        //Blume-Capel ordinal MRF variable ---------------------------------------
        denominator_prop = 0.0;
        denominator_curr = 0.0;
        for(int category = 0; category < no_categories[variable2] + 1; category++) {
          exponent = (thresholds(variable2, 0) -
            .5 * main_difference(variable2, 0)) *
          category;
          exponent += (thresholds(variable2, 1) -
            .5 * main_difference(variable2, 1)) *
          (category - reference_category[variable2]) *
          (category - reference_category[variable2]);
          exponent+=  category * rest_score - bound;
          denominator_prop +=
            std::exp(exponent + .25 * category * obs_score2 * proposed_state);
          denominator_curr +=
            std::exp(exponent + .25 * category * obs_score2 * current_state);
        }
      }
      pseudolikelihood_ratio -= std::log(denominator_prop);
      pseudolikelihood_ratio += std::log(denominator_curr);
    }
  }
  return pseudolikelihood_ratio;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of a cross-lagged
//  interaction parameter (nuisance)
// ----------------------------------------------------------------------------|
void compare_metropolis_cross_lagged(NumericMatrix thresholds,
                                     NumericMatrix cross_lagged,
                                     NumericMatrix main_difference,
                                     IntegerMatrix observations_gr1,
                                     IntegerMatrix observations_gr2,
                                     IntegerVector no_categories,
                                     NumericMatrix proposal_sd_cross_lagged,
                                     double cross_lagged_scale,
                                     int no_persons_gr1,
                                     int no_persons_gr2,
                                     int no_variables,
                                     NumericMatrix rest_matrix_gr1,
                                     NumericMatrix rest_matrix_gr2,
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

  for(int variable1 = 0; variable1 <  no_variables; variable1++) {
    for(int variable2 = variable1; variable2 <  no_variables; variable2++) {
      current_state = cross_lagged(variable1, variable2);
      proposed_state = R::rnorm(current_state,
                                proposal_sd_cross_lagged(variable1, variable2));

      log_prob = compare_log_pseudolikelihood_ratio_cross_lagged(thresholds,
                                                                 main_difference,
                                                                 observations_gr1,
                                                                 observations_gr2,
                                                                 no_categories,
                                                                 no_persons_gr1,
                                                                 no_persons_gr2,
                                                                 variable1,
                                                                 variable2,
                                                                 proposed_state,
                                                                 current_state,
                                                                 rest_matrix_gr1,
                                                                 rest_matrix_gr2,
                                                                 variable_bool,
                                                                 reference_category);

      log_prob += R::dcauchy(proposed_state, 0.0, cross_lagged_scale, true);
      log_prob -= R::dcauchy(current_state, 0.0, cross_lagged_scale, true);

      U = R::unif_rand();
      if(std::log(U) < log_prob) {
        double state_difference = proposed_state - current_state;
        cross_lagged(variable1, variable2) = proposed_state;
        cross_lagged(variable2, variable1) = proposed_state;

        //Update the rest score matrices
        if(variable1 != variable2) {
          for(int person = 0; person < no_persons_gr1; person++) {
            rest_matrix_gr1(person, variable1) +=
              .25 *
              observations_gr2(person, variable2) *
              state_difference;

            rest_matrix_gr1(person, variable2) +=
              .25 *
              observations_gr2(person, variable1) *
              state_difference;

            rest_matrix_gr2(person, variable2) +=
              .25 *
              observations_gr1(person, variable1) *
              state_difference;

            rest_matrix_gr2(person, variable1) +=
              .25 *
              observations_gr1(person, variable2) *
              state_difference;
          }
        } else {
          for(int person = 0; person < no_persons_gr1; person++) {
            rest_matrix_gr1(person, variable1) +=
              .25 *
              observations_gr2(person, variable2) *
              state_difference;

            rest_matrix_gr2(person, variable2) +=
              .25 *
              observations_gr1(person, variable1) *
              state_difference;
          }
        }
      }

      if(log_prob > 0) {
        log_prob = 1.0;
      } else {
        log_prob = std::exp(log_prob);
      }

      double update_proposal_sd =
        proposal_sd_cross_lagged(variable1, variable2) +
        (log_prob - target_ar) * std::exp(-log(t) * phi);

      if(std::isnan(update_proposal_sd) == true) {
        update_proposal_sd = 1.0;
      }

      if(update_proposal_sd < epsilon_lo) {
        update_proposal_sd = epsilon_lo;
      } else if (update_proposal_sd > epsilon_hi) {
        update_proposal_sd = epsilon_hi;
      }

      proposal_sd_cross_lagged(variable1, variable2) = update_proposal_sd;
    }
  }
}

// ----------------------------------------------------------------------------|
// A Gibbs step for graphical model parameters for Bayesian parameter comparison
// ----------------------------------------------------------------------------|
List compare_gibbs_step_gm(NumericMatrix interactions,
                           NumericMatrix thresholds,
                           NumericMatrix cross_lagged,
                           NumericMatrix pairwise_difference,
                           NumericMatrix main_difference,
                           IntegerMatrix observations_gr1,
                           IntegerMatrix observations_gr2,
                           IntegerVector no_categories,
                           NumericMatrix proposal_sd_interaction,
                           NumericMatrix proposal_sd_cross_lagged,
                           NumericMatrix proposal_sd_pairwise_difference,
                           NumericMatrix proposal_sd_blumecapel,
                           NumericMatrix proposal_sd_main_difference,
                           double interaction_scale,
                           double cross_lagged_scale,
                           int no_persons_gr1,
                           int no_persons_gr2,
                           int no_variables,
                           NumericMatrix rest_matrix_gr1,
                           NumericMatrix rest_matrix_gr2,
                           double phi,
                           double target_ar,
                           int t,
                           double epsilon_lo,
                           double epsilon_hi,
                           LogicalVector variable_bool,
                           IntegerVector reference_category,
                           IntegerMatrix indicator,
                           double pairwise_difference_scale,
                           double main_difference_scale,
                           NumericMatrix inclusion_probability_difference,
                           int no_interactions,
                           IntegerMatrix index,
                           IntegerMatrix n_cat_obs_gr1,
                           IntegerMatrix n_cat_obs_gr2,
                           IntegerMatrix sufficient_blume_capel_gr1,
                           IntegerMatrix sufficient_blume_capel_gr2,
                           double threshold_alpha,
                           double threshold_beta,
                           bool paired,
                           bool difference_selection) {

  //Between model move for differences in interaction parameters
  if(difference_selection == true) {
    compare_metropolis_pairwise_difference_between_model(thresholds,
                                                         pairwise_difference,
                                                         main_difference,
                                                         observations_gr1,
                                                         observations_gr2,
                                                         no_categories,
                                                         indicator,
                                                         proposal_sd_pairwise_difference,
                                                         pairwise_difference_scale,
                                                         inclusion_probability_difference,
                                                         no_persons_gr1,
                                                         no_persons_gr2,
                                                         no_interactions,
                                                         index,
                                                         rest_matrix_gr1,
                                                         rest_matrix_gr2,
                                                         phi,
                                                         target_ar,
                                                         t,
                                                         epsilon_lo,
                                                         epsilon_hi,
                                                         variable_bool,
                                                         reference_category);
  }

  // //Within model move for differences in interaction parameters
  compare_metropolis_pairwise_difference(thresholds,
                                         pairwise_difference,
                                         main_difference,
                                         observations_gr1,
                                         observations_gr2,
                                         no_categories,
                                         indicator,
                                         proposal_sd_pairwise_difference,
                                         pairwise_difference_scale,
                                         no_persons_gr1,
                                         no_persons_gr2,
                                         no_variables,
                                         rest_matrix_gr1,
                                         rest_matrix_gr2,
                                         phi,
                                         target_ar,
                                         t,
                                         epsilon_lo,
                                         epsilon_hi,
                                         variable_bool,
                                         reference_category);

  //Within model move for the main interaction parameters
  compare_metropolis_interaction(interactions,
                                 thresholds,
                                 main_difference,
                                 observations_gr1,
                                 observations_gr2,
                                 no_categories,
                                 proposal_sd_interaction,
                                 interaction_scale,
                                 no_persons_gr1,
                                 no_persons_gr2,
                                 no_variables,
                                 rest_matrix_gr1,
                                 rest_matrix_gr2,
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
      if(difference_selection == true) {
        //Between model move for the differences in category thresholds
        compare_metropolis_main_difference_regular_between_model(thresholds,
                                                                 main_difference,
                                                                 n_cat_obs_gr1,
                                                                 n_cat_obs_gr2,
                                                                 no_categories,
                                                                 indicator,
                                                                 proposal_sd_main_difference,
                                                                 main_difference_scale,
                                                                 no_persons_gr1,
                                                                 no_persons_gr2,
                                                                 variable,
                                                                 rest_matrix_gr1,
                                                                 rest_matrix_gr2,
                                                                 inclusion_probability_difference);
      }

      //Within model move for the differences in category thresholds
      compare_metropolis_main_difference_regular(thresholds,
                                                 main_difference,
                                                 n_cat_obs_gr1,
                                                 n_cat_obs_gr2,
                                                 no_categories,
                                                 indicator,
                                                 proposal_sd_main_difference,
                                                 main_difference_scale,
                                                 no_persons_gr1,
                                                 no_persons_gr2,
                                                 variable,
                                                 rest_matrix_gr1,
                                                 rest_matrix_gr2,
                                                 phi,
                                                 target_ar,
                                                 t,
                                                 epsilon_lo,
                                                 epsilon_hi);

      //Within model move for the main category thresholds
      compare_metropolis_threshold_regular(thresholds,
                                           main_difference,
                                           no_categories,
                                           n_cat_obs_gr1,
                                           n_cat_obs_gr2,
                                           no_persons_gr1,
                                           no_persons_gr2,
                                           variable,
                                           threshold_alpha,
                                           threshold_beta,
                                           rest_matrix_gr1,
                                           rest_matrix_gr2);

    } else {
      //Between model move for the differences in category thresholds
      if(difference_selection == true) {
        compare_metropolis_main_difference_blumecapel_between_model(thresholds,
                                                                    main_difference,
                                                                    sufficient_blume_capel_gr1,
                                                                    sufficient_blume_capel_gr1,
                                                                    no_categories,
                                                                    indicator,
                                                                    proposal_sd_main_difference,
                                                                    main_difference_scale,
                                                                    no_persons_gr1,
                                                                    no_persons_gr2,
                                                                    variable,
                                                                    rest_matrix_gr1,
                                                                    rest_matrix_gr2,
                                                                    inclusion_probability_difference,
                                                                    reference_category);
      }

      //Within model move for the differences in category thresholds
      compare_metropolis_main_difference_blumecapel(thresholds,
                                                    main_difference,
                                                    no_categories,
                                                    sufficient_blume_capel_gr1,
                                                    sufficient_blume_capel_gr2,
                                                    no_persons_gr1,
                                                    no_persons_gr2,
                                                    variable,
                                                    reference_category,
                                                    main_difference_scale,
                                                    rest_matrix_gr1,
                                                    rest_matrix_gr2,
                                                    proposal_sd_main_difference,
                                                    phi,
                                                    target_ar,
                                                    t,
                                                    epsilon_lo,
                                                    epsilon_hi);

      //Within model move for the main category thresholds
      compare_metropolis_threshold_blumecapel(thresholds,
                                              main_difference,
                                              no_categories,
                                              sufficient_blume_capel_gr1,
                                              sufficient_blume_capel_gr2,
                                              no_persons_gr1,
                                              no_persons_gr2,
                                              variable,
                                              reference_category,
                                              threshold_alpha,
                                              threshold_beta,
                                              rest_matrix_gr1,
                                              rest_matrix_gr2,
                                              proposal_sd_blumecapel,
                                              phi,
                                              target_ar,
                                              t,
                                              epsilon_lo,
                                              epsilon_hi);

    }
  }

  if(paired == true) {
    compare_metropolis_cross_lagged(thresholds,
                                    cross_lagged,
                                    main_difference,
                                    observations_gr1,
                                    observations_gr2,
                                    no_categories,
                                    proposal_sd_cross_lagged,
                                    cross_lagged_scale,
                                    no_persons_gr1,
                                    no_persons_gr2,
                                    no_variables,
                                    rest_matrix_gr1,
                                    rest_matrix_gr2,
                                    phi,
                                    target_ar,
                                    t,
                                    epsilon_lo,
                                    epsilon_hi,
                                    variable_bool,
                                    reference_category);
  }

  return List::create(Named("indicator") = indicator,
                      Named("interactions") = interactions,
                      Named("cross_lagged") = cross_lagged,
                      Named("thresholds") = thresholds,
                      Named("pairwise_difference") = pairwise_difference,
                      Named("main_difference") = main_difference,
                      Named("rest_matrix_gr1") = rest_matrix_gr1,
                      Named("rest_matrix_gr2") = rest_matrix_gr2,
                      Named("proposal_sd_interaction") = proposal_sd_interaction,
                      Named("proposal_sd_interaction_difference") = proposal_sd_pairwise_difference,
                      Named("proposal_sd_main_difference") = proposal_sd_main_difference);
}

// ----------------------------------------------------------------------------|
// The Gibbs sampler for Bayesian parameter comparisons
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List compare_gibbs_sampler(IntegerMatrix observations_gr1,
                           IntegerMatrix observations_gr2,
                           IntegerVector no_categories,
                           double interaction_scale,
                           double cross_lagged_scale,
                           double pairwise_difference_scale,
                           double main_difference_scale,
                           String pairwise_difference_prior,
                           String main_difference_prior,
                           NumericMatrix inclusion_probability_difference,
                           double pairwise_beta_bernoulli_alpha,
                           double pairwise_beta_bernoulli_beta,
                           double main_beta_bernoulli_alpha,
                           double main_beta_bernoulli_beta,
                           IntegerMatrix Index,
                           int iter,
                           int burnin,
                           IntegerMatrix n_cat_obs_gr1,
                           IntegerMatrix n_cat_obs_gr2,
                           IntegerMatrix sufficient_blume_capel_gr1,
                           IntegerMatrix sufficient_blume_capel_gr2,
                           double threshold_alpha,
                           double threshold_beta,
                           bool na_impute,
                           IntegerMatrix missing_index_gr1,
                           IntegerMatrix missing_index_gr2,
                           LogicalVector variable_bool,
                           IntegerVector reference_category,
                           bool paired,
                           bool save = false,
                           bool display_progress = false,
                           bool difference_selection = true) {
  int cntr;
  int no_variables = observations_gr1.ncol();
  int no_persons_gr1 = observations_gr1.nrow();
  int no_persons_gr2 = observations_gr2.nrow();
  int no_interactions = Index.nrow();
  int no_thresholds = sum(no_categories);
  int max_no_categories = max(no_categories);

  // Matrices with model parameters --------------------------------------------
  NumericMatrix interactions(no_variables, no_variables);
  NumericMatrix cross_lagged(no_variables, no_variables);
  NumericMatrix pairwise_difference(no_variables, no_variables);
  NumericMatrix thresholds(no_variables, max_no_categories);
  NumericMatrix main_difference(no_variables, max_no_categories);
  IntegerMatrix indicator(no_variables, no_variables);

  // Matrices with standard deviations for adaptive Metropolis -----------------
  NumericMatrix proposal_sd_interaction(no_variables, no_variables);
  NumericMatrix proposal_sd_cross_lagged(no_variables, no_variables);
  NumericMatrix proposal_sd_pairwise_difference(no_variables, no_variables);
  NumericMatrix proposal_sd_main_difference(no_variables, max_no_categories);
  NumericMatrix proposal_sd_blumecapel(no_variables, max_no_categories);

  //Initialize them at 1.0
  for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
      indicator(variable1, variable2) = 1;
      indicator(variable2, variable1) = 1;
      proposal_sd_interaction(variable1, variable2) = 1.0;
      proposal_sd_interaction(variable2, variable1) = 1.0;
      proposal_sd_pairwise_difference(variable1, variable2) = 1.0;
      proposal_sd_pairwise_difference(variable2, variable1) = 1.0;
    }
    indicator(variable1, variable1) = 1;
    for(int category = 0; category < no_categories[variable1]; category++) {
      proposal_sd_main_difference(variable1, category) = 1.0;
    }
  }
  indicator(no_variables - 1, no_variables - 1) = 1;
  for(int category = 0; category < no_categories[no_variables-1]; category++) {
    proposal_sd_main_difference(no_variables-1, category) = 1.0;
  }
  if(paired == true) {
    for(int variable1 = 0; variable1 < no_variables; variable1++) {
      for(int variable2 = variable1; variable2 < no_variables; variable2++) {
        proposal_sd_cross_lagged(variable1, variable2) = 1.0;
        proposal_sd_cross_lagged(variable2, variable1) = 1.0;
      }
      proposal_sd_blumecapel(variable1, 0) = 1.0;
      proposal_sd_blumecapel(variable1, 1) = 1.0;
    }
  }

  //Parameters for the Robbins-Monro approach for adaptive Metropolis ----------
  double phi = .75;
  double target_ar = 0.234;
  double epsilon_lo;
  if(no_persons_gr1 > no_persons_gr2) {
    epsilon_lo = 1 / no_persons_gr2;
  } else {
    epsilon_lo = 1 / no_persons_gr1;
  }
  double epsilon_hi = 2.0;

  //Randomized index for the pairwise updates ----------------------------------
  IntegerVector v = seq(0, no_interactions - 1);
  IntegerVector order(no_interactions);
  IntegerMatrix index(no_interactions, 3);

  //Output matrices (resizing based on ``save'') -------------------------------
  int nrow = no_variables;
  int ncol_edges = no_variables;
  int ncol_thresholds = max_no_categories;
  int ncol_cross_lagged = no_variables;

  if(save == true) {
    nrow = iter;
    ncol_edges = no_interactions;
    ncol_thresholds = no_thresholds;
    if(paired == true) {
      ncol_cross_lagged = no_variables + no_interactions;
    }
  }

  NumericMatrix out_interactions(nrow, ncol_edges);
  NumericMatrix out_cross_lagged(nrow, ncol_cross_lagged);
  NumericMatrix out_thresholds(nrow, ncol_thresholds);
  NumericMatrix out_pairwise_difference(nrow, ncol_edges);
  NumericMatrix out_main_difference(nrow, ncol_thresholds);

  NumericMatrix out_indicator_pairwise_difference(nrow, ncol_edges);
  if(save == false) {
    nrow = 1;
  }
  NumericMatrix out_indicator_main_difference(nrow, no_variables);

  //These matrices will contain the rest scores in the pseudolikelihoods -------
  NumericMatrix rest_matrix_gr1(no_persons_gr1, no_variables);
  NumericMatrix rest_matrix_gr2(no_persons_gr2, no_variables);

  //Progress bar ---------------------------------------------------------------
  Progress p(iter + burnin, display_progress);

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------
  for(int iteration = 0; iteration < burnin; iteration++) {
    if (Progress::check_abort()) {
      if(paired == true) {
        return List::create(Named("pairwise_difference_indicator") = out_indicator_pairwise_difference,
                            Named("main_difference_indicator") = out_indicator_main_difference,
                            Named("interactions") = out_interactions,
                            Named("cross_lagged") = out_cross_lagged,
                            Named("thresholds") = out_thresholds,
                            Named("pairwise_difference") = out_pairwise_difference,
                            Named("main_difference") = out_main_difference);
      } else {
        return List::create(Named("pairwise_difference_indicator") = out_indicator_pairwise_difference,
                            Named("main_difference_indicator") = out_indicator_main_difference,
                            Named("interactions") = out_interactions,
                            Named("thresholds") = out_thresholds,
                            Named("pairwise_difference") = out_pairwise_difference,
                            Named("main_difference") = out_main_difference);
      }
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
      List out = compare_impute_missing_data(interactions,
                                             thresholds,
                                             cross_lagged,
                                             pairwise_difference,
                                             main_difference,
                                             observations_gr1,
                                             observations_gr2,
                                             n_cat_obs_gr1,
                                             n_cat_obs_gr2,
                                             sufficient_blume_capel_gr1,
                                             sufficient_blume_capel_gr2,
                                             no_categories,
                                             rest_matrix_gr1,
                                             rest_matrix_gr2,
                                             missing_index_gr1,
                                             missing_index_gr2,
                                             variable_bool,
                                             reference_category,
                                             paired);

      IntegerMatrix observations_gr1 = out["observations_gr1"];
      IntegerMatrix observations_gr2 = out["observations_gr2"];
      IntegerMatrix n_cat_obs_gr1 = out["n_cat_obs_gr1"];
      IntegerMatrix n_cat_obs_gr2 = out["n_cat_obs_gr2"];
      IntegerMatrix sufficient_blume_capel_gr1 = out["sufficient_blume_capel_gr1"];
      IntegerMatrix sufficient_blume_capel_gr2 = out["sufficient_blume_capel_gr2"];
      NumericMatrix rest_matrix_gr1 = out["rest_matrix_gr1"];
      NumericMatrix rest_matrix_gr2 = out["rest_matrix_gr2"];
    }

    List out = compare_gibbs_step_gm(interactions,
                                     thresholds,
                                     cross_lagged,
                                     pairwise_difference,
                                     main_difference,
                                     observations_gr1,
                                     observations_gr2,
                                     no_categories,
                                     proposal_sd_interaction,
                                     proposal_sd_cross_lagged,
                                     proposal_sd_pairwise_difference,
                                     proposal_sd_blumecapel,
                                     proposal_sd_main_difference,
                                     interaction_scale,
                                     cross_lagged_scale,
                                     no_persons_gr1,
                                     no_persons_gr2,
                                     no_variables,
                                     rest_matrix_gr1,
                                     rest_matrix_gr2,
                                     phi,
                                     target_ar,
                                     iteration,
                                     epsilon_lo,
                                     epsilon_hi,
                                     variable_bool,
                                     reference_category,
                                     indicator,
                                     pairwise_difference_scale,
                                     main_difference_scale,
                                     inclusion_probability_difference,
                                     no_interactions,
                                     index,
                                     n_cat_obs_gr1,
                                     n_cat_obs_gr2,
                                     sufficient_blume_capel_gr1,
                                     sufficient_blume_capel_gr2,
                                     threshold_alpha,
                                     threshold_beta,
                                     paired,
                                     difference_selection);

    IntegerMatrix indicator = out["indicator"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix cross_lagged = out["cross_lagged"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix pairwise_difference = out["pairwise_difference"];
    NumericMatrix main_difference = out["main_difference"];
    NumericMatrix rest_matrix_gr1 = out["rest_matrix_gr1"];
    NumericMatrix rest_matrix_gr2 = out["rest_matrix_gr2"];
    NumericMatrix proposal_sd_interaction = out["proposal_sd_interaction"];
    NumericMatrix proposal_sd_pairwise_difference = out["proposal_sd_interaction_difference"];
    NumericMatrix proposal_sd_main_difference = out["proposal_sd_main_difference"];

    if(difference_selection == true) {
      if(pairwise_difference_prior == "Beta-Bernoulli") {
        int sumG = 0;
        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            sumG += indicator(i, j);
          }
        }
        double probability = R::rbeta(pairwise_beta_bernoulli_alpha + sumG,
                                      pairwise_beta_bernoulli_beta + no_interactions - sumG);

        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            inclusion_probability_difference(i, j) = probability;
            inclusion_probability_difference(j, i) = probability;
          }
        }
      }

      if(main_difference_prior == "Beta-Bernoulli") {
        int sumG = 0;
        for(int i = 0; i < no_variables; i++) {
          sumG += indicator(i, i);
        }
        double probability = R::rbeta(main_beta_bernoulli_alpha + sumG,
                                      main_beta_bernoulli_beta + no_variables - sumG);

        for(int i = 0; i < no_variables; i++) {
          inclusion_probability_difference(i, i) = probability;
        }
      }
    }
  }

  //The post burn-in iterations ------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    if (Progress::check_abort()) {
      if(paired == true) {
        return List::create(Named("pairwise_difference_indicator") = out_indicator_pairwise_difference,
                            Named("main_difference_indicator") = out_indicator_main_difference,
                            Named("interactions") = out_interactions,
                            Named("cross_lagged") = out_cross_lagged,
                            Named("thresholds") = out_thresholds,
                            Named("pairwise_difference") = out_pairwise_difference,
                            Named("main_difference") = out_main_difference);
      } else {
        return List::create(Named("pairwise_difference_indicator") = out_indicator_pairwise_difference,
                            Named("main_difference_indicator") = out_indicator_main_difference,
                            Named("interactions") = out_interactions,
                            Named("thresholds") = out_thresholds,
                            Named("pairwise_difference") = out_pairwise_difference,
                            Named("main_difference") = out_main_difference);
      }
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
      List out = compare_impute_missing_data(interactions,
                                             thresholds,
                                             cross_lagged,
                                             pairwise_difference,
                                             main_difference,
                                             observations_gr1,
                                             observations_gr2,
                                             n_cat_obs_gr1,
                                             n_cat_obs_gr2,
                                             sufficient_blume_capel_gr1,
                                             sufficient_blume_capel_gr2,
                                             no_categories,
                                             rest_matrix_gr1,
                                             rest_matrix_gr2,
                                             missing_index_gr1,
                                             missing_index_gr2,
                                             variable_bool,
                                             reference_category,
                                             paired);

      IntegerMatrix observations_gr1 = out["observations_gr1"];
      IntegerMatrix observations_gr2 = out["observations_gr2"];
      IntegerMatrix n_cat_obs_gr1 = out["n_cat_obs_gr1"];
      IntegerMatrix n_cat_obs_gr2 = out["n_cat_obs_gr2"];
      IntegerMatrix sufficient_blume_capel_gr1 = out["sufficient_blume_capel_gr1"];
      IntegerMatrix sufficient_blume_capel_gr2 = out["sufficient_blume_capel_gr2"];
      NumericMatrix rest_matrix_gr1 = out["rest_matrix_gr1"];
      NumericMatrix rest_matrix_gr2 = out["rest_matrix_gr2"];
    }

    List out = compare_gibbs_step_gm(interactions,
                                     thresholds,
                                     cross_lagged,
                                     pairwise_difference,
                                     main_difference,
                                     observations_gr1,
                                     observations_gr2,
                                     no_categories,
                                     proposal_sd_interaction,
                                     proposal_sd_cross_lagged,
                                     proposal_sd_pairwise_difference,
                                     proposal_sd_blumecapel,
                                     proposal_sd_main_difference,
                                     interaction_scale,
                                     cross_lagged_scale,
                                     no_persons_gr1,
                                     no_persons_gr2,
                                     no_variables,
                                     rest_matrix_gr1,
                                     rest_matrix_gr2,
                                     phi,
                                     target_ar,
                                     iteration,
                                     epsilon_lo,
                                     epsilon_hi,
                                     variable_bool,
                                     reference_category,
                                     indicator,
                                     pairwise_difference_scale,
                                     main_difference_scale,
                                     inclusion_probability_difference,
                                     no_interactions,
                                     index,
                                     n_cat_obs_gr1,
                                     n_cat_obs_gr2,
                                     sufficient_blume_capel_gr1,
                                     sufficient_blume_capel_gr2,
                                     threshold_alpha,
                                     threshold_beta,
                                     paired,
                                     difference_selection);

    IntegerMatrix indicator = out["indicator"];
    NumericMatrix interactions = out["interactions"];
    NumericMatrix cross_lagged = out["cross_lagged"];
    NumericMatrix thresholds = out["thresholds"];
    NumericMatrix pairwise_difference = out["pairwise_difference"];
    NumericMatrix main_difference = out["main_difference"];
    NumericMatrix rest_matrix_gr1 = out["rest_matrix_gr1"];
    NumericMatrix rest_matrix_gr2 = out["rest_matrix_gr2"];
    NumericMatrix proposal_sd_interaction = out["proposal_sd_interaction"];
    NumericMatrix proposal_sd_pairwise_difference = out["proposal_sd_interaction_difference"];
    NumericMatrix proposal_sd_main_difference = out["proposal_sd_main_difference"];

    if(difference_selection == true) {
      if(pairwise_difference_prior == "Beta-Bernoulli") {
        int sumG = 0;
        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            sumG += indicator(i, j);
          }
        }
        double probability = R::rbeta(pairwise_beta_bernoulli_alpha + sumG,
                                      pairwise_beta_bernoulli_beta + no_interactions - sumG);

        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            inclusion_probability_difference(i, j) = probability;
            inclusion_probability_difference(j, i) = probability;
          }
        }
      }

      if(main_difference_prior == "Beta-Bernoulli") {
        int sumG = 0;
        for(int i = 0; i < no_variables; i++) {
          sumG += indicator(i, i);
        }
        double probability = R::rbeta(main_beta_bernoulli_alpha + sumG,
                                      main_beta_bernoulli_beta + no_variables - sumG);

        for(int i = 0; i < no_variables; i++) {
          inclusion_probability_difference(i, i) = probability;
        }
      }
    }

    //Output -------------------------------------------------------------------
    if(save == true) {
      //Save raw samples -------------------------------------------------------
      cntr = 0;
      for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
        for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
          out_indicator_pairwise_difference(iteration, cntr) = indicator(variable1, variable2);
          out_interactions(iteration, cntr) = interactions(variable1, variable2);
          out_pairwise_difference(iteration, cntr) = pairwise_difference(variable1, variable2);
          cntr++;
        }
      }

      if(paired == true) {
        cntr = 0;
        for(int variable1 = 0; variable1 < no_variables; variable1++) {
          for(int variable2 = variable1; variable2 < no_variables; variable2++) {
            out_cross_lagged(iteration, cntr) = cross_lagged(variable1, variable2);
            cntr++;
          }
        }
      }

      cntr = 0;
      for(int variable = 0; variable < no_variables; variable++) {
        if(variable_bool[variable] == true) {
          out_indicator_main_difference(iteration, variable) = indicator(variable, variable);
          for(int category = 0; category < no_categories[variable]; category++) {
            out_thresholds(iteration, cntr) = thresholds(variable, category);
            out_main_difference(iteration, cntr) = main_difference(variable, category);
            cntr++;
          }
        } else {
          out_thresholds(iteration, cntr) = thresholds(variable, 0);
          out_main_difference(iteration, cntr) = main_difference(variable, 0);
          cntr++;
          out_thresholds(iteration, cntr) = thresholds(variable, 1);
          out_main_difference(iteration, cntr) = main_difference(variable, 1);
          cntr++;
        }
      }
    } else {
      //Compute running averages -----------------------------------------------
      for(int variable1 = 0; variable1 < no_variables; variable1++) {
        for(int variable2 = variable1; variable2 < no_variables; variable2++) {
          out_cross_lagged(variable1, variable2) *= iteration;
          out_cross_lagged(variable1, variable2) += cross_lagged(variable1, variable2);
          out_cross_lagged(variable1, variable2) /= iteration + 1;
          out_cross_lagged(variable2, variable1) = out_cross_lagged(variable1, variable2);
        }
      }

      for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
        for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
          out_indicator_pairwise_difference(variable1, variable2) *= iteration;
          out_indicator_pairwise_difference(variable1, variable2) += indicator(variable1, variable2);
          out_indicator_pairwise_difference(variable1, variable2) /= iteration + 1;
          out_indicator_pairwise_difference(variable2, variable1) = out_indicator_pairwise_difference(variable1, variable2);

          out_interactions(variable1, variable2) *= iteration;
          out_interactions(variable1, variable2) += interactions(variable1, variable2);
          out_interactions(variable1, variable2) /= iteration + 1;
          out_interactions(variable2, variable1) = out_interactions(variable1, variable2);

          out_pairwise_difference(variable1, variable2) *= iteration;
          out_pairwise_difference(variable1, variable2) += pairwise_difference(variable1, variable2);
          out_pairwise_difference(variable1, variable2) /= iteration + 1;
          out_pairwise_difference(variable2, variable1) = out_pairwise_difference(variable1, variable2);
        }

        if(variable_bool[variable1] == true) {
          for(int category = 0; category < no_categories[variable1]; category++) {
            out_thresholds(variable1, category) *= iteration;
            out_thresholds(variable1, category) += thresholds(variable1, category);
            out_thresholds(variable1, category) /= iteration + 1;

            out_main_difference(variable1, category) *= iteration;
            out_main_difference(variable1, category) += main_difference(variable1, category);
            out_main_difference(variable1, category) /= iteration + 1;

            out_indicator_main_difference(0, variable1) *= iteration;
            out_indicator_main_difference(0, variable1) += indicator(variable1, variable1);
            out_indicator_main_difference(0, variable1) /= iteration + 1;
          }
        } else {
          out_thresholds(variable1, 0) *= iteration;
          out_thresholds(variable1, 0) += thresholds(variable1, 0);
          out_thresholds(variable1, 0) /= iteration + 1;
          out_thresholds(variable1, 1) *= iteration;
          out_thresholds(variable1, 1) += thresholds(variable1, 1);
          out_thresholds(variable1, 1) /= iteration + 1;

          out_main_difference(variable1, 0) *= iteration;
          out_main_difference(variable1, 0) += main_difference(variable1, 0);
          out_main_difference(variable1, 0) /= iteration + 1;
          out_main_difference(variable1, 1) *= iteration;
          out_main_difference(variable1, 1) += main_difference(variable1, 1);
          out_main_difference(variable1, 1) /= iteration + 1;

          out_indicator_main_difference(0, variable1) *= iteration;
          out_indicator_main_difference(0, variable1) += indicator(variable1, variable1);
          out_indicator_main_difference(0, variable1) /= iteration + 1;
        }
      }
      if(variable_bool[no_variables - 1] == true) {
        for(int category = 0; category < no_categories[no_variables - 1]; category++) {
          out_thresholds(no_variables - 1, category) *= iteration;
          out_thresholds(no_variables - 1, category) += thresholds(no_variables - 1, category);
          out_thresholds(no_variables - 1, category) /= iteration + 1;

          out_main_difference(no_variables - 1, category) *= iteration;
          out_main_difference(no_variables - 1, category) += main_difference(no_variables - 1, category);
          out_main_difference(no_variables - 1, category) /= iteration + 1;
        }
      } else {
        out_thresholds(no_variables - 1, 0) *= iteration;
        out_thresholds(no_variables - 1, 0) += thresholds(no_variables - 1, 0);
        out_thresholds(no_variables - 1, 0) /= iteration + 1;
        out_thresholds(no_variables - 1, 1) *= iteration;
        out_thresholds(no_variables - 1, 1) += thresholds(no_variables - 1, 1);
        out_thresholds(no_variables - 1, 1) /= iteration + 1;

        out_main_difference(no_variables - 1, 0) *= iteration;
        out_main_difference(no_variables - 1, 0) += main_difference(no_variables - 1, 0);
        out_main_difference(no_variables - 1, 0) /= iteration + 1;
        out_main_difference(no_variables - 1, 1) *= iteration;
        out_main_difference(no_variables - 1, 1) += main_difference(no_variables - 1, 1);
        out_main_difference(no_variables - 1, 1) /= iteration + 1;
      }
      out_indicator_main_difference(0, no_variables - 1) *= iteration;
      out_indicator_main_difference(0, no_variables - 1) += indicator(no_variables - 1, no_variables - 1);
      out_indicator_main_difference(0, no_variables - 1) /= iteration + 1;
    }
  }

  if(paired == true) {
    return List::create(Named("pairwise_difference_indicator") = out_indicator_pairwise_difference,
                        Named("main_difference_indicator") = out_indicator_main_difference,
                        Named("interactions") = out_interactions,
                        Named("cross_lagged") = out_cross_lagged,
                        Named("thresholds") = out_thresholds,
                        Named("pairwise_difference") = out_pairwise_difference,
                        Named("main_difference") = out_main_difference);
  } else {
    return List::create(Named("pairwise_difference_indicator") = out_indicator_pairwise_difference,
                        Named("main_difference_indicator") = out_indicator_main_difference,
                        Named("interactions") = out_interactions,
                        Named("thresholds") = out_thresholds,
                        Named("pairwise_difference") = out_pairwise_difference,
                        Named("main_difference") = out_main_difference);
  }
}