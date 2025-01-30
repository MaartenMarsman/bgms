// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "gibbs_functions.h"
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// Impute missing data for two-samples designs
// ----------------------------------------------------------------------------|
List compare_ttest_impute_missing_data(NumericMatrix thresholds_gr1,
                                       NumericMatrix thresholds_gr2,
                                       NumericMatrix interactions_gr1,
                                       NumericMatrix interactions_gr2,
                                       IntegerMatrix observations_gr1,
                                       IntegerMatrix observations_gr2,
                                       IntegerMatrix n_cat_obs_gr1,
                                       IntegerMatrix n_cat_obs_gr2,
                                       IntegerMatrix sufficient_blume_capel_gr1,
                                       IntegerMatrix sufficient_blume_capel_gr2,
                                       IntegerVector no_categories_gr1,
                                       IntegerVector no_categories_gr2,
                                       NumericMatrix rest_matrix_gr1,
                                       NumericMatrix rest_matrix_gr2,
                                       IntegerMatrix missing_index_gr1,
                                       IntegerMatrix missing_index_gr2,
                                       LogicalVector ordinal_variable,
                                       IntegerVector reference_category) {
  int no_variables = observations_gr1.ncol();
  int no_missings_gr1 = missing_index_gr1.nrow();
  int no_missings_gr2 = missing_index_gr2.nrow();
  int max_no_categories = max(no_categories_gr1);
  if(max(no_categories_gr2) > max(no_categories_gr1)) {
    max_no_categories = max(no_categories_gr2);
  }
  NumericVector probabilities(max_no_categories + 1);
  double exponent, rest_score, cumsum, u;
  int score, person, variable, new_observation, old_observation;

  //Impute missing data (if there are any) for group 1 -------------------------
  if(no_missings_gr1 > 1) {
    for(int missing = 0; missing < no_missings_gr1; missing++) {
      //Which observation to impute? -------------------------------------------
      person = missing_index_gr1(missing, 0);
      variable = missing_index_gr1(missing, 1);

      //Generate a new observation from the ordinal MRF ------------------------
      rest_score = rest_matrix_gr1(person, variable);
      if(ordinal_variable[variable] == true) {
        //Regular binary or ordinal variable -----------------------------------
        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < no_categories_gr1[variable]; category++) {
          exponent = thresholds_gr1(variable, category);
          exponent += (category + 1) * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }
      } else {
        //Blume-Capel variable -------------------------------------------------
        cumsum = 0.0;
        for(int category = 0; category < no_categories_gr1[variable] + 1; category++) {
          exponent = thresholds_gr1(variable, 0) * category;
          exponent += thresholds_gr1(variable, 1) *
            (category - reference_category[variable]) *
            (category - reference_category[variable]);
          exponent += category * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category] = cumsum;
        }
      }
      u = cumsum * R::unif_rand();
      score = 0;
      while (u > probabilities[score]) {
        score++;
      }

      //Update current data with newly generated observation -------------------
      new_observation = score;
      old_observation = observations_gr1(person, variable);
      if(old_observation != new_observation) {
        //Update raw data ------------------------------------------------------
        observations_gr1(person, variable) = new_observation;

        //Update pre-computed statistics ---------------------------------------
        if(ordinal_variable[variable] == true) {
          //Regular binary or ordinal variable ---------------------------------
          n_cat_obs_gr1(old_observation, variable)--;
          n_cat_obs_gr1(new_observation, variable)++;
        } else {
          //Blume-Capel variable -----------------------------------------------
          sufficient_blume_capel_gr1(0, variable) -= old_observation;
          sufficient_blume_capel_gr1(0, variable) += new_observation;
          sufficient_blume_capel_gr1(1, variable) -=
            (old_observation - reference_category[variable]) *
            (old_observation - reference_category[variable]);
          sufficient_blume_capel_gr1(1, variable) +=
            (new_observation - reference_category[variable]) *
            (new_observation - reference_category[variable]);
        }

        //Update rest score ----------------------------------------------------
        for(int vertex = 0; vertex < no_variables; vertex++) {
          rest_matrix_gr1(person, vertex) -= old_observation *
            interactions_gr1(vertex, variable);
          rest_matrix_gr1(person, vertex) += new_observation *
            interactions_gr1(vertex, variable);
        }
      }
    }
  }

  //Impute missing data (if there are any) for group 1 -------------------------
  if(no_missings_gr2 > 1) {
    for(int missing = 0; missing < no_missings_gr2; missing++) {
      //Which observation to impute? -------------------------------------------
      person = missing_index_gr2(missing, 0) - 1; //R to C++ indexing
      variable = missing_index_gr2(missing, 1) - 1; //R to C++ indexing

      //Generate a new observation from the ordinal MRF ------------------------
      rest_score = rest_matrix_gr2(person, variable);
      if(ordinal_variable[variable] == true) {
        //Regular binary or ordinal variable -----------------------------------
        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < no_categories_gr2[variable]; category++) {
          exponent = thresholds_gr2(variable, category);
          exponent += (category + 1) * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }
      } else {
        //Blume-Capel variable -------------------------------------------------
        cumsum = 0.0;
        for(int category = 0; category <= no_categories_gr2[variable]; category++) {
          exponent = thresholds_gr2(variable, 0) * category;
          exponent += thresholds_gr2(variable, 1) *
            (category - reference_category[variable]) *
            (category - reference_category[variable]);
          exponent += category * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category] = cumsum;
        }
      }
      u = cumsum * R::unif_rand();
      score = 0;
      while (u > probabilities[score]) {
        score++;
      }

      //Update current data with newly generated observation -------------------
      new_observation = score;
      old_observation = observations_gr2(person, variable);
      if(old_observation != new_observation) {
        //Update raw data ------------------------------------------------------
        observations_gr2(person, variable) = new_observation;

        //Update pre-computed statistics ---------------------------------------
        if(ordinal_variable[variable] == true) {
          //Regular binary or ordinal variable ---------------------------------
          n_cat_obs_gr2(old_observation, variable)--;
          n_cat_obs_gr2(new_observation, variable)++;
        } else {
          //Blume-Capel variable -----------------------------------------------
          sufficient_blume_capel_gr2(0, variable) -= old_observation;
          sufficient_blume_capel_gr2(0, variable) += new_observation;
          sufficient_blume_capel_gr2(1, variable) -=
            (old_observation - reference_category[variable]) *
            (old_observation - reference_category[variable]);
          sufficient_blume_capel_gr2(1, variable) +=
            (new_observation - reference_category[variable]) *
            (new_observation - reference_category[variable]);
        }

        //Update rest score ----------------------------------------------------
        for(int vertex = 0; vertex < no_variables; vertex++) {
          rest_matrix_gr2(person, vertex) -= old_observation *
            interactions_gr2(vertex, variable);
          rest_matrix_gr2(person, vertex) += new_observation *
            interactions_gr2(vertex, variable);
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
double compare_ttest_log_pseudolikelihood_ratio_interaction(NumericMatrix thresholds_gr1,
                                                            NumericMatrix thresholds_gr2,
                                                            IntegerMatrix observations_gr1,
                                                            IntegerMatrix observations_gr2,
                                                            IntegerVector no_categories_gr1,
                                                            IntegerVector no_categories_gr2,
                                                            int no_persons_gr1,
                                                            int no_persons_gr2,
                                                            int variable1,
                                                            int variable2,
                                                            double proposed_state,
                                                            double current_state,
                                                            NumericMatrix rest_matrix_gr1,
                                                            NumericMatrix rest_matrix_gr2,
                                                            LogicalVector ordinal_variable,
                                                            IntegerVector reference_category) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  double delta_state = proposed_state - current_state;

  // --------------------------------------------------------------------------|
  // Compute log pseudolikelihood ratio for group 1 ----------------------------
  // --------------------------------------------------------------------------|
  for(int person = 0; person < no_persons_gr1; person++) {
    obs_score1 = observations_gr1(person, variable1);
    obs_score2 = observations_gr1(person, variable2);

    pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * delta_state;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix_gr1(person, variable1) -
      obs_score2 * current_state;

    if(rest_score > 0) {
      bound = no_categories_gr1[variable1] * rest_score;
    } else {
      bound = 0.0;
    }

    if(ordinal_variable[variable1] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories_gr1[variable1]; category++) {
        score = category + 1;
        exponent = thresholds_gr1(variable1, category) +
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
      for(int category = 0; category < no_categories_gr1[variable1] + 1; category++) {
        exponent = thresholds_gr1(variable1, 0) * category;
        exponent += thresholds_gr1(variable1, 1) *
          (category - reference_category[variable1]) *
          (category - reference_category[variable1]);
        exponent += category * rest_score - bound;
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
      bound = no_categories_gr1[variable2] * rest_score;
    } else {
      bound = 0.0;
    }

    if(ordinal_variable[variable2] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories_gr1[variable2]; category++) {
        score = category + 1;
        exponent = thresholds_gr1(variable2, category) +
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
      for(int category = 0; category < no_categories_gr1[variable2] + 1; category++) {
        exponent = thresholds_gr1(variable2, 0) * category;
        exponent += thresholds_gr1(variable2, 1) *
          (category - reference_category[variable2]) *
          (category - reference_category[variable2]);
        exponent +=  category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * obs_score1 * proposed_state);
        denominator_curr +=
          std::exp(exponent + category * obs_score1 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);
  }

  // --------------------------------------------------------------------------|
  // Compute log pseudolikelihood ratio for group 2 ----------------------------
  // --------------------------------------------------------------------------|
  for(int person = 0; person < no_persons_gr2; person++) {
    obs_score1 = observations_gr2(person, variable1);
    obs_score2 = observations_gr2(person, variable2);

    pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * delta_state;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix_gr2(person, variable1) - obs_score2 * current_state;

    if(rest_score > 0) {
      bound = no_categories_gr2[variable1] * rest_score;
    } else {
      bound = 0.0;
    }

    if(ordinal_variable[variable1] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories_gr2[variable1]; category++) {
        score = category + 1;
        exponent = thresholds_gr2(variable1, category) +
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
      for(int category = 0; category < no_categories_gr2[variable1] + 1; category++) {
        exponent = thresholds_gr2(variable1, 0) * category;
        exponent += thresholds_gr2(variable1, 1) *
          (category - reference_category[variable1]) *
          (category - reference_category[variable1]);
        exponent += category * rest_score - bound;
        denominator_prop +=
          std::exp(exponent + category * obs_score2 * proposed_state);
        denominator_curr +=
          std::exp(exponent + category * obs_score2 * current_state);
      }
    }
    pseudolikelihood_ratio -= std::log(denominator_prop);
    pseudolikelihood_ratio += std::log(denominator_curr);

    //variable 2 log pseudolikelihood ratio
    rest_score = rest_matrix_gr2(person, variable2) - obs_score1 * current_state;

    if(rest_score > 0) {
      bound = no_categories_gr2[variable2] * rest_score;
    } else {
      bound = 0.0;
    }

    if(ordinal_variable[variable2] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories_gr2[variable2]; category++) {
        score = category + 1;
        exponent = thresholds_gr2(variable2, category) +
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
      for(int category = 0; category < no_categories_gr2[variable2] + 1; category++) {
        exponent = thresholds_gr2(variable2, 0) * category;
        exponent += thresholds_gr2(variable2, 1) *
          (category - reference_category[variable2]) *
          (category - reference_category[variable2]);
        exponent +=  category * rest_score - bound;
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
void compare_ttest_metropolis_interaction(NumericMatrix interactions,
                                          NumericMatrix thresholds_gr1,
                                          NumericMatrix thresholds_gr2,
                                          IntegerMatrix observations_gr1,
                                          IntegerMatrix observations_gr2,
                                          IntegerVector no_categories_gr1,
                                          IntegerVector no_categories_gr2,
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
                                          LogicalVector ordinal_variable,
                                          IntegerVector reference_category) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  for(int variable1 = 0; variable1 <  no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  no_variables; variable2++) {
      current_state = interactions(variable1, variable2);
      proposed_state = R::rnorm(current_state, proposal_sd_interaction(variable1, variable2));

      log_prob = compare_ttest_log_pseudolikelihood_ratio_interaction(thresholds_gr1,
                                                                      thresholds_gr2,
                                                                      observations_gr1,
                                                                      observations_gr2,
                                                                      no_categories_gr1,
                                                                      no_categories_gr2,
                                                                      no_persons_gr1,
                                                                      no_persons_gr2,
                                                                      variable1,
                                                                      variable2,
                                                                      proposed_state,
                                                                      current_state,
                                                                      rest_matrix_gr1,
                                                                      rest_matrix_gr2,
                                                                      ordinal_variable,
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

      update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

      proposal_sd_interaction(variable1, variable2) = update_proposal_sd;
      proposal_sd_interaction(variable2, variable1) = update_proposal_sd;
    }
  }
}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for the difference
//  in a pairwise interaction for a two independent samples design
// ----------------------------------------------------------------------------|
double compare_ttest_log_pseudolikelihood_ratio_pairwise_difference(NumericMatrix thresholds_gr1,
                                                                    NumericMatrix thresholds_gr2,
                                                                    IntegerMatrix observations_gr1,
                                                                    IntegerMatrix observations_gr2,
                                                                    IntegerVector no_categories_gr1,
                                                                    IntegerVector no_categories_gr2,
                                                                    int no_persons_gr1,
                                                                    int no_persons_gr2,
                                                                    int variable1,
                                                                    int variable2,
                                                                    double proposed_state,
                                                                    double current_state,
                                                                    NumericMatrix rest_matrix_gr1,
                                                                    NumericMatrix rest_matrix_gr2,
                                                                    LogicalVector ordinal_variable,
                                                                    IntegerVector reference_category) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  double delta_state = proposed_state - current_state;

  // --------------------------------------------------------------------------|
  // Compute log pseudolikelihood ratio for group 1 ----------------------------
  // --------------------------------------------------------------------------|
  for(int person = 0; person < no_persons_gr1; person++) {
    obs_score1 = observations_gr1(person, variable1);
    obs_score2 = observations_gr1(person, variable2);

    pseudolikelihood_ratio -= obs_score1 * obs_score2 * delta_state;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix_gr1(person, variable1) +
      .5 * obs_score2 * current_state;

    if(rest_score > 0) {
      bound = no_categories_gr1[variable1] * rest_score;
    } else {
      bound = 0.0;
    }

    if(ordinal_variable[variable1] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories_gr1[variable1]; category++) {
        score = category + 1;
        exponent = thresholds_gr1(variable1, category) +
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
      for(int category = 0; category < no_categories_gr1[variable1] + 1; category++) {
        exponent = thresholds_gr1(variable1, 0) * category;
        exponent += thresholds_gr1(variable1, 1) *
          (category - reference_category[variable1]) *
          (category - reference_category[variable1]);
        exponent += category * rest_score - bound;
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
      bound = no_categories_gr1[variable2] * rest_score;
    } else {
      bound = 0.0;
    }

    if(ordinal_variable[variable2] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories_gr1[variable2]; category++) {
        score = category + 1;
        exponent = thresholds_gr1(variable2, category) +
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
      for(int category = 0; category < no_categories_gr1[variable2] + 1; category++) {
        exponent = thresholds_gr1(variable2, 0) * category;
        exponent += thresholds_gr1(variable2, 1) *
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

  // --------------------------------------------------------------------------|
  // Compute log pseudolikelihood ratio for group 2 ----------------------------
  // --------------------------------------------------------------------------|
  for(int person = 0; person < no_persons_gr2; person++) {
    obs_score1 = observations_gr2(person, variable1);
    obs_score2 = observations_gr2(person, variable2);

    pseudolikelihood_ratio += obs_score1 * obs_score2 * delta_state;

    //variable 1 log pseudolikelihood ratio
    rest_score = rest_matrix_gr2(person, variable1) -
      .5 * obs_score2 * current_state;

    if(rest_score > 0) {
      bound = no_categories_gr2[variable1] * rest_score;
    } else {
      bound = 0.0;
    }

    if(ordinal_variable[variable1] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories_gr2[variable1]; category++) {
        score = category + 1;
        exponent = thresholds_gr2(variable1, category) +
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
      for(int category = 0; category < no_categories_gr2[variable1] + 1; category++) {
        exponent = thresholds_gr2(variable1, 0) * category;
        exponent += thresholds_gr2(variable1, 1) *
          (category - reference_category[variable1]) *
          (category - reference_category[variable1]);
        exponent += category * rest_score - bound;
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
      bound = no_categories_gr2[variable2] * rest_score;
    } else {
      bound = 0.0;
    }

    if(ordinal_variable[variable2] == true) {
      //Regular binary or ordinal MRF variable ---------------------------------
      denominator_prop = std::exp(-bound);
      denominator_curr = std::exp(-bound);
      for(int category = 0; category < no_categories_gr2[variable2]; category++) {
        score = category + 1;
        exponent = thresholds_gr2(variable2, category) +
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
      for(int category = 0; category < no_categories_gr2[variable2] + 1; category++) {
        exponent = thresholds_gr2(variable2, 0) * category;
        exponent += thresholds_gr2(variable2, 1) *
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
void compare_ttest_metropolis_pairwise_difference(NumericMatrix pairwise_difference,
                                                  NumericMatrix thresholds_gr1,
                                                  NumericMatrix thresholds_gr2,
                                                  IntegerMatrix observations_gr1,
                                                  IntegerMatrix observations_gr2,
                                                  IntegerVector no_categories_gr1,
                                                  IntegerVector no_categories_gr2,
                                                  IntegerMatrix indicator,
                                                  int no_persons_gr1,
                                                  int no_persons_gr2,
                                                  int no_variables,
                                                  NumericMatrix rest_matrix_gr1,
                                                  NumericMatrix rest_matrix_gr2,
                                                  NumericMatrix proposal_sd_pairwise_difference,
                                                  double pairwise_difference_scale,
                                                  double phi,
                                                  double target_ar,
                                                  int t,
                                                  double epsilon_lo,
                                                  double epsilon_hi,
                                                  LogicalVector ordinal_variable,
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

        log_prob = compare_ttest_log_pseudolikelihood_ratio_pairwise_difference(thresholds_gr1,
                                                                                thresholds_gr2,
                                                                                observations_gr1,
                                                                                observations_gr2,
                                                                                no_categories_gr1,
                                                                                no_categories_gr2,
                                                                                no_persons_gr1,
                                                                                no_persons_gr2,
                                                                                variable1,
                                                                                variable2,
                                                                                proposed_state,
                                                                                current_state,
                                                                                rest_matrix_gr1,
                                                                                rest_matrix_gr2,
                                                                                ordinal_variable,
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

        update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

        proposal_sd_pairwise_difference(variable1, variable2) = update_proposal_sd;
        proposal_sd_pairwise_difference(variable2, variable1) = update_proposal_sd;
      }
    }
  }
}

// ----------------------------------------------------------------------------|
// Between model MH algorithm for the pairwise differences
// ----------------------------------------------------------------------------|
void compare_ttest_metropolis_pairwise_difference_between_model(IntegerVector indicator,
                                                                NumericMatrix inclusion_probability_difference,
                                                                NumericMatrix pairwise_difference,
                                                                NumericMatrix thresholds_gr1,
                                                                NumericMatrix thresholds_gr2,
                                                                IntegerMatrix observations_gr1,
                                                                IntegerMatrix observations_gr2,
                                                                IntegerVector no_categories_gr1,
                                                                IntegerVector no_categories_gr2,
                                                                int no_persons_gr1,
                                                                int no_persons_gr2,
                                                                int no_interactions,
                                                                IntegerMatrix index,
                                                                NumericMatrix rest_matrix_gr1,
                                                                NumericMatrix rest_matrix_gr2,
                                                                NumericMatrix proposal_sd_pairwise_difference,
                                                                double pairwise_difference_scale,
                                                                double phi,
                                                                double target_ar,
                                                                int t,
                                                                double epsilon_lo,
                                                                double epsilon_hi,
                                                                LogicalVector ordinal_variable,
                                                                IntegerVector reference_category) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  int variable1;
  int variable2;

  for(int cntr = 0; cntr < no_interactions; cntr ++) {
    variable1 = index(cntr, 1);
    variable2 = index(cntr, 2);

    current_state = pairwise_difference(variable1, variable2);

    if(indicator(variable1, variable2) == 0) {
      proposed_state = R::rnorm(current_state, proposal_sd_pairwise_difference(variable1, variable2));
    } else {
      proposed_state = 0.0;
    }

    log_prob = compare_ttest_log_pseudolikelihood_ratio_pairwise_difference(thresholds_gr1,
                                                                            thresholds_gr2,
                                                                            observations_gr1,
                                                                            observations_gr2,
                                                                            no_categories_gr1,
                                                                            no_categories_gr2,
                                                                            no_persons_gr1,
                                                                            no_persons_gr2,
                                                                            variable1,
                                                                            variable2,
                                                                            proposed_state,
                                                                            current_state,
                                                                            rest_matrix_gr1,
                                                                            rest_matrix_gr2,
                                                                            ordinal_variable,
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
      indicator(variable2, variable1) = indicator(variable1, variable2);

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
void compare_ttest_metropolis_threshold_regular(NumericMatrix thresholds,
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
  NumericVector g1(no_persons_gr1);
  NumericVector q1(no_persons_gr1);
  NumericVector g2(no_persons_gr2);
  NumericVector q2(no_persons_gr2);

  double log_prob, rest_score;
  double a, b, c;
  double tmp;
  double current_state, proposed_state;
  double U;
  double exp_current, exp_proposed;
  double exponent;

  for(int category = 0; category < no_categories[variable]; category++) {
    current_state = thresholds(variable, category);
    exp_current = std::exp(current_state);
    c = (threshold_alpha + threshold_beta) / (1 + exp_current);

    for(int person = 0; person < no_persons_gr1; person++) {
      g1[person] = 1.0;
      q1[person] = 1.0;
      rest_score = rest_matrix_gr1(person, variable);
      for(int cat = 0; cat < no_categories[variable]; cat++) {
        if(cat != category) {
          exponent = thresholds(variable, cat);
          exponent -= .5 * main_difference(variable, cat);
          exponent += (cat + 1) * rest_score;
          g1[person] += std::exp(exponent);
        }
      }
      exponent = (category + 1) * rest_score;
      exponent -= .5 * main_difference(variable, category);
      q1[person] = std::exp(exponent);
      c +=  q1[person] / (g1[person] + q1[person] * exp_current);
    }
    for(int person = 0; person < no_persons_gr2; person++) {
      g2[person] = 1.0;
      q2[person] = 1.0;
      rest_score = rest_matrix_gr2(person, variable);
      for(int cat = 0; cat < no_categories[variable]; cat++) {
        if(cat != category) {
          exponent = thresholds(variable, cat);
          exponent += .5 * main_difference(variable, cat);
          exponent += (cat + 1) * rest_score;
          g2[person] += std::exp(exponent);
        }
      }
      exponent = (category + 1) * rest_score;
      exponent += .5 * main_difference(variable, category);
      q2[person] = std::exp(exponent);
      c +=  q2[person] / (g2[person] + q2[person] * exp_current);
    }
    tmp = no_persons_gr1;
    tmp += no_persons_gr2;
    tmp += threshold_alpha;
    tmp += threshold_beta;
    tmp -= exp_current * c;
    c = c / tmp;

    //Proposal is generalized beta-prime.
    a = n_cat_obs_gr1(category + 1, variable);
    a += n_cat_obs_gr2(category + 1, variable);
    a += threshold_alpha;
    b = no_persons_gr1;
    b += no_persons_gr2;
    b += threshold_beta;
    b -= n_cat_obs_gr1(category + 1, variable);
    b -= n_cat_obs_gr2(category + 1, variable);
    tmp = R::rbeta(a, b);
    proposed_state = std::log(tmp / (1  - tmp) / c);
    exp_proposed = exp(proposed_state);

    //Compute log_acceptance probability for Metropolis.
    //First, we use g and q above to compute the ratio of pseudolikelihoods
    log_prob = 0;
    for(int person = 0; person < no_persons_gr1; person++) {
      log_prob += std::log(g1[person] + q1[person] * exp_current);
      log_prob -= std::log(g1[person] + q1[person] * exp_proposed);
    }
    for(int person = 0; person < no_persons_gr2; person++) {
      log_prob += std::log(g2[person] + q2[person] * exp_current);
      log_prob -= std::log(g2[person] + q2[person] * exp_proposed);
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
double compare_ttest_log_pseudolikelihood_ratio_main_difference(NumericMatrix thresholds,
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
    current_thresholds[cat] = thresholds(variable, cat);
    current_thresholds[cat] -= .5 * main_difference(variable, cat);
    proposed_thresholds[cat] = current_thresholds[cat];
  }
  proposed_thresholds[category] = thresholds(variable, category);
  proposed_thresholds[category] -= .5 * proposed_state;

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
    current_thresholds[cat] = thresholds(variable, cat);
    current_thresholds[cat] += .5 * main_difference(variable, cat);
    proposed_thresholds[cat] = current_thresholds[cat];
  }
  proposed_thresholds[category] = thresholds(variable, category);
  proposed_thresholds[category] += .5 * proposed_state;

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
void compare_ttest_metropolis_main_difference_regular(NumericMatrix thresholds,
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
      if(n_cat_obs_gr1(category + 1, variable) * n_cat_obs_gr2(category + 1, variable) > 0) {

        current_state = main_difference(variable, category);
        proposed_state = R::rnorm(current_state,
                                  proposal_sd_main_difference(variable, category));

        log_prob = compare_ttest_log_pseudolikelihood_ratio_main_difference(thresholds,
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

        update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

        proposal_sd_main_difference(variable, category) = update_proposal_sd;
      } else {
        main_difference(variable, category) = 0.0;
      }
    }
  }
}

// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for a vector of
// category threshold differences for the two independent samples design
// ----------------------------------------------------------------------------|
double compare_ttest_log_pseudolikelihood_ratio_main_differences(NumericMatrix thresholds,
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
void compare_ttest_metropolis_main_difference_regular_between_model(NumericMatrix thresholds,
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
    if(n_cat_obs_gr1(category + 1, variable) * n_cat_obs_gr2(category + 1, variable) > 0) {

      current_state = main_difference(variable, category);
      current_states[category] = current_state;

      if(indicator(variable, variable) == 0) {
        proposed_state = R::rnorm(current_state, proposal_sd_main_difference(variable, category));
        proposed_states[category] = proposed_state;
        log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
        log_prob -= R::dnorm(proposed_state,
                             current_state,
                             proposal_sd_main_difference(variable, category),
                             true);
      } else {
        proposed_state = 0.0;
        proposed_states[category] = proposed_state;
        log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);
        log_prob += R::dnorm(current_state,
                             proposed_state,
                             proposal_sd_main_difference(variable, category),
                             true);
      }
    } else {
      current_states[category] = 0.0;
      proposed_states[category] = 0.0;
    }
  }

  if(indicator(variable, variable) == 0) {
    log_prob += std::log(inclusion_probability_difference(variable, variable));
    log_prob -= std::log(1 - inclusion_probability_difference(variable, variable));
  } else {
    log_prob -= log(inclusion_probability_difference(variable, variable));
    log_prob += log(1 - inclusion_probability_difference(variable, variable));
  }

  log_prob += compare_ttest_log_pseudolikelihood_ratio_main_differences(thresholds,
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
// The log pseudolikelihood ratio [proposed against current] for the two
// category threshold parameters of the Blume-Capel model
// ----------------------------------------------------------------------------|
double compare_ttest_log_pseudolikelihood_ratio_thresholds_blumecapel(double linear_current,
                                                                      double quadratic_current,
                                                                      double linear_proposed,
                                                                      double quadratic_proposed,
                                                                      int variable,
                                                                      IntegerVector reference_category,
                                                                      NumericMatrix main_difference,
                                                                      IntegerMatrix sufficient_blume_capel_gr1,
                                                                      IntegerMatrix sufficient_blume_capel_gr2,
                                                                      int no_persons_gr1,
                                                                      int no_persons_gr2,
                                                                      NumericMatrix rest_matrix_gr1,
                                                                      NumericMatrix rest_matrix_gr2,
                                                                      IntegerVector no_categories) {
  NumericVector constant_numerator (no_categories[variable] + 1);
  NumericVector constant_denominator (no_categories[variable] + 1);
  double lbound, bound;
  double log_prob, rest_score, numerator, denominator, exponent;
  int linear_score, quadratic_score;

  //----------------------------------------------------------------------------
  //Compute the log acceptance probability -------------------------------------
  //----------------------------------------------------------------------------

  //Precompute common terms for group 1 for computational efficiency -----------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    linear_score = category;
    quadratic_score =
      (category - reference_category[variable]) *
      (category - reference_category[variable]);

    constant_numerator[category] = linear_current * linear_score;
    constant_numerator[category] -= .5 * main_difference(variable, 0) * linear_score;
    constant_numerator[category] += quadratic_current * quadratic_score;
    constant_numerator[category] -= .5 * main_difference(variable, 1) * quadratic_score;

    constant_denominator[category] = linear_proposed * linear_score ;
    constant_denominator[category] -= .5 * main_difference(variable, 0) * linear_score;
    constant_denominator[category] += quadratic_proposed * quadratic_score;
    constant_denominator[category] -= .5 * main_difference(variable, 1) * quadratic_score;
  }

  //Precompute bounds for group 1 for numeric stability ------------------------
  double tmp_num = max(constant_numerator);
  double tmp_den = max(constant_denominator);
  if(tmp_num > 0) {
    if(tmp_num > tmp_den) {
      lbound = tmp_num;
    } else {
      lbound = tmp_den;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log pseudolikelihood ratio for group 1--------------------------
  log_prob = sufficient_blume_capel_gr1(0, variable) * linear_proposed;
  log_prob += sufficient_blume_capel_gr1(1, variable) * quadratic_proposed;
  log_prob -= sufficient_blume_capel_gr1(0, variable) * linear_current;
  log_prob -= sufficient_blume_capel_gr1(1, variable) * quadratic_current;

  for(int person = 0; person < no_persons_gr1; person++) {
    rest_score = rest_matrix_gr1(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = 0.0;
    denominator = 0.0;
    for(int category = 0; category < no_categories[variable] + 1; category ++) {
      exponent = category * rest_score - bound;
      numerator += std::exp(constant_numerator[category] + exponent);
      denominator += std::exp(constant_denominator[category] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  //Precompute common terms for group 2 for computational efficiency -----------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    linear_score = category;
    quadratic_score =
      (category - reference_category[variable]) *
      (category - reference_category[variable]);

    constant_numerator[category] = linear_current * linear_score ;
    constant_numerator[category] += .5 * main_difference(variable, 0) * linear_score;
    constant_numerator[category] += quadratic_current * quadratic_score;
    constant_numerator[category] += .5 * main_difference(variable, 1) * quadratic_score;

    constant_denominator[category] = linear_proposed * linear_score;
    constant_denominator[category] += .5 * main_difference(variable, 0) * linear_score;
    constant_denominator[category] += quadratic_proposed * quadratic_score;
    constant_denominator[category] += .5 * main_difference(variable, 1) * quadratic_score;
  }

  //Precompute bounds for group 2 for numeric stability ------------------------
  tmp_num = max(constant_numerator);
  tmp_den = max(constant_denominator);
  if(tmp_num > 0) {
    if(tmp_num > tmp_den) {
      lbound = tmp_num;
    } else {
      lbound = tmp_den;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log pseudolikelihood ratio for group 2--------------------------
  log_prob += sufficient_blume_capel_gr2(0, variable) * linear_proposed;
  log_prob += sufficient_blume_capel_gr2(1, variable) * quadratic_proposed;
  log_prob -= sufficient_blume_capel_gr2(0, variable) * linear_current;
  log_prob -= sufficient_blume_capel_gr2(1, variable) * quadratic_current;

  for(int person = 0; person < no_persons_gr2; person++) {
    rest_score = rest_matrix_gr2(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = 0.0;
    denominator = 0.0;
    for(int category = 0; category < no_categories[variable] + 1; category ++) {
      exponent = category * rest_score - bound;
      numerator += std::exp(constant_numerator[category] + exponent);
      denominator += std::exp(constant_denominator[category] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  return log_prob;
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the overall category
//  threshold parameters (nuisance) -- Blume-Capel ordinal variable
// ----------------------------------------------------------------------------|
void compare_ttest_metropolis_threshold_blumecapel(NumericMatrix thresholds,
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
  NumericVector constant_numerator (no_categories[variable] + 1);
  NumericVector constant_denominator (no_categories[variable] + 1);

  //----------------------------------------------------------------------------
  // Adaptive Metropolis for the linear Blume-Capel parameter
  //----------------------------------------------------------------------------
  current_state = thresholds(variable, 0);
  proposed_state = R::rnorm(current_state, proposal_sd_blumecapel(variable, 0));

  //----------------------------------------------------------------------------
  //Compute the log acceptance probability -------------------------------------
  //----------------------------------------------------------------------------
  double linear_current = current_state;
  double quadratic_current = thresholds(variable, 1);
  double linear_proposed = proposed_state;
  double quadratic_proposed = thresholds(variable, 1);

  log_prob = compare_ttest_log_pseudolikelihood_ratio_thresholds_blumecapel(linear_current,
                                                                            quadratic_current,
                                                                            linear_proposed,
                                                                            quadratic_proposed,
                                                                            variable,
                                                                            reference_category,
                                                                            main_difference,
                                                                            sufficient_blume_capel_gr1,
                                                                            sufficient_blume_capel_gr2,
                                                                            no_persons_gr1,
                                                                            no_persons_gr2,
                                                                            rest_matrix_gr1,
                                                                            rest_matrix_gr2,
                                                                            no_categories);

  //Compute the prior ratio ----------------------------------------------------
  log_prob += threshold_alpha * (proposed_state - current_state);
  log_prob += (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(current_state));
  log_prob -= (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(proposed_state));

  //Metropolis step ------------------------------------------------------------
  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    thresholds(variable, 0) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ------------------------------
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

  update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

  proposal_sd_blumecapel(variable, 0) = update_proposal_sd;

  //---------------------------------------------------------------------------|
  // Adaptive Metropolis for the quadratic Blume-Capel parameter
  //---------------------------------------------------------------------------|
  current_state = thresholds(variable, 1);
  proposed_state = R::rnorm(current_state, proposal_sd_blumecapel(variable, 1));

  //----------------------------------------------------------------------------
  //Compute the log acceptance probability -------------------------------------
  //----------------------------------------------------------------------------

  linear_current = thresholds(variable, 0);
  quadratic_current = current_state;
  linear_proposed = thresholds(variable, 0);
  quadratic_proposed =  proposed_state;

  log_prob = compare_ttest_log_pseudolikelihood_ratio_thresholds_blumecapel(linear_current,
                                                                            quadratic_current,
                                                                            linear_proposed,
                                                                            quadratic_proposed,
                                                                            variable,
                                                                            reference_category,
                                                                            main_difference,
                                                                            sufficient_blume_capel_gr1,
                                                                            sufficient_blume_capel_gr2,
                                                                            no_persons_gr1,
                                                                            no_persons_gr2,
                                                                            rest_matrix_gr1,
                                                                            rest_matrix_gr2,
                                                                            no_categories);

  //Compute the prior ratio ----------------------------------------------------
  log_prob += threshold_alpha * (proposed_state - current_state);
  log_prob += (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(current_state));
  log_prob -= (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(proposed_state));

  //Metropolis step ------------------------------------------------------------
  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    thresholds(variable, 1) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ------------------------------
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

  update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

  proposal_sd_blumecapel(variable, 1) = update_proposal_sd;
}


// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for the two
// category threshold differences of the Blume-Capel model
// ----------------------------------------------------------------------------|
double compare_ttest_log_pseudolikelihood_ratio_main_difference_blumecapel(double linear_current,
                                                                           double quadratic_current,
                                                                           double linear_proposed,
                                                                           double quadratic_proposed,
                                                                           int variable,
                                                                           IntegerVector reference_category,
                                                                           NumericMatrix thresholds,
                                                                           IntegerMatrix sufficient_blume_capel_gr1,
                                                                           IntegerMatrix sufficient_blume_capel_gr2,
                                                                           int no_persons_gr1,
                                                                           int no_persons_gr2,
                                                                           NumericMatrix rest_matrix_gr1,
                                                                           NumericMatrix rest_matrix_gr2,
                                                                           IntegerVector no_categories) {
  NumericVector constant_numerator (no_categories[variable] + 1);
  NumericVector constant_denominator (no_categories[variable] + 1);
  double lbound, bound;
  double log_prob, rest_score, numerator, denominator, exponent;
  int linear_score, quadratic_score;

  //----------------------------------------------------------------------------
  //Compute the log acceptance probability -------------------------------------
  //----------------------------------------------------------------------------

  double linear_sufficient = .5 * (sufficient_blume_capel_gr2(0, variable) -
                                sufficient_blume_capel_gr1(0, variable));
  double quadratic_sufficient = .5 * (sufficient_blume_capel_gr2(1, variable) -
                                   sufficient_blume_capel_gr1(1, variable));

  log_prob = linear_proposed * linear_sufficient;
  log_prob -= linear_current * linear_sufficient;
  log_prob += quadratic_proposed * quadratic_sufficient;
  log_prob -= quadratic_current * quadratic_sufficient;

  //Precompute common terms for group 1 for computational efficiency -----------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    linear_score = category;
    quadratic_score =
      (category - reference_category[variable]) *
      (category - reference_category[variable]);

    constant_numerator[category] = thresholds(variable, 0) * linear_score ;
    constant_numerator[category] -= .5 * linear_current * linear_score;
    constant_numerator[category] += thresholds(variable, 1) * quadratic_score;
    constant_numerator[category] -= .5 * quadratic_current * quadratic_score;

    constant_denominator[category] = thresholds(variable, 0) * linear_score ;
    constant_denominator[category] -= .5 * linear_proposed * linear_score;
    constant_denominator[category] += thresholds(variable, 1) * quadratic_score;
    constant_denominator[category] -= .5 * quadratic_proposed * quadratic_score;
  }

  //Precompute bounds for group 1 for numeric stability ------------------------
  double tmp_num = max(constant_numerator);
  double tmp_den = max(constant_denominator);
  if(tmp_num > 0) {
    if(tmp_num > tmp_den) {
      lbound = tmp_num;
    } else {
      lbound = tmp_den;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log pseudolikelihood ratio for group 1--------------------------
  for(int person = 0; person < no_persons_gr1; person++) {
    rest_score = rest_matrix_gr1(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = 0.0;
    denominator = 0.0;
    for(int category = 0; category < no_categories[variable] + 1; category ++) {
      exponent = category * rest_score - bound;
      numerator += std::exp(constant_numerator[category] + exponent);
      denominator += std::exp(constant_denominator[category] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  //Precompute common terms for group 2 for computational efficiency -----------
  for(int category = 0; category < no_categories[variable] + 1; category ++) {
    linear_score = category;
    quadratic_score =
      (category - reference_category[variable]) *
      (category - reference_category[variable]);

    constant_numerator[category] = thresholds(variable, 0) * linear_score ;
    constant_numerator[category] += .5 * linear_current * linear_score;
    constant_numerator[category] += thresholds(variable, 1) * quadratic_score;
    constant_numerator[category] += .5 * quadratic_current * quadratic_score;

    constant_denominator[category] = thresholds(variable, 0) * linear_score ;
    constant_denominator[category] += .5 * linear_proposed * linear_score;
    constant_denominator[category] += thresholds(variable, 1) * quadratic_score;
    constant_denominator[category] += .5 * quadratic_proposed * quadratic_score;
  }

  //Precompute bounds for group 2 for numeric stability ------------------------
  tmp_num = max(constant_numerator);
  tmp_den = max(constant_denominator);
  if(tmp_num > 0) {
    if(tmp_num > tmp_den) {
      lbound = tmp_num;
    } else {
      lbound = tmp_den;
    }
  } else {
    lbound = 0.0;
  }

  //Compute the log pseudolikelihood ratio for group 2--------------------------
  for(int person = 0; person < no_persons_gr2; person++) {
    rest_score = rest_matrix_gr2(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = 0.0;
    denominator = 0.0;
    for(int category = 0; category < no_categories[variable] + 1; category ++) {
      exponent = category * rest_score - bound;
      numerator += std::exp(constant_numerator[category] + exponent);
      denominator += std::exp(constant_denominator[category] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  return log_prob;
}



// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the category threshold
//  difference parameters-- Blume-Capel ordinal variable
// ----------------------------------------------------------------------------|
void compare_ttest_metropolis_main_difference_blumecapel(NumericMatrix thresholds,
                                                         NumericMatrix main_difference,
                                                         IntegerMatrix indicator,
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
  NumericVector constant_numerator (no_categories[variable] + 1);
  NumericVector constant_denominator (no_categories[variable] + 1);

  // Check if the variable is active for sampling
  if (indicator(variable, variable) == 0) {
    return; // Skip if variable is inactive
  }

  //---------------------------------------------------------------------------|
  // Adaptive Metropolis for the difference in the linear Blume-Capel parameter
  //---------------------------------------------------------------------------|
  current_state = main_difference(variable, 0);
  proposed_state = R::rnorm(current_state,
                            proposal_sd_main_difference(variable, 0));

  //----------------------------------------------------------------------------
  //Compute the log acceptance probability -------------------------------------
  //----------------------------------------------------------------------------

  double linear_current = current_state;
  double quadratic_current = main_difference(variable, 1);
  double linear_proposed = proposed_state;
  double quadratic_proposed = main_difference(variable, 1);

  log_prob = compare_ttest_log_pseudolikelihood_ratio_main_difference_blumecapel(linear_current,
                                                                                 quadratic_current,
                                                                                 linear_proposed,
                                                                                 quadratic_proposed,
                                                                                 variable,
                                                                                 reference_category,
                                                                                 thresholds,
                                                                                 sufficient_blume_capel_gr1,
                                                                                 sufficient_blume_capel_gr2,
                                                                                 no_persons_gr1,
                                                                                 no_persons_gr2,
                                                                                 rest_matrix_gr1,
                                                                                 rest_matrix_gr2,
                                                                                 no_categories);

  //Compute the prior ratio ---------------------------------------------------
  log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
  log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);

  //Metropolis step ------------------------------------------------------------
  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    main_difference(variable, 0) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ------------------------------
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

  update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

  proposal_sd_main_difference(variable, 0) = update_proposal_sd;


  //----------------------------------------------------------------------------
  // Adaptive Metropolis for the difference in quadratic Blume-Capel parameter
  //----------------------------------------------------------------------------
  current_state = main_difference(variable, 1);
  proposed_state = R::rnorm(current_state,
                            proposal_sd_main_difference(variable, 1));

  //----------------------------------------------------------------------------
  // Compute the log acceptance probability ------------------------------------
  //----------------------------------------------------------------------------

  linear_current = main_difference(variable, 0);
  quadratic_current = current_state;
  linear_proposed = main_difference(variable, 0);
  quadratic_proposed = proposed_state;

  log_prob = compare_ttest_log_pseudolikelihood_ratio_main_difference_blumecapel(linear_current,
                                                                                 quadratic_current,
                                                                                 linear_proposed,
                                                                                 quadratic_proposed,
                                                                                 variable,
                                                                                 reference_category,
                                                                                 thresholds,
                                                                                 sufficient_blume_capel_gr1,
                                                                                 sufficient_blume_capel_gr2,
                                                                                 no_persons_gr1,
                                                                                 no_persons_gr2,
                                                                                 rest_matrix_gr1,
                                                                                 rest_matrix_gr2,
                                                                                 no_categories);

  //Compute the prior ratio -------------------------------------------------
  log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
  log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);

  //Metropolis step ------------------------------------------------------------
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

  update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

  proposal_sd_main_difference(variable, 1) = update_proposal_sd;
}


// ----------------------------------------------------------------------------|
// Between model MH algorithm for the threshold differences
//    -- blume capel ordinal variable
// ----------------------------------------------------------------------------|
void compare_ttest_metropolis_main_difference_blumecapel_between_model(NumericMatrix thresholds,
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
  double log_prob;
  double U;
  NumericVector proposed_states(2);
  NumericVector current_states(2);

  //--------------------------------------------------------------------------
  // Adaptive Metropolis for the difference in Blume-Capel parameters
  //--------------------------------------------------------------------------
  current_states[0] = main_difference(variable, 0);
  current_states[1] = main_difference(variable, 1);

  if(indicator(variable, variable) == 0) {
    proposed_states[0] = R::rnorm(current_states[0], proposal_sd_main_difference(variable, 0));
    proposed_states[1] = R::rnorm(current_states[1], proposal_sd_main_difference(variable, 1));
  } else {
    proposed_states[0] = 0.0;
    proposed_states[1] = 0.0;
  }

  //--------------------------------------------------------------------------
  // Compute the log acceptance probability ----------------------------------
  //--------------------------------------------------------------------------

  double linear_current = current_states[0];
  double quadratic_current = current_states[1];
  double linear_proposed = proposed_states[0];
  double quadratic_proposed = proposed_states[1];

  log_prob = compare_ttest_log_pseudolikelihood_ratio_main_difference_blumecapel(linear_current,
                                                                                 quadratic_current,
                                                                                 linear_proposed,
                                                                                 quadratic_proposed,
                                                                                 variable,
                                                                                 reference_category,
                                                                                 thresholds,
                                                                                 sufficient_blume_capel_gr1,
                                                                                 sufficient_blume_capel_gr2,
                                                                                 no_persons_gr1,
                                                                                 no_persons_gr2,
                                                                                 rest_matrix_gr1,
                                                                                 rest_matrix_gr2,
                                                                                 no_categories);

  //Compute the parameter prior and proposal ratios ----------------------------
  if(indicator(variable, variable) == 0) {
    log_prob += R::dcauchy(proposed_states[0], 0.0, main_difference_scale, true);
    log_prob += R::dcauchy(proposed_states[1], 0.0, main_difference_scale, true);
    log_prob -= R::dnorm(proposed_states[0],
                         current_states[0],
                                       proposal_sd_main_difference(variable, 0),
                                       true);
    log_prob -= R::dnorm(proposed_states[1],
                         current_states[1],
                                       proposal_sd_main_difference(variable, 1),
                                       true);
  } else {
    log_prob -= R::dcauchy(current_states[0], 0.0, main_difference_scale, true);
    log_prob -= R::dcauchy(current_states[1], 0.0, main_difference_scale, true);
    log_prob += R::dnorm(current_states[0],
                         proposed_states[0],
                                        proposal_sd_main_difference(variable, 0),
                                        true);
    log_prob += R::dnorm(current_states[1],
                         proposed_states[1],
                                        proposal_sd_main_difference(variable, 1),
                                        true);
  }

  //Compute the prior inclusion ratios -----------------------------------------
  if(indicator(variable, variable) == 0) {
    log_prob += std::log(inclusion_probability_difference(variable, variable));
    log_prob -= std::log(1 - inclusion_probability_difference(variable, variable));
  } else {
    log_prob -= std::log(inclusion_probability_difference(variable, variable));
    log_prob += std::log(1 - inclusion_probability_difference(variable, variable));
  }

  //Metropolis step ------------------------------------------------------------
  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    indicator(variable, variable) = 1 - indicator(variable, variable);
    main_difference(variable, 0) = proposed_states[0];
    main_difference(variable, 1) = proposed_states[1];
  }
}

// ----------------------------------------------------------------------------|
// A Gibbs step for graphical model parameters for Bayesian parameter comparison
// ----------------------------------------------------------------------------|
List compare_ttest_gibbs_step_gm(NumericMatrix thresholds_gr1,
                                 NumericMatrix thresholds_gr2,
                                 NumericMatrix interactions_gr1,
                                 NumericMatrix interactions_gr2,
                                 IntegerMatrix indicator,
                                 NumericMatrix inclusion_probability_difference,
                                 IntegerMatrix observations_gr1,
                                 IntegerMatrix observations_gr2,
                                 IntegerVector no_categories_gr1,
                                 IntegerVector no_categories_gr2,
                                 int no_persons_gr1,
                                 int no_persons_gr2,
                                 int no_variables,
                                 IntegerMatrix n_cat_obs_gr1,
                                 IntegerMatrix n_cat_obs_gr2,
                                 IntegerMatrix sufficient_blume_capel_gr1,
                                 IntegerMatrix sufficient_blume_capel_gr2,
                                 NumericMatrix rest_matrix_gr1,
                                 NumericMatrix rest_matrix_gr2,
                                 NumericMatrix proposal_sd_interaction,
                                 NumericMatrix proposal_sd_main_difference,
                                 NumericMatrix proposal_sd_pairwise_difference,
                                 NumericMatrix proposal_sd_blumecapel_gr1,
                                 NumericMatrix proposal_sd_blumecapel_gr2,
                                 double interaction_scale,
                                 double main_difference_scale,
                                 double pairwise_difference_scale,
                                 double threshold_alpha,
                                 double threshold_beta,
                                 double phi,
                                 double target_ar,
                                 int t,
                                 double epsilon_lo,
                                 double epsilon_hi,
                                 LogicalVector ordinal_variable,
                                 IntegerVector reference_category,
                                 bool independent_thresholds,
                                 bool difference_selection,
                                 int no_interactions,
                                 IntegerMatrix index) {

  NumericMatrix thresholds(no_variables, max(no_categories_gr1));
  NumericMatrix main_difference(no_variables, max(no_categories_gr1));
  NumericMatrix interactions(no_interactions, no_interactions);
  NumericMatrix pairwise_difference(no_interactions, no_interactions);
  if(independent_thresholds == false) {
    for(int variable = 0; variable < no_variables; variable++) {
      if(ordinal_variable[variable] == true) {
        for(int category = 0; category < no_categories_gr1[variable]; category++) {
          main_difference(variable, category) =
            thresholds_gr2(variable, category) -
            thresholds_gr1(variable, category);

          thresholds(variable, category) =
            .5 * thresholds_gr2(variable, category) +
            .5 * thresholds_gr1(variable, category);
        }
      } else {
        main_difference(variable, 0) =
          thresholds_gr2(variable, 0) -
          thresholds_gr1(variable, 0);

        thresholds(variable, 0) =
          .5 * thresholds_gr2(variable, 0) +
          .5 * thresholds_gr1(variable, 0);

        main_difference(variable, 1) =
        thresholds_gr2(variable, 1) -
        thresholds_gr1(variable, 1);

        thresholds(variable, 1) =
          .5 * thresholds_gr2(variable, 1) +
          .5 * thresholds_gr1(variable, 1);
      }
    }
  }
  for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
      pairwise_difference(variable1, variable2) =
        interactions_gr2(variable1, variable2) -
        interactions_gr1(variable1, variable2);

      interactions(variable1, variable2) =
        .5 * interactions_gr1(variable1, variable2) +
        .5 * interactions_gr2(variable1, variable2);

      pairwise_difference(variable2, variable1) =
      pairwise_difference(variable1, variable2);

      interactions(variable2, variable1) =
        interactions(variable1, variable2);
    }
  }

  //Between model move for differences in interaction parameters
  if(difference_selection == true) {
    compare_ttest_metropolis_pairwise_difference_between_model(indicator,
                                                               inclusion_probability_difference,
                                                               pairwise_difference,
                                                               thresholds_gr1,
                                                               thresholds_gr2,
                                                               observations_gr1,
                                                               observations_gr2,
                                                               no_categories_gr1,
                                                               no_categories_gr2,
                                                               no_persons_gr1,
                                                               no_persons_gr2,
                                                               no_interactions,
                                                               index,
                                                               rest_matrix_gr1,
                                                               rest_matrix_gr2,
                                                               proposal_sd_pairwise_difference,
                                                               pairwise_difference_scale,
                                                               phi,
                                                               target_ar,
                                                               t,
                                                               epsilon_lo,
                                                               epsilon_hi,
                                                               ordinal_variable,
                                                               reference_category);
  }

  //Within model move for differences in interaction parameters
  compare_ttest_metropolis_pairwise_difference(pairwise_difference,
                                               thresholds_gr1,
                                               thresholds_gr2,
                                               observations_gr1,
                                               observations_gr2,
                                               no_categories_gr1,
                                               no_categories_gr2,
                                               indicator,
                                               no_persons_gr1,
                                               no_persons_gr2,
                                               no_variables,
                                               rest_matrix_gr1,
                                               rest_matrix_gr2,
                                               proposal_sd_pairwise_difference,
                                               pairwise_difference_scale,
                                               phi,
                                               target_ar,
                                               t,
                                               epsilon_lo,
                                               epsilon_hi,
                                               ordinal_variable,
                                               reference_category);

  //Within model move for the pairwise interaction parameters
  compare_ttest_metropolis_interaction(interactions,
                                       thresholds_gr1,
                                       thresholds_gr2,
                                       observations_gr1,
                                       observations_gr2,
                                       no_categories_gr1,
                                       no_categories_gr2,
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
                                       ordinal_variable,
                                       reference_category);

  //Update threshold parameters
  if(independent_thresholds == false) {
    for(int variable = 0; variable < no_variables; variable++) {
      if(ordinal_variable[variable] == true) {
        if(difference_selection == true) {
          //Between model move for the differences in category thresholds
          compare_ttest_metropolis_main_difference_regular_between_model(thresholds,
                                                                         main_difference,
                                                                         n_cat_obs_gr1,
                                                                         n_cat_obs_gr2,
                                                                         no_categories_gr1,
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
        compare_ttest_metropolis_main_difference_regular(thresholds,
                                                         main_difference,
                                                         n_cat_obs_gr1,
                                                         n_cat_obs_gr2,
                                                         no_categories_gr1,
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
        compare_ttest_metropolis_threshold_regular(thresholds,
                                                   main_difference,
                                                   no_categories_gr1,
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
          compare_ttest_metropolis_main_difference_blumecapel_between_model(thresholds,
                                                                            main_difference,
                                                                            sufficient_blume_capel_gr1,
                                                                            sufficient_blume_capel_gr2,
                                                                            no_categories_gr1,
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
        compare_ttest_metropolis_main_difference_blumecapel(thresholds,
                                                            main_difference,
                                                            indicator,
                                                            no_categories_gr1,
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
        compare_ttest_metropolis_threshold_blumecapel(thresholds,
                                                      main_difference,
                                                      no_categories_gr1,
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
                                                      proposal_sd_blumecapel_gr1,
                                                      phi,
                                                      target_ar,
                                                      t,
                                                      epsilon_lo,
                                                      epsilon_hi);
      }
    }
  } else {
    for(int variable = 0; variable < no_variables; variable++) {
      if(ordinal_variable[variable] == true) {
        //Within model move for the main category thresholds
        metropolis_thresholds_regular(thresholds_gr1,
                                      observations_gr1,
                                      no_categories_gr1,
                                      n_cat_obs_gr1,
                                      no_persons_gr1,
                                      variable,
                                      threshold_alpha,
                                      threshold_beta,
                                      rest_matrix_gr1);

        //Within model move for the main category thresholds
        metropolis_thresholds_regular(thresholds_gr2,
                                      observations_gr2,
                                      no_categories_gr2,
                                      n_cat_obs_gr2,
                                      no_persons_gr2,
                                      variable,
                                      threshold_alpha,
                                      threshold_beta,
                                      rest_matrix_gr2);
      } else {
        metropolis_thresholds_blumecapel(thresholds_gr1,
                                         observations_gr1,
                                         no_categories_gr1,
                                         sufficient_blume_capel_gr1,
                                         no_persons_gr1,
                                         variable,
                                         reference_category,
                                         threshold_alpha,
                                         threshold_beta,
                                         rest_matrix_gr1,
                                         proposal_sd_blumecapel_gr1,
                                         phi,
                                         target_ar,
                                         t,
                                         epsilon_lo,
                                         epsilon_hi);

        metropolis_thresholds_blumecapel(thresholds_gr2,
                                         observations_gr2,
                                         no_categories_gr2,
                                         sufficient_blume_capel_gr2,
                                         no_persons_gr2,
                                         variable,
                                         reference_category,
                                         threshold_alpha,
                                         threshold_beta,
                                         rest_matrix_gr2,
                                         proposal_sd_blumecapel_gr2,
                                         phi,
                                         target_ar,
                                         t,
                                         epsilon_lo,
                                         epsilon_hi);
      }
    }
  }

  if(independent_thresholds == false) {
    for(int variable = 0; variable < no_variables; variable++) {
      if(ordinal_variable[variable] == true) {
        for(int category = 0; category < no_categories_gr1[variable]; category++) {
          thresholds_gr1(variable, category) =
            thresholds(variable, category) -
            .5 * main_difference(variable, category);

            thresholds_gr2(variable, category) =
            thresholds(variable, category) +
            .5 * main_difference(variable, category);
        }
      } else {
        thresholds_gr1(variable, 0) =
          thresholds(variable, 0) -
          .5 * main_difference(variable, 0);

          thresholds_gr2(variable, 0) =
          thresholds(variable, 0) +
          .5 * main_difference(variable, 0);

          thresholds_gr1(variable, 1) =
          thresholds(variable, 1) -
          .5 * main_difference(variable, 1);

          thresholds_gr2(variable, 1) =
          thresholds(variable, 1) +
          .5 * main_difference(variable, 1);
      }
    }
  }

  for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
      interactions_gr1(variable1, variable2) =
        interactions(variable1, variable2) -
        .5 * pairwise_difference(variable1, variable2);

        interactions_gr2(variable1, variable2) =
        interactions(variable1, variable2) +
        .5 * pairwise_difference(variable1, variable2);

        interactions_gr1(variable2, variable1) =
        interactions_gr1(variable1, variable2);

        interactions_gr2(variable2, variable1) =
          interactions_gr2(variable1, variable2);
    }
  }

  return List::create(Named("indicator") = indicator,
                      Named("interactions_gr1") = interactions_gr1,
                      Named("interactions_gr2") = interactions_gr2,
                      Named("thresholds_gr1") = thresholds_gr1,
                      Named("thresholds_gr2") = thresholds_gr2,
                      Named("rest_matrix_gr1") = rest_matrix_gr1,
                      Named("rest_matrix_gr2") = rest_matrix_gr2,
                      Named("proposal_sd_interaction") = proposal_sd_interaction,
                      Named("proposal_sd_pairwise_difference") = proposal_sd_pairwise_difference,
                      Named("proposal_sd_main_difference") = proposal_sd_main_difference,
                      Named("proposal_sd_blumecapel_gr1") = proposal_sd_blumecapel_gr1,
                      Named("proposal_sd_blumecapel_gr2") = proposal_sd_blumecapel_gr2);
}

// ----------------------------------------------------------------------------|
// The Gibbs sampler for Bayesian parameter comparisons
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List compare_ttest_gibbs_sampler(IntegerMatrix observations_gr1,
                                 IntegerMatrix observations_gr2,
                                 IntegerVector no_categories_gr1,
                                 IntegerVector no_categories_gr2,
                                 double interaction_scale,
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
                                 LogicalVector ordinal_variable,
                                 IntegerVector reference_category,
                                 bool independent_thresholds,
                                 bool save = false,
                                 bool display_progress = false,
                                 bool difference_selection = true) {
  int cntr;
  int no_variables = observations_gr1.ncol();
  int no_persons_gr1 = observations_gr1.nrow();
  int no_persons_gr2 = observations_gr2.nrow();
  int no_interactions = Index.nrow();
  int no_thresholds_gr1 = 0;
  int no_thresholds_gr2 = 0;
  for(int v = 0; v < no_variables; v++) {
    if(ordinal_variable[v]) {
      no_thresholds_gr1 += no_categories_gr1[v];
      no_thresholds_gr2 += no_categories_gr2[v];
    } else {
      no_thresholds_gr1 += 2;
      no_thresholds_gr2 += 2;
    }
  }

  int max_no_categories_gr1 = max(no_categories_gr1);
  int max_no_categories_gr2 = max(no_categories_gr2);

  // Matrices with model parameters --------------------------------------------
  NumericMatrix interactions_gr1(no_variables, no_variables);
  NumericMatrix interactions_gr2(no_variables, no_variables);
  NumericMatrix thresholds_gr1(no_variables, max_no_categories_gr1);
  NumericMatrix thresholds_gr2(no_variables, max_no_categories_gr2);
  IntegerMatrix indicator(no_variables, no_variables);
  std::fill(indicator.begin(), indicator.end(), 1.0);

  // Matrices with standard deviations for adaptive Metropolis -----------------
  NumericMatrix proposal_sd_main_difference(no_variables, max_no_categories_gr1);
  NumericMatrix proposal_sd_interaction(no_variables, no_variables);
  NumericMatrix proposal_sd_pairwise_difference(no_variables, no_variables);
  NumericMatrix proposal_sd_blumecapel_gr1(no_variables, 2);
  NumericMatrix proposal_sd_blumecapel_gr2(no_variables, 2);

  std::fill(proposal_sd_main_difference.begin(), proposal_sd_main_difference.end(), 1.0);
  std::fill(proposal_sd_interaction.begin(), proposal_sd_interaction.end(), 1.0);
  std::fill(proposal_sd_pairwise_difference.begin(), proposal_sd_pairwise_difference.end(), 1.0);
  std::fill(proposal_sd_blumecapel_gr1.begin(), proposal_sd_blumecapel_gr1.end(), 1.0);
  std::fill(proposal_sd_blumecapel_gr2.begin(), proposal_sd_blumecapel_gr2.end(), 1.0);

  //Parameters for the Robbins-Monro approach for adaptive Metropolis ----------
  double phi = 0.75;
  double target_ar = 0.234;
  double epsilon_lo = 1.0 / static_cast<double>(no_persons_gr1);
  if(no_persons_gr1 > no_persons_gr2) {
    epsilon_lo = 1.0 / static_cast<double>(no_persons_gr2);
  }
  double epsilon_hi = 2.0;

  //Randomized index for the pairwise updates ----------------------------------
  IntegerVector v = seq(0, no_interactions - 1);
  IntegerVector order(no_interactions);
  IntegerMatrix index(no_interactions, 3);

  //Output matrices (resizing based on ``save'' and ``independent_thresholds'')
  int nrow = no_variables;
  int ncol_edges = no_variables;
  int ncol_thresholds_gr1 = max_no_categories_gr1;
  int ncol_thresholds_gr2 = max_no_categories_gr2;
  int ncol_main_difference = max_no_categories_gr1;

  if(save == true) {
    nrow = iter;
    ncol_edges = no_interactions;
    ncol_thresholds_gr1 = no_thresholds_gr1;
    ncol_thresholds_gr2 = no_thresholds_gr2;
    ncol_main_difference = no_thresholds_gr1;
  }
  NumericMatrix out_interactions(nrow, ncol_edges);
  NumericMatrix out_pairwise_difference(nrow, ncol_edges);
  NumericMatrix out_indicator_pairwise_difference(nrow, ncol_edges);

  int nrow_main_difference = nrow;
  int nrow_thresholds_gr2 = nrow;
  if(independent_thresholds == true) {
    ncol_main_difference = 1;
    nrow_main_difference = 1;
  } else {
    ncol_thresholds_gr2 = 1;
    nrow_thresholds_gr2 = 1;
  }
  NumericMatrix out_thresholds(nrow, ncol_thresholds_gr1);
  NumericMatrix out_main_difference(nrow_main_difference, ncol_main_difference);
  NumericMatrix out_thresholds_gr2(nrow_thresholds_gr2, ncol_thresholds_gr2);

  int ncol = no_variables;
  if(save == false) {
    nrow = 1;
  }
  if(independent_thresholds == true) {
    nrow = 1;
    ncol = 1;
  }
  NumericMatrix out_indicator_main_difference(nrow, ncol);

  //These matrices will contain the rest scores in the pseudolikelihoods -------
  NumericMatrix rest_matrix_gr1(no_persons_gr1, no_variables);
  NumericMatrix rest_matrix_gr2(no_persons_gr2, no_variables);

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------

  //When difference_selection = true we do 2 * burnin iterations. The first
  // burnin iterations without selection to ensure good starting values, and
  // proposal calibration. The second burnin iterations with selection.

  int first_burnin = burnin;
  int second_burnin = 0;
  if(difference_selection == true)
    second_burnin = burnin;
  bool input_difference_selection = difference_selection;
  difference_selection = false;

  //Progress bar ---------------------------------------------------------------
  Progress p(iter + first_burnin + second_burnin, display_progress);

  for(int iteration = 0; iteration < first_burnin + second_burnin; iteration++) {
    if(iteration >= first_burnin) {
      difference_selection = input_difference_selection;
    }
    if (Progress::check_abort()) {
      if(independent_thresholds == true) {
        return List::create(Named("pairwise_difference_indicator") = out_indicator_pairwise_difference,
                            Named("interactions") = out_interactions,
                            Named("thresholds_gr1") = out_thresholds,
                            Named("thresholds_gr2") = out_thresholds_gr2,
                            Named("pairwise_difference") = out_pairwise_difference);
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
      List out = compare_ttest_impute_missing_data(thresholds_gr1,
                                                   thresholds_gr2,
                                                   interactions_gr1,
                                                   interactions_gr2,
                                                   observations_gr1,
                                                   observations_gr2,
                                                   n_cat_obs_gr1,
                                                   n_cat_obs_gr2,
                                                   sufficient_blume_capel_gr1,
                                                   sufficient_blume_capel_gr2,
                                                   no_categories_gr1,
                                                   no_categories_gr2,
                                                   rest_matrix_gr1,
                                                   rest_matrix_gr2,
                                                   missing_index_gr1,
                                                   missing_index_gr2,
                                                   ordinal_variable,
                                                   reference_category);

      IntegerMatrix observations_gr1 = out["observations_gr1"];
      IntegerMatrix observations_gr2 = out["observations_gr2"];
      IntegerMatrix n_cat_obs_gr1 = out["n_cat_obs_gr1"];
      IntegerMatrix n_cat_obs_gr2 = out["n_cat_obs_gr2"];
      IntegerMatrix sufficient_blume_capel_gr1 = out["sufficient_blume_capel_gr1"];
      IntegerMatrix sufficient_blume_capel_gr2 = out["sufficient_blume_capel_gr2"];
      NumericMatrix rest_matrix_gr1 = out["rest_matrix_gr1"];
      NumericMatrix rest_matrix_gr2 = out["rest_matrix_gr2"];
    }

    List out = compare_ttest_gibbs_step_gm(thresholds_gr1,
                                           thresholds_gr2,
                                           interactions_gr1,
                                           interactions_gr2,
                                           indicator,
                                           inclusion_probability_difference,
                                           observations_gr1,
                                           observations_gr2,
                                           no_categories_gr1,
                                           no_categories_gr2,
                                           no_persons_gr1,
                                           no_persons_gr2,
                                           no_variables,
                                           n_cat_obs_gr1,
                                           n_cat_obs_gr2,
                                           sufficient_blume_capel_gr1,
                                           sufficient_blume_capel_gr2,
                                           rest_matrix_gr1,
                                           rest_matrix_gr2,
                                           proposal_sd_interaction,
                                           proposal_sd_main_difference,
                                           proposal_sd_pairwise_difference,
                                           proposal_sd_blumecapel_gr1,
                                           proposal_sd_blumecapel_gr2,
                                           interaction_scale,
                                           main_difference_scale,
                                           pairwise_difference_scale,
                                           threshold_alpha,
                                           threshold_beta,
                                           phi,
                                           target_ar,
                                           iteration,
                                           epsilon_lo,
                                           epsilon_hi,
                                           ordinal_variable,
                                           reference_category,
                                           independent_thresholds,
                                           difference_selection,
                                           no_interactions,
                                           index);

    IntegerMatrix indicator = out["indicator"];
    NumericMatrix interactions_gr1 = out["interactions_gr1"];
    NumericMatrix interactions_gr2 = out["interactions_gr2"];
    NumericMatrix thresholds_gr1 = out["thresholds_gr1"];
    NumericMatrix thresholds_gr2 = out["thresholds_gr2"];
    NumericMatrix rest_matrix_gr1 = out["rest_matrix_gr1"];
    NumericMatrix rest_matrix_gr2 = out["rest_matrix_gr2"];
    NumericMatrix proposal_sd_interaction = out["proposal_sd_interaction"];
    NumericMatrix proposal_sd_pairwise_difference = out["proposal_sd_pairwise_difference"];
    NumericMatrix proposal_sd_main_difference = out["proposal_sd_main_difference"];
    NumericMatrix proposal_sd_blumecapel_gr1 = out["proposal_sd_blumecapel_gr1"];
    NumericMatrix proposal_sd_blumecapel_gr2 = out["proposal_sd_blumecapel_gr2"];

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
  //To ensure that difference_selection is reinstated to the input value -------
  difference_selection = input_difference_selection;

  //The post burn-in iterations ------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    if (Progress::check_abort()) {
      if(independent_thresholds == true) {
        return List::create(Named("pairwise_difference_indicator") = out_indicator_pairwise_difference,
                            Named("interactions") = out_interactions,
                            Named("thresholds_gr1") = out_thresholds,
                            Named("thresholds_gr2") = out_thresholds_gr2,
                            Named("pairwise_difference") = out_pairwise_difference);
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
      List out = compare_ttest_impute_missing_data(thresholds_gr1,
                                                   thresholds_gr2,
                                                   interactions_gr1,
                                                   interactions_gr2,
                                                   observations_gr1,
                                                   observations_gr2,
                                                   n_cat_obs_gr1,
                                                   n_cat_obs_gr2,
                                                   sufficient_blume_capel_gr1,
                                                   sufficient_blume_capel_gr2,
                                                   no_categories_gr1,
                                                   no_categories_gr2,
                                                   rest_matrix_gr1,
                                                   rest_matrix_gr2,
                                                   missing_index_gr1,
                                                   missing_index_gr2,
                                                   ordinal_variable,
                                                   reference_category);

      IntegerMatrix observations_gr1 = out["observations_gr1"];
      IntegerMatrix observations_gr2 = out["observations_gr2"];
      IntegerMatrix n_cat_obs_gr1 = out["n_cat_obs_gr1"];
      IntegerMatrix n_cat_obs_gr2 = out["n_cat_obs_gr2"];
      IntegerMatrix sufficient_blume_capel_gr1 = out["sufficient_blume_capel_gr1"];
      IntegerMatrix sufficient_blume_capel_gr2 = out["sufficient_blume_capel_gr2"];
      NumericMatrix rest_matrix_gr1 = out["rest_matrix_gr1"];
      NumericMatrix rest_matrix_gr2 = out["rest_matrix_gr2"];
    }

    List out = compare_ttest_gibbs_step_gm(thresholds_gr1,
                                           thresholds_gr2,
                                           interactions_gr1,
                                           interactions_gr2,
                                           indicator,
                                           inclusion_probability_difference,
                                           observations_gr1,
                                           observations_gr2,
                                           no_categories_gr1,
                                           no_categories_gr2,
                                           no_persons_gr1,
                                           no_persons_gr2,
                                           no_variables,
                                           n_cat_obs_gr1,
                                           n_cat_obs_gr2,
                                           sufficient_blume_capel_gr1,
                                           sufficient_blume_capel_gr2,
                                           rest_matrix_gr1,
                                           rest_matrix_gr2,
                                           proposal_sd_interaction,
                                           proposal_sd_main_difference,
                                           proposal_sd_pairwise_difference,
                                           proposal_sd_blumecapel_gr1,
                                           proposal_sd_blumecapel_gr2,
                                           interaction_scale,
                                           main_difference_scale,
                                           pairwise_difference_scale,
                                           threshold_alpha,
                                           threshold_beta,
                                           phi,
                                           target_ar,
                                           iteration,
                                           epsilon_lo,
                                           epsilon_hi,
                                           ordinal_variable,
                                           reference_category,
                                           independent_thresholds,
                                           difference_selection,
                                           no_interactions,
                                           index);

    IntegerMatrix indicator = out["indicator"];
    NumericMatrix interactions_gr1 = out["interactions_gr1"];
    NumericMatrix interactions_gr2 = out["interactions_gr2"];
    NumericMatrix thresholds_gr1 = out["thresholds_gr1"];
    NumericMatrix thresholds_gr2 = out["thresholds_gr2"];
    NumericMatrix rest_matrix_gr1 = out["rest_matrix_gr1"];
    NumericMatrix rest_matrix_gr2 = out["rest_matrix_gr2"];
    NumericMatrix proposal_sd_interaction = out["proposal_sd_interaction"];
    NumericMatrix proposal_sd_pairwise_difference = out["proposal_sd_pairwise_difference"];
    NumericMatrix proposal_sd_main_difference = out["proposal_sd_main_difference"];
    NumericMatrix proposal_sd_blumecapel_gr1 = out["proposal_sd_blumecapel_gr1"];
    NumericMatrix proposal_sd_blumecapel_gr2 = out["proposal_sd_blumecapel_gr2"];


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
          double pairwise_difference =
            interactions_gr2(variable1, variable2) -
            interactions_gr1(variable1, variable2);
          double interaction =
            interactions_gr1(variable1, variable2) +
            .5 * pairwise_difference;

            out_indicator_pairwise_difference(iteration, cntr) =
            indicator(variable1, variable2);
            out_interactions(iteration, cntr) = interaction;
            out_pairwise_difference(iteration, cntr) = pairwise_difference;
            cntr++;
        }
      }

      if(independent_thresholds == true) {
        int cntr_gr1 = 0;
        int cntr_gr2 = 0;
        for(int variable = 0; variable < no_variables; variable++) {
          if(ordinal_variable[variable] == true) {
            for(int category = 0; category < no_categories_gr1[variable]; category++) {
              out_thresholds(iteration, cntr_gr1) = thresholds_gr1(variable, category);
              cntr_gr1++;
            }
          } else {
            out_thresholds(iteration, cntr_gr1) = thresholds_gr1(variable, 0);
            cntr_gr1++;
            out_thresholds(iteration, cntr_gr1) = thresholds_gr1(variable, 1);
            cntr_gr1++;
          }
        }
        for(int variable = 0; variable < no_variables; variable++) {
          if(ordinal_variable[variable] == true) {
            for(int category = 0; category < no_categories_gr2[variable]; category++) {
              out_thresholds_gr2(iteration, cntr_gr2) = thresholds_gr2(variable, category);
              cntr_gr2++;
            }
          } else {
            out_thresholds_gr2(iteration, cntr_gr2) = thresholds_gr2(variable, 0);
            cntr_gr2++;
            out_thresholds_gr2(iteration, cntr_gr2) = thresholds_gr2(variable, 1);
            cntr_gr2++;
          }
        }
      } else {
        cntr = 0;
        for(int variable = 0; variable < no_variables; variable++) {
          if(ordinal_variable[variable] == true) {
            out_indicator_main_difference(iteration, variable) =
              indicator(variable, variable);
            for(int category = 0; category < no_categories_gr1[variable]; category++) {
              double main_difference =
                thresholds_gr2(variable, category) -
                thresholds_gr1(variable, category);
              double threshold =
                .5 * thresholds_gr2(variable, category) +
                .5 * thresholds_gr1(variable, category);

              out_thresholds(iteration, cntr) = threshold;
              out_main_difference(iteration, cntr) = main_difference;
              cntr++;
            }
          } else {
            double main_difference =
              thresholds_gr2(variable, 0) -
              thresholds_gr1(variable, 0);
            double threshold =
              .5 * thresholds_gr2(variable, 0) +
              .5 * thresholds_gr1(variable, 0);

            out_thresholds(iteration, cntr) = threshold;
            out_main_difference(iteration, cntr) = main_difference;
            cntr++;

            main_difference =
              thresholds_gr2(variable, 1) -
              thresholds_gr1(variable, 1);
            threshold =
              .5 * thresholds_gr2(variable, 1) +
              .5 * thresholds_gr1(variable, 1);

            out_thresholds(iteration, cntr) = threshold;
            out_main_difference(iteration, cntr) = main_difference;
            cntr++;
          }
        }
      }

    } else {
      //Compute running averages -----------------------------------------------
      for(int variable1 = 0; variable1 < no_variables - 1; variable1++) {
        for(int variable2 = variable1 + 1; variable2 < no_variables; variable2++) {
          double pairwise_difference =
            interactions_gr2(variable1, variable2) -
            interactions_gr1(variable1, variable2);
          double interaction =
            interactions_gr1(variable1, variable2) +
            .5 * pairwise_difference;

            out_indicator_pairwise_difference(variable1, variable2) *= iteration;
            out_indicator_pairwise_difference(variable1, variable2) += indicator(variable1, variable2);
            out_indicator_pairwise_difference(variable1, variable2) /= iteration + 1;
            out_indicator_pairwise_difference(variable2, variable1) = out_indicator_pairwise_difference(variable1, variable2);

            out_interactions(variable1, variable2) *= iteration;
            out_interactions(variable1, variable2) += interaction;
            out_interactions(variable1, variable2) /= iteration + 1;
            out_interactions(variable2, variable1) = out_interactions(variable1, variable2);

            out_pairwise_difference(variable1, variable2) *= iteration;
            out_pairwise_difference(variable1, variable2) += pairwise_difference;
            out_pairwise_difference(variable1, variable2) /= iteration + 1;
            out_pairwise_difference(variable2, variable1) = out_pairwise_difference(variable1, variable2);
        }
      }

      if(independent_thresholds == true) {
        for(int variable = 0; variable < no_variables; variable++) {
          if(ordinal_variable[variable] == true) {
            for(int category = 0; category < no_categories_gr1[variable]; category++) {
              out_thresholds(variable, category) *= iteration;
              out_thresholds(variable, category) += thresholds_gr1(variable, category);
              out_thresholds(variable, category) /= iteration + 1;
            }
            for(int category = 0; category < no_categories_gr2[variable]; category++) {
              out_thresholds_gr2(variable, category) *= iteration;
              out_thresholds_gr2(variable, category) += thresholds_gr2(variable, category);
              out_thresholds_gr2(variable, category) /= iteration + 1;
            }
          } else {
            out_thresholds(variable, 0) *= iteration;
            out_thresholds(variable, 0) += thresholds_gr1(variable, 0);
            out_thresholds(variable, 0) /= iteration + 1;
            out_thresholds(variable, 1) *= iteration;
            out_thresholds(variable, 1) += thresholds_gr1(variable, 1);
            out_thresholds(variable, 1) /= iteration + 1;

            out_thresholds_gr2(variable, 0) *= iteration;
            out_thresholds_gr2(variable, 0) += thresholds_gr2(variable, 0);
            out_thresholds_gr2(variable, 0) /= iteration + 1;
            out_thresholds_gr2(variable, 1) *= iteration;
            out_thresholds_gr2(variable, 1) += thresholds_gr2(variable, 1);
            out_thresholds_gr2(variable, 1) /= iteration + 1;
          }
        }
      } else {
        for(int variable = 0; variable < no_variables; variable++) {
          if(ordinal_variable[variable] == true) {
            for(int category = 0; category < no_categories_gr1[variable]; category++) {
              double main_difference =
                thresholds_gr2(variable, category) -
                thresholds_gr1(variable, category);
              double threshold =
                .5 * thresholds_gr1(variable, category) +
                .5 * thresholds_gr2(variable, category);

              out_thresholds(variable, category) *= iteration;
              out_thresholds(variable, category) += threshold;
              out_thresholds(variable, category) /= iteration + 1;

              out_main_difference(variable, category) *= iteration;
              out_main_difference(variable, category) += main_difference;
              out_main_difference(variable, category) /= iteration + 1;

              out_indicator_main_difference(0, variable) *= iteration;
              out_indicator_main_difference(0, variable) += indicator(variable, variable);
              out_indicator_main_difference(0, variable) /= iteration + 1;
            }
          } else {
            double main_difference =
              thresholds_gr2(variable, 0) -
              thresholds_gr1(variable, 0);
            double threshold =
              .5 * thresholds_gr1(variable, 0) +
              .5 * thresholds_gr2(variable, 0);

            out_thresholds(variable, 0) *= iteration;
            out_thresholds(variable, 0) += threshold;
            out_thresholds(variable, 0) /= iteration + 1;
            out_main_difference(variable, 0) *= iteration;
            out_main_difference(variable, 0) += main_difference;
            out_main_difference(variable, 0) /= iteration + 1;

            main_difference =
              thresholds_gr2(variable, 1) -
              thresholds_gr1(variable, 1);
            threshold =
              .5 * thresholds_gr1(variable, 1) +
              .5 * thresholds_gr2(variable, 1);

            out_thresholds(variable, 1) *= iteration;
            out_thresholds(variable, 1) += threshold;
            out_thresholds(variable, 1) /= iteration + 1;

            out_main_difference(variable, 1) *= iteration;
            out_main_difference(variable, 1) += main_difference;
            out_main_difference(variable, 1) /= iteration + 1;

            out_indicator_main_difference(0, variable) *= iteration;
            out_indicator_main_difference(0, variable) += indicator(variable, variable);
            out_indicator_main_difference(0, variable) /= iteration + 1;
          }
        }
      }
    }
  }

  if(independent_thresholds == true) {
    return List::create(Named("pairwise_difference_indicator") = out_indicator_pairwise_difference,
                        Named("interactions") = out_interactions,
                        Named("thresholds_gr1") = out_thresholds,
                        Named("thresholds_gr2") = out_thresholds_gr2,
                        Named("pairwise_difference") = out_pairwise_difference);
  } else {
    return List::create(Named("pairwise_difference_indicator") = out_indicator_pairwise_difference,
                        Named("main_difference_indicator") = out_indicator_main_difference,
                        Named("interactions") = out_interactions,
                        Named("thresholds") = out_thresholds,
                        Named("pairwise_difference") = out_pairwise_difference,
                        Named("main_difference") = out_main_difference);
  }
}