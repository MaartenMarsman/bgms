// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include "gibbs_functions.h"
using namespace Rcpp;


// ----------------------------------------------------------------------------|
// List of to do's
// 1. Add between model moves.
// 2. Use "index" to vary between-model updates to pairwise effects.
// 3. Add g prior.
// 4. Output matrices (resizing based on ``save'' and ``independent_thresholds'')
// ----------------------------------------------------------------------------|


// ----------------------------------------------------------------------------|
// Create the group-specific thresholds for one variable
// ----------------------------------------------------------------------------|
NumericVector group_thresholds_for_variable(int variable,
                                            LogicalVector ordinal_variable,
                                            int group,
                                            int no_groups,
                                            IntegerMatrix no_categories,
                                            NumericMatrix main_effects,
                                            IntegerMatrix main_index,
                                            NumericMatrix projection,
                                            bool independent_thresholds) {
  int category_index;
  int vector_length = 2;
  if(ordinal_variable[variable] == true) {
    vector_length = no_categories(variable, group);
  }
  NumericVector GroupThresholds(vector_length);

  if(independent_thresholds == false) {
    if(ordinal_variable[variable] == true) {
      for(int category = 0; category < no_categories(variable, group); category++) {
        category_index = main_index(variable, 0) + category;
        GroupThresholds[category] = main_effects(category_index, 0);
        for(int h = 0; h < no_groups - 1; h++) {
          GroupThresholds[category] += projection(group, h) *
            main_effects(category_index, h + 1);
        }
      }
    } else {
      for(int category = 0; category < 2; category++) {
        category_index = main_index(variable, 0) + category;
        GroupThresholds[category] = main_effects(category_index, 0);
        for(int h = 0; h < no_groups - 1; h++) {
          GroupThresholds[category] += projection(group, h) *
            main_effects(category_index, h + 1);
        }
      }
    }
  } else {
    if(ordinal_variable[variable] == true) {
      for(int category = 0; category < no_categories(variable, group); category++) {
        category_index = main_index(variable, 0) + category;
        GroupThresholds[category] = main_effects(category_index, group);
      }
    } else {
      category_index = main_index(variable, 0);
      GroupThresholds[0] = main_effects(category_index, group);
      GroupThresholds[1] = main_effects(category_index + 1, group);
    }
  }

  return GroupThresholds;
}


// ----------------------------------------------------------------------------|
// Impute missing data for independent samples designs
// ----------------------------------------------------------------------------|
List compare_anova_impute_missing_data(NumericMatrix main_effects,
                                       NumericMatrix pairwise_effects,
                                       IntegerMatrix main_index,
                                       IntegerMatrix pairwise_index,
                                       NumericMatrix projection,
                                       IntegerMatrix observations,
                                       int no_groups,
                                       IntegerVector person_group_indicator,
                                       List n_cat_obs,
                                       List sufficient_blume_capel,
                                       IntegerMatrix no_categories,
                                       NumericMatrix rest_matrix,
                                       IntegerMatrix missing_index,
                                       LogicalVector ordinal_variable,
                                       IntegerVector reference_category,
                                       bool independent_thresholds) {

  int no_variables = observations.ncol();
  int no_missings = missing_index.nrow();
  int max_no_categories = max(no_categories);

  NumericVector probabilities(max_no_categories + 1);
  double exponent, rest_score, cumsum, u;
  int score, person, variable, new_observation, old_observation;

  //Impute missing data --------------------------------------------------------
  for(int missing = 0; missing < no_missings; missing++) {
    //Which observation to impute? -------------------------------------------
    person = missing_index(missing, 0) - 1; //R to C++ indexing
    variable = missing_index(missing, 1) - 1; //R to C++ indexing
    int gr = person_group_indicator[person];

    NumericVector GroupThresholds = group_thresholds_for_variable(variable,
                                                                  ordinal_variable,
                                                                  gr,
                                                                  no_groups,
                                                                  no_categories,
                                                                  main_effects,
                                                                  main_index,
                                                                  projection,
                                                                  independent_thresholds);

    //Generate a new observation from the ordinal MRF ------------------------
    rest_score = rest_matrix(person, variable);
    if(ordinal_variable[variable] == true) {
      //Regular binary or ordinal variable -----------------------------------
      cumsum = 1.0;
      probabilities[0] = 1.0;
      for(int category = 0; category < no_categories(variable, gr); category++) {
        exponent = GroupThresholds(category);
        exponent += (category + 1) * rest_score;
        cumsum += std::exp(exponent);
        probabilities[category + 1] = cumsum;
      }
    } else {
      //Blume-Capel variable -------------------------------------------------
      cumsum = 0.0;
      for(int category = 0; category < no_categories(variable, gr) + 1; category++) {
        exponent = GroupThresholds[0] * category;
        exponent += GroupThresholds[1] *
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
    old_observation = observations(person, variable);
    if(old_observation != new_observation) {
      //Update raw data ------------------------------------------------------
      observations(person, variable) = new_observation;

      //Update pre-computed statistics ---------------------------------------
      if(ordinal_variable[variable] == true) {
        //Regular binary or ordinal variable ---------------------------------
        IntegerMatrix n_cat_obs_gr = n_cat_obs[gr];
        n_cat_obs_gr(old_observation, variable)--;
        n_cat_obs_gr(new_observation, variable)++;
        n_cat_obs[gr] = n_cat_obs_gr;
      } else {
        //Blume-Capel variable -----------------------------------------------
        IntegerMatrix sufficient_blume_capel_gr = sufficient_blume_capel[gr];
        sufficient_blume_capel_gr(0, variable) -= old_observation;
        sufficient_blume_capel_gr(0, variable) += new_observation;
        sufficient_blume_capel_gr(1, variable) -=
          (old_observation - reference_category[variable]) *
          (old_observation - reference_category[variable]);
        sufficient_blume_capel_gr(1, variable) +=
          (new_observation - reference_category[variable]) *
          (new_observation - reference_category[variable]);
        sufficient_blume_capel[gr] = sufficient_blume_capel_gr;
      }

      //Update rest score ----------------------------------------------------
      for(int vertex = 0; vertex < no_variables; vertex++) {
        int index = pairwise_index(vertex, variable);
        double GroupInteraction = pairwise_effects(index, 0);
        for(int h = 0; h < no_groups - 1; h++) {
          GroupInteraction += projection(gr, h) * pairwise_effects(index, h + 1);
        }
        rest_matrix(person, vertex) -= old_observation * GroupInteraction;
        rest_matrix(person, vertex) += new_observation * GroupInteraction;
      }
    }
  }

  return List::create(Named("observations") = observations,
                      Named("n_cat_obs") = n_cat_obs,
                      Named("sufficient_blume_capel") = sufficient_blume_capel,
                      Named("rest_matrix") = rest_matrix);
}


// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for an interaction
//  for an independent samples design
// ----------------------------------------------------------------------------|
double compare_anova_log_pseudolikelihood_ratio_interaction(NumericMatrix main_effects,
                                                            IntegerMatrix main_index,
                                                            NumericMatrix projection,
                                                            IntegerMatrix observations,
                                                            int no_groups,
                                                            IntegerMatrix group_index,
                                                            IntegerMatrix no_categories,
                                                            bool independent_thresholds,
                                                            int no_persons,
                                                            int variable1,
                                                            int variable2,
                                                            double proposed_state,
                                                            double current_state,
                                                            NumericMatrix rest_matrix,
                                                            LogicalVector ordinal_variable,
                                                            IntegerVector reference_category) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  double delta_state = proposed_state - current_state;

  // --------------------------------------------------------------------------|
  // Compute log pseudo-likelihood ratio ---------------------------------------
  // --------------------------------------------------------------------------|
  for(int gr = 0; gr < no_groups; gr++) {
    NumericVector GroupThresholds_v1 = group_thresholds_for_variable(variable1,
                                                                     ordinal_variable,
                                                                     gr,
                                                                     no_groups,
                                                                     no_categories,
                                                                     main_effects,
                                                                     main_index,
                                                                     projection,
                                                                     independent_thresholds);

    NumericVector GroupThresholds_v2 = group_thresholds_for_variable(variable2,
                                                                     ordinal_variable,
                                                                     gr,
                                                                     no_groups,
                                                                     no_categories,
                                                                     main_effects,
                                                                     main_index,
                                                                     projection,
                                                                     independent_thresholds);

    //Loop over the pseudo-likelihoods for persons in group gr -----------------
    for(int person = group_index(gr, 0); person < group_index(gr, 1) + 1; person++) {
      obs_score1 = observations(person, variable1);
      obs_score2 = observations(person, variable2);

      pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * delta_state;

      //Log pseudo-likelihood for variable 1 -----------------------------------
      rest_score = rest_matrix(person, variable1) - obs_score2 * current_state; //This should contains the pairwise difference effects.

      if(rest_score > 0) {
        bound = no_categories(variable1, gr) * rest_score;
      } else {
        bound = 0.0;
      }

      if(ordinal_variable[variable1] == true) {
        //Regular binary or ordinal MRF variable -------------------------------
        denominator_prop = std::exp(-bound);
        denominator_curr = std::exp(-bound);
        for(int cat = 0; cat < no_categories(variable1, gr); cat++) {
          score = cat + 1;
          exponent = GroupThresholds_v1[cat] +
            score * rest_score -
            bound;
          denominator_prop +=
            std::exp(exponent + score * obs_score2 * proposed_state);
          denominator_curr +=
            std::exp(exponent + score * obs_score2 * current_state);
        }
      } else {
        //Blume-Capel ordinal MRF variable -------------------------------------
        denominator_prop = 0.0;
        denominator_curr = 0.0;
        for(int cat = 0; cat < no_categories(variable1, gr) + 1; cat++) {
          exponent = GroupThresholds_v1[0] * cat;
          exponent += GroupThresholds_v1[1] *
            (cat - reference_category[variable1]) *
            (cat - reference_category[variable1]);
          exponent += cat * rest_score - bound;
          denominator_prop +=
            std::exp(exponent + cat * obs_score2 * proposed_state);
          denominator_curr +=
            std::exp(exponent + cat * obs_score2 * current_state);
        }
      }
      pseudolikelihood_ratio -= std::log(denominator_prop);
      pseudolikelihood_ratio += std::log(denominator_curr);

      //Log pseudo-likelihood for variable 2 -----------------------------------
      rest_score = rest_matrix(person, variable2) - obs_score1 * current_state;

      if(rest_score > 0) {
        bound = no_categories(variable2, gr) * rest_score;
      } else {
        bound = 0.0;
      }

      if(ordinal_variable[variable2] == true) {
        //Regular binary or ordinal MRF variable ---------------------------------
        denominator_prop = std::exp(-bound);
        denominator_curr = std::exp(-bound);
        for(int cat = 0; cat < no_categories(variable2, gr); cat++) {
          score = cat + 1;
          exponent = GroupThresholds_v2[cat] +
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
        for(int cat = 0; cat < no_categories(variable2, gr) + 1; cat++) {
          exponent = GroupThresholds_v2[0] * cat;
          exponent += GroupThresholds_v2[1] *
            (cat - reference_category[variable2]) *
            (cat - reference_category[variable2]);
          exponent +=  cat * rest_score - bound;
          denominator_prop +=
            std::exp(exponent + cat * obs_score1 * proposed_state);
          denominator_curr +=
            std::exp(exponent + cat * obs_score1 * current_state);
        }
      }
      pseudolikelihood_ratio -= std::log(denominator_prop);
      pseudolikelihood_ratio += std::log(denominator_curr);
    }
  }
  return pseudolikelihood_ratio;
}


// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of an overall pairwise
//  interaction parameter (nuisance)
// ----------------------------------------------------------------------------|
void compare_anova_metropolis_interaction(NumericMatrix main_effects,
                                          NumericMatrix pairwise_effects,
                                          IntegerMatrix main_index,
                                          IntegerMatrix pairwise_index,
                                          NumericMatrix projection,
                                          IntegerMatrix observations,
                                          int no_groups,
                                          IntegerMatrix group_index,
                                          IntegerMatrix no_categories,
                                          bool independent_thresholds,
                                          int no_persons,
                                          NumericMatrix rest_matrix,
                                          LogicalVector ordinal_variable,
                                          IntegerVector reference_category,
                                          NumericMatrix proposal_sd_pairwise,
                                          double interaction_scale,
                                          int no_variables,
                                          double phi,
                                          double target_ar,
                                          int t,
                                          double epsilon_lo,
                                          double epsilon_hi) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  for(int variable1 = 0; variable1 <  no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  no_variables; variable2++) {

      int int_index = pairwise_index(variable1, variable2);

      current_state = pairwise_effects(int_index, 0);
      proposed_state = R::rnorm(current_state, proposal_sd_pairwise(int_index, 0));

      log_prob = compare_anova_log_pseudolikelihood_ratio_interaction(main_effects,
                                                                      main_index,
                                                                      projection,
                                                                      observations,
                                                                      no_groups,
                                                                      group_index,
                                                                      no_categories,
                                                                      independent_thresholds,
                                                                      no_persons,
                                                                      variable1,
                                                                      variable2,
                                                                      proposed_state,
                                                                      current_state,
                                                                      rest_matrix,
                                                                      ordinal_variable,
                                                                      reference_category);

      log_prob += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_prob -= R::dcauchy(current_state, 0.0, interaction_scale, true);

      U = R::unif_rand();
      if(std::log(U) < log_prob) {
        double state_difference = proposed_state - current_state;
        pairwise_effects(int_index, 0) = proposed_state;

        //Update the matrix of rest scores
        for(int person = 0; person < no_persons; person++) {
          rest_matrix(person, variable1) += observations(person, variable2) *
            state_difference;
          rest_matrix(person, variable2) += observations(person, variable1) *
            state_difference;
        }
      }

      if(log_prob > 0) {
        log_prob = 1.0;
      } else {
        log_prob = std::exp(log_prob);
      }

      double update_proposal_sd =
        proposal_sd_pairwise(int_index, 0) +
        (log_prob - target_ar) * std::exp(-log(t) * phi);

      if(std::isnan(update_proposal_sd) == true) {
        update_proposal_sd = 1.0;
      }

      update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

      proposal_sd_pairwise(int_index, 0) = update_proposal_sd;
    }
  }
}


// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for the difference
//  in a pairwise interaction for an independent samples design
// ----------------------------------------------------------------------------|
double compare_anova_log_pseudolikelihood_ratio_pairwise_difference(NumericMatrix main_effects,
                                                                    IntegerMatrix main_index,
                                                                    NumericMatrix projection,
                                                                    IntegerMatrix observations,
                                                                    int no_groups,
                                                                    IntegerMatrix group_index,
                                                                    IntegerMatrix no_categories,
                                                                    bool independent_thresholds,
                                                                    int no_persons,
                                                                    int variable1,
                                                                    int variable2,
                                                                    int h,
                                                                    double proposed_state,
                                                                    double current_state,
                                                                    NumericMatrix rest_matrix,
                                                                    LogicalVector ordinal_variable,
                                                                    IntegerVector reference_category) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  double delta_state = proposed_state - current_state;

  // --------------------------------------------------------------------------|
  // Compute log pseudo-likelihood ratio ---------------------------------------
  // --------------------------------------------------------------------------|
  for(int gr = 0; gr < no_groups; gr++) {

    NumericVector GroupThresholds_v1 = group_thresholds_for_variable(variable1,
                                                                     ordinal_variable,
                                                                     gr,
                                                                     no_groups,
                                                                     no_categories,
                                                                     main_effects,
                                                                     main_index,
                                                                     projection,
                                                                     independent_thresholds);

    NumericVector GroupThresholds_v2 = group_thresholds_for_variable(variable2,
                                                                     ordinal_variable,
                                                                     gr,
                                                                     no_groups,
                                                                     no_categories,
                                                                     main_effects,
                                                                     main_index,
                                                                     projection,
                                                                     independent_thresholds);

    double P = projection(gr, h);

    //Loop over the pseudo-likelihoods for persons in group gr -----------------
    for(int person = group_index(gr, 0); person < group_index(gr, 1) + 1; person++) {
      obs_score1 = observations(person, variable1);
      obs_score2 = observations(person, variable2);

      pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * delta_state * P;

      //Log pseudo-likelihood for variable 1 -----------------------------------
      rest_score = rest_matrix(person, variable1) - obs_score2 * current_state * P;

      if(rest_score > 0) {
        bound = no_categories(variable1, gr) * rest_score;
      } else {
        bound = 0.0;
      }

      if(ordinal_variable[variable1] == true) {
        //Regular binary or ordinal MRF variable -------------------------------
        denominator_prop = std::exp(-bound);
        denominator_curr = std::exp(-bound);
        for(int cat = 0; cat < no_categories(variable1, gr); cat++) {
          score = cat + 1;
          exponent = GroupThresholds_v1[cat] +
            score * rest_score -
            bound;
          denominator_prop +=
            std::exp(exponent + score * obs_score2 * proposed_state * P);
          denominator_curr +=
            std::exp(exponent + score * obs_score2 * current_state * P);
        }
      } else {
        //Blume-Capel ordinal MRF variable -------------------------------------
        denominator_prop = 0.0;
        denominator_curr = 0.0;
        for(int cat = 0; cat < no_categories(variable1, gr) + 1; cat++) {
          exponent = GroupThresholds_v1[0] * cat;
          exponent += GroupThresholds_v1[1] *
            (cat - reference_category[variable1]) *
            (cat - reference_category[variable1]);
          exponent += cat * rest_score - bound;
          denominator_prop +=
            std::exp(exponent + cat * obs_score2 * proposed_state * P);
          denominator_curr +=
            std::exp(exponent + cat * obs_score2 * current_state * P);
        }
      }
      pseudolikelihood_ratio -= std::log(denominator_prop);
      pseudolikelihood_ratio += std::log(denominator_curr);

      //Log pseudo-likelihood for variable 2 -----------------------------------
      rest_score = rest_matrix(person, variable2) - obs_score1 * current_state * P;

      if(rest_score > 0) {
        bound = no_categories(variable2, gr) * rest_score;
      } else {
        bound = 0.0;
      }

      if(ordinal_variable[variable2] == true) {
        //Regular binary or ordinal MRF variable ---------------------------------
        denominator_prop = std::exp(-bound);
        denominator_curr = std::exp(-bound);
        for(int cat = 0; cat < no_categories(variable2, gr); cat++) {
          score = cat + 1;
          exponent = GroupThresholds_v2[cat] +
            score * rest_score -
            bound;
          denominator_prop +=
            std::exp(exponent + score * obs_score1 * proposed_state * P);
          denominator_curr +=
            std::exp(exponent + score * obs_score1 * current_state * P);
        }
      } else {
        //Blume-Capel ordinal MRF variable ---------------------------------------
        denominator_prop = 0.0;
        denominator_curr = 0.0;
        for(int cat = 0; cat < no_categories(variable2, gr) + 1; cat++) {
          exponent = GroupThresholds_v2[0] * cat;
          exponent += GroupThresholds_v2[1] *
            (cat - reference_category[variable2]) *
            (cat - reference_category[variable2]);
          exponent +=  cat * rest_score - bound;
          denominator_prop +=
            std::exp(exponent + cat * obs_score1 * proposed_state * P);
          denominator_curr +=
            std::exp(exponent + cat * obs_score1 * current_state * P);
        }
      }
      pseudolikelihood_ratio -= std::log(denominator_prop);
      pseudolikelihood_ratio += std::log(denominator_curr);
    }
  }
  return pseudolikelihood_ratio;
}


// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of pairwise differences
// ----------------------------------------------------------------------------|
void compare_anova_metropolis_pairwise_difference(NumericMatrix main_effects,
                                                  NumericMatrix pairwise_effects,
                                                  IntegerMatrix main_index,
                                                  IntegerMatrix pairwise_index,
                                                  NumericMatrix projection,
                                                  IntegerMatrix observations,
                                                  int no_groups,
                                                  IntegerMatrix group_index,
                                                  IntegerMatrix no_categories,
                                                  bool independent_thresholds,
                                                  IntegerMatrix indicator,
                                                  int no_persons,
                                                  NumericMatrix rest_matrix,
                                                  LogicalVector ordinal_variable,
                                                  IntegerVector reference_category,
                                                  NumericMatrix proposal_sd_pairwise,
                                                  double interaction_scale,
                                                  int no_variables,
                                                  double phi,
                                                  double target_ar,
                                                  int t,
                                                  double epsilon_lo,
                                                  double epsilon_hi) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  for(int variable1 = 0; variable1 <  no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  no_variables; variable2++) {
      if(indicator(variable1, variable2) == 1) {
        int int_index = pairwise_index(variable1, variable2);

        for(int h = 0; h < no_groups - 1; h++) {
          current_state = pairwise_effects(int_index, h + 1);
          proposed_state = R::rnorm(current_state, proposal_sd_pairwise(int_index, h + 1));

          log_prob = compare_anova_log_pseudolikelihood_ratio_pairwise_difference(main_effects,
                                                                                  main_index,
                                                                                  projection,
                                                                                  observations,
                                                                                  no_groups,
                                                                                  group_index,
                                                                                  no_categories,
                                                                                  independent_thresholds,
                                                                                  no_persons,
                                                                                  variable1,
                                                                                  variable2,
                                                                                  h,
                                                                                  proposed_state,
                                                                                  current_state,
                                                                                  rest_matrix,
                                                                                  ordinal_variable,
                                                                                  reference_category);

          log_prob += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
          log_prob -= R::dcauchy(current_state, 0.0, interaction_scale, true);

          U = R::unif_rand();
          if(std::log(U) < log_prob) {
            pairwise_effects(int_index, h + 1) = proposed_state;

            //Update the matrix of rest scores
            for(int gr = 0; gr < no_groups; gr++) {
              double state_difference = proposed_state - current_state;
              state_difference *= projection(gr, h);
              for(int person = group_index(gr, 0); person < group_index(gr, 1) + 1; person++) {
                rest_matrix(person, variable1) += observations(person, variable2) *
                  state_difference;
                rest_matrix(person, variable2) += observations(person, variable1) *
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
            proposal_sd_pairwise(int_index, h + 1) +
            (log_prob - target_ar) * std::exp(-log(t) * phi);

          if(std::isnan(update_proposal_sd) == true) {
            update_proposal_sd = 1.0;
          }

          update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

          proposal_sd_pairwise(int_index, h + 1) = update_proposal_sd;
        }
      }
    }
  }
}


// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the overall category
//  threshold parameters (nuisance) -- regular ordinal variable
// ----------------------------------------------------------------------------|
void compare_anova_metropolis_threshold_regular(NumericMatrix main_effects,
                                                IntegerMatrix main_index,
                                                NumericMatrix projection,
                                                IntegerMatrix observations,
                                                int no_groups,
                                                IntegerMatrix group_index,
                                                IntegerMatrix no_categories,
                                                int no_persons,
                                                NumericMatrix rest_matrix,
                                                List n_cat_obs,
                                                double threshold_alpha,
                                                double threshold_beta,
                                                int variable) {
  NumericVector q(no_persons);
  NumericVector r(no_persons);
  NumericVector GroupThresholds(no_categories(variable, 0));

  double log_prob, rest_score;
  double a, b, c;
  double tmp;
  double current_state, proposed_state;
  double U;
  double exp_current, exp_proposed;
  double exponent;

  for(int category = 0; category < no_categories(variable, 0); category++) {
    int cat_index = main_index(variable, 0);
    current_state = main_effects(cat_index + category, 0);

    exp_current = std::exp(current_state);
    c = (threshold_alpha + threshold_beta) / (1 + exp_current);

    for(int gr = 0; gr < no_groups; gr++) {

      // Group threshold parameters --------------------------------------------
      for(int cat = 0; cat < no_categories(variable, 0); cat++) {
        GroupThresholds[cat] = main_effects(cat_index + cat, 0);
        for(int h = 0; h < no_groups-1; h++) {
          GroupThresholds[cat] += projection(gr, h) *
            main_effects(cat_index + cat, h + 1);
        }
      }
      GroupThresholds[category] -= main_effects(cat_index + category, 0);       //Store only pairwise difference in element "category"

      for(int person = group_index(gr, 0); person < group_index(gr, 1) + 1; person++) {
        q[person] = 1.0;
        r[person] = 1.0;
        rest_score = rest_matrix(person, variable);
        for(int cat = 0; cat < no_categories(variable, 0); cat++) {
          if(cat != category) {
            exponent = GroupThresholds[cat];
            exponent += (cat + 1) * rest_score;
            q[person] += std::exp(exponent);
          }
        }
        exponent = (category + 1) * rest_score;
        exponent += GroupThresholds[category];
        r[person] = std::exp(exponent);
        c +=  r[person] / (q[person] + r[person] * exp_current);
      }
    }

    tmp = no_persons;
    tmp += threshold_alpha;
    tmp += threshold_beta;
    tmp -= exp_current * c;
    c = c / tmp;

    //Proposal is generalized beta-prime.
    a = threshold_alpha;
    b = no_persons;
    b += threshold_beta;
    for(int gr = 0; gr < no_groups; gr++) {
      IntegerMatrix n_cat_obs_gr = n_cat_obs[gr];
      a += n_cat_obs_gr(category + 1, variable);
      b -= n_cat_obs_gr(category + 1, variable);
    }

    tmp = R::rbeta(a, b);
    proposed_state = std::log(tmp / (1  - tmp) / c);
    exp_proposed = exp(proposed_state);

    //Compute log_acceptance probability for Metropolis.
    //First, we use q and r above to compute the ratio of pseudo-likelihoods
    log_prob = 0;
    for(int gr = 0; gr < no_groups; gr++) {
      for(int person = group_index(gr, 0); person < group_index(gr, 1) + 1; person++) {
        log_prob += std::log(q[person] + r[person] * exp_current);
        log_prob -= std::log(q[person] + r[person] * exp_proposed);
      }
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
      main_effects(cat_index + category, 0) = proposed_state;
    }
  }
}


// ----------------------------------------------------------------------------|
// The log pseudo-likelihood ratio [proposed against current] for a category
//  threshold difference for an independent samples design
// ----------------------------------------------------------------------------|
double compare_anova_log_pseudolikelihood_ratio_main_difference(NumericMatrix main_effects,
                                                                IntegerMatrix main_index,
                                                                NumericMatrix projection,
                                                                IntegerMatrix observations,
                                                                int no_groups,
                                                                IntegerMatrix group_index,
                                                                IntegerMatrix no_categories,
                                                                int no_persons,
                                                                NumericMatrix rest_matrix,
                                                                List n_cat_obs,
                                                                int variable,
                                                                int category,
                                                                int h,
                                                                double proposed_state,
                                                                double current_state) {

  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_proposed, denominator_current, exponent;
  int score;
  double delta_state = proposed_state - current_state;

  for(int gr = 0; gr < no_groups; gr++) {
    NumericVector proposed_thresholds(no_categories(variable, gr));
    NumericVector current_thresholds(no_categories(variable, gr));

    for(int cat = 0; cat < no_categories(variable, gr); cat++) {
      int cat_index = main_index(variable, 0) + cat;
      current_thresholds[cat] = main_effects(cat_index, 0);
      for(int hh = 0; hh < no_groups - 1; hh++) {
        current_thresholds[cat] += projection(gr, hh) * main_effects(cat_index, hh + 1);
      }
      proposed_thresholds[cat] = current_thresholds[cat];
    }

    proposed_thresholds[category] -= projection(gr, h) * current_state;
    proposed_thresholds[category] += projection(gr, h) * proposed_state;

    //Compute the pseudo-ikelihood ratio
    IntegerMatrix n_cat_obs_gr = n_cat_obs[gr];
    pseudolikelihood_ratio += delta_state *
      projection(gr, h) *
      n_cat_obs_gr(category + 1, variable);

    for(int person = group_index(gr, 0); person < group_index(gr, 1) + 1; person++) {
      rest_score = rest_matrix(person, variable);
      if(rest_score > 0) {
        bound = no_categories(variable, gr) * rest_score;
      } else {
        bound = 0.0;
      }

      denominator_proposed = std::exp(-bound);
      denominator_current = std::exp(-bound);
      for(int cat = 0; cat < no_categories(variable, gr); cat++) {
        score = cat + 1;
        exponent = score * rest_score - bound;
        denominator_proposed += std::exp(exponent + proposed_thresholds[cat]);
        denominator_current += std::exp(exponent + current_thresholds[cat]);
      }
      pseudolikelihood_ratio -= std::log(denominator_proposed);
      pseudolikelihood_ratio += std::log(denominator_current);
    }
  }

  return pseudolikelihood_ratio;
}


// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the category threshold
//  difference parameter -- regular ordinal variable
// ----------------------------------------------------------------------------|
void compare_anova_metropolis_main_difference_regular(NumericMatrix main_effects,
                                                      IntegerMatrix main_index,
                                                      NumericMatrix projection,
                                                      IntegerMatrix observations,
                                                      int no_groups,
                                                      IntegerMatrix group_index,
                                                      IntegerMatrix no_categories,
                                                      int no_persons,
                                                      NumericMatrix rest_matrix,
                                                      List n_cat_obs,
                                                      int variable,
                                                      IntegerMatrix indicator,
                                                      NumericMatrix proposal_sd_main,
                                                      double main_difference_scale,
                                                      double phi,
                                                      double target_ar,
                                                      int t,
                                                      double epsilon_lo,
                                                      double epsilon_hi) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;
  int cat_index;

  if(indicator(variable, variable) == 1) {
    for(int category = 0; category < no_categories(variable, 0); category++) {  //Should be the same across groups for this model choice.
      cat_index = main_index(variable, 0) + category;
      for(int h = 0; h < no_groups - 1; h++) {
        current_state = main_effects(cat_index, h + 1);
        proposed_state = R::rnorm(current_state,
                                  proposal_sd_main(cat_index, h));              //No sd needed for overall threshold, so start at zero.

        log_prob = compare_anova_log_pseudolikelihood_ratio_main_difference(main_effects,
                                                                            main_index,
                                                                            projection,
                                                                            observations,
                                                                            no_groups,
                                                                            group_index,
                                                                            no_categories,
                                                                            no_persons,
                                                                            rest_matrix,
                                                                            n_cat_obs,
                                                                            variable,
                                                                            category,
                                                                            h,
                                                                            proposed_state,
                                                                            current_state);

        log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
        log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);

        U = R::unif_rand();
        if(std::log(U) < log_prob) {
          main_effects(cat_index, h + 1) = proposed_state;
        }

        // Update the Robbins-Monro parameters
        if(log_prob > 0) {
          log_prob = 1;
        } else {
          log_prob = std::exp(log_prob);
        }

        double update_proposal_sd =
          proposal_sd_main(cat_index, h) +
          (log_prob - target_ar) * std::exp(-log(t) * phi);

        if(std::isnan(update_proposal_sd) == true) {
          update_proposal_sd = 1.0;
        }

        update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

        proposal_sd_main(cat_index, h) = update_proposal_sd;
      }
    }
  }
}


// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for the two
// category threshold parameters of the Blume-Capel model
// ----------------------------------------------------------------------------|
double compare_anova_log_pseudolikelihood_ratio_thresholds_blumecapel(double linear_current,
                                                                      double quadratic_current,
                                                                      double linear_proposed,
                                                                      double quadratic_proposed,
                                                                      int variable,
                                                                      IntegerVector reference_category,
                                                                      NumericMatrix main_effects,
                                                                      IntegerMatrix main_index,
                                                                      NumericMatrix projection,
                                                                      List sufficient_blume_capel,
                                                                      int no_persons,
                                                                      int no_groups,
                                                                      IntegerMatrix group_index,
                                                                      NumericMatrix rest_matrix,
                                                                      IntegerMatrix no_categories) {
  double lbound, bound;
  double rest_score, numerator, denominator, exponent;
  double pseudolikelihood_ratio = 0.0;
  int linear_score, quadratic_score;
  int cat_index = main_index(variable, 0);

  //----------------------------------------------------------------------------
  //Compute the log pseudo-likelihood ratio ------------------------------------
  //----------------------------------------------------------------------------
  for(int gr = 0; gr < no_groups; gr++) {
    NumericVector constant_numerator (no_categories(variable, gr) + 1);
    NumericVector constant_denominator (no_categories(variable, gr) + 1);

    //Pre-compute common terms
    for(int category = 0; category < no_categories(variable, gr) + 1; category++) {
      linear_score = category;
      quadratic_score =
        (category - reference_category[variable]) *
        (category - reference_category[variable]);

      constant_numerator[category] = linear_current * linear_score;
      constant_numerator[category] += quadratic_current * quadratic_score;
      constant_denominator[category] = linear_proposed * linear_score ;
      constant_denominator[category] += quadratic_proposed * quadratic_score;
      for(int h = 0; h < no_groups-1; h++) {
        constant_numerator[category] += projection(gr, h) *
          main_effects(cat_index, h + 1) *
          linear_score;
        constant_numerator[category] += projection(gr, h) *
          main_effects(cat_index + 1, h + 1) *
          quadratic_score;
        constant_denominator[category] += projection(gr, h) *
          main_effects(cat_index, h + 1) *
          linear_score;
        constant_denominator[category] += projection(gr, h) *
          main_effects(cat_index + 1, h + 1) *
          quadratic_score;
      }
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

    IntegerMatrix sufficient_blume_capel_gr = sufficient_blume_capel[gr];

    pseudolikelihood_ratio += sufficient_blume_capel_gr(0, variable) * linear_proposed;
    pseudolikelihood_ratio += sufficient_blume_capel_gr(1, variable) * quadratic_proposed;
    pseudolikelihood_ratio -= sufficient_blume_capel_gr(0, variable) * linear_current;
    pseudolikelihood_ratio -= sufficient_blume_capel_gr(1, variable) * quadratic_current;

    for(int person = group_index(gr, 0); person < group_index(gr, 1) + 1; person++) {
      rest_score = rest_matrix(person, variable);
      if(rest_score > 0) {
        bound = no_categories(variable, gr) * rest_score + lbound;
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
      pseudolikelihood_ratio += std::log(numerator);
      pseudolikelihood_ratio -= std::log(denominator);
    }
  }

  return pseudolikelihood_ratio;
}


// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the overall category
//  threshold parameters (nuisance) -- Blume-Capel ordinal variable
// ----------------------------------------------------------------------------|
void compare_anova_metropolis_threshold_blumecapel(NumericMatrix main_effects,
                                                   IntegerMatrix main_index,
                                                   NumericMatrix projection,
                                                   IntegerMatrix no_categories,
                                                   List sufficient_blume_capel,
                                                   int no_persons,
                                                   int no_groups,
                                                   IntegerMatrix group_index,
                                                   int variable,
                                                   IntegerVector reference_category,
                                                   double threshold_alpha,
                                                   double threshold_beta,
                                                   NumericMatrix rest_matrix,
                                                   NumericMatrix proposal_sd_main,
                                                   double phi,
                                                   double target_ar,
                                                   int t,
                                                   double epsilon_lo,
                                                   double epsilon_hi) {
  double log_prob, U;
  double current_state, proposed_state;

  //----------------------------------------------------------------------------
  // Adaptive Metropolis for the linear Blume-Capel parameter
  //----------------------------------------------------------------------------
  int cat_index = main_index(variable, 0);
  current_state = main_effects(cat_index, 0);
  proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index, 0));

  //----------------------------------------------------------------------------
  //Compute the log acceptance probability -------------------------------------
  //----------------------------------------------------------------------------
  double linear_current = current_state;
  double quadratic_current = main_effects(cat_index + 1, 0);
  double linear_proposed = proposed_state;
  double quadratic_proposed = main_effects(cat_index + 1, 0);

  log_prob = compare_anova_log_pseudolikelihood_ratio_thresholds_blumecapel(linear_current,
                                                                            quadratic_current,
                                                                            linear_proposed,
                                                                            quadratic_proposed,
                                                                            variable,
                                                                            reference_category,
                                                                            main_effects,
                                                                            main_index,
                                                                            projection,
                                                                            sufficient_blume_capel,
                                                                            no_persons,
                                                                            no_groups,
                                                                            group_index,
                                                                            rest_matrix,
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
    main_effects(cat_index, 0) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ------------------------------
  if(log_prob > 0) {
    log_prob = 1;
  } else {
    log_prob = std::exp(log_prob);
  }

  double update_proposal_sd = proposal_sd_main(cat_index, 0) +
    (log_prob - target_ar) * std::exp(-log(t) * phi);

  if(std::isnan(update_proposal_sd) == true) {
    update_proposal_sd = 1.0;
  }

  update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

  proposal_sd_main(cat_index, 0) = update_proposal_sd;

  //---------------------------------------------------------------------------|
  // Adaptive Metropolis for the quadratic Blume-Capel parameter
  //---------------------------------------------------------------------------|
  current_state = main_effects(cat_index + 1, 0);
  proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index + 1, 0));

  //----------------------------------------------------------------------------
  //Compute the log acceptance probability -------------------------------------
  //----------------------------------------------------------------------------

  linear_current = main_effects(cat_index, 0);
  quadratic_current = current_state;
  linear_proposed =  main_effects(cat_index, 0);
  quadratic_proposed =  proposed_state;

  log_prob = compare_anova_log_pseudolikelihood_ratio_thresholds_blumecapel(linear_current,
                                                                            quadratic_current,
                                                                            linear_proposed,
                                                                            quadratic_proposed,
                                                                            variable,
                                                                            reference_category,
                                                                            main_effects,
                                                                            main_index,
                                                                            projection,
                                                                            sufficient_blume_capel,
                                                                            no_persons,
                                                                            no_groups,
                                                                            group_index,
                                                                            rest_matrix,
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
    main_effects(cat_index + 1, 0) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ------------------------------
  if(log_prob > 0) {
    log_prob = 1;
  } else {
    log_prob = std::exp(log_prob);
  }

  update_proposal_sd = proposal_sd_main(cat_index + 1, 0) +
    (log_prob - target_ar) * std::exp(-log(t) * phi);

  if(std::isnan(update_proposal_sd) == true) {
    update_proposal_sd = 1.0;
  }

  update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

  proposal_sd_main(cat_index + 1, 0) = update_proposal_sd;
}


// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for the two
// category threshold differences of the Blume-Capel model
// ----------------------------------------------------------------------------|
double compare_anova_log_pseudolikelihood_ratio_main_difference_blumecapel(double linear_current,
                                                                           double quadratic_current,
                                                                           double linear_proposed,
                                                                           double quadratic_proposed,
                                                                           int variable,
                                                                           int h,
                                                                           IntegerVector reference_category,
                                                                           NumericMatrix main_effects,
                                                                           IntegerMatrix main_index,
                                                                           NumericMatrix projection,
                                                                           List sufficient_blume_capel,
                                                                           int no_persons,
                                                                           int no_groups,
                                                                           IntegerMatrix group_index,
                                                                           NumericMatrix rest_matrix,
                                                                           IntegerMatrix no_categories) {
  double lbound, bound;
  double rest_score, numerator, denominator, exponent;
  double pseudolikelihood_ratio = 0.0;
  int linear_score, quadratic_score;
  int cat_index = main_index(variable, 0);

  //----------------------------------------------------------------------------
  //Compute the log pseudo-likelihood ratio ------------------------------------
  //----------------------------------------------------------------------------
  for(int gr = 0; gr < no_groups; gr++) {
    NumericVector constant_numerator (no_categories(variable, gr) + 1);
    NumericVector constant_denominator (no_categories(variable, gr) + 1);
    double P = projection(gr, h);

    //Pre-compute common terms
    for(int category = 0; category < no_categories(variable, gr) + 1; category++) {
      linear_score = category;
      quadratic_score =
        (category - reference_category[variable]) *
        (category - reference_category[variable]);

      constant_numerator[category] = main_effects(cat_index, 0) * linear_score;
      constant_numerator[category] +=  main_effects(cat_index + 1, 0) * quadratic_score;
      constant_denominator[category] =  main_effects(cat_index, 0) * linear_score ;
      constant_denominator[category] += main_effects(cat_index + 1, 0) * quadratic_score;
      for(int hh = 0; hh < no_groups-1; hh++) {
        constant_numerator[category] += projection(gr, hh) *
          main_effects(cat_index, hh + 1) *
          linear_score;
        constant_numerator[category] += projection(gr, hh) *
          main_effects(cat_index + 1, hh + 1) *
          quadratic_score;
        constant_denominator[category] += projection(gr, hh) *
          main_effects(cat_index, hh + 1) *
          linear_score;
        constant_denominator[category] += projection(gr, hh) *
          main_effects(cat_index + 1, hh + 1) *
          quadratic_score;
      }

      constant_denominator[category] -= P *
        main_effects(cat_index, h + 1) *
        linear_score;
      constant_denominator[category] -= P *
        main_effects(cat_index + 1, h + 1) *
        quadratic_score;
      constant_denominator[category] += P *
        linear_proposed *
        linear_score;
      constant_denominator[category] += P *
        quadratic_proposed *
        quadratic_score;
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

    IntegerMatrix sufficient_blume_capel_gr = sufficient_blume_capel[gr];

    pseudolikelihood_ratio += sufficient_blume_capel_gr(0, variable) * P * linear_proposed;
    pseudolikelihood_ratio += sufficient_blume_capel_gr(1, variable) * P * quadratic_proposed;
    pseudolikelihood_ratio -= sufficient_blume_capel_gr(0, variable) * P * linear_current;
    pseudolikelihood_ratio -= sufficient_blume_capel_gr(1, variable) * P * quadratic_current;

    for(int person = group_index(gr, 0); person < group_index(gr, 1) + 1; person++) {
      rest_score = rest_matrix(person, variable);
      if(rest_score > 0) {
        bound = no_categories(variable, gr) * rest_score + lbound;
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
      pseudolikelihood_ratio += std::log(numerator);
      pseudolikelihood_ratio -= std::log(denominator);
    }
  }

  return pseudolikelihood_ratio;
}


// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the category threshold
//  difference parameters-- Blume-Capel ordinal variable
// ----------------------------------------------------------------------------|
void compare_anova_metropolis_main_difference_blumecapel(NumericMatrix main_effects,
                                                         IntegerMatrix main_index,
                                                         NumericMatrix projection,
                                                         IntegerMatrix no_categories,
                                                         List sufficient_blume_capel,
                                                         int no_persons,
                                                         int no_groups,
                                                         IntegerMatrix group_index,
                                                         int variable,
                                                         IntegerVector reference_category,
                                                         double threshold_alpha,
                                                         double threshold_beta,
                                                         NumericMatrix rest_matrix,
                                                         IntegerMatrix indicator,
                                                         NumericMatrix proposal_sd_main,
                                                         double phi,
                                                         double target_ar,
                                                         int t,
                                                         double epsilon_lo,
                                                         double epsilon_hi) {
  double log_prob, U;
  double current_state, proposed_state;

  if(indicator(variable, variable) == 1) {
    //Loop over the group difference effects--------------------------------------
    for(int h = 0; h < no_groups - 1; h++) {
      //----------------------------------------------------------------------------
      // Adaptive Metropolis for the linear Blume-Capel parameter
      //----------------------------------------------------------------------------
      int cat_index = main_index(variable, 0);
      current_state = main_effects(cat_index, h + 1);
      proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index, h + 1));

      //----------------------------------------------------------------------------
      //Compute the log acceptance probability -------------------------------------
      //----------------------------------------------------------------------------
      double linear_current = current_state;
      double quadratic_current = main_effects(cat_index + 1, h + 1);
      double linear_proposed = proposed_state;
      double quadratic_proposed = main_effects(cat_index + 1, h + 1);

      log_prob = compare_anova_log_pseudolikelihood_ratio_main_difference_blumecapel(linear_current,
                                                                                     quadratic_current,
                                                                                     linear_proposed,
                                                                                     quadratic_proposed,
                                                                                     variable,
                                                                                     h,
                                                                                     reference_category,
                                                                                     main_effects,
                                                                                     main_index,
                                                                                     projection,
                                                                                     sufficient_blume_capel,
                                                                                     no_persons,
                                                                                     no_groups,
                                                                                     group_index,
                                                                                     rest_matrix,
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
        main_effects(cat_index, h + 1) = proposed_state;
      }

      //Robbins-Monro update of the proposal variance ------------------------------
      if(log_prob > 0) {
        log_prob = 1;
      } else {
        log_prob = std::exp(log_prob);
      }

      double update_proposal_sd = proposal_sd_main(cat_index, h + 1) +
        (log_prob - target_ar) * std::exp(-log(t) * phi);

      if(std::isnan(update_proposal_sd) == true) {
        update_proposal_sd = 1.0;
      }

      update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

      proposal_sd_main(cat_index, h + 1) = update_proposal_sd;

      //---------------------------------------------------------------------------|
      // Adaptive Metropolis for the quadratic Blume-Capel parameter
      //---------------------------------------------------------------------------|
      current_state = main_effects(cat_index + 1, h + 1);
      proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index + 1, h + 1));

      //----------------------------------------------------------------------------
      //Compute the log acceptance probability -------------------------------------
      //----------------------------------------------------------------------------

      linear_current = main_effects(cat_index, h + 1);
      quadratic_current = current_state;
      linear_proposed =  main_effects(cat_index, h + 1);
      quadratic_proposed =  proposed_state;

      log_prob = compare_anova_log_pseudolikelihood_ratio_main_difference_blumecapel(linear_current,
                                                                                     quadratic_current,
                                                                                     linear_proposed,
                                                                                     quadratic_proposed,
                                                                                     variable,
                                                                                     h,
                                                                                     reference_category,
                                                                                     main_effects,
                                                                                     main_index,
                                                                                     projection,
                                                                                     sufficient_blume_capel,
                                                                                     no_persons,
                                                                                     no_groups,
                                                                                     group_index,
                                                                                     rest_matrix,
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
        main_effects(cat_index + 1, h + 1) = proposed_state;
      }

      //Robbins-Monro update of the proposal variance ------------------------------
      if(log_prob > 0) {
        log_prob = 1;
      } else {
        log_prob = std::exp(log_prob);
      }

      update_proposal_sd = proposal_sd_main(cat_index + 1, h + 1) +
        (log_prob - target_ar) * std::exp(-log(t) * phi);

      if(std::isnan(update_proposal_sd) == true) {
        update_proposal_sd = 1.0;
      }

      update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);

      proposal_sd_main(cat_index + 1, h + 1) = update_proposal_sd;
    }
  }
}


// ----------------------------------------------------------------------------|
// The log pseudolikelihood ratio [proposed against current] for the difference
//  in a pairwise interaction for an independent samples design
// ----------------------------------------------------------------------------|
double compare_anova_log_pseudolikelihood_ratio_pairwise_difference_between_model(NumericMatrix main_effects,
                                                                                  IntegerMatrix main_index,
                                                                                  NumericMatrix projection,
                                                                                  IntegerMatrix observations,
                                                                                  int no_groups,
                                                                                  IntegerMatrix group_index,
                                                                                  IntegerMatrix no_categories,
                                                                                  bool independent_thresholds,
                                                                                  int no_persons,
                                                                                  int variable1,
                                                                                  int variable2,
                                                                                  NumericVector proposed_states,
                                                                                  NumericVector current_states,
                                                                                  NumericMatrix rest_matrix,
                                                                                  LogicalVector ordinal_variable,
                                                                                  IntegerVector reference_category) {
  double rest_score, bound;
  double pseudolikelihood_ratio = 0.0;
  double denominator_prop, denominator_curr, exponent;
  int score, obs_score1, obs_score2;

  // --------------------------------------------------------------------------|
  // Compute log pseudo-likelihood ratio ---------------------------------------
  // --------------------------------------------------------------------------|
  for(int gr = 0; gr < no_groups; gr++) {

    NumericVector GroupThresholds_v1 = group_thresholds_for_variable(variable1,
                                                                     ordinal_variable,
                                                                     gr,
                                                                     no_groups,
                                                                     no_categories,
                                                                     main_effects,
                                                                     main_index,
                                                                     projection,
                                                                     independent_thresholds);

    NumericVector GroupThresholds_v2 = group_thresholds_for_variable(variable2,
                                                                     ordinal_variable,
                                                                     gr,
                                                                     no_groups,
                                                                     no_categories,
                                                                     main_effects,
                                                                     main_index,
                                                                     projection,
                                                                     independent_thresholds);

    double proposed_pairwise_group_effect = 0.0;
    double current_pairwise_group_effect = 0.0;
    for(int h = 0; h < no_groups - 1; h++) {
      double P = projection(gr, h);
      proposed_pairwise_group_effect += P * proposed_states[h];
      current_pairwise_group_effect += P * current_states[h];
    }

    //Loop over the pseudo-likelihoods for persons in group gr -----------------
    for(int person = group_index(gr, 0); person < group_index(gr, 1) + 1; person++) {
      obs_score1 = observations(person, variable1);
      obs_score2 = observations(person, variable2);

      pseudolikelihood_ratio += 2 * obs_score1 * obs_score2 * proposed_pairwise_group_effect;
      pseudolikelihood_ratio -= 2 * obs_score1 * obs_score2 * current_pairwise_group_effect;

      //Log pseudo-likelihood for variable 1 -----------------------------------
      rest_score = rest_matrix(person, variable1) - obs_score2 * current_pairwise_group_effect;

      if(rest_score > 0) {
        bound = no_categories(variable1, gr) * rest_score;
      } else {
        bound = 0.0;
      }

      if(ordinal_variable[variable1] == true) {
        //Regular binary or ordinal MRF variable -------------------------------
        denominator_prop = std::exp(-bound);
        denominator_curr = std::exp(-bound);
        for(int cat = 0; cat < no_categories(variable1, gr); cat++) {
          score = cat + 1;
          exponent = GroupThresholds_v1[cat] +
            score * rest_score -
            bound;
          denominator_prop +=
            std::exp(exponent + score * obs_score2 * proposed_pairwise_group_effect);
          denominator_curr +=
            std::exp(exponent + score * obs_score2 * current_pairwise_group_effect);
        }
      } else {
        //Blume-Capel ordinal MRF variable -------------------------------------
        denominator_prop = 0.0;
        denominator_curr = 0.0;
        for(int cat = 0; cat < no_categories(variable1, gr) + 1; cat++) {
          exponent = GroupThresholds_v1[0] * cat;
          exponent += GroupThresholds_v1[1] *
            (cat - reference_category[variable1]) *
            (cat - reference_category[variable1]);
          exponent += cat * rest_score - bound;
          denominator_prop +=
            std::exp(exponent + cat * obs_score2 * proposed_pairwise_group_effect);
          denominator_curr +=
            std::exp(exponent + cat * obs_score2 * current_pairwise_group_effect);
        }
      }
      pseudolikelihood_ratio -= std::log(denominator_prop);
      pseudolikelihood_ratio += std::log(denominator_curr);

      //Log pseudo-likelihood for variable 2 -----------------------------------
      rest_score = rest_matrix(person, variable2) - obs_score1 * current_pairwise_group_effect;

      if(rest_score > 0) {
        bound = no_categories(variable2, gr) * rest_score;
      } else {
        bound = 0.0;
      }

      if(ordinal_variable[variable2] == true) {
        //Regular binary or ordinal MRF variable ---------------------------------
        denominator_prop = std::exp(-bound);
        denominator_curr = std::exp(-bound);
        for(int cat = 0; cat < no_categories(variable2, gr); cat++) {
          score = cat + 1;
          exponent = GroupThresholds_v2[cat] +
            score * rest_score -
            bound;
          denominator_prop +=
            std::exp(exponent + score * obs_score1 * proposed_pairwise_group_effect);
          denominator_curr +=
            std::exp(exponent + score * obs_score1 * current_pairwise_group_effect);
        }
      } else {
        //Blume-Capel ordinal MRF variable ---------------------------------------
        denominator_prop = 0.0;
        denominator_curr = 0.0;
        for(int cat = 0; cat < no_categories(variable2, gr) + 1; cat++) {
          exponent = GroupThresholds_v2[0] * cat;
          exponent += GroupThresholds_v2[1] *
            (cat - reference_category[variable2]) *
            (cat - reference_category[variable2]);
          exponent +=  cat * rest_score - bound;
          denominator_prop +=
            std::exp(exponent + cat * obs_score1 * proposed_pairwise_group_effect);
          denominator_curr +=
            std::exp(exponent + cat * obs_score1 * current_pairwise_group_effect);
        }
      }
      pseudolikelihood_ratio -= std::log(denominator_prop);
      pseudolikelihood_ratio += std::log(denominator_curr);
    }
  }
  return pseudolikelihood_ratio;
}



// ----------------------------------------------------------------------------|
// MH algorithm to sample from the full-conditional of the threshold parameters
//   for a regular binary or ordinal variable
// ----------------------------------------------------------------------------|
void compare_anova_metropolis_thresholds_regular_free(NumericMatrix main_effects,
                                                      IntegerMatrix main_index,
                                                      IntegerMatrix observations,
                                                      int no_groups,
                                                      IntegerMatrix group_index,
                                                      IntegerMatrix no_categories,
                                                      NumericMatrix rest_matrix,
                                                      List n_cat_obs,
                                                      double threshold_alpha,
                                                      double threshold_beta,
                                                      int variable,
                                                      int group) {

  int no_persons = 1 + group_index(group, 1) - group_index(group, 0);

  NumericVector q(no_persons);
  NumericVector r(no_persons);

  double log_prob, rest_score;
  double a, b, c;
  double tmp;
  double current_state, proposed_state;
  double U;
  double exp_current, exp_proposed;
  IntegerMatrix n_cat_obs_gr = n_cat_obs[group];


  for(int category = 0; category < no_categories(variable, group); category++) {
    int cat_index = main_index(variable, 0);
    current_state = main_effects(cat_index + category, group);
    exp_current = std::exp(current_state);
    c = (threshold_alpha + threshold_beta) / (1 + exp_current);

    for(int person = 0; person < no_persons; person++) {
      q[person] = 1.0;
      r[person] = 1.0;
      rest_score = rest_matrix(group_index(group, 0) + person, variable);
      for(int cat = 0; cat < no_categories[variable]; cat++) {
        if(cat != category) {
          q[person] += std::exp(main_effects(cat_index + cat, group) +
            (cat + 1) * rest_score);
        }
      }
      r[person] = std::exp((category + 1) * rest_score);
      c +=  r[person] / (q[person] + r[person] * exp_current);
    }

    c = c / ((no_persons + threshold_alpha + threshold_beta) -
      exp_current * c);

    //Proposal is generalized beta-prime.
    a = n_cat_obs_gr(category + 1, variable) + threshold_alpha;
    b = no_persons + threshold_beta - n_cat_obs_gr(category + 1, variable);
    tmp = R::rbeta(a, b);
    proposed_state = std::log(tmp / (1  - tmp) / c);
    exp_proposed = exp(proposed_state);

    //Compute log_acceptance probability for Metropolis.
    //First, we use g and q above to compute the ratio of pseudolikelihoods
    log_prob = 0;
    for(int person = 0; person < no_persons; person++) {
      log_prob += std::log(q[person] + r[person] * exp_current);
      log_prob -= std::log(q[person] + r[person] * exp_proposed);
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
      main_effects(cat_index + category, group) = proposed_state;
    }
  }
}

// ----------------------------------------------------------------------------|
// Adaptive Metropolis algorithm to sample from the full-conditional of the
//   threshold parameters for a Blume-Capel ordinal variable
// ----------------------------------------------------------------------------|
void compare_anova_metropolis_thresholds_blumecapel_free(NumericMatrix main_effects,
                                                         IntegerMatrix main_index,
                                                         IntegerMatrix observations,
                                                         int no_groups,
                                                         IntegerMatrix group_index,
                                                         IntegerVector reference_category,
                                                         IntegerMatrix no_categories,
                                                         List sufficient_blume_capel,
                                                         int no_persons,
                                                         NumericMatrix rest_matrix,
                                                         List n_cat_obs,
                                                         double threshold_alpha,
                                                         double threshold_beta,
                                                         int variable,
                                                         int group,
                                                         NumericMatrix proposal_sd_main,
                                                         double phi,
                                                         double target_ar,
                                                         int t,
                                                         double epsilon_lo,
                                                         double epsilon_hi) {

  double log_prob, U;
  double current_state, proposed_state, difference;
  double numerator, denominator;
  double lbound, bound, exponent, rest_score;

  NumericMatrix sufficient_blume_capel_group = sufficient_blume_capel[group];

  NumericVector constant_numerator (no_categories(variable, group) + 1);
  NumericVector constant_denominator (no_categories(variable, group) + 1);

  //----------------------------------------------------------------------------
  //Adaptive Metropolis for the linear Blume-Capel parameter
  //----------------------------------------------------------------------------
  int cat_index = main_index(variable, 0);
  current_state = main_effects(cat_index, group);
  proposed_state = R::rnorm(current_state,
                            proposal_sd_main(cat_index, group));

  //Precompute terms for the log acceptance probability ------------------------
  difference = proposed_state - current_state;

  for(int category = 0; category < no_categories(variable, group) + 1; category ++) {
    exponent = main_effects(cat_index + 1, group) *
      (category - reference_category[variable]) *
      (category - reference_category[variable]);
    constant_numerator[category] = current_state * category + exponent;
    constant_denominator[category] = proposed_state * category + exponent;
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

  //Compute the log acceptance probability -------------------------------------
  log_prob = threshold_alpha * difference;
  log_prob += sufficient_blume_capel_group(0, variable) * difference;

  for(int person = group_index(group, 0); person < group_index(group, 1) + 1; person++) {
    rest_score = rest_matrix(person, variable);
    if(rest_score > 0) {
      bound = no_categories(variable, group) * rest_score + lbound;
    } else {
      bound = lbound;
    }
    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int category = 0; category < no_categories(variable, group); category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }
    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }

  log_prob += (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(current_state));
  log_prob -= (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(proposed_state));

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    main_effects(cat_index, group) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ------------------------------
  if(log_prob > 0) {
    log_prob = 1;
  } else {
    log_prob = std::exp(log_prob);
  }

  double update_proposal_sd = proposal_sd_main(cat_index, group) +
    (log_prob - target_ar) * std::exp(-log(t) * phi);

  if(std::isnan(update_proposal_sd) == true) {
    update_proposal_sd = 1.0;
  }

  update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);
  proposal_sd_main(cat_index, group) = update_proposal_sd;

  //----------------------------------------------------------------------------
  //Adaptive Metropolis for the quadratic Blume-Capel parameter
  //----------------------------------------------------------------------------
  current_state = main_effects(cat_index + 1, group);
  proposed_state = R::rnorm(current_state,
                            proposal_sd_main(cat_index + 1, group));

  //Precompute terms for the log acceptance probability ------------------------
  difference = proposed_state - current_state;

  for(int category = 0; category < no_categories(variable, group) + 1; category ++) {
    exponent = main_effects(cat_index, group) * category;
    int score = (category - reference_category[variable]) *
      (category - reference_category[variable]);

    constant_numerator[category] = current_state * score + exponent;
    constant_denominator[category] = proposed_state * score + exponent;
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

  //Compute the log acceptance probability -------------------------------------
  log_prob = threshold_alpha * difference;
  log_prob += sufficient_blume_capel_group(1, variable) * difference;

  for(int person = group_index(group, 0); person < group_index(group, 1) + 1; person++) {
    rest_score = rest_matrix(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);

    for(int category = 0; category < no_categories(variable, group); category ++) {
      exponent = (category + 1) * rest_score - bound;
      numerator += std::exp(constant_numerator[category + 1] + exponent);
      denominator += std::exp(constant_denominator[category + 1] + exponent);
    }

    log_prob += std::log(numerator);
    log_prob -= std::log(denominator);
  }
  log_prob += (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(current_state));
  log_prob -= (threshold_alpha + threshold_beta) *
    std::log(1 + std::exp(proposed_state));

  U = R::unif_rand();
  if(std::log(U) < log_prob) {
    main_effects(cat_index + 1, group) = proposed_state;
  }

  //Robbins-Monro update of the proposal variance ------------------------------
  if(log_prob > 0) {
    log_prob = 1;
  } else {
    log_prob = std::exp(log_prob);
  }

  update_proposal_sd = proposal_sd_main(cat_index + 1, group) +
    (log_prob - target_ar) * std::exp(-log(t) * phi);

  if(std::isnan(update_proposal_sd) == true) {
    update_proposal_sd = 1.0;
  }

  update_proposal_sd = std::clamp(update_proposal_sd, epsilon_lo, epsilon_hi);
  proposal_sd_main(cat_index + 1, group) = update_proposal_sd;

}


// ----------------------------------------------------------------------------|
// A Gibbs step for graphical model parameters for Bayesian parameter comparison
// ----------------------------------------------------------------------------|
List compare_anova_gibbs_step_gm(NumericMatrix main_effects,
                                 IntegerMatrix main_index,
                                 NumericMatrix pairwise_effects,
                                 IntegerMatrix pairwise_index,
                                 NumericMatrix projection,
                                 IntegerMatrix no_categories,
                                 IntegerMatrix observations,
                                 int no_persons,
                                 int no_groups,
                                 IntegerMatrix group_index,
                                 List n_cat_obs,
                                 List sufficient_blume_capel,
                                 NumericMatrix rest_matrix,
                                 bool independent_thresholds,
                                 LogicalVector ordinal_variable,
                                 IntegerVector reference_category,
                                 IntegerMatrix indicator,
                                 NumericMatrix inclusion_probability_difference,
                                 NumericMatrix proposal_sd_main,
                                 NumericMatrix proposal_sd_pairwise,
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
                                 bool difference_selection,
                                 int no_pairwise,
                                 int no_variables,
                                 IntegerMatrix index) {

  NumericMatrix proposal_sd_blumecapel_group(no_variables, 2);


  //Within model move for the pairwise interaction parameters
  compare_anova_metropolis_interaction(main_effects,
                                       pairwise_effects,
                                       main_index,
                                       pairwise_index,
                                       projection,
                                       observations,
                                       no_groups,
                                       group_index,
                                       no_categories,
                                       independent_thresholds,
                                       no_persons,
                                       rest_matrix,
                                       ordinal_variable,
                                       reference_category,
                                       proposal_sd_pairwise,
                                       interaction_scale,
                                       no_variables,
                                       phi,
                                       target_ar,
                                       t,
                                       epsilon_lo,
                                       epsilon_hi);

  //Between model move for differences in interaction parameters
  // if(difference_selection == true) {
  //
  // }

  //Within model move for differences in interaction parameters
  compare_anova_metropolis_pairwise_difference(main_effects,
                                               pairwise_effects,
                                               main_index,
                                               pairwise_index,
                                               projection,
                                               observations,
                                               no_groups,
                                               group_index,
                                               no_categories,
                                               independent_thresholds,
                                               indicator,
                                               no_persons,
                                               rest_matrix,
                                               ordinal_variable,
                                               reference_category,
                                               proposal_sd_pairwise,
                                               interaction_scale,
                                               no_variables,
                                               phi,
                                               target_ar,
                                               t,
                                               epsilon_lo,
                                               epsilon_hi);

  //Update threshold parameters
  if(independent_thresholds == false) {
    for(int variable = 0; variable < no_variables; variable++) {
      if(ordinal_variable[variable] == true) {
        //Within model move for the main category thresholds
        compare_anova_metropolis_threshold_regular(main_effects,
                                                   main_index,
                                                   projection,
                                                   observations,
                                                   no_groups,
                                                   group_index,
                                                   no_categories,
                                                   no_persons,
                                                   rest_matrix,
                                                   n_cat_obs,
                                                   threshold_alpha,
                                                   threshold_beta,
                                                   variable);

        // if(difference_selection == true) {
        //   //Between model move for the differences in category thresholds
        // }

        //Within model move for the differences in category thresholds
        compare_anova_metropolis_main_difference_regular(main_effects,
                                                         main_index,
                                                         projection,
                                                         observations,
                                                         no_groups,
                                                         group_index,
                                                         no_categories,
                                                         no_persons,
                                                         rest_matrix,
                                                         n_cat_obs,
                                                         variable,
                                                         indicator,
                                                         proposal_sd_main,
                                                         main_difference_scale,
                                                         phi,
                                                         target_ar,
                                                         t,
                                                         epsilon_lo,
                                                         epsilon_hi);
      } else {
        //Within model move for the main category thresholds
        compare_anova_metropolis_threshold_blumecapel(main_effects,
                                                      main_index,
                                                      projection,
                                                      no_categories,
                                                      sufficient_blume_capel,
                                                      no_persons,
                                                      no_groups,
                                                      group_index,
                                                      variable,
                                                      reference_category,
                                                      threshold_alpha,
                                                      threshold_beta,
                                                      rest_matrix,
                                                      proposal_sd_main,
                                                      phi,
                                                      target_ar,
                                                      t,
                                                      epsilon_lo,
                                                      epsilon_hi);

        //Between model move for the differences in category thresholds
        // if(difference_selection == true) {
        //
        // }

        //Within model move for the differences in category thresholds
        compare_anova_metropolis_main_difference_blumecapel(main_effects,
                                                            main_index,
                                                            projection,
                                                            no_categories,
                                                            sufficient_blume_capel,
                                                            no_persons,
                                                            no_groups,
                                                            group_index,
                                                            variable,
                                                            reference_category,
                                                            threshold_alpha,
                                                            threshold_beta,
                                                            rest_matrix,
                                                            indicator,
                                                            proposal_sd_main,
                                                            phi,
                                                            target_ar,
                                                            t,
                                                            epsilon_lo,
                                                            epsilon_hi);
      }
    }
  } else {
    for(int variable = 0; variable < no_variables; variable++) {
      for(int group = 0; group < no_groups; group++) {
        if(ordinal_variable[variable] == true) {
          compare_anova_metropolis_thresholds_regular_free(main_effects,
                                                           main_index,
                                                           observations,
                                                           no_groups,
                                                           group_index,
                                                           no_categories,
                                                           rest_matrix,
                                                           n_cat_obs,
                                                           threshold_alpha,
                                                           threshold_beta,
                                                           variable,
                                                           group);
        } else {
          compare_anova_metropolis_thresholds_blumecapel_free(main_effects,
                                                              main_index,
                                                              observations,
                                                              no_groups,
                                                              group_index,
                                                              reference_category,
                                                              no_categories,
                                                              sufficient_blume_capel,
                                                              no_persons,
                                                              rest_matrix,
                                                              n_cat_obs,
                                                              threshold_alpha,
                                                              threshold_beta,
                                                              variable,
                                                              group,
                                                              proposal_sd_main,
                                                              phi,
                                                              target_ar,
                                                              t,
                                                              epsilon_lo,
                                                              epsilon_hi);
        }
      }
    }
  }

  return List::create(Named("indicator") = indicator,
                      Named("main_effects") = main_effects,
                      Named("pairwise_effects") = pairwise_effects,
                      Named("rest_matrix") = rest_matrix,
                      Named("proposal_sd_main") = proposal_sd_main,
                      Named("proposal_sd_pairwise") = proposal_sd_pairwise);
}


// ----------------------------------------------------------------------------|
// The Gibbs sampler for Bayesian parameter comparisons
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
List compare_anova_gibbs_sampler(IntegerMatrix observations,
                                 IntegerMatrix main_index,
                                 IntegerMatrix pairwise_index,
                                 NumericMatrix projection,
                                 IntegerMatrix no_categories,
                                 int no_groups,
                                 IntegerMatrix group_index,
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
                                 List n_cat_obs,
                                 List sufficient_blume_capel,
                                 double threshold_alpha,
                                 double threshold_beta,
                                 bool na_impute,
                                 IntegerMatrix missing_index,
                                 LogicalVector ordinal_variable,
                                 IntegerVector reference_category,
                                 bool independent_thresholds,
                                 bool save = false,
                                 bool display_progress = false,
                                 bool difference_selection = true) {

  int no_variables = observations.ncol();
  int no_persons = observations.nrow();
  int no_main = 1 + main_index(no_variables - 1, 1);
  int no_pairwise = no_variables * (no_variables - 1) / 2;

  IntegerVector person_group_indicator (no_persons);
  for(int group = 0; group < no_groups; group++) {
    for(int person = group_index(group, 0); person < group_index(group, 1) + 1; person++) {
      person_group_indicator[person] = group;
    }
  }

  // Matrices with model parameters --------------------------------------------
  NumericMatrix main_effects(no_main, no_groups);
  NumericMatrix pairwise_effects(no_pairwise, no_groups);
  IntegerMatrix indicator(no_variables, no_variables);
  std::fill(indicator.begin(), indicator.end(), 1.0);

  // Matrices with standard deviations for adaptive Metropolis -----------------
  NumericMatrix proposal_sd_main(no_main, no_groups);
  NumericMatrix proposal_sd_pairwise(no_pairwise, no_groups);

  std::fill(proposal_sd_main.begin(), proposal_sd_main.end(), 1.0);
  std::fill(proposal_sd_pairwise.begin(), proposal_sd_pairwise.end(), 1.0);

  //Parameters for the Robbins-Monro approach for adaptive Metropolis ----------
  double phi = 0.75;
  double target_ar = 0.234;
  double epsilon_lo = 1.0 / static_cast<double>(no_persons);
  double epsilon_hi = 2.0;

  //Randomized index for the pairwise updates ----------------------------------
  IntegerVector v = seq(0, no_pairwise - 1);
  IntegerVector order(no_pairwise);
  IntegerMatrix index(no_pairwise, 3);

  NumericMatrix out_main(no_main, no_groups);
  NumericMatrix out_pairwise(no_pairwise, no_groups);
  NumericMatrix out_indicator(no_variables, no_variables);

  //These matrices will contain the rest scores in the pseudo-likelihoods ------
  NumericMatrix rest_matrix(no_persons, no_variables);

  //The Gibbs sampler ----------------------------------------------------------
  //First, we do burn-in iterations---------------------------------------------

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
      return List::create(Named("main") = out_main,
                          Named("pairwise") = out_pairwise,
                          Named("indicator") = out_indicator);

    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_pairwise,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_pairwise; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    if(na_impute == true) {
      List out = compare_anova_impute_missing_data(main_effects,
                                                   pairwise_effects,
                                                   main_index,
                                                   pairwise_index,
                                                   projection,
                                                   observations,
                                                   no_groups,
                                                   person_group_indicator,
                                                   n_cat_obs,
                                                   sufficient_blume_capel,
                                                   no_categories,
                                                   rest_matrix,
                                                   missing_index,
                                                   ordinal_variable,
                                                   reference_category,
                                                   independent_thresholds);

      IntegerMatrix observations = out["observations"];
      List n_cat_obs = out["n_cat_obs"];
      List sufficient_blume_capel = out["sufficient_blume_capel"];
      NumericMatrix rest_matrix = out["rest_matrix"];
    }

    List out = compare_anova_gibbs_step_gm(main_effects,
                                           main_index,
                                           pairwise_effects,
                                           pairwise_index,
                                           projection,
                                           no_categories,
                                           observations,
                                           no_persons,
                                           no_groups,
                                           group_index,
                                           n_cat_obs,
                                           sufficient_blume_capel,
                                           rest_matrix,
                                           independent_thresholds,
                                           ordinal_variable,
                                           reference_category,
                                           indicator,
                                           inclusion_probability_difference,
                                           proposal_sd_main,
                                           proposal_sd_pairwise,
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
                                           difference_selection,
                                           no_pairwise,
                                           no_variables,
                                           index);

    IntegerMatrix indicator = out["indicator"];
    NumericMatrix main_effects = out["main_effects"];
    NumericMatrix pairwise_effects = out["pairwise_effects"];
    NumericMatrix rest_matrix = out["rest_matrix"];
    NumericMatrix proposal_sd_pairwise = out["proposal_sd_pairwise"];
    NumericMatrix proposal_sd_main = out["proposal_sd_main"];

    if(difference_selection == true) {
      if(pairwise_difference_prior == "Beta-Bernoulli") {
        int sumG = 0;
        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            sumG += indicator(i, j);
          }
        }
        double probability = R::rbeta(pairwise_beta_bernoulli_alpha + sumG,
                                      pairwise_beta_bernoulli_beta + no_pairwise - sumG);

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
      return List::create(Named("main") = out_main,
                          Named("pairwise") = out_pairwise,
                          Named("indicator") = out_indicator);
    }
    p.increment();

    //Update interactions and model (between model move) -----------------------
    //Use a random order to update the edge - interaction pairs ----------------
    order = sample(v,
                   no_pairwise,
                   false,
                   R_NilValue);

    for(int cntr = 0; cntr < no_pairwise; cntr++) {
      index(cntr, 0) = Index(order[cntr], 0);
      index(cntr, 1) = Index(order[cntr], 1);
      index(cntr, 2) = Index(order[cntr], 2);
    }

    if(na_impute == true) {
      List out = compare_anova_impute_missing_data(main_effects,
                                                   pairwise_effects,
                                                   main_index,
                                                   pairwise_index,
                                                   projection,
                                                   observations,
                                                   no_groups,
                                                   person_group_indicator,
                                                   n_cat_obs,
                                                   sufficient_blume_capel,
                                                   no_categories,
                                                   rest_matrix,
                                                   missing_index,
                                                   ordinal_variable,
                                                   reference_category,
                                                   independent_thresholds);

      IntegerMatrix observations = out["observations"];
      List n_cat_obs = out["n_cat_obs"];
      List sufficient_blume_capel = out["sufficient_blume_capel"];
      NumericMatrix rest_matrix = out["rest_matrix"];
    }

    List out = compare_anova_gibbs_step_gm(main_effects,
                                           main_index,
                                           pairwise_effects,
                                           pairwise_index,
                                           projection,
                                           no_categories,
                                           observations,
                                           no_persons,
                                           no_groups,
                                           group_index,
                                           n_cat_obs,
                                           sufficient_blume_capel,
                                           rest_matrix,
                                           independent_thresholds,
                                           ordinal_variable,
                                           reference_category,
                                           indicator,
                                           inclusion_probability_difference,
                                           proposal_sd_main,
                                           proposal_sd_pairwise,
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
                                           difference_selection,
                                           no_pairwise,
                                           no_variables,
                                           index);

    IntegerMatrix indicator = out["indicator"];
    NumericMatrix main_effects = out["main_effects"];
    NumericMatrix pairwise_effects = out["pairwise_effects"];
    NumericMatrix rest_matrix = out["rest_matrix"];
    NumericMatrix proposal_sd_pairwise = out["proposal_sd_pairwise"];
    NumericMatrix proposal_sd_main = out["proposal_sd_main"];

    if(difference_selection == true) {
      if(pairwise_difference_prior == "Beta-Bernoulli") {
        int sumG = 0;
        for(int i = 0; i < no_variables - 1; i++) {
          for(int j = i + 1; j < no_variables; j++) {
            sumG += indicator(i, j);
          }
        }
        double probability = R::rbeta(pairwise_beta_bernoulli_alpha + sumG,
                                      pairwise_beta_bernoulli_beta + no_pairwise - sumG);

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
    //Compute running averages -----------------------------------------------
    for(int row = 0; row < no_main; row++) {
      for(int column = 0; column < no_groups; column++) {

        out_main(row, column) *= iteration;
        out_main(row, column) += main_effects(row, column);
        out_main(row, column) /= iteration + 1;

      }
    }

    for(int row = 0; row < no_pairwise; row++) {
      for(int column = 0; column < no_groups; column++) {

        out_pairwise(row, column) *= iteration;
        out_pairwise(row, column) += pairwise_effects(row, column);
        out_pairwise(row, column) /= iteration + 1;

      }
    }

    for(int i = 0; i < no_variables - 1; i++) {
      for(int j = i + 1; j < no_variables; j++) {
        out_indicator(i, j) *= iteration;
        out_indicator(i, j) += static_cast<double>(indicator(i, j));
        out_indicator(i, j) /= iteration + 1;
        out_indicator(j, i) = out_indicator(i, j);
      }
    }
  }

  return List::create(Named("main") = out_main,
                      Named("pairwise") = out_pairwise,
                      Named("indicator") = out_indicator);
}