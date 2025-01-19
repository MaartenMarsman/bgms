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
// 5. Add "enum" to handle threshold scenarios
// ----------------------------------------------------------------------------|

double robbins_monro_update(double current_sd,
                            double log_ar,
                            double target_ar,
                            double epsilon_lo,
                            double epsilon_hi,
                            double rm_weight) {

  // Normalize the acceptance probability
  double obs_ar = 1.0;
  if(log_ar < 0) {
    obs_ar = std::exp(log_ar);
  }

  // Update the proposal standard deviation
  double update = current_sd + (obs_ar - target_ar) * rm_weight;

  // Handle NaN cases by resetting to default value
  if (std::isnan(update)) {
    update = 1.0; // Default proposal standard deviation
  }

  // Clamp the updated standard deviation within bounds
  return std::clamp(update, epsilon_lo, epsilon_hi);
}


NumericVector group_thresholds_for_variable (int variable,
                                             LogicalVector ordinal_variable,
                                             int group,
                                             int no_groups,
                                             IntegerMatrix no_categories,
                                             NumericMatrix main_effects,
                                             IntegerMatrix main_index,
                                             NumericMatrix projection,
                                             bool independent_thresholds) {
  int vector_length = (ordinal_variable[variable]) ? no_categories(variable, group) : 2;
  NumericVector GroupThresholds(vector_length);
  int base_category_index = main_index(variable, 0);

  if (!independent_thresholds) {
    int n_cats = no_categories(variable, group);
    for (int category = 0; category < n_cats; category++) {
      int category_index = base_category_index + category;
      double threshold = main_effects(category_index, 0);

      for (int h = 0; h < no_groups - 1; h++) {
        threshold += projection(group, h) * main_effects(category_index, h + 1);
      }

      GroupThresholds[category] = threshold;
    }
  } else {
    if (ordinal_variable[variable]) {
      int n_cats = no_categories(variable, group);
      for (int category = 0; category < n_cats; category++) {
        int category_index = base_category_index + category;
        GroupThresholds[category] = main_effects(category_index, group);
      }
    } else {
      GroupThresholds[0] = main_effects(base_category_index, group);
      GroupThresholds[1] = main_effects(base_category_index + 1, group);
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
    person = missing_index(missing, 0);
    variable = missing_index(missing, 1);
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
      for(int category = 1; category <= no_categories(variable, gr); category++) {
        exponent = GroupThresholds(category - 1);
        exponent += category * rest_score;
        cumsum += std::exp(exponent);
        probabilities[category] = cumsum;
      }
    } else {
      //Blume-Capel variable -------------------------------------------------
      cumsum = 0.0;
      for(int category = 0; category <= no_categories(variable, gr); category++) {
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
  double delta_state = 2.0 * (proposed_state - current_state);
  double pseudolikelihood_ratio = 0.0;

  // --------------------------------------------------------------------------|
  // Compute log pseudo-likelihood ratio ---------------------------------------
  // --------------------------------------------------------------------------|
  for(int gr = 0; gr < no_groups; gr++) {
    NumericVector GroupThresholds_v1 = group_thresholds_for_variable (
      variable1, ordinal_variable, gr, no_groups, no_categories, main_effects,
      main_index, projection, independent_thresholds);

    NumericVector GroupThresholds_v2 = group_thresholds_for_variable (
      variable2, ordinal_variable, gr, no_groups, no_categories, main_effects,
      main_index, projection, independent_thresholds);

    int n_cats_v1 = no_categories(variable1, gr);
    int n_cats_v2 = no_categories(variable2, gr);

    //Loop over the pseudo-likelihoods for persons in group gr -----------------
    for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
      int obs_score1 = observations(person, variable1);
      int obs_score2 = observations(person, variable2);

      pseudolikelihood_ratio += obs_score1 * obs_score2 * delta_state;

      for (int variable = 1; variable <= 2; variable++) {
        int var = (variable == 1) ? variable1 : variable2;
        int obs = (variable == 1) ? obs_score2 : obs_score1;
        double obs_current = obs * current_state;
        double obs_proposed = obs * proposed_state;

        int n_cats = (variable == 1) ? n_cats_v1 : n_cats_v2;
        NumericVector& GroupThresholds = (variable == 1) ? GroupThresholds_v1 : GroupThresholds_v2;

        double rest_score = rest_matrix(person, var) - obs_current;
        double bound = (rest_score > 0) ? n_cats * rest_score : 0.0;
        double denominator_prop = 0.0;
        double denominator_curr = 0.0;

        if(ordinal_variable[var]) {
          // Regular binary or ordinal MRF variable
          denominator_prop += std::exp(-bound);
          denominator_curr += std::exp(-bound);
          for(int cat = 0; cat < n_cats; cat++) {
            int score = cat + 1;
            double exponent = GroupThresholds[cat] + score * rest_score - bound;
            denominator_prop += std::exp(exponent + score * obs_proposed);
            denominator_curr += std::exp(exponent + score * obs_current);
          }
        } else {
          // Blume-Capel ordinal MRF variable
          for(int cat = 0; cat <= n_cats; cat++) {
            double exponent = GroupThresholds[0] * cat;
            exponent += GroupThresholds[1] *
              (cat - reference_category[var]) *
              (cat - reference_category[var]);
            exponent += cat * rest_score - bound;
            denominator_prop += std::exp(exponent + cat * obs_proposed);
            denominator_curr += std::exp(exponent + cat * obs_current);
          }
        }
        pseudolikelihood_ratio -= std::log(denominator_prop);
        pseudolikelihood_ratio += std::log(denominator_curr);
      }
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
  double log_prob;
  double exp_neg_log_t_phi = std::exp(-std::log(t) * phi);

  for(int variable1 = 0; variable1 <  no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  no_variables; variable2++) {
      int int_index = pairwise_index(variable1, variable2);

      double current_state = pairwise_effects(int_index, 0);
      double proposed_state = R::rnorm(current_state, proposal_sd_pairwise(int_index, 0));

      log_prob = compare_anova_log_pseudolikelihood_ratio_interaction(
        main_effects, main_index, projection, observations, no_groups,
        group_index, no_categories, independent_thresholds, no_persons,
        variable1, variable2, proposed_state, current_state, rest_matrix,
        ordinal_variable, reference_category);

      log_prob += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_prob -= R::dcauchy(current_state, 0.0, interaction_scale, true);

      double U = R::unif_rand();
      if(std::log(U) < log_prob) {
        double state_difference = proposed_state - current_state;
        pairwise_effects(int_index, 0) = proposed_state;

        //Update the matrix of rest scores
        for(int person = 0; person < no_persons; person++) {
          double obs1 = observations(person, variable1);
          double obs2 = observations(person, variable2);
          rest_matrix(person, variable1) += obs2 * state_difference;
          rest_matrix(person, variable2) += obs1 * state_difference;
        }
      }

      // Robbins-Monro update to the proposal sd
      proposal_sd_pairwise(int_index, 0) = robbins_monro_update(
        proposal_sd_pairwise(int_index, 0), log_prob, target_ar, epsilon_lo,
        epsilon_hi, exp_neg_log_t_phi);
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
  double pseudolikelihood_ratio = 0.0;
  double delta_state = 2.0 * (proposed_state - current_state);

  // Loop over groups
  for(int gr = 0; gr < no_groups; gr++) {

    NumericVector GroupThresholds_v1 = group_thresholds_for_variable(
      variable1, ordinal_variable, gr, no_groups, no_categories, main_effects,
      main_index, projection, independent_thresholds);

    NumericVector GroupThresholds_v2 = group_thresholds_for_variable(
      variable2, ordinal_variable, gr, no_groups, no_categories, main_effects,
      main_index, projection, independent_thresholds);

    double P = projection(gr, h);
    double delta_state_group = delta_state * P;

    // Cache the number of categories
    int n_cats_v1 = no_categories(variable1, gr);
    int n_cats_v2 = no_categories(variable2, gr);

    //Loop over the pseudo-likelihoods for persons in group gr -----------------
    for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
      // Cache observation scores and scaled terms
      int obs_score1 = observations(person, variable1);
      int obs_score2 = observations(person, variable2);
      double obs_proposed_p1 = obs_score2 * proposed_state * P;
      double obs_current_p1 = obs_score2 * current_state * P;
      double obs_proposed_p2 = obs_score1 * proposed_state * P;
      double obs_current_p2 = obs_score1 * current_state * P;

      pseudolikelihood_ratio +=  obs_score1 * obs_score2 * delta_state_group;

      // Process each variable
      for (int variable = 1; variable <= 2; variable++) {
        int var = (variable == 1) ? variable1 : variable2;
        int n_cats = (variable == 1) ? n_cats_v1 : n_cats_v2;
        NumericVector& GroupThresholds = (variable == 1) ? GroupThresholds_v1 : GroupThresholds_v2;
        double obs_proposed_p = (variable == 1) ? obs_proposed_p1 : obs_proposed_p2;
        double obs_current_p = (variable == 1) ? obs_current_p1 : obs_current_p2;

        // Calculate rest_score
        double rest_score = rest_matrix(person, var) - obs_current_p;
        double bound = (rest_score > 0) ? n_cats * rest_score : 0.0;

        // Initialize denominators
        double denominator_prop = 0.0;
        double denominator_curr = 0.0;

        // Compute denominators
        if (ordinal_variable[var]) {
          // Binary or ordinal MRF variable
          denominator_prop += std::exp(-bound);
          denominator_curr += std::exp(-bound);

          for (int cat = 0; cat < n_cats; cat++) {
            int score = cat + 1;
            double exponent = GroupThresholds[cat] + score * rest_score - bound;
            denominator_prop += std::exp(exponent + score * obs_proposed_p);
            denominator_curr += std::exp(exponent + score * obs_current_p);
          }
        } else {
          // Blume-Capel ordinal MRF variable
          for (int cat = 0; cat <= n_cats; cat++) {
            double exponent = GroupThresholds[0] * cat +
              GroupThresholds[1] * (cat - reference_category[var]) *
              (cat - reference_category[var]) +
              cat * rest_score - bound;
            denominator_prop += std::exp(exponent + cat * obs_proposed_p);
            denominator_curr += std::exp(exponent + cat * obs_current_p);
          }
        }

        // Update pseudolikelihood ratio
        pseudolikelihood_ratio -= std::log(denominator_prop);
        pseudolikelihood_ratio += std::log(denominator_curr);
      }
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
  double exp_neg_log_t_phi = std::exp(-std::log(t) * phi);

  // Loop over variables pairs
  for(int variable1 = 0; variable1 <  no_variables - 1; variable1++) {
    for(int variable2 = variable1 + 1; variable2 <  no_variables; variable2++) {
      if (indicator(variable1, variable2) == 1) {
        int int_index = pairwise_index(variable1, variable2);

        // Loop over groups
        for(int h = 1; h < no_groups; h++) {
          double current_state = pairwise_effects(int_index, h);
          double proposed_state = R::rnorm(current_state, proposal_sd_pairwise(int_index, h));

          // Compute log pseudo-likelihood ratio
          double log_prob = compare_anova_log_pseudolikelihood_ratio_pairwise_difference(
            main_effects, main_index, projection, observations, no_groups,
            group_index, no_categories, independent_thresholds, no_persons,
            variable1, variable2, h - 1, proposed_state, current_state,
            rest_matrix, ordinal_variable, reference_category);

          log_prob += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
          log_prob -= R::dcauchy(current_state, 0.0, interaction_scale, true);

          double U = R::unif_rand();
          if(std::log(U) < log_prob) {
            pairwise_effects(int_index, h) = proposed_state;

            // Update rest matrix
            for(int gr = 0; gr < no_groups; gr++) {
              double state_difference = (proposed_state - current_state) * projection(gr, h - 1);
              for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
                double obs1 = observations(person, variable1);
                double obs2 = observations(person, variable2);
                rest_matrix(person, variable1) += obs2 * state_difference;
                rest_matrix(person, variable2) += obs1 * state_difference;
              }
            }
          }

          // Robbins-Monro update to the proposal sd
          proposal_sd_pairwise(int_index, h) = robbins_monro_update(
            proposal_sd_pairwise(int_index, h), log_prob, target_ar, epsilon_lo,
            epsilon_hi, exp_neg_log_t_phi);
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

  // Cache the number of categories for this variable
  int n_cats = no_categories(variable, 0);
  int cat_index = main_index(variable, 0);

  NumericVector GroupThresholds(n_cats);

  // Iterate over categories
  for(int category = 0; category < n_cats; category++) {
    double current_state = main_effects(cat_index + category, 0);
    double exp_current = std::exp(current_state);
    double c = (threshold_alpha + threshold_beta) / (1 + exp_current);

    // Update group thresholds
    for (int gr = 0; gr < no_groups; gr++) {
      for (int cat = 0; cat < n_cats; cat++) {
        double threshold = main_effects(cat_index + cat, 0);
        for (int h = 1; h < no_groups; h++) {
          threshold += projection(gr, h - 1) * main_effects(cat_index + cat, h);
        }
        GroupThresholds[cat] = threshold;
      }
      GroupThresholds[category] -= main_effects(cat_index + category, 0);

      // Compute `q` and `r` for each person in the group
      for (int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
        double rest_score = rest_matrix(person, variable);
        double q_person = 1.0;
        for (int cat = 0; cat < n_cats; ++cat) {
          if (cat != category) {
            double exponent = GroupThresholds[cat] + (cat + 1) * rest_score;
            q_person += std::exp(exponent);
          }
        }
        double exponent_r = GroupThresholds[category] + (category + 1) * rest_score;
        double r_person = std::exp(exponent_r);
        q[person] = q_person;
        r[person] = r_person;
        c += r_person / (q_person + r_person * exp_current);
      }
    }

    // Update c
    double tmp = no_persons + threshold_alpha + threshold_beta - exp_current * c;
    c /= tmp;

    // Generalized beta-prime proposal
    double a = threshold_alpha;
    double b = no_persons + threshold_beta;
    for(int gr = 0; gr < no_groups; gr++) {
      IntegerMatrix n_cat_obs_gr = n_cat_obs[gr];
      a += n_cat_obs_gr(category + 1, variable);
      b -= n_cat_obs_gr(category + 1, variable);
    }

    double tmp_beta = R::rbeta(a, b);
    double proposed_state = std::log(tmp_beta / (1  - tmp_beta) / c);
    double exp_proposed = std::exp(proposed_state);

    // Compute log acceptance probability
    double log_prob = 0.0;

    // Pseudo-likelihood ratio
    for (int gr = 0; gr < no_groups; ++gr) {
      for (int person = group_index(gr, 0); person <= group_index(gr, 1); ++person) {
        log_prob += std::log(q[person] + r[person] * exp_current);
        log_prob -= std::log(q[person] + r[person] * exp_proposed);
      }
    }

    // Prior ratio
    log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + exp_proposed);
    log_prob += (threshold_alpha + threshold_beta) * std::log(1 + exp_current);

    // Proposal ratio
    log_prob -= (a + b) * std::log(1 + c * exp_current);
    log_prob += (a + b) * std::log(1 + c * exp_proposed);

    // Metropolis-Hastings acceptance step
    double U = std::log(R::unif_rand());
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

  double pseudolikelihood_ratio = 0.0;
  double delta_state = proposed_state - current_state;
  int cat_index = main_index(variable, 0);

  // Loop over groups
  for(int gr = 0; gr < no_groups; gr++) {
    int n_cats = no_categories(variable, gr);
    NumericVector current_thresholds(n_cats);
    NumericVector proposed_thresholds(n_cats);
    double P = projection(gr, h);

    // Compute current and proposed thresholds
    for(int cat = 0; cat < n_cats; cat++) {
      int full_cat_index = cat_index + cat;
      double threshold = main_effects(full_cat_index, 0);
      for (int hh = 1; hh < no_groups; ++hh) {
        threshold += projection(gr, hh - 1) * main_effects(full_cat_index, hh);
      }
      current_thresholds[cat] = threshold;
      proposed_thresholds[cat] = threshold;
    }

    // Adjust vector of proposals for the specific category
    proposed_thresholds[category] -= P * current_state;
    proposed_thresholds[category] += P * proposed_state;

    // Contribution from delta_state
    IntegerMatrix n_cat_obs_gr = n_cat_obs[gr];
    pseudolikelihood_ratio += delta_state * P * n_cat_obs_gr(category + 1, variable);

    // Loop over persons in the group

    for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
      double rest_score = rest_matrix(person, variable);
      double bound = (rest_score > 0) ? n_cats * rest_score : 0.0;

      double denominator_proposed = std::exp(-bound);
      double denominator_current = std::exp(-bound);
      for(int cat = 0; cat < n_cats; cat++) {
        double exponent = (cat + 1) * rest_score - bound;
        denominator_proposed += std::exp(exponent + proposed_thresholds[cat]);
        denominator_current += std::exp(exponent + current_thresholds[cat]);
      }
      // Update pseudolikelihood ratio
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
  // Precompute Robbins-Monro term
  double exp_neg_log_t_phi = std::exp(-log(t) * phi);

  // Check if the variable is active
  if (indicator(variable, variable) != 1) {
    return;
  }

  // Look up and cache base category index
  int base_cat_index = main_index(variable, 0);

  // Loop over categories
  int n_cats = no_categories(variable, 0);
  for(int category = 0; category < n_cats; category++) {
    int cat_index = base_cat_index + category;

    // Loop over groups (starting from h = 1)
    for(int h = 1; h < no_groups; h++) {
      double current_state = main_effects(cat_index, h);
      double proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index, h));

      // Compute log pseudo-likelihood ratio
      double log_prob = compare_anova_log_pseudolikelihood_ratio_main_difference(
        main_effects, main_index, projection, observations, no_groups,
        group_index, no_categories, no_persons, rest_matrix, n_cat_obs,
        variable, category, h - 1, proposed_state, current_state);

      // Add prior contributions
      log_prob += R::dcauchy(proposed_state, 0.0, main_difference_scale, true);
      log_prob -= R::dcauchy(current_state, 0.0, main_difference_scale, true);

      // Metropolis-Hastings acceptance step
      double U = R::unif_rand();
      if(std::log(U) < log_prob) {
        main_effects(cat_index, h) = proposed_state;
      }

      // Robbins-Monro update to the proposal sd
      proposal_sd_main(cat_index, h) = robbins_monro_update(
        proposal_sd_main(cat_index, h), log_prob, target_ar, epsilon_lo,
        epsilon_hi, exp_neg_log_t_phi);
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
    for(int category = 0; category <= no_categories(variable, gr); category++) {
      linear_score = category;
      quadratic_score = (category - reference_category[variable]) *
        (category - reference_category[variable]);

      constant_numerator[category] = linear_current * linear_score;
      constant_numerator[category] += quadratic_current * quadratic_score;
      constant_denominator[category] = linear_proposed * linear_score ;
      constant_denominator[category] += quadratic_proposed * quadratic_score;
      for(int h = 1; h < no_groups; h++) {
        double P = projection(gr, h - 1);
        double m1 =  main_effects(cat_index, h);
        double m2 =  main_effects(cat_index + 1, h);
        constant_numerator[category] += P * m1 * linear_score;
        constant_numerator[category] += P * m2 * quadratic_score;
        constant_denominator[category] += P * m1 * linear_score;
        constant_denominator[category] += P * m2 * quadratic_score;
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

    for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
      rest_score = rest_matrix(person, variable);
      if(rest_score > 0) {
        bound = no_categories(variable, gr) * rest_score + lbound;
      } else {
        bound = lbound;
      }

      numerator = 0.0;
      denominator = 0.0;
      for(int category = 0; category <= no_categories[variable]; category ++) {
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
  double exp_neg_log_t_phi = std::exp(-log(t) * phi);
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

  // Robbins-Monro update to the proposal sd
  proposal_sd_main(cat_index, 0) = robbins_monro_update(
    proposal_sd_main(cat_index, 0), log_prob, target_ar, epsilon_lo,
    epsilon_hi, exp_neg_log_t_phi);

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

  // Robbins-Monro update to the proposal sd
  proposal_sd_main(cat_index + 1, 0) = robbins_monro_update(
    proposal_sd_main(cat_index + 1, 0), log_prob, target_ar, epsilon_lo,
    epsilon_hi, exp_neg_log_t_phi);

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
    for(int category = 0; category <= no_categories(variable, gr); category++) {
      linear_score = category;
      quadratic_score = (category - reference_category[variable]) *
        (category - reference_category[variable]);

      constant_numerator[category] = main_effects(cat_index, 0) * linear_score;
      constant_numerator[category] +=  main_effects(cat_index + 1, 0) * quadratic_score;
      constant_denominator[category] =  main_effects(cat_index, 0) * linear_score ;
      constant_denominator[category] += main_effects(cat_index + 1, 0) * quadratic_score;
      for(int hh = 1; hh < no_groups; hh++) {
        double P = projection(gr, hh - 1);
        double m1 =  main_effects(cat_index, hh);
        double m2 =  main_effects(cat_index + 1, hh);
        constant_numerator[category] += P * m1 * linear_score;
        constant_numerator[category] += P * m2 * quadratic_score;
        constant_denominator[category] += P * m1 * linear_score;
        constant_denominator[category] += P * m2 * quadratic_score;
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

    double sp0 = sufficient_blume_capel_gr(0, variable) * P;
    double sp1 = sufficient_blume_capel_gr(1, variable) * P;
    pseudolikelihood_ratio += sp0 * linear_proposed;
    pseudolikelihood_ratio += sp1 * quadratic_proposed;
    pseudolikelihood_ratio -= sp0 * linear_current;
    pseudolikelihood_ratio -= sp1 * quadratic_current;

    for(int person = group_index(gr, 0); person <= group_index(gr, 1); person++) {
      rest_score = rest_matrix(person, variable);
      if(rest_score > 0) {
        bound = no_categories(variable, gr) * rest_score + lbound;
      } else {
        bound = lbound;
      }

      numerator = 0.0;
      denominator = 0.0;
      for(int category = 0; category <= no_categories[variable]; category ++) {
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
  double exp_neg_log_t_phi = std::exp(-log(t) * phi);

  if(indicator(variable, variable) == 1) {
    //Loop over the group difference effects--------------------------------------
    for(int h = 1; h < no_groups; h++) {
      //----------------------------------------------------------------------------
      // Adaptive Metropolis for the linear Blume-Capel parameter
      //----------------------------------------------------------------------------
      int cat_index = main_index(variable, 0);
      current_state = main_effects(cat_index, h);
      proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index, h));

      //----------------------------------------------------------------------------
      //Compute the log acceptance probability -------------------------------------
      //----------------------------------------------------------------------------
      double linear_current = current_state;
      double quadratic_current = main_effects(cat_index + 1, h);
      double linear_proposed = proposed_state;
      double quadratic_proposed = main_effects(cat_index + 1, h);

      log_prob = compare_anova_log_pseudolikelihood_ratio_main_difference_blumecapel(linear_current,
                                                                                     quadratic_current,
                                                                                     linear_proposed,
                                                                                     quadratic_proposed,
                                                                                     variable,
                                                                                     h - 1,
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
        main_effects(cat_index, h) = proposed_state;
      }

      // Robbins-Monro update to the proposal sd
      proposal_sd_main(cat_index, h) = robbins_monro_update(
        proposal_sd_main(cat_index, h), log_prob, target_ar, epsilon_lo,
        epsilon_hi, exp_neg_log_t_phi);

      //---------------------------------------------------------------------------|
      // Adaptive Metropolis for the quadratic Blume-Capel parameter
      //---------------------------------------------------------------------------|
      current_state = main_effects(cat_index + 1, h);
      proposed_state = R::rnorm(current_state, proposal_sd_main(cat_index + 1, h));

      //----------------------------------------------------------------------------
      //Compute the log acceptance probability -------------------------------------
      //----------------------------------------------------------------------------

      linear_current = main_effects(cat_index, h);
      quadratic_current = current_state;
      linear_proposed =  main_effects(cat_index, h);
      quadratic_proposed =  proposed_state;

      log_prob = compare_anova_log_pseudolikelihood_ratio_main_difference_blumecapel(linear_current,
                                                                                     quadratic_current,
                                                                                     linear_proposed,
                                                                                     quadratic_proposed,
                                                                                     variable,
                                                                                     h - 1,
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
        main_effects(cat_index + 1, h) = proposed_state;
      }

      // Robbins-Monro update to the proposal sd
      proposal_sd_main(cat_index + 1, h) = robbins_monro_update(
        proposal_sd_main(cat_index + 1, h), log_prob, target_ar, epsilon_lo,
        epsilon_hi, exp_neg_log_t_phi);

    }
  }
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

  // Number of persons in the group
  int no_persons = group_index(group, 1) - group_index(group, 0) + 1;

  // Cache group-specific category observations
  IntegerMatrix n_cat_obs_gr = n_cat_obs[group];

  // Base category index for the variable
  int base_cat_index = main_index(variable, 0);
  int n_cats = no_categories(variable, group);

  // Temporary storage for pseudo-likelihood elements
  NumericVector q(no_persons);
  NumericVector r(no_persons);

  // Loop over categories
  for(int category = 0; category < n_cats; category++) {
    double current_state = main_effects(base_cat_index + category, group);
    double exp_current = std::exp(current_state);

    // Initialize `c`
    double c = (threshold_alpha + threshold_beta) / (1 + exp_current);

    // Loop over persons
    for(int person = 0; person < no_persons; person++) {
      double rest_score = rest_matrix(group_index(group, 0) + person, variable);

      // Compute `q` and `r` for this person
      double q_person = 1.0;
      for(int cat = 0; cat < n_cats; cat++) {
        if(cat != category) {
          q_person += std::exp(main_effects(base_cat_index + cat, group) + (cat + 1) * rest_score);
        }
      }
      double r_person = std::exp((category + 1) * rest_score);

      // Update `q` and `r` vectors
      q[person] = q_person;
      r[person] = r_person;

      // Update `c`
      c += r_person / (q_person + r_person * exp_current);
    }

    // Finalize `c`
    c /= (no_persons + threshold_alpha + threshold_beta - exp_current * c);

    // Generalized beta-prime proposal
    double a = n_cat_obs_gr(category + 1, variable) + threshold_alpha;
    double b = no_persons + threshold_beta - n_cat_obs_gr(category + 1, variable);
    double tmp = R::rbeta(a, b);
    double proposed_state = std::log(tmp / ((1 - tmp) * c));
    double exp_proposed = std::exp(proposed_state);


    // Compute log acceptance probability
    double log_prob = 0.0;

    // Pseudo-likelihood ratio
    for (int person = 0; person < no_persons; ++person) {
      log_prob += std::log(q[person] + r[person] * exp_current);
      log_prob -= std::log(q[person] + r[person] * exp_proposed);
    }

    // Prior ratio
    log_prob -= (threshold_alpha + threshold_beta) * std::log(1 + exp_proposed);
    log_prob += (threshold_alpha + threshold_beta) * std::log(1 + exp_current);

    // Proposal ratio
    log_prob -= (a + b) * std::log(1 + c * exp_current);
    log_prob += (a + b) * std::log(1 + c * exp_proposed);

    // Metropolis-Hastings acceptance step
    double U = std::log(R::unif_rand());
    if(U < log_prob) {
      main_effects(base_cat_index + category, group) = proposed_state;
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
  double exp_neg_log_t_phi = std::exp(-log(t) * phi);

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

  for(int category = 0; category <= no_categories(variable, group); category ++) {
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

  for(int person = group_index(group, 0); person <= group_index(group, 1); person++) {
    rest_score = rest_matrix(person, variable);
    if(rest_score > 0) {
      bound = no_categories(variable, group) * rest_score + lbound;
    } else {
      bound = lbound;
    }
    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);
    for(int score = 1; score <= no_categories(variable, group); score++) {
      exponent = score * rest_score - bound;
      numerator += std::exp(constant_numerator[score] + exponent);
      denominator += std::exp(constant_denominator[score] + exponent);
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

  // Robbins-Monro update to the proposal sd
  proposal_sd_main(cat_index, group) = robbins_monro_update(
    proposal_sd_main(cat_index, group), log_prob, target_ar, epsilon_lo,
    epsilon_hi, exp_neg_log_t_phi);


  //----------------------------------------------------------------------------
  //Adaptive Metropolis for the quadratic Blume-Capel parameter
  //----------------------------------------------------------------------------
  current_state = main_effects(cat_index + 1, group);
  proposed_state = R::rnorm(current_state,
                            proposal_sd_main(cat_index + 1, group));

  //Precompute terms for the log acceptance probability ------------------------
  difference = proposed_state - current_state;

  for(int category = 0; category <= no_categories(variable, group); category ++) {
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

  for(int person = group_index(group, 0); person <= group_index(group, 1); person++) {
    rest_score = rest_matrix(person, variable);
    if(rest_score > 0) {
      bound = no_categories[variable] * rest_score + lbound;
    } else {
      bound = lbound;
    }

    numerator = std::exp(constant_numerator[0] - bound);
    denominator = std::exp(constant_denominator[0] - bound);

    for(int score = 1; score <= no_categories(variable, group); score ++) {
      exponent = score * rest_score - bound;
      numerator += std::exp(constant_numerator[score] + exponent);
      denominator += std::exp(constant_denominator[score] + exponent);
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

  // Robbins-Monro update to the proposal sd
  proposal_sd_main(cat_index + 1, group) = robbins_monro_update(
    proposal_sd_main(cat_index + 1, group), log_prob, target_ar, epsilon_lo,
    epsilon_hi, exp_neg_log_t_phi);

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
                                 int no_variables) {

  NumericMatrix proposal_sd_blumecapel_group(no_variables, 2);


  // Update pairwise interaction parameters
  compare_anova_metropolis_interaction(
    main_effects, pairwise_effects, main_index,
    pairwise_index, projection, observations, no_groups, group_index,
    no_categories, independent_thresholds, no_persons, rest_matrix,
    ordinal_variable, reference_category, proposal_sd_pairwise,
    interaction_scale, no_variables, phi, target_ar, t, epsilon_lo, epsilon_hi);

  // Update pairwise differences
  compare_anova_metropolis_pairwise_difference(
    main_effects, pairwise_effects, main_index, pairwise_index, projection,
    observations, no_groups, group_index, no_categories, independent_thresholds,
    indicator, no_persons, rest_matrix, ordinal_variable, reference_category,
    proposal_sd_pairwise, interaction_scale, no_variables, phi, target_ar, t, epsilon_lo,
    epsilon_hi);

  // Update thresholds based on main_model input
  for (int variable = 0; variable < no_variables; variable++) {
    if (independent_thresholds) {
      // Group-specific thresholds (main_model = "Free")
      for (int group = 0; group < no_groups; ++group) {
        if (ordinal_variable[variable]) {
          compare_anova_metropolis_thresholds_regular_free(
            main_effects, main_index, observations, no_groups, group_index,
            no_categories, rest_matrix, n_cat_obs, threshold_alpha, threshold_beta,
            variable, group);
        } else {
          compare_anova_metropolis_thresholds_blumecapel_free(
            main_effects, main_index, observations, no_groups, group_index,
            reference_category, no_categories, sufficient_blume_capel, no_persons,
            rest_matrix, n_cat_obs, threshold_alpha, threshold_beta, variable, group,
            proposal_sd_main, phi, target_ar, t, epsilon_lo, epsilon_hi);
        }
      }
    } else {
      // Model group differences in thresholds (main_model != "Free")
      if (ordinal_variable[variable]) {
        compare_anova_metropolis_threshold_regular(
          main_effects, main_index, projection, observations, no_groups,
          group_index, no_categories, no_persons, rest_matrix, n_cat_obs,
          threshold_alpha, threshold_beta, variable);

        compare_anova_metropolis_main_difference_regular(
          main_effects, main_index, projection, observations, no_groups,
          group_index, no_categories, no_persons, rest_matrix, n_cat_obs, variable,
          indicator, proposal_sd_main, main_difference_scale, phi, target_ar, t,
          epsilon_lo, epsilon_hi);
      } else {
        compare_anova_metropolis_threshold_blumecapel(
          main_effects, main_index, projection, no_categories, sufficient_blume_capel,
          no_persons, no_groups, group_index, variable, reference_category,
          threshold_alpha, threshold_beta, rest_matrix, proposal_sd_main, phi, target_ar,
          t, epsilon_lo, epsilon_hi);

        compare_anova_metropolis_main_difference_blumecapel(
          main_effects, main_index, projection, no_categories, sufficient_blume_capel,
          no_persons, no_groups, group_index, variable, reference_category,
          threshold_alpha, threshold_beta, rest_matrix, indicator, proposal_sd_main,
          phi, target_ar, t, epsilon_lo, epsilon_hi);
      }
    }
  }

  // Return updated parameters
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

  // Dimensions and parameters
  int no_variables = observations.ncol();
  int no_persons = observations.nrow();
  int no_main = 1 + main_index(no_variables - 1, 1);
  int no_pairwise = no_variables * (no_variables - 1) / 2;

  // Initialize person-group indicator
  IntegerVector person_group_indicator (no_persons);
  for(int group = 0; group < no_groups; group++) {
    for(int person = group_index(group, 0); person < group_index(group, 1) + 1; person++) {
      person_group_indicator[person] = group;
    }
  }

  // Initialize model parameters
  NumericMatrix main_effects(no_main, no_groups);
  NumericMatrix pairwise_effects(no_pairwise, no_groups);
  IntegerMatrix indicator(no_variables, no_variables);
  std::fill(indicator.begin(), indicator.end(), 1);

  // Adaptive Metropolis proposal standard deviations
  NumericMatrix proposal_sd_main(no_main, no_groups);
  NumericMatrix proposal_sd_pairwise(no_pairwise, no_groups);
  std::fill(proposal_sd_main.begin(), proposal_sd_main.end(), 1.0);
  std::fill(proposal_sd_pairwise.begin(), proposal_sd_pairwise.end(), 1.0);

  // Robbins-Monro parameters
  double phi = 0.75, target_ar = 0.234;
  double epsilon_lo = 1.0 / static_cast<double>(no_persons);
  double epsilon_hi = 2.0;

  // Rest matrix for pseudo-likelihoods
  NumericMatrix rest_matrix(no_persons, no_variables);

  // Output matrices
  NumericMatrix out_main(no_main, no_groups);
  NumericMatrix out_pairwise(no_pairwise, no_groups);
  NumericMatrix out_indicator(no_variables, no_variables);

  // Progress bar
  Progress p(iter + burnin, display_progress);

  // Gibbs sampling iterations
  for (int iteration = 0; iteration < iter + burnin; ++iteration) {
    if (Progress::check_abort()) {
      return List::create(Named("main") = out_main,
                          Named("pairwise") = out_pairwise,
                          Named("indicator") = out_indicator);

    }
    p.increment();

    // Handle missing data if required
    if (na_impute) {
      List impute_out = compare_anova_impute_missing_data(
        main_effects, pairwise_effects, main_index, pairwise_index,
        projection, observations, no_groups, person_group_indicator,
        n_cat_obs, sufficient_blume_capel, no_categories, rest_matrix,
        missing_index, ordinal_variable, reference_category, independent_thresholds);

      IntegerMatrix observations_tmp = impute_out["observations"];
      List n_cat_obs_tmp = impute_out["n_cat_obs"];
      List sufficient_blume_capel_tmp = impute_out["sufficient_blume_capel"];
      NumericMatrix rest_matrix_tmp = impute_out["rest_matrix"];

      // Reassign to original variables
      observations = observations_tmp;
      n_cat_obs = n_cat_obs_tmp;
      sufficient_blume_capel = sufficient_blume_capel_tmp;
      rest_matrix = rest_matrix_tmp;
    }

    // Update parameters using Gibbs step
    List step_out = compare_anova_gibbs_step_gm(
      main_effects, main_index, pairwise_effects, pairwise_index, projection,
      no_categories, observations, no_persons, no_groups, group_index, n_cat_obs,
      sufficient_blume_capel, rest_matrix, independent_thresholds, ordinal_variable,
      reference_category, indicator, inclusion_probability_difference,
      proposal_sd_main, proposal_sd_pairwise, interaction_scale,
      main_difference_scale, pairwise_difference_scale, threshold_alpha,
      threshold_beta, phi, target_ar, iteration, epsilon_lo, epsilon_hi,
      difference_selection, no_pairwise, no_variables);

    IntegerMatrix indicator_tmp = step_out["indicator"];
    NumericMatrix main_effects_tmp = step_out["main_effects"];
    NumericMatrix pairwise_effects_tmp = step_out["pairwise_effects"];
    NumericMatrix rest_matrix_tmp = step_out["rest_matrix"];
    NumericMatrix proposal_sd_pairwise_tmp = step_out["proposal_sd_pairwise"];
    NumericMatrix proposal_sd_main_tmp = step_out["proposal_sd_main"];

    // Update matrices from Gibbs step output
    indicator = indicator_tmp;
    main_effects = main_effects_tmp;
    pairwise_effects = pairwise_effects_tmp;
    rest_matrix = rest_matrix_tmp;
    proposal_sd_pairwise = proposal_sd_pairwise_tmp;
    proposal_sd_main = proposal_sd_main_tmp;


    // Update inclusion probabilities if difference_selection is enabled
    if (difference_selection) {
      int sumG = 0;

      // Pairwise differences
      if (pairwise_difference_prior == "Beta-Bernoulli") {
        for (int i = 0; i < no_variables - 1; ++i) {
          for (int j = i + 1; j < no_variables; ++j) {
            sumG += indicator(i, j);
          }
        }
        double prob = R::rbeta(pairwise_beta_bernoulli_alpha + sumG,
                               pairwise_beta_bernoulli_beta + no_pairwise - sumG);
        std::fill(inclusion_probability_difference.begin(),
                  inclusion_probability_difference.end(), prob);
      }

      // Main differences
      if (main_difference_prior == "Beta-Bernoulli") {
        sumG = 0;
        for (int i = 0; i < no_variables; ++i) {
          sumG += indicator(i, i);
        }
        double prob = R::rbeta(main_beta_bernoulli_alpha + sumG,
                               main_beta_bernoulli_beta + no_variables - sumG);
        for (int i = 0; i < no_variables; ++i) {
          inclusion_probability_difference(i, i) = prob;
        }
      }
    }

    // Save results after burn-in
    if (iteration >= burnin) {
      int iter_adj = iteration - burnin;
      for (int row = 0; row < no_main; ++row) {
        for (int col = 0; col < no_groups; ++col) {
          out_main(row, col) = (out_main(row, col) * iter_adj + main_effects(row, col)) /
            (iter_adj + 1);
        }
      }
      for (int row = 0; row < no_pairwise; ++row) {
        for (int col = 0; col < no_groups; ++col) {
          out_pairwise(row, col) = (out_pairwise(row, col) * iter_adj + pairwise_effects(row, col)) /
            (iter_adj + 1);
        }
      }
      for (int i = 0; i < no_variables - 1; ++i) {
        for (int j = i + 1; j < no_variables; ++j) {
          out_indicator(i, j) = (out_indicator(i, j) * iter_adj + indicator(i, j)) /
            (iter_adj + 1);
          out_indicator(j, i) = out_indicator(i, j);
        }
      }
    }
  }

  return List::create(Named("main") = out_main,
                      Named("pairwise") = out_pairwise,
                      Named("indicator") = out_indicator);
}