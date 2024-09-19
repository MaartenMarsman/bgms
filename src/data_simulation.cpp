#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix sample_omrf_gibbs(int no_states,
                                int no_variables,
                                IntegerVector no_categories,
                                NumericMatrix interactions,
                                NumericMatrix thresholds,
                                int iter) {

  IntegerMatrix observations(no_states, no_variables);
  int max_no_categories = max(no_categories);
  NumericVector probabilities(max_no_categories + 1);
  double exponent = 0.0;
  double rest_score = 0.0;
  double cumsum = 0.0;
  double u = 0.0;
  int score = 0;

  //Random (uniform) starting values -------------------------------------------
  for(int variable = 0; variable < no_variables; variable++) {
    for(int person =  0; person < no_states; person++) {
      cumsum = 1.0;
      probabilities[0] = 1.0;
      for(int category = 0; category < no_categories[variable]; category++) {
        cumsum += 1;
        probabilities[category + 1] = cumsum;
      }

      u = cumsum * R::unif_rand();

      score = 0;
      while (u > probabilities[score]) {
        score++;
      }
      observations(person, variable) = score;
    }
  }

  //The Gibbs sampler ----------------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    for(int variable = 0; variable < no_variables; variable++) {
      for(int person =  0; person < no_states; person++) {
        rest_score = 0.0;
        for(int vertex = 0; vertex < no_variables; vertex++) {
          rest_score += observations(person, vertex) *
            interactions(vertex, variable);
        }

        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < no_categories[variable]; category++) {
          exponent = thresholds(variable, category);
          exponent += (category + 1) * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }

        u = cumsum * R::unif_rand();

        score = 0;
        while (u > probabilities[score]) {
          score++;
        }
        observations(person, variable) = score;
      }
    }
    Rcpp::checkUserInterrupt();
  }

  return observations;
}

// [[Rcpp::export]]
IntegerMatrix sample_bcomrf_gibbs(int no_states,
                                  int no_variables,
                                  IntegerVector no_categories,
                                  NumericMatrix interactions,
                                  NumericMatrix thresholds,
                                  StringVector variable_type,
                                  IntegerVector reference_category,
                                  int iter) {

  IntegerMatrix observations(no_states, no_variables);
  int max_no_categories = max(no_categories);
  NumericVector probabilities(max_no_categories + 1);
  double exponent = 0.0;
  double rest_score = 0.0;
  double cumsum = 0.0;
  double u = 0.0;
  int score = 0;

  //Random (uniform) starting values -------------------------------------------
  for(int variable = 0; variable < no_variables; variable++) {
    for(int person =  0; person < no_states; person++) {
      cumsum = 1.0;
      probabilities[0] = 1.0;
      for(int category = 0; category < no_categories[variable]; category++) {
        cumsum += 1;
        probabilities[category + 1] = cumsum;
      }

      u = cumsum * R::unif_rand();

      score = 0;
      while (u > probabilities[score]) {
        score++;
      }
      observations(person, variable) = score;
    }
  }

  //The Gibbs sampler ----------------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    for(int variable = 0; variable < no_variables; variable++) {
      for(int person =  0; person < no_states; person++) {
        rest_score = 0.0;
        for(int vertex = 0; vertex < no_variables; vertex++) {
          rest_score += observations(person, vertex) *
            interactions(vertex, variable);
        }

        if(variable_type[variable] == "blume-capel") {
          cumsum = 0.0;
          for(int category = 0; category < no_categories[variable] + 1; category++) {
            //The linear term of the Blume-Capel variable
            exponent = thresholds(variable, 0) * category;
            //The quadratic term of the Blume-Capel variable
            exponent += thresholds(variable, 1) *
              (category - reference_category[variable]) *
              (category - reference_category[variable]);
            //The pairwise interactions
            exponent += category * rest_score;
            cumsum += std::exp(exponent);
            probabilities[category] = cumsum;
          }
        } else {
          cumsum = 1.0;
          probabilities[0] = cumsum;
          for(int category = 0; category < no_categories[variable]; category++) {
            exponent = thresholds(variable, category);
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
        observations(person, variable) = score;
      }
    }
    Rcpp::checkUserInterrupt();
  }

  return observations;
}