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

// [[Rcpp::export]]
List sample_twosample_omrf_gibbs(int no_states_gr1,
                                 int no_states_gr2,
                                 int no_variables,
                                 IntegerVector no_categories,
                                 NumericMatrix interactions,
                                 NumericMatrix thresholds,
                                 NumericMatrix main_difference,
                                 NumericMatrix pairwise_difference,
                                 NumericMatrix cross_lagged,
                                 bool paired,
                                 int iter) {

  IntegerMatrix observations_gr1(no_states_gr1, no_variables);
  IntegerMatrix observations_gr2(no_states_gr2, no_variables);
  int max_no_categories = max(no_categories);

  NumericVector probabilities(max_no_categories + 1);
  double exponent = 0.0;
  double rest_score = 0.0;
  double cumsum = 0.0;
  double u = 0.0;
  int score = 0;

  //Random (uniform) starting values -------------------------------------------
  for(int variable = 0; variable < no_variables; variable++) {
    for(int person =  0; person < no_states_gr1; person++) {
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
      observations_gr1(person, variable) = score;
    }
    for(int person =  0; person < no_states_gr2; person++) {
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
      observations_gr2(person, variable) = score;
    }
  }

  //The Gibbs sampler ----------------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    for(int variable = 0; variable < no_variables; variable++) {

      //Group 1 ----------------------------------------------------------------
      for(int person = 0; person < no_states_gr1; person++) {
        rest_score = 0.0;

        for(int vertex = 0; vertex < no_variables; vertex++) {
          rest_score +=
            observations_gr1(person, vertex) *
            interactions(vertex, variable);

          rest_score -=
            .5 *
            observations_gr1(person, vertex) *
            pairwise_difference(vertex, variable);
        }

        if(paired == true) {
          for(int vertex = 0; vertex < no_variables; vertex++) {
            rest_score +=
              .25 *
              observations_gr2(person, vertex) *
              cross_lagged (vertex, variable);
          }
        }

        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < no_categories[variable]; category++) {
          exponent = thresholds(variable, category);
          exponent -= .5 * main_difference(variable, category);
          exponent += (category + 1) * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }

        u = cumsum * R::unif_rand();

        score = 0;
        while (u > probabilities[score]) {
          score++;
        }
        observations_gr1(person, variable) = score;
      }

      //Group 2 ----------------------------------------------------------------
      for(int person =  0; person < no_states_gr2; person++) {
        rest_score = 0.0;

        for(int vertex = 0; vertex < no_variables; vertex++) {
          rest_score +=
            observations_gr2(person, vertex) *
            interactions(vertex, variable);

          rest_score +=
            .5 *
            observations_gr2(person, vertex) *
            pairwise_difference(vertex, variable);
        }

        if(paired == true) {
          for(int vertex = 0; vertex < no_variables; vertex++) {
            rest_score +=
              .25 *
              observations_gr1(person, vertex) *
              cross_lagged(vertex, variable);
          }
        }

        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < no_categories[variable]; category++) {
          exponent = thresholds(variable, category);
          exponent += .5 * main_difference(variable, category);
          exponent += (category + 1) * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }

        u = cumsum * R::unif_rand();

        score = 0;
        while (u > probabilities[score]) {
          score++;
        }
        observations_gr2(person, variable) = score;
      }
    }
    Rcpp::checkUserInterrupt();
  }

  return List::create(Named("x") = observations_gr1,
                      Named("y") = observations_gr2);
}

// [[Rcpp::export]]
List sample_twosample_bcomrf_gibbs(int no_states_gr1,
                                   int no_states_gr2,
                                   int no_variables,
                                   IntegerVector no_categories,
                                   NumericMatrix interactions,
                                   NumericMatrix thresholds,
                                   NumericMatrix main_difference,
                                   NumericMatrix pairwise_difference,
                                   NumericMatrix cross_lagged,
                                   bool paired,
                                   StringVector variable_type,
                                   IntegerVector reference_category,
                                   int iter) {

  IntegerMatrix observations_gr1(no_states_gr1, no_variables);
  IntegerMatrix observations_gr2(no_states_gr2, no_variables);
  int max_no_categories = max(no_categories);
  NumericVector probabilities(max_no_categories + 1);
  double exponent = 0.0;
  double rest_score = 0.0;
  double cumsum = 0.0;
  double u = 0.0;
  int score = 0;

  //Random (uniform) starting values -------------------------------------------
  for(int variable = 0; variable < no_variables; variable++) {
    for(int person =  0; person < no_states_gr1; person++) {
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
      observations_gr1(person, variable) = score;
    }
    for(int person =  0; person < no_states_gr2; person++) {
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
      observations_gr2(person, variable) = score;
    }
  }

  //The Gibbs sampler ----------------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    for(int variable = 0; variable < no_variables; variable++) {
      //Group 1-----------------------------------------------------------------
      for(int person =  0; person < no_states_gr1; person++) {
        rest_score = 0.0;
        for(int vertex = 0; vertex < no_variables; vertex++) {
          rest_score +=
            observations_gr1(person, vertex) *
            interactions(vertex, variable);

          rest_score -=
            .5 *
            observations_gr1(person, vertex) *
            pairwise_difference(vertex, variable);
        }

        if(paired == true) {
          for(int vertex = 0; vertex < no_variables; vertex++) {
            rest_score +=
              .25 *
              observations_gr2(person, vertex) *
              cross_lagged(vertex, variable);
          }
        }

        if(variable_type[variable] == "blume-capel") {
          cumsum = 0.0;
          for(int category = 0; category < no_categories[variable] + 1; category++) {
            //The linear term of the Blume-Capel variable
            exponent = thresholds(variable, 0) * category;
            exponent -= .5 * main_difference(variable, 0) * category;
            //The quadratic term of the Blume-Capel variable
            exponent += thresholds(variable, 1) *
              (category - reference_category[variable]) *
              (category - reference_category[variable]);
            exponent -=
              .5 *
              main_difference(variable, 1) *
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
            exponent -= .5 * main_difference(variable, category);
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
        observations_gr1(person, variable) = score;
      }

      //Group 2 ----------------------------------------------------------------
      for(int person =  0; person < no_states_gr2; person++) {
        rest_score = 0.0;
        for(int vertex = 0; vertex < no_variables; vertex++) {
          rest_score +=
            observations_gr2(person, vertex) *
            interactions(vertex, variable);
          rest_score +=
            .5 *
            observations_gr2(person, vertex) *
            pairwise_difference(vertex, variable);
        }

        if(paired == true) {
          for(int vertex = 0; vertex < no_variables; vertex++) {
            rest_score +=
              .25 *
              observations_gr1(person, vertex) *
              cross_lagged(vertex, variable);
          }
        }

        if(variable_type[variable] == "blume-capel") {
          cumsum = 0.0;
          for(int category = 0; category < no_categories[variable] + 1; category++) {
            //The linear term of the Blume-Capel variable
            exponent = thresholds(variable, 0) * category;
            exponent += .5 * main_difference(variable, 0) * category;

            //The quadratic term of the Blume-Capel variable
            exponent +=
              thresholds(variable, 1) *
              (category - reference_category[variable]) *
              (category - reference_category[variable]);
            exponent +=
              .5 *
              main_difference(variable, 1) *
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
            exponent += .5 * main_difference(variable, category);
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
        observations_gr2(person, variable) = score;
      }
    }
    Rcpp::checkUserInterrupt();
  }

  return List::create(Named("x") = observations_gr1,
                      Named("y") = observations_gr2);
}