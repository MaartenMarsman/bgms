#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix sample_mrf_gibbs(int no_states,
                               int no_nodes,
                               IntegerVector no_categories,
                               NumericMatrix interactions,
                               NumericMatrix thresholds,
                               int iter) {

  IntegerMatrix observations(no_states, no_nodes);
  int max_no_categories = 0;
  for(int node = 0; node < no_nodes; node++) {
    if(no_categories[node] > max_no_categories) {
      max_no_categories = no_categories[node];
    }
  }
  NumericVector probabilities(max_no_categories + 1);
  double exponent = 0.0;
  double rest_score = 0.0;
  double cumsum = 0.0;
  double bound;
  double u = 0.0;
  int score = 0;

  //Random (uniform) starting values -------------------------------------------
  for(int node = 0; node < no_nodes; node++) {
    for(int person =  0; person < no_states; person++) {

      //Compute the full-conditional's probabilities per category --------------
      cumsum = 1.0;
      probabilities[0] = 1.0;
      for(int category = 0; category < no_categories[node]; category++) {
        cumsum += 1;
        probabilities[category + 1] = cumsum;
      }


      //Sample a category or score -------------------------------------------
      u = cumsum * R::unif_rand();
      score = 0;
      while (u > probabilities[score]) {
        score++;
      }

      observations(person, node) = score;
    }
  }

  //The Gibbs sampler ----------------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    for(int node = 0; node < no_nodes; node++) {
      for(int person =  0; person < no_states; person++) {

        //Compute the rest score constant in the generating probabilities ------
        rest_score = 0.0;
        for(int vertex = 0; vertex < no_nodes; vertex++) {
          rest_score += observations(person, vertex) *
            interactions(vertex, node);
        }

        //Compute a bound to keep exponents from overflowing -------------------
        if(rest_score > 0) {
          bound = no_categories[node] * rest_score;
        } else {
          bound = 0.0;
        }

        //Compute the full-conditional's probabilities per category ------------
        cumsum = std::exp(-bound);
        probabilities[0] = cumsum;
        for(int category = 0; category < no_categories[node]; category++) {
          exponent = thresholds(node, category);
          exponent += (category + 1) * rest_score;
          exponent -= bound;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }

        //Scaling the probabilities such that they sum to one ------------------
        for(int category = 0; category <= no_categories[node]; category++) {
          probabilities[category] /= cumsum;
        }

        //Sample a category or score -------------------------------------------
        u = R::unif_rand();
        score = 0;
        while (u > probabilities[score]) {
          score++;
        }

        observations(person, node) = score;
      }
    }
    Rcpp::checkUserInterrupt();
  }

  return observations;
}

// [[Rcpp::export]]
IntegerMatrix sample_panel_mrf_gibbs(int no_states,
                                     int no_nodes,
                                     int no_timepoints,
                                     IntegerVector no_categories,
                                     NumericMatrix cross_sectional_interactions,
                                     NumericMatrix cross_lagged_interactions,
                                     NumericMatrix thresholds,
                                     NumericMatrix null_interactions,
                                     NumericMatrix null_thresholds,
                                     int iter) {

  IntegerVector start(no_timepoints + 1);
  for(int t = 0; t <= no_timepoints; t++) {
    start[t] = t * no_nodes;
  }

  IntegerMatrix observations(no_states, no_nodes * (no_timepoints + 1));
  int max_no_categories = 0;
  for(int node = 0; node < no_nodes; node++) {
    if(no_categories[node] > max_no_categories) {
      max_no_categories = no_categories[node];
    }
  }
  NumericVector probabilities(max_no_categories + 1);
  double exponent = 0.0;
  double rest_score = 0.0;
  double cumsum = 0.0;
  double bound;
  double u = 0.0;
  int score = 0;

  //Random (uniform) starting values -------------------------------------------
  for(int t = 0; t <= no_timepoints; t++) {
    for(int node = 0; node < no_nodes; node++) {
      for(int person =  0; person < no_states; person++) {

        //Compute the full-conditional's probabilities per category ------------
        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < no_categories[node]; category++) {
          cumsum += 1.0;
          probabilities[category + 1] = cumsum;
        }

        //Sample a category or score -------------------------------------------
        u = cumsum * R::unif_rand();
        score = 0;
        while (u > probabilities[score]) {
          score++;
        }

        observations(person, start[t] + node) = score;
      }
    }
  }

  //The Gibbs sampler ----------------------------------------------------------
  for(int iteration = 0; iteration < iter; iteration++) {
    //We start with sampling states at t = 0 -----------------------------------
    for(int node = 0; node < no_nodes; node++) {
      for(int person =  0; person < no_states; person++) {

        //Compute the rest score constant in the generating probabilities ------
        rest_score = 0.0;
        for(int vertex = 0; vertex < no_nodes; vertex++) {
          rest_score += observations(person, vertex) *
            null_interactions(vertex, node);
        }

        //Compute a bound to keep exponents from overflowing -------------------
        if(rest_score > 0) {
          bound = no_categories[node] * rest_score;
        } else {
          bound = 0.0;
        }

        //Compute the full-conditional's probabilities per category ------------
        cumsum = std::exp(-bound);
        probabilities[0] = cumsum;
        for(int category = 0; category < no_categories[node]; category++) {
          exponent = null_thresholds(node, category);
          exponent += (category + 1) * rest_score;
          exponent -= bound;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }

        //Scaling the probabilities such that they sum to one ------------------
        for(int category = 0; category <= no_categories[node]; category++) {
          probabilities[category] /= cumsum;
        }

        //Sample a category or score -------------------------------------------
        u = R::unif_rand();
        score = 0;
        while (u > probabilities[score]) {
          score++;
        }

        observations(person, node) = score;
      }
    }
    Rcpp::checkUserInterrupt();
  }

  //Now we sample the rest of the time points in turn --------------------------
  for(int t = 1; t <= no_timepoints; t++) {
    for(int node = 0; node < no_nodes; node++) {
      for(int person =  0; person < no_states; person++) {

        //Compute the rest score constant in the generating probabilities ------
        rest_score = 0.0;
        for(int node2 = 0; node2 < no_nodes; node2++) {
          // Cross-sectional elements ------------------------------------------
          rest_score += cross_sectional_interactions(node, node2) *
            observations(person, start[t] + node2);

          // Cross-lagged elements ---------------------------------------------
          rest_score += cross_lagged_interactions(node, node2) *
            observations(person, start[t - 1] + node2);
        }

        //Compute a bound to keep exponents from overflowing -------------------
        if(rest_score > 0) {
          bound = no_categories[node] * rest_score;
        } else {
          bound = 0.0;
        }

        //Compute the full-conditional's probabilities per category ------------
        cumsum = std::exp(-bound);
        probabilities[0] = cumsum;
        for(int category = 0; category < no_categories[node]; category++) {
          exponent = thresholds(node + start[t-1], category);
          exponent += (category + 1) * rest_score;
          exponent -= bound;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }

        //Scaling the probabilities such that they sum to one ------------------
        for(int category = 0; category <= no_categories[node]; category++) {
          probabilities[category] /= cumsum;
        }

        //Sample a category or score -------------------------------------------
        u = R::unif_rand();
        score = 0;
        while (u > probabilities[score]) {
          score++;
        }

        observations(person, start[t] + node) = score;
      }
    }
    Rcpp::checkUserInterrupt();
  }

  return observations;
}
