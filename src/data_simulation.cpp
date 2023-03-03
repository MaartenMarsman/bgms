#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix sample_omrf_gibbs(int no_states,
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
  double u = 0.0;
  int score = 0;
  
  //Random (uniform) starting values -------------------------------------------
  for(int node = 0; node < no_nodes; node++) {
    for(int person =  0; person < no_states; person++) {
      cumsum = 1.0;
      probabilities[0] = 1.0;
      for(int category = 0; category < no_categories[node]; category++) {
        cumsum += 1;
        probabilities[category + 1] = cumsum;
      }
      
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
        rest_score = 0.0;
        for(int vertex = 0; vertex < no_nodes; vertex++) {
          rest_score += observations(person, vertex) *
            interactions(vertex, node);
        }
        
        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < no_categories[node]; category++) {
          exponent = thresholds(node, category); 
          exponent += (category + 1) * rest_score;
          cumsum += std::exp(exponent);
          probabilities[category + 1] = cumsum;
        }

        u = cumsum * R::unif_rand();
        
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