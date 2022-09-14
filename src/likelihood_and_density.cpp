#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double joint_log_density(NumericMatrix interactions,
                         NumericMatrix thresholds,
                         IntegerMatrix observations,
                         IntegerVector no_categories,
                         NumericMatrix interaction_var,
                         NumericMatrix threshold_var) {
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  double rest_score = 0.0;
  double bound =  0.0;
  double density = 0.0;
  double denominator = 0.0;
  double exponent = 0.0;
  int score = 0;
  
  //Contributions of the full-conditional of nodes (pseudolikelihoods) --------
  for(int s = 0; s <  no_nodes; s++) {
    //Numerator of full-conditional of node s (pseudolikelihood) --------------
    for(int person = 0; person < no_persons; person++) {
      rest_score = 0.0;
      for(int node = 0; node < no_nodes; node++) {
        rest_score += (no_categories[node] - observations(person, node)) *
          interactions(node, s);  
      }
      density +=  (no_categories[s] - observations(person,s)) * 
        rest_score;
      bound = no_categories[s] * rest_score;
      density -= bound;
      denominator = std::exp(-bound);
      for(int category = 0; category < no_categories[s]; category++) {
        if(observations(person, s) == category) {
          density += thresholds(s, category);
        }
        score = no_categories[s] - category;
        exponent = thresholds(s, category) + 
          score * rest_score - 
          bound;
        denominator += std::exp(exponent);
      }
      //Denominator of full-conditional of node s (pseudolikelihood) ----------
      density -= log(denominator);
    }
  }
  //Contribution of the prior densities (interactions) ------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      //If Var == Inf, density is pseudolikelihood ----------------------------
      if(std::isfinite(interaction_var(s, t)))
        density += R::dnorm(interactions(s, t), 
                            0.0, 
                            std::sqrt(interaction_var(s, t)), 
                            true);
    }
  }
  
  //Contribution of the prior densities (thresholds) --------------------------
  for(int s = 0; s < no_nodes; s++) {
    for(int category = 0; category < no_categories[s]; category++) {
      //If Var == Inf, density is pseudolikelihood ----------------------------
      if(std::isfinite(threshold_var(s, category)))
        density += R::dnorm(thresholds(s, category), 
                            0.0, 
                            std::sqrt(threshold_var(s, category)), 
                            true);
    }
  }
  return density;
}

// [[Rcpp::export]]
double joint_log_density_cauchy(NumericMatrix interactions,
                                NumericMatrix thresholds,
                                IntegerMatrix observations,
                                NumericMatrix threshold_var,
                                double cauchy_scale,
                                IntegerVector no_categories) {
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  double rest_score = 0.0;
  double bound =  0.0;
  double density = 0.0;
  double denominator = 0.0;
  double exponent = 0.0;
  int score = 0;
  
  //Contributions of the full-conditional of nodes (pseudolikelihoods) --------
  for(int s = 0; s <  no_nodes; s++) {
    //Numerator of full-conditional of node s (pseudolikelihood) --------------
    for(int person = 0; person < no_persons; person++) {
      rest_score = 0.0;
      for(int node = 0; node < no_nodes; node++) {
        rest_score += (no_categories[node] - observations(person, node)) *
          interactions(node, s);  
      }
      density +=  (no_categories[s] - observations(person,s)) * 
        rest_score;
      bound = no_categories[s] * rest_score;
      density -= bound;
      denominator = std::exp(-bound);
      for(int category = 0; category < no_categories[s]; category++) {
        if(observations(person, s) == category) {
          density += thresholds(s, category);
        }
        score = no_categories[s] - category;
        exponent = thresholds(s, category) + 
          score * rest_score - 
          bound;
        denominator += std::exp(exponent);
      }
      //Denominator of full-conditional of node s (pseudolikelihood) ----------
      density -= log(denominator);
    }
  }
  //Contribution of the prior densities (interactions) ------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      //If Var == Inf, density is pseudolikelihood ----------------------------
      density += R::dcauchy(interactions(s, t), 
                            0.0, 
                            cauchy_scale, 
                            true);
    }
  }
  
  //Contribution of the prior densities (thresholds) --------------------------
  for(int s = 0; s < no_nodes; s++) {
    for(int category = 0; category < no_categories[s]; category++) {
      //If Var == Inf, density is pseudolikelihood ----------------------------
      if(std::isfinite(threshold_var(s, category)))
        density += R::dnorm(thresholds(s, category), 
                            0.0, 
                            std::sqrt(threshold_var(s, category)), 
                            true);
    }
  }
  return density;
}


// [[Rcpp::export]]
double emvs_joint_log_density(NumericMatrix interactions,
                              NumericMatrix thresholds,
                              IntegerMatrix observations,
                              IntegerVector no_categories,
                              double xi,
                              NumericMatrix slab_var,
                              NumericMatrix threshold_var,
                              double theta = 0.5,
                              bool hierarchical = false,
                              double alpha = 1.0,
                              double beta = 1.0) {
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  double rest_score = 0.0;
  double bound =  0.0;
  double density = 0.0;
  double denominator = 0.0;
  double exponent = 0.0;
  int score = 0;
  
  //Contributions of the full-conditional of nodes (pseudolikelihoods) --------
  for(int s = 0; s <  no_nodes; s++) {
    //Numerator of full-conditional of node s (pseudolikelihood) --------------
    for(int person = 0; person < no_persons; person++) {
      rest_score = 0.0;
      for(int node = 0; node < no_nodes; node++) {
        rest_score += (no_categories[node] - observations(person, node)) *
          interactions(node, s);  
      }
      density +=  (no_categories[s] - observations(person,s)) * 
        rest_score;
      bound = no_categories[s] * rest_score;
      density -= bound;
      denominator = std::exp(-bound);
      for(int category = 0; category < no_categories[s]; category++) {
        if(observations(person, s) == category) {
          density += thresholds(s, category);
        }
        score = no_categories[s] - category;
        exponent = thresholds(s, category) + 
          score * rest_score - 
          bound;
        denominator += std::exp(exponent);
      }
      //Denominator of full-conditional of node s (pseudolikelihood) ----------
      density -= log(denominator);
    }
  }
  //Contribution of the prior densities (interactions) ------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      //If Var == Inf, density is pseudolikelihood ----------------------------
      density += std::log(
        theta * R::dnorm(interactions(s, t), 
                         0.0, 
                         std::sqrt(slab_var(s, t)), 
                         false) + 
                           (1 - theta) * R::dnorm(interactions(s, t), 
                            0.0, 
                            std::sqrt(slab_var(s, t) * xi / no_persons), 
                            false));
    }
  }
  
  //Contribution of the prior densities (thresholds) --------------------------
  for(int s = 0; s < no_nodes; s++) {
    for(int category = 0; category < no_categories[s]; category++) {
      //If Var == Inf, density is pseudolikelihood ----------------------------
      if(std::isfinite(threshold_var(s, category)))
        density += R::dnorm(thresholds(s, category), 
                            0.0, 
                            std::sqrt(threshold_var(s, category)), 
                            true);
    }
  }
  
  if(hierarchical == true)
    density += R::dbeta(theta, alpha, beta, true);
  
  return density;
}
