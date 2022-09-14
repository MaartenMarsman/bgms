#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gradient_interactions(NumericMatrix interactions,
                                    NumericMatrix thresholds,
                                    IntegerMatrix observations,
                                    IntegerVector no_categories,
                                    NumericMatrix interaction_var) {
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  int no_parameters = no_nodes * (no_nodes - 1) / 2;
  
  NumericVector gradient (no_parameters);
  double bound = 0.0;
  double bound_s = 0.0;
  double bound_t = 0.0;
  double denominator = 0.0;
  double numerator = 0.0;
  double rest_score = 0.0;
  double exponent = 0.0;
  int score = 0;
  int counter = -1;
  
  //Gradient of interactions parameters ---------------------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      counter += 1;
      gradient[counter] = 0.0;
      
      bound_s = thresholds(s, 0);
      for(int category = 1; category < no_categories[s]; category++) {
        if(thresholds(s, category) > bound_s) {
          bound_s = thresholds(s, category);
        }
      }
      bound_t = thresholds(t, 0);
      for(int category = 1; category < no_categories[t]; category++) {
        if(thresholds(t, category) > bound_t) {
          bound_t = thresholds(t, category);
        }
      }
      
      for(int i = 0; i < no_persons; i ++) {
        //Sufficient statistic ------------------------------------------------
        gradient[counter] += 2 * 
          (no_categories[s] - observations(i, s)) * 
          (no_categories[t] - observations(i, t));
        
        //Contribution of the full-conditional of node s (pseudolikelihood) ---
        rest_score = 0.0;
        for(int node = 0; node < no_nodes; node++) {
          //interactions(s, s) = 0.0;
          rest_score += (no_categories[node] - observations(i, node)) * 
            interactions(node, s);
        }
        if(rest_score > 0) {
          bound = bound_s + no_categories[s] * rest_score;
        } else {
          bound = bound_s;
        }
        
        denominator = std::exp(-bound);
        numerator = 0.0;
        for(int category = 0; category < no_categories[s]; category++) {
          score = no_categories[s] - category;
          exponent = thresholds(s, category) + 
            score * rest_score - 
            bound;
          numerator += score * 
            (no_categories[t] - observations(i, t)) * 
            std::exp(exponent); 
          denominator += std::exp(exponent);
        }
        gradient[counter] -= numerator / denominator;
        
        //Contribution of the full-conditional of node t (pseudolikelihood) ---
        rest_score = 0.0;
        for(int node = 0; node < no_nodes; node++) {
          //interactions(s, s) = 0.0;
          rest_score += (no_categories[node] - observations(i, node)) * 
            interactions(node, t);
        }
        if(rest_score > 0) {
          bound = bound_t + no_categories[t] * rest_score;
        } else {
          bound = bound_t;
        }
        
        denominator = std::exp(-bound);
        numerator = 0.0;
        for(int category = 0; category < no_categories[t]; category++) {
          score = no_categories[t] - category;
          exponent = thresholds(t, category) + 
            score * rest_score - 
            bound;
          numerator += score * 
            (no_categories[s] - observations(i, s)) * 
            std::exp(exponent); 
          denominator += std::exp(exponent);
        }
        gradient[counter] -= numerator / denominator;
      }
      
      //Contribution of the prior density --------------------------------------
      if(std::isfinite(interaction_var(s, t)))
        gradient[counter] -= interactions(s, t) / interaction_var(s, t);
    }
  }
  return gradient;
}

// [[Rcpp::export]]
NumericVector gradient_interactions_cauchy(NumericMatrix interactions,
                                           NumericMatrix thresholds,
                                           IntegerMatrix observations,
                                           IntegerVector no_categories,
                                           double cauchy_scale) {
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  int no_parameters = no_nodes * (no_nodes - 1) / 2;
  
  NumericVector gradient (no_parameters);
  double bound = 0.0;
  double bound_s = 0.0;
  double bound_t = 0.0;
  double denominator = 0.0;
  double numerator = 0.0;
  double rest_score = 0.0;
  double exponent = 0.0;
  int score = 0;
  int counter = -1;
  
  //Gradient of interactions parameters ---------------------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      counter += 1;
      gradient[counter] = 0.0;
      
      bound_s = thresholds(s, 0);
      for(int category = 1; category < no_categories[s]; category++) {
        if(thresholds(s, category) > bound_s) {
          bound_s = thresholds(s, category);
        }
      }
      bound_t = thresholds(t, 0);
      for(int category = 1; category < no_categories[t]; category++) {
        if(thresholds(t, category) > bound_t) {
          bound_t = thresholds(t, category);
        }
      }
      
      for(int i = 0; i < no_persons; i ++) {
        //Sufficient statistic ------------------------------------------------
        gradient[counter] += 2 * 
          (no_categories[s] - observations(i, s)) * 
          (no_categories[t] - observations(i, t));
        
        //Contribution of the full-conditional of node s (pseudolikelihood) ---
        rest_score = 0.0;
        for(int node = 0; node < no_nodes; node++) {
          //interactions(s, s) = 0.0;
          rest_score += (no_categories[node] - observations(i, node)) * 
            interactions(node, s);
        }
        if(rest_score > 0) {
          bound = bound_s + no_categories[s] * rest_score;
        } else {
          bound = bound_s;
        }
        
        denominator = std::exp(-bound);
        numerator = 0.0;
        for(int category = 0; category < no_categories[s]; category++) {
          score = no_categories[s] - category;
          exponent = thresholds(s, category) + 
            score * rest_score - 
            bound;
          numerator += score * 
            (no_categories[t] - observations(i, t)) * 
            std::exp(exponent); 
          denominator += std::exp(exponent);
        }
        gradient[counter] -= numerator / denominator;
        
        //Contribution of the full-conditional of node t (pseudolikelihood) ---
        rest_score = 0.0;
        for(int node = 0; node < no_nodes; node++) {
          //interactions(s, s) = 0.0;
          rest_score += (no_categories[node] - observations(i, node)) * 
            interactions(node, t);
        }
        if(rest_score > 0) {
          bound = bound_t + no_categories[t] * rest_score;
        } else {
          bound = bound_t;
        }
        
        denominator = std::exp(-bound);
        numerator = 0.0;
        for(int category = 0; category < no_categories[t]; category++) {
          score = no_categories[t] - category;
          exponent = thresholds(t, category) + 
            score * rest_score - 
            bound;
          numerator += score * 
            (no_categories[s] - observations(i, s)) * 
            std::exp(exponent); 
          denominator += std::exp(exponent);
        }
        gradient[counter] -= numerator / denominator;
      }
      
      //Contribution of the Cauchy prior density -------------------------------
      gradient[counter] -= 2 * interactions(s, t) / 
        (interactions(s, t) * interactions(s, t) + cauchy_scale * cauchy_scale);
    }
  }
  return gradient;
}


// [[Rcpp::export]]
NumericVector gradient_thresholds(NumericMatrix interactions,
                                  NumericMatrix thresholds,
                                  IntegerMatrix observations,
                                  IntegerVector no_categories,
                                  NumericMatrix threshold_var) {
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  int no_parameters = 0;
  for(int node = 0; node < no_nodes; node++)
    no_parameters += no_categories[node];
  
  NumericVector gradient (no_parameters);
  double bound = 0.0;
  double bound_s = 0.0;
  double denominator = 0.0;
  double numerator = 0.0;
  double rest_score = 0.0;
  double exponent = 0.0;
  int score = 0;
  int counter = -1;
  int observed_category = 0;
  
  //Gradient of thresholds parameters -----------------------------------------
  for(int s = 0; s < no_nodes; s++) {
    counter += 1;
    //Contribution of the prior distribution ----------------------------------
    for(int category = 0; category < no_categories[s]; category++) {
      if(std::isfinite(threshold_var(s, category)))
        gradient[counter + category] = 
          -thresholds(s, category) / threshold_var(s, category);
    }
    //Contribution of the full-conditional of node s (pseudolikelihood) -------
    bound_s = thresholds(s, 0);
    for(int category = 1; category < no_categories[s]; category++) {
      if(thresholds(s, category) > bound_s) {
        bound_s = thresholds(s, category);
      }
    }
    for(int person = 0; person < no_persons; person++) {
      observed_category = observations(person, s);
      if(observed_category < no_categories[s]) {
        gradient[counter + observed_category] += 1;//observed_category
      }
      rest_score = 0.0;
      for(int node = 0; node < no_nodes; node++) {
        rest_score += (no_categories[node] - observations(person, node)) * 
          interactions(node, s);  
      }
      bound = bound_s + no_categories[s] * rest_score;
      
      denominator = std::exp(-bound);
      for(int category = 0; category < no_categories[s]; category++) {
        score = no_categories[s] - category;
        exponent = thresholds(s, category) + 
          score * rest_score - 
          bound;
        denominator += std::exp(exponent);
      }
      for(int category = 0; category < no_categories[s]; category++) {
        score = no_categories[s] - category;
        exponent = thresholds(s, category) +  
          score * rest_score - 
          bound;
        numerator = std::exp(exponent); 
        gradient[counter + category] -= numerator / denominator;  
      }
    }
    counter += no_categories[s] - 1;
  }
  return gradient;
}

// [[Rcpp::export]]
NumericMatrix hessian_interactions(NumericMatrix interactions,
                                   NumericMatrix thresholds,
                                   IntegerMatrix observations,
                                   IntegerVector no_categories,
                                   NumericMatrix interaction_var) {
  int no_nodes = observations.ncol();                   
  int no_persons = observations.nrow();                 
  int no_parameters = no_nodes * (no_nodes - 1) / 2;
  
  NumericVector node_bound (no_nodes);
  NumericMatrix hessian (no_parameters, no_parameters);
  IntegerMatrix index (no_parameters, no_parameters);
  int score = 0;
  int counter = 0;
  int score_weight = 0;
  double rest_score = 0.0;
  double bound = 0.0;
  double denominator = 0.0;
  double exponent = 0.0;
  double expected_square = 0.0;
  double expected_value = 0.0;
  
  //create indexing matrix for partial derivatives
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      index(s, t) = counter;
      index(t, s) = counter;
      counter++;
    }
  }
  
  //initial bounds
  for(int node = 0; node < no_nodes; node++) {
    node_bound(node) = thresholds(node, 0);
    for(int category = 1; category < no_categories[node]; category++) {
      if(thresholds(node, category) > node_bound(node)) {
        node_bound(node) = thresholds(node, category);
      }
    }
  }
  
  // compute hessian ----------------------------------------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      //hessian(par, par)
      for(int person = 0; person < no_persons; person++) {
        //Contribution of the full-conditional of s (pseudolikelihood) --------
        bound = node_bound[s];
        rest_score = 0.0;
        for(int node = 0; node < no_nodes; node++) {
          //interactions(s, s) = 0.0;
          rest_score += (no_categories[node] - observations(person, node)) * 
            interactions(s, node);
        }
        if(rest_score > 0) 
          bound += no_categories[s] * rest_score;
        
        denominator = std::exp(-bound);
        expected_square = 0.0;
        expected_value = 0.0;
        for(int category = 0; category < no_categories[s]; category++) {
          score = no_categories[s] - category;
          exponent = thresholds(s, category) + 
            score * rest_score - 
            bound;
          denominator += std::exp(exponent);
          
          score *= (no_categories[t] - observations(person, t));
          expected_square += score * score * std::exp(exponent); 
          expected_value += score * std::exp(exponent); 
        }
        
        expected_square /= denominator; //E([SCORE * (J - obs[t])] ^ 2)
        expected_value /= denominator; //E(SCORE * (J - obs[t]))
        
        hessian(index(s, t), index(s, t)) -= expected_square;
        hessian(index(s, t), index(s, t)) += expected_value * expected_value;
        
        //Contribution of the full-conditional of t (pseudolikelihood) ---
        bound = node_bound[t];
        rest_score = 0.0;
        for(int node = 0; node < no_nodes; node++) {
          //interactions(s, s) = 0.0;
          rest_score += (no_categories[node] - observations(person, node)) * 
            interactions(t, node);
        }
        if(rest_score > 0)
          bound += no_categories[t] * rest_score;
        
        denominator = std::exp(-bound);
        expected_square = 0.0;
        expected_value = 0.0;
        for(int category = 0; category < no_categories[t]; category++) {
          score = no_categories[t] - category;
          exponent = thresholds(t, category) + 
            score * rest_score - 
            bound;
          denominator += std::exp(exponent);
          
          score *= (no_categories[s] - observations(person, s));
          expected_square += score * score * std::exp(exponent); 
          
          expected_value += score * std::exp(exponent); 
        }
        expected_square /= denominator; //E([SCORE * (J - obs[s])] ^ 2)
        expected_value /= denominator;  //E([SCORE * (J - obs[s])])
        
        hessian(index(s, t), index(s, t)) -= expected_square;
        hessian(index(s, t), index(s, t)) += expected_value * expected_value;
        
        //Partial derivatives ---------------------------------------------------     
        if(t < no_nodes - 1) {
          for(int h = t + 1; h < no_nodes; h++) {
            //Contribution to d^2 / d[s,t]d[s,h] = d^2 / d[s,t] d[h,s]
            bound = node_bound[s];
            rest_score = 0.0;
            for(int node = 0; node < no_nodes; node++) {
              //interactions(s, s) = 0.0;
              rest_score += (no_categories[node] - observations(person, node)) * 
                interactions(s, node);
            }
            if(rest_score > 0) 
              bound += no_categories[s] * rest_score;
            
            denominator = std::exp(-bound);
            score_weight = (no_categories[t] - observations(person, t)) *
              (no_categories[h] - observations(person, h));
            expected_square = 0.0;
            expected_value = 0.0;
            for(int category = 0; category < no_categories[s]; category++) {
              score = no_categories[s] - category;
              exponent = thresholds(s, category) + 
                score * rest_score - 
                bound;
              denominator += std::exp(exponent);
              expected_square += score * score * std::exp(exponent); 
              expected_value += score * std::exp(exponent); 
            }
            
            expected_square /= denominator; //E([SCORE * (J - obs[t])] ^ 2)
            expected_value /= denominator; //E(SCORE * (J - obs[t]))
            
            //Contribution to d^2 / d[s,t]d[s,h] = d^2 / d[s,t] d[h,s]
            hessian(index(s, t), index(s, h)) += 
              score_weight * expected_value * expected_value;
            hessian(index(s, t), index(s, h)) -= 
              score_weight * expected_square;
            hessian(index(s, h), index(s, t)) += 
              score_weight * expected_value * expected_value;
            hessian(index(s, h), index(s, t)) -= 
              score_weight * expected_square;
            
            //Contribution to d^2 / d[s,t]d[t,h] = d^2 / d[s,t] d[h,t]
            bound = node_bound[t];
            rest_score = 0.0;
            for(int node = 0; node < no_nodes; node++) {
              //interactions(t, t) = 0.0;
              rest_score += (no_categories[node] - observations(person, node)) * 
                interactions(t, node);
            }
            if(rest_score > 0) 
              bound += no_categories[t] * rest_score;
            
            denominator = std::exp(-bound);
            score_weight = (no_categories[s] - observations(person, s)) *
              (no_categories[h] - observations(person, h));
            expected_square = 0.0;
            expected_value = 0.0;
            for(int category = 0; category < no_categories[t]; category++) {
              score = no_categories[t] - category;
              exponent = thresholds(t, category) + 
                score * rest_score - 
                bound;
              denominator += std::exp(exponent);
              expected_square += score * score * std::exp(exponent); 
              expected_value += score * std::exp(exponent); 
            }
            
            expected_square /= denominator; //E([SCORE] ^ 2)
            expected_value /= denominator; //E(SCORE)
            
            //Contribution to d^2 / d[s,t]d[t,h] = d^2 / d[s,t] d[h,t]
            
            hessian(index(s, t), index(t, h)) += 
              score_weight * expected_value * expected_value;
            hessian(index(s, t), index(t, h)) -= 
              score_weight * expected_square; 
            hessian(index(t, h), index(s, t)) += 
              score_weight * expected_value * expected_value;
            hessian(index(t, h), index(s, t)) -= 
              score_weight * expected_square; 
          }
        }
        if(t > s + 1) {
          for(int h = s + 1; h < t; h++) {
            //Contribution to d^2 / d[s,t]d[h,t] = d^2 / d[h,t] d[s,t]
            bound = node_bound[t];
            rest_score = 0.0;
            for(int node = 0; node < no_nodes; node++) {
              //interactions(s, s) = 0.0;
              rest_score += (no_categories[node] - observations(person, node)) * 
                interactions(t, node);
            }
            if(rest_score > 0) 
              bound += no_categories[t] * rest_score;
            
            denominator = std::exp(-bound);
            score_weight = (no_categories[s] - observations(person, s)) *
              (no_categories[h] - observations(person, h));
            expected_square = 0.0;
            expected_value = 0.0;
            for(int category = 0; category < no_categories[t]; category++) {
              score = no_categories[t] - category;
              exponent = thresholds(t, category) + 
                score * rest_score - 
                bound;
              denominator += std::exp(exponent);
              expected_square += score * score * std::exp(exponent); 
              expected_value += score * std::exp(exponent); 
            }
            
            expected_square /= denominator; //E([SCORE * (J - obs[t])] ^ 2)
            expected_value /= denominator; //E(SCORE * (J - obs[t]))
            
            //Contribution to d^2 / d[s,t]d[h,t] = d^2 / d[h,t] d[s,t]
            hessian(index(s, t), index(h, t)) += 
              score_weight * expected_value * expected_value;
            hessian(index(s, t), index(h, t)) -= 
              score_weight * expected_square;
            hessian(index(h, t), index(s, t)) += 
              score_weight * expected_value * expected_value;
            hessian(index(h, t), index(s, t)) -= 
              score_weight * expected_square;
          }
        } 
      }
      
      //Contribution of the prior density -------------------------------------
      if(std::isfinite(interaction_var(s, t)))
        hessian(index(s, t), index(s, t)) -= 1 / interaction_var(s, t);
    }
  }
  
  return hessian;  
}

// [[Rcpp::export]]
NumericMatrix hessian_interactions_cauchy(NumericMatrix interactions,
                                          NumericMatrix thresholds,
                                          IntegerMatrix observations,
                                          IntegerVector no_categories,
                                          double cauchy_scale) {
  int no_nodes = observations.ncol();                   
  int no_persons = observations.nrow();                 
  int no_parameters = no_nodes * (no_nodes - 1) / 2;
  
  NumericVector node_bound (no_nodes);
  NumericMatrix hessian (no_parameters, no_parameters);
  IntegerMatrix index (no_parameters, no_parameters);
  int score = 0;
  int counter = 0;
  int score_weight = 0;
  double rest_score = 0.0;
  double bound = 0.0;
  double denominator = 0.0;
  double exponent = 0.0;
  double expected_square = 0.0;
  double expected_value = 0.0;
  
  //create indexing matrix for partial derivatives
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      index(s, t) = counter;
      index(t, s) = counter;
      counter++;
    }
  }
  
  //initial bounds
  for(int node = 0; node < no_nodes; node++) {
    node_bound(node) = thresholds(node, 0);
    for(int category = 1; category < no_categories[node]; category++) {
      if(thresholds(node, category) > node_bound(node)) {
        node_bound(node) = thresholds(node, category);
      }
    }
  }
  
  // compute hessian ----------------------------------------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      //hessian(par, par)
      for(int person = 0; person < no_persons; person++) {
        //Contribution of the full-conditional of s (pseudolikelihood) --------
        bound = node_bound[s];
        rest_score = 0.0;
        for(int node = 0; node < no_nodes; node++) {
          //interactions(s, s) = 0.0;
          rest_score += (no_categories[node] - observations(person, node)) * 
            interactions(s, node);
        }
        if(rest_score > 0) 
          bound += no_categories[s] * rest_score;
        
        denominator = std::exp(-bound);
        expected_square = 0.0;
        expected_value = 0.0;
        for(int category = 0; category < no_categories[s]; category++) {
          score = no_categories[s] - category;
          exponent = thresholds(s, category) + 
            score * rest_score - 
            bound;
          denominator += std::exp(exponent);
          
          score *= (no_categories[t] - observations(person, t));
          expected_square += score * score * std::exp(exponent); 
          expected_value += score * std::exp(exponent); 
        }
        
        expected_square /= denominator; //E([SCORE * (J - obs[t])] ^ 2)
        expected_value /= denominator; //E(SCORE * (J - obs[t]))
        
        hessian(index(s, t), index(s, t)) -= expected_square;
        hessian(index(s, t), index(s, t)) += expected_value * expected_value;
        
        //Contribution of the full-conditional of t (pseudolikelihood) ---
        bound = node_bound[t];
        rest_score = 0.0;
        for(int node = 0; node < no_nodes; node++) {
          //interactions(s, s) = 0.0;
          rest_score += (no_categories[node] - observations(person, node)) * 
            interactions(t, node);
        }
        if(rest_score > 0)
          bound += no_categories[t] * rest_score;
        
        denominator = std::exp(-bound);
        expected_square = 0.0;
        expected_value = 0.0;
        for(int category = 0; category < no_categories[t]; category++) {
          score = no_categories[t] - category;
          exponent = thresholds(t, category) + 
            score * rest_score - 
            bound;
          denominator += std::exp(exponent);
          
          score *= (no_categories[s] - observations(person, s));
          expected_square += score * score * std::exp(exponent); 
          
          expected_value += score * std::exp(exponent); 
        }
        expected_square /= denominator; //E([SCORE * (J - obs[s])] ^ 2)
        expected_value /= denominator;  //E([SCORE * (J - obs[s])])
        
        hessian(index(s, t), index(s, t)) -= expected_square;
        hessian(index(s, t), index(s, t)) += expected_value * expected_value;
        
        //Partial derivatives ---------------------------------------------------     
        if(t < no_nodes - 1) {
          for(int h = t + 1; h < no_nodes; h++) {
            //Contribution to d^2 / d[s,t]d[s,h] = d^2 / d[s,t] d[h,s]
            bound = node_bound[s];
            rest_score = 0.0;
            for(int node = 0; node < no_nodes; node++) {
              //interactions(s, s) = 0.0;
              rest_score += (no_categories[node] - observations(person, node)) * 
                interactions(s, node);
            }
            if(rest_score > 0) 
              bound += no_categories[s] * rest_score;
            
            denominator = std::exp(-bound);
            score_weight = (no_categories[t] - observations(person, t)) *
              (no_categories[h] - observations(person, h));
            expected_square = 0.0;
            expected_value = 0.0;
            for(int category = 0; category < no_categories[s]; category++) {
              score = no_categories[s] - category;
              exponent = thresholds(s, category) + 
                score * rest_score - 
                bound;
              denominator += std::exp(exponent);
              expected_square += score * score * std::exp(exponent); 
              expected_value += score * std::exp(exponent); 
            }
            
            expected_square /= denominator; //E([SCORE * (J - obs[t])] ^ 2)
            expected_value /= denominator; //E(SCORE * (J - obs[t]))
            
            //Contribution to d^2 / d[s,t]d[s,h] = d^2 / d[s,t] d[h,s]
            hessian(index(s, t), index(s, h)) += 
              score_weight * expected_value * expected_value;
            hessian(index(s, t), index(s, h)) -= 
              score_weight * expected_square;
            hessian(index(s, h), index(s, t)) += 
              score_weight * expected_value * expected_value;
            hessian(index(s, h), index(s, t)) -= 
              score_weight * expected_square;
            
            //Contribution to d^2 / d[s,t]d[t,h] = d^2 / d[s,t] d[h,t]
            bound = node_bound[t];
            rest_score = 0.0;
            for(int node = 0; node < no_nodes; node++) {
              //interactions(t, t) = 0.0;
              rest_score += (no_categories[node] - observations(person, node)) * 
                interactions(t, node);
            }
            if(rest_score > 0) 
              bound += no_categories[t] * rest_score;
            
            denominator = std::exp(-bound);
            score_weight = (no_categories[s] - observations(person, s)) *
              (no_categories[h] - observations(person, h));
            expected_square = 0.0;
            expected_value = 0.0;
            for(int category = 0; category < no_categories[t]; category++) {
              score = no_categories[t] - category;
              exponent = thresholds(t, category) + 
                score * rest_score - 
                bound;
              denominator += std::exp(exponent);
              expected_square += score * score * std::exp(exponent); 
              expected_value += score * std::exp(exponent); 
            }
            
            expected_square /= denominator; //E([SCORE] ^ 2)
            expected_value /= denominator; //E(SCORE)
            
            //Contribution to d^2 / d[s,t]d[t,h] = d^2 / d[s,t] d[h,t]
            
            hessian(index(s, t), index(t, h)) += 
              score_weight * expected_value * expected_value;
            hessian(index(s, t), index(t, h)) -= 
              score_weight * expected_square; 
            hessian(index(t, h), index(s, t)) += 
              score_weight * expected_value * expected_value;
            hessian(index(t, h), index(s, t)) -= 
              score_weight * expected_square; 
          }
        }
        if(t > s + 1) {
          for(int h = s + 1; h < t; h++) {
            //Contribution to d^2 / d[s,t]d[h,t] = d^2 / d[h,t] d[s,t]
            bound = node_bound[t];
            rest_score = 0.0;
            for(int node = 0; node < no_nodes; node++) {
              //interactions(s, s) = 0.0;
              rest_score += (no_categories[node] - observations(person, node)) * 
                interactions(t, node);
            }
            if(rest_score > 0) 
              bound += no_categories[t] * rest_score;
            
            denominator = std::exp(-bound);
            score_weight = (no_categories[s] - observations(person, s)) *
              (no_categories[h] - observations(person, h));
            expected_square = 0.0;
            expected_value = 0.0;
            for(int category = 0; category < no_categories[t]; category++) {
              score = no_categories[t] - category;
              exponent = thresholds(t, category) + 
                score * rest_score - 
                bound;
              denominator += std::exp(exponent);
              expected_square += score * score * std::exp(exponent); 
              expected_value += score * std::exp(exponent); 
            }
            
            expected_square /= denominator; //E([SCORE * (J - obs[t])] ^ 2)
            expected_value /= denominator; //E(SCORE * (J - obs[t]))
            
            //Contribution to d^2 / d[s,t]d[h,t] = d^2 / d[h,t] d[s,t]
            hessian(index(s, t), index(h, t)) += 
              score_weight * expected_value * expected_value;
            hessian(index(s, t), index(h, t)) -= 
              score_weight * expected_square;
            hessian(index(h, t), index(s, t)) += 
              score_weight * expected_value * expected_value;
            hessian(index(h, t), index(s, t)) -= 
              score_weight * expected_square;
          }
        } 
      }
      
      //Contribution of the prior density --------------------------------------
      double tmp_m =  cauchy_scale * cauchy_scale - 
        interactions(s, t) * interactions(s, t);
      double tmp_p = interactions(s, t) * interactions(s, t) + 
        cauchy_scale * cauchy_scale;
      hessian(index(s, t), index(s, t)) += 2 * tmp_m / (tmp_p * tmp_p);
    }
  }
  
  return hessian;  
}


// [[Rcpp::export]]
NumericMatrix hessian_thresholds(NumericMatrix interactions,
                                 NumericMatrix thresholds,
                                 IntegerMatrix observations,
                                 IntegerVector no_categories,
                                 NumericMatrix threshold_var) {
  int no_nodes = observations.ncol();                   
  int no_persons = observations.nrow();                 
  int no_parameters = 0;
  for(int node = 0; node < no_nodes; node++)
    no_parameters += no_categories[node];
  
  NumericVector node_bound (no_nodes);
  NumericMatrix hessian (no_parameters, no_parameters);
  int score = 0;
  int counter = -1;
  double rest_score = 0.0;
  double bound = 0.0;
  double bound_s = 0.0;
  double denominator = 0.0;
  double exponent = 0.0;
  double numerator = 0.0;
  double prob = 0.0;
  
  //Gradient of thresholds parameters -----------------------------------------
  for(int s = 0; s < no_nodes; s++) {
    counter += 1;
    
    //Contribution of the full-conditional of node s (pseudolikelihood) -------
    bound_s = thresholds(s, 0);
    for(int category = 1; category < no_categories[s]; category++) {
      if(thresholds(s, category) > bound_s) {
        bound_s = thresholds(s, category);
      }
    }
    
    for(int person = 0; person < no_persons; person++) {
      bound = bound_s;
      rest_score = 0.0;
      for(int node = 0; node < no_nodes; node++) {
        rest_score += (no_categories[node] - observations(person, node)) * 
          interactions(node, s);  
      }
      if(rest_score > 0) 
        bound += no_categories[s] * rest_score;
      
      denominator = std::exp(-bound);
      for(int category = 0; category < no_categories[s]; category++) {
        score = no_categories[s] - category;
        exponent = thresholds(s, category) + 
          score * rest_score - 
          bound;
        denominator += std::exp(exponent);
      }
      for(int category = 0; category < no_categories[s]; category++) {
        score = no_categories[s] - category;
        exponent = thresholds(s, category) +  
          score * rest_score - 
          bound;
        numerator = std::exp(exponent); 
        prob = numerator / denominator;
        hessian(counter + category, counter + category) -= prob;  
        hessian(counter + category, counter + category) += prob * prob;  
      }
      
      for(int j = 0; j < no_categories[s] - 1; j++) {
        for(int h = j + 1; h < no_categories[s]; h++) {
          score = no_categories[s] - j;
          exponent = thresholds(s, j) +  
            score * rest_score - 
            bound;
          
          numerator = std::exp(exponent);
          
          score = no_categories[s] - h;
          exponent = thresholds(s, h) +  
            score * rest_score - 
            bound;
          
          numerator *= std::exp(exponent);
          
          hessian(counter + j, counter + h) += numerator / 
            (denominator * denominator);  
          hessian(counter + h, counter + j) += numerator / 
            (denominator * denominator);  
        }
      }
    }
    for(int category = 0; category < no_categories[s]; category++) {
      if(std::isfinite(threshold_var(s, category))) {
        hessian(counter + category, counter + category) -= 
          1 / threshold_var(s, category);  
      }
    }
    counter += no_categories[s] - 1;
  }
  return hessian;  
}

// [[Rcpp::export]]
NumericMatrix hessian_crossparameters(NumericMatrix interactions,
                                      NumericMatrix thresholds,
                                      IntegerMatrix observations,
                                      IntegerVector no_categories) {
  int no_nodes = observations.ncol();                   
  int no_persons = observations.nrow();
  int no_thresholds = 0;
  for(int node = 0; node < no_nodes; node++)
    no_thresholds += no_categories[node];
  int no_interactions = no_nodes * (no_nodes - 1) / 2;
  NumericVector node_bound (no_nodes);
  NumericMatrix hessian (no_interactions, no_thresholds);
  IntegerMatrix index_interactions (no_interactions, no_interactions);
  IntegerMatrix index_thresholds (no_thresholds, no_thresholds);
  int counter = 0;
  int score = 0;
  double rest_score = 0.0;
  double bound = 0.0;
  double denominator = 0.0;
  double exponent = 0.0;
  double probability = 0.0;
  double expected_value = 0.0;
  
  //create indexing matrix for interactions
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      index_interactions(s, t) = counter;
      index_interactions(t, s) = counter;
      counter++;
    }
  }
  //create indexing matrix for thresholds
  counter = 0;
  for(int s = 0; s < no_nodes; s++) {
    for(int cat = 0; cat < no_categories[s]; cat++) {
      index_thresholds(s, cat) = counter;
      counter++;
    }
  }
  
  //initial bounds
  for(int node = 0; node < no_nodes; node++) {
    node_bound(node) = thresholds(node, 0);
    for(int category = 1; category < no_categories[node]; category++) {
      if(thresholds(node, category) > node_bound(node)) {
        node_bound(node) = thresholds(node, category);
      }
    }
  }
  
  //
  for(int s = 0; s < no_nodes; s++) {
    for(int cat = 0; cat < no_categories[s]; cat++) {
      for(int t = 0; t < no_nodes; t++) {
        if(t != s) {
          for(int person = 0; person < no_persons; person++) {
            //Contribution of the full-conditional of s (pseudolikelihood) ----
            bound = node_bound[s];
            rest_score = 0.0;
            for(int node = 0; node < no_nodes; node++) {
              //interactions(s, s) = 0.0;
              rest_score += (no_categories[node] - observations(person, node)) * 
                interactions(s, node);
            }
            if(rest_score > 0) 
              bound += no_categories[s] * rest_score;
            
            denominator = std::exp(-bound);
            probability = 0.0;
            expected_value = 0.0;
            for(int category = 0; category < no_categories[s]; category++) {
              score = no_categories[s] - category;
              exponent = thresholds(s, category) + 
                score * rest_score - 
                bound;
              denominator += std::exp(exponent);
              if(category == cat) {
                probability = std::exp(exponent);
              }
              expected_value += score * std::exp(exponent); 
            }
            probability /= denominator; //p(Score  = j)
            expected_value /= denominator; //E(SCORE)
            
            hessian(index_interactions(s, t), 
                    index_thresholds(s, cat)) -= 
                      (no_categories[t] - observations(person, t)) *
                      probability * ((no_categories[s] - cat)  - expected_value);
          }
        }
      }
    }
  }
  
  return hessian;  
}