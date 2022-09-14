#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector em_gamma(NumericMatrix interactions, 
                       NumericMatrix slab_var,
                       double theta,
                       double xi,
                       int no_persons) {
  int no_nodes = interactions.ncol();
  NumericMatrix gamma(no_nodes, no_nodes);
  double tmp1;
  double tmp2;
  
  for(int node1 = 0; node1 < no_nodes - 1; node1++) {
    for(int node2 = node1 + 1; node2 < no_nodes; node2++) {
      tmp1 = theta * R::dnorm(interactions(node1, node2), 
                              0.0,
                              std::sqrt(slab_var(node1, node2)),
                              false);
      tmp2 = (1 - theta) * R::dnorm(interactions(node1, node2), 
              0.0,
              std::sqrt(slab_var(node1, node2) * xi / no_persons),
              false);
      gamma(node1, node2) = tmp1 / (tmp1 + tmp2);
      gamma(node2, node1) = tmp1 / (tmp1 + tmp2);
    }
  }
  
  return gamma;
}


// [[Rcpp::export]]
NumericVector em_interaction_var(NumericMatrix gamma, 
                                 NumericMatrix slab_var,
                                 double theta,
                                 double xi,
                                 int no_persons) {
  int no_nodes = gamma.ncol();
  NumericMatrix interaction_var(no_nodes, no_nodes);
  double tmp1;
  double tmp2;
  
  for(int node1 = 0; node1 < no_nodes - 1; node1++) {
    for(int node2 = node1 + 1; node2 < no_nodes; node2++) {
      tmp1 = gamma(node1, node2) / slab_var(node1, node2);
      tmp2 = (1 - gamma(node1, node2)) / 
        (slab_var(node1, node2) * xi / no_persons);
      interaction_var(node1, node2) = 1 / (tmp1 + tmp2);
      interaction_var(node2, node1) = 1 / (tmp1 + tmp2);
    }
  }
  
  return interaction_var;
}



