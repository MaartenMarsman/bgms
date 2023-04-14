#include <Rcpp.h>
#include "utility_functions.h"
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// The c++ code below is based on the R code accompanying the paper:
//  Geng, J., Bhattacharya, A., & Pati, D. (2019). Probabilistic Community
//  Detection With Unknown Number of Communities, Journal of the American
//  Statistical Association, 114:526, 893-905, DOI:10.1080/01621459.2018.1458618
// ----------------------------------------------------------------------------|

// ----------------------------------------------------------------------------|
// Compute partition coefficient for the MFM - SBM
// ----------------------------------------------------------------------------|
NumericVector compute_Vn_mfm_sbm(int no_nodes,
                                 double dirichlet_gamma,
                                 int t_max) {
  NumericVector log_Vn(t_max);
  double r;
  double tmp;

  for(int t = 1; t <= t_max; t++) {
    r = -INFINITY;
    for(int k = t; k < 500; k++) {
      tmp = 0.0;
      for(int tt = 1 - t; tt < 1; tt ++) {
        tmp += std::log(dirichlet_gamma * k + tt);
      }
      for(int n = 0; n < no_nodes; n++) {
        tmp -= std::log(dirichlet_gamma * k + n);
      }
      tmp -= std::lgamma(k + 1);

      // Compute the maximum between r and tmp
      if (tmp > r) {
        r = std::log(std::exp(r - tmp) + 1) + tmp;
      } else {
        r = std::log(1 + std::exp(tmp - r)) + r;
      }
    }
    log_Vn(t-1) = r - log(exp(1) - 1);
  }
  return log_Vn;
}

// ----------------------------------------------------------------------------|
// Compute log-likelihood for the MFM - SBM
// ----------------------------------------------------------------------------|
double log_likelihood_mfm_sbm(IntegerVector cluster_assign,
                              NumericMatrix cluster_probs,
                              IntegerMatrix gamma,
                              int node,
                              int no_nodes) {
  double output = 0;

  for(int j = 0; j < no_nodes; j++) {
    if(j != node) {
      if(j < node) {
        output += gamma(j, node) *
          std::log(cluster_probs(cluster_assign[j], cluster_assign[node]));
        output += (1 - gamma(j, node)) *
          std::log(1 - cluster_probs(cluster_assign[j], cluster_assign[node]));
      } else {
        output += gamma(node, j) *
          std::log(cluster_probs(cluster_assign[node], cluster_assign[j]));
        output += (1 - gamma(node, j)) *
          std::log(1 - cluster_probs(cluster_assign[node], cluster_assign[j]));
      }
    }
  }

  return output;
}

// ----------------------------------------------------------------------------|
// Compute log-marginal for the MFM - SBM
// ----------------------------------------------------------------------------|
double log_marginal_mfm_sbm(IntegerVector cluster_assign,
                            IntegerMatrix gamma,
                            int node,
                            int no_nodes,
                            double beta_alpha,
                            double beta_beta) {

  int no_clusters_excl_node = max(cluster_assign);
  //*std::max_element(cluster_assign.begin(), cluster_assign.end());
  double output = 0;
  int sumG;
  int sumN;

  output -= no_clusters_excl_node * R::lbeta(beta_alpha, beta_beta);

  for(int c = 0; c < no_clusters_excl_node; c++) {
    sumG = 0;
    sumN = 0;
    for(int i = 0; i < no_nodes; i++) {
      if(cluster_assign[i] == c) {
        sumG += gamma(node, i);
        sumN++;
      }
    }
    output += R::lbeta(sumG + beta_alpha, sumN - sumG + beta_beta);
  }
  return output;
}

// ----------------------------------------------------------------------------|
// Helper function to update sumG in sample_block_probs_mfm_sbm()
// ----------------------------------------------------------------------------|
inline void update_sumG(int &sumG,
                        const IntegerVector &cluster_assign,
                        const IntegerMatrix &gamma,
                        int r,
                        int s,
                        int no_nodes) {
  for(int node1 = 0; node1 < no_nodes - 1; node1++) {
    if(cluster_assign[node1] == r) {
      for(int node2 = node1 + 1; node2 < no_nodes; node2++) {
        if(cluster_assign[node2] == s) {
          sumG += gamma(node1, node2);
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------|
// Sample the cluster assignment in sample_block_allocations_mfm_sbm()
// ----------------------------------------------------------------------------|
int sample_cluster(NumericVector cluster_prob) {
  NumericVector cum_prob = cumsum(cluster_prob);
  double u = R::runif(0, max(cum_prob));

  for (int i = 0; i < cum_prob.size(); i++) {
    if (u <= cum_prob[i]) {
      return i;
    }
  }
  return cum_prob.size();
}

// ----------------------------------------------------------------------------|
// Sample the block allocations for the MFM - SBM
// ----------------------------------------------------------------------------|
IntegerVector block_allocations_mfm_sbm(IntegerVector cluster_assign,
                                        int no_nodes,
                                        NumericVector log_Vn,
                                        NumericMatrix block_probs,
                                        IntegerMatrix gamma,
                                        int dirichlet_gamma,
                                        double beta_alpha,
                                        double beta_beta) {
  int old;
  int cluster;
  int no_clusters;
  double prob;
  double loglike;
  double logmarg;

  for(int node = 0; node < no_nodes; node++) {
    old = cluster_assign[node];

    IntegerVector cluster_size = table_cpp(cluster_assign);
    no_clusters = cluster_size.size();

    if (cluster_size[old] == 1) {
      // Singleton cluster.

      // Cluster sizes without node
      IntegerVector cluster_size_node = clone(cluster_size);
      cluster_size_node[old] -= (1 + dirichlet_gamma);

      // Compute probabilities for sampling process
      NumericVector cluster_prob(no_clusters + 1);
      for (int c = 0; c <= no_clusters; c++) {
        IntegerVector cluster_assign_tmp = clone(cluster_assign);
        cluster_assign_tmp[node] = c;

        if (c < no_clusters) {
          loglike = log_likelihood_mfm_sbm(cluster_assign_tmp,
                                           block_probs,
                                           gamma,
                                           node,
                                           no_nodes);

          prob = (dirichlet_gamma + cluster_size_node[c]) *
            std::exp(loglike);
        } else {
          logmarg = log_marginal_mfm_sbm(cluster_assign_tmp,
                                         gamma,
                                         node,
                                         no_nodes,
                                         beta_alpha,
                                         beta_beta);

          prob = dirichlet_gamma *
            std::exp(logmarg) *
            std::exp(log_Vn[no_clusters - 1] - log_Vn[no_clusters - 2]);
        }

        cluster_prob[c] = prob;
      }


      //Choose the cluster number for node
      cluster = sample_cluster(cluster_prob);

      // Remove the empty cluster
      if (cluster == no_clusters) {
        cluster_assign[node] = old;
      } else {
        cluster_assign[node] = cluster;
        for (int i = 0; i < no_nodes; i++) {
          if (cluster_assign[i] > old) {
            cluster_assign[i] -= 1;
          }
        }
        block_probs = remove_row_col_matrix(block_probs, old);
      }
    } else {
      // Cluster sizes without node
      IntegerVector cluster_size_node = clone(cluster_size);
      cluster_size_node[old] -= 1;

      // Compute probabilities for sampling process
      NumericVector cluster_prob(no_clusters + 1);
      for (int c = 0; c <= no_clusters; c++) {
        IntegerVector cluster_assign_tmp = clone(cluster_assign);
        cluster_assign_tmp[node] = c;

        if (c < no_clusters) {
          loglike = log_likelihood_mfm_sbm(cluster_assign_tmp,
                                           block_probs,
                                           gamma,
                                           node,
                                           no_nodes);

          prob = (dirichlet_gamma + cluster_size_node[c]) *
            std::exp(loglike);
        } else {
          logmarg = log_marginal_mfm_sbm(cluster_assign_tmp,
                                         gamma,
                                         node,
                                         no_nodes,
                                         beta_alpha,
                                         beta_beta);

          prob = dirichlet_gamma *
            std::exp(logmarg) *
            std::exp(log_Vn[no_clusters] - log_Vn[no_clusters - 1]);
        }

        cluster_prob[c] = prob;
      }


      //Choose the cluster number for node
      cluster = sample_cluster(cluster_prob);

      cluster_assign[node] = cluster;

      if (cluster == no_clusters) {
        block_probs = add_row_col_block_prob_matrix(block_probs,
                                                    beta_alpha,
                                                    beta_beta);
      }
    }
  }
  return cluster_assign;
}

// ----------------------------------------------------------------------------|
// Sample the block parameters for the MFM - SBM
// ----------------------------------------------------------------------------|
NumericMatrix block_probs_mfm_sbm(IntegerVector cluster_assign,
                                  IntegerMatrix gamma,
                                  int no_nodes,
                                  double beta_alpha,
                                  double beta_beta) {

  IntegerVector cluster_size = table_cpp(cluster_assign);
  int no_clusters = cluster_size.size();

  NumericMatrix block_probs(no_clusters, no_clusters);

  int sumG;
  int size;

  for(int r = 0; r < no_clusters; r++) {
    for(int s = r; s < no_clusters; s++) {
      sumG = 0;

      if(r == s) {
        update_sumG(sumG, cluster_assign, gamma, r, r, no_nodes);
        size = cluster_size[r] * (cluster_size[r] - 1) / 2;
      } else {
        update_sumG(sumG, cluster_assign, gamma, r, s, no_nodes);
        update_sumG(sumG, cluster_assign, gamma, s, r, no_nodes);
        size = cluster_size[s] * cluster_size[r];
      }

      block_probs(r, s) = R::rbeta(sumG + beta_alpha, size - sumG + beta_beta);
      block_probs(s, r) = block_probs(r, s);
    }
  }
  return block_probs;
}