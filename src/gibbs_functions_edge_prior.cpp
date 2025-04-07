#include <RcppArmadillo.h>

#include <Rcpp.h>
using namespace Rcpp;

// ----------------------------------------------------------------------------|
// The c++ code below is based on the R code accompanying the paper:
//  Geng, J., Bhattacharya, A., & Pati, D. (2019). Probabilistic Community
//  Detection With Unknown Number of Communities, Journal of the American
//  Statistical Association, 114:526, 893-905, DOI:10.1080/01621459.2018.1458618
// ----------------------------------------------------------------------------|

// ----------------------------------------------------------------------------|
// A c++ version of table
// ----------------------------------------------------------------------------|
arma::uvec table_cpp(arma::uvec x) {
  arma::uword n = x.n_elem;
  arma::uword m = arma::max(x);
  arma::uvec counts(m+1,arma::fill::zeros);

  for (arma::uword i = 0; i < n; i++) {
    counts(x(i))++;
  }

  return counts;
}

// ----------------------------------------------------------------------------|
// Remove row i and column i from a matrix
// ----------------------------------------------------------------------------|
arma::mat remove_row_col_matrix(arma::mat X, arma::uword i) {
  arma::uword n = X.n_rows;

  // Handle special case of matrix with only two rows and columns
  if (n == 2) {
    if(i == 0) {
      X = X(arma::span(1, 1), arma::span(1, 1));
    } else {
      X = X(arma::span(0, 0), arma::span(0, 0));
    }
    return X;
  }

  // Remove row i
  for (arma::uword j = i; j < n - 1; j++) {
    X.row(j) = X(j+1);  // Shift all rows below i up by one
  }
  X = X.rows(arma::span(0, n - 2));  // Remove last row

  // Remove column i
  for (arma::uword j = i; j < n - 1; j++) {
    X.col(j) = X.col(j+1);  // Shift all columns right of i to the left by one
  }
  X = X.cols(arma::span(0,n-2));  // Remove last column

  return X;
}

// ----------------------------------------------------------------------------|
// Add a row and column to a matrix (and fill with beta variables)
// ----------------------------------------------------------------------------|
arma::mat add_row_col_block_prob_matrix(arma::mat X,
                                            double beta_alpha,
                                            double beta_beta) {
  arma::uword dim = X.n_rows;
  arma::mat Y(dim+1,dim+1,arma::fill::zeros);

  for(arma::uword r = 0; r < dim; r++) {
    for(arma::uword c = 0; c < dim; c++) {
      Y(r, c) = X(r, c);
    }
  }

  for(arma::uword i = 0; i < dim; i++) {
    Y(dim, i) = R::rbeta(beta_alpha, beta_beta);
    Y(i, dim) = Y(dim, i);
  }
  Y(dim, dim) = R::rbeta(beta_alpha, beta_beta);

  return Y;
}


// ----------------------------------------------------------------------------|
// Compute partition coefficient for the MFM - SBM
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
arma::vec compute_Vn_mfm_sbm(arma::uword no_variables,
                                 double dirichlet_alpha,
                                 arma::uword t_max,
                                 double lambda) {
  arma::vec log_Vn(t_max);
  double r;
  double tmp;

  for(arma::uword t = 1; t <= t_max; t++) {
    r = -INFINITY;
    for(arma::uword k = t; k < 500; k++) {
      tmp = 0.0;
      for(arma::uword tt = 1 - t; tt < 1; tt ++) {
        tmp += std::log(dirichlet_alpha * k + tt);
      }
      for(arma::uword n = 0; n < no_variables; n++) {
        tmp -= std::log(dirichlet_alpha * k + n);
      }
      tmp -= std::lgamma(k + 1);

      // Add the poisson term
      double log_norm_factor = log(1.0 - exp(R::dpois(0, lambda, true)));
      tmp += R::dpois(k-1, lambda, true) - log_norm_factor;

      // Compute the maximum between r and tmp
      if (tmp > r) {
        r = std::log(std::exp(r - tmp) + 1) + tmp;
      } else {
        r = std::log(1 + std::exp(tmp - r)) + r;
      }
    }
    log_Vn(t-1) = r - std::log(std::exp(1) - 1);
  }
  return log_Vn;
}

// ----------------------------------------------------------------------------|
// Compute log-likelihood for the MFM - SBM
// ----------------------------------------------------------------------------|
double log_likelihood_mfm_sbm(arma::uvec cluster_assign,
                              arma::mat cluster_probs,
                              arma::umat indicator,
                              arma::uword node,
                              arma::uword no_variables) {
  double output = 0;

  for(arma::uword j = 0; j < no_variables; j++) {
    if(j != node) {
      if(j < node) {
        output += indicator(j, node) *
          std::log(cluster_probs(cluster_assign(j), cluster_assign(node)));
        output += (1 - indicator(j, node)) *
          std::log(1 - cluster_probs(cluster_assign(j), cluster_assign(node)));
      } else {
        output += indicator(node, j) *
          std::log(cluster_probs(cluster_assign(node), cluster_assign(j)));
        output += (1 - indicator(node, j)) *
          std::log(1 - cluster_probs(cluster_assign(node), cluster_assign(j)));
      }
    }
  }

  return output;
}

// ----------------------------------------------------------------------------|
// Compute log-marginal for the MFM - SBM
// ----------------------------------------------------------------------------|
double log_marginal_mfm_sbm(arma::uvec cluster_assign,
                            arma::umat indicator,
                            arma::uword node,
                            arma::uword no_variables,
                            double beta_bernoulli_alpha,
                            double beta_bernoulli_beta) {

  arma::uword no_clusters_excl_node = arma::max(cluster_assign);
  //*std::max_element(cluster_assign.begin(), cluster_assign.end());
  double output = 0;
  double sumG;
  double sumN;

  output -= no_clusters_excl_node * R::lbeta(beta_bernoulli_alpha, beta_bernoulli_beta);

  for(arma::uword c = 0; c < no_clusters_excl_node; c++) {
    sumG = 0;
    sumN = 0;
    for(arma::uword i = 0; i < no_variables; i++) {
      if(cluster_assign(i) == c) {
        sumG += static_cast<double>(indicator(node, i));
        sumN += 1.0;
      }
    }
    output += R::lbeta(sumG + beta_bernoulli_alpha, sumN - sumG + beta_bernoulli_beta);
  }
  return output;
}

// ----------------------------------------------------------------------------|
// Helper function to update sumG in sample_block_probs_mfm_sbm()
// ----------------------------------------------------------------------------|
inline void update_sumG(double &sumG,
                        const arma::uvec &cluster_assign,
                        const arma::umat &indicator,
                        arma::uword r,
                        arma::uword s,
                        arma::uword no_variables) {
  for(arma::uword node1 = 0; node1 < no_variables - 1; node1++) {
    if(cluster_assign(node1) == r) {
      for(arma::uword node2 = node1 + 1; node2 < no_variables; node2++) {
        if(cluster_assign(node2) == s) {
          sumG += static_cast<double>(indicator(node1, node2));
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------|
// Sample the cluster assignment in sample_block_allocations_mfm_sbm()
// ----------------------------------------------------------------------------|
arma::uword sample_cluster(arma::vec cluster_prob) {
  arma::vec cum_prob = arma::cumsum(cluster_prob);
  double u = R::runif(0, arma::max(cum_prob));

  for (arma::uword i = 0; i < cum_prob.n_elem; i++) {
    if (u <= cum_prob(i)) {
      return i;
    }
  }
  return cum_prob.n_elem;
}

// ----------------------------------------------------------------------------|
// Sample the block allocations for the MFM - SBM
// ----------------------------------------------------------------------------|
arma::uvec block_allocations_mfm_sbm(arma::uvec cluster_assign,
                                        arma::uword no_variables,
                                        arma::vec log_Vn,
                                        arma::mat block_probs,
                                        arma::umat indicator,
                                        arma::uword dirichlet_alpha,
                                        double beta_bernoulli_alpha,
                                        double beta_bernoulli_beta) {
  arma::uword old;
  arma::uword cluster;
  arma::uword no_clusters;
  double prob;
  double loglike;
  double logmarg;

  // Generate a randomized order using Rcpp's sample function
  arma::ivec indices = Rcpp::sample(no_variables, no_variables);

  for (arma::uword idx = 0; idx < no_variables; ++idx) {
    arma::uword node = indices(idx) - 1; // Convert to zero-based index
    old = cluster_assign(node);

    arma::uvec cluster_size = table_cpp(cluster_assign);
    no_clusters = cluster_size.n_elem;

    if (cluster_size(old) == 1) {
      // Singleton cluster.

      // Cluster sizes without node
      arma::uvec cluster_size_node = cluster_size;
      cluster_size_node(old) -= (1 + dirichlet_alpha);

      // Compute probabilities for sampling process
      arma::vec cluster_prob(no_clusters + 1);
      for (arma::uword c = 0; c <= no_clusters; c++) {
        arma::uvec cluster_assign_tmp = cluster_assign;
        cluster_assign_tmp(node) = c;

        if (c < no_clusters) {
          loglike = log_likelihood_mfm_sbm(cluster_assign_tmp,
                                           block_probs,
                                           indicator,
                                           node,
                                           no_variables);

          prob = (static_cast<double>(dirichlet_alpha) + static_cast<double>(cluster_size_node(c))) *
            std::exp(loglike);
        } else {
          logmarg = log_marginal_mfm_sbm(cluster_assign_tmp,
                                         indicator,
                                         node,
                                         no_variables,
                                         beta_bernoulli_alpha,
                                         beta_bernoulli_beta);

          prob = static_cast<double>(dirichlet_alpha) *
            std::exp(logmarg) *
            std::exp(log_Vn(no_clusters - 1) - log_Vn(no_clusters - 2));
        }

        cluster_prob(c) = prob;
      }


      //Choose the cluster number for node
      cluster = sample_cluster(cluster_prob);

      // Remove the empty cluster
      if (cluster == no_clusters) {
        cluster_assign(node) = old;
      } else {
        cluster_assign(node) = cluster;
        for (arma::uword i = 0; i < no_variables; i++) {
          if (cluster_assign(i) > old) {
            cluster_assign(i) -= 1;
          }
        }
        block_probs = remove_row_col_matrix(block_probs, old);
      }
    } else {
      // Cluster sizes without node
      arma::uvec cluster_size_node = cluster_size;
      cluster_size_node(old) -= 1;

      // Compute probabilities for sampling process
      arma::vec cluster_prob(no_clusters + 1);
      for (arma::uword c = 0; c <= no_clusters; c++) {
        arma::uvec cluster_assign_tmp = cluster_assign;
        cluster_assign_tmp(node) = c;

        if (c < no_clusters) {
          loglike = log_likelihood_mfm_sbm(cluster_assign_tmp,
                                           block_probs,
                                           indicator,
                                           node,
                                           no_variables);

          prob = (static_cast<double>(dirichlet_alpha) + static_cast<double>(cluster_size_node(c))) *
            std::exp(loglike);
        } else {
          logmarg = log_marginal_mfm_sbm(cluster_assign_tmp,
                                         indicator,
                                         node,
                                         no_variables,
                                         beta_bernoulli_alpha,
                                         beta_bernoulli_beta);

          prob = static_cast<double>(dirichlet_alpha) *
            std::exp(logmarg) *
            std::exp(log_Vn(no_clusters) - log_Vn(no_clusters - 1));
        }

        cluster_prob(c) = prob;
      }


      //Choose the cluster number for node
      cluster = sample_cluster(cluster_prob);

      cluster_assign(node) = cluster;

      if (cluster == no_clusters) {
        block_probs = add_row_col_block_prob_matrix(block_probs,
                                                    beta_bernoulli_alpha,
                                                    beta_bernoulli_beta);
      }
    }
  }
  return cluster_assign;

}

// ----------------------------------------------------------------------------|
// Sample the block parameters for the MFM - SBM
// ----------------------------------------------------------------------------|
arma::mat block_probs_mfm_sbm(arma::uvec cluster_assign,
                                  arma::umat indicator,
                                  arma::uword no_variables,
                                  double beta_bernoulli_alpha,
                                  double beta_bernoulli_beta) {

  arma::uvec cluster_size = table_cpp(cluster_assign);
  arma::uword no_clusters = cluster_size.n_elem;

  arma::mat block_probs(no_clusters, no_clusters);
  arma::mat theta(no_variables, no_variables);

  double sumG;
  double size;

  for(arma::uword r = 0; r < no_clusters; r++) {
    for(arma::uword s = r; s < no_clusters; s++) {
      sumG = 0;
      if(r == s) {
        update_sumG(sumG, cluster_assign, indicator, r, r, no_variables);
        size = static_cast<double>(cluster_size(r)) * (static_cast<double>(cluster_size(r)) - 1) / 2;
      } else {
        update_sumG(sumG, cluster_assign, indicator, r, s, no_variables);
        update_sumG(sumG, cluster_assign, indicator, s, r, no_variables);
        size = static_cast<double>(cluster_size(s)) * static_cast<double>(cluster_size(r));
      }
      block_probs(r, s) = R::rbeta(sumG + beta_bernoulli_alpha, size - sumG + beta_bernoulli_beta);
      block_probs(s, r) = block_probs(r, s);
    }
  }

  return block_probs;
}
