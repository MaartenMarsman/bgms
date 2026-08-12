#include <RcppArmadillo.h>
#include <stdexcept>
#include "rng/rng_utils.h"
#include "math/explog_macros.h"
#include "priors/edge_prior_correction.h"

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
// Add a row and column to a matrix (and fill with beta variables)
// Modified to support separate within/between cluster hyperparameters
// ----------------------------------------------------------------------------|
arma::mat add_row_col_block_prob_matrix(arma::mat X,
                                        double beta_alpha,
                                        double beta_beta,
                                        SafeRNG& rng,
                                        double beta_bernoulli_alpha_between,
                                        double beta_bernoulli_beta_between) {
  arma::uword dim = X.n_rows;
  arma::mat Y(dim+1,dim+1,arma::fill::zeros);

  for(arma::uword r = 0; r < dim; r++) {
    for(arma::uword c = 0; c < dim; c++) {
      Y(r, c) = X(r, c);
    }
  }

  // Add new row and column for the new cluster
  for(arma::uword i = 0; i < dim; i++) {
    // Between-cluster edge probabilities (new cluster to existing clusters)
    Y(dim, i) = rbeta(rng, beta_bernoulli_alpha_between, beta_bernoulli_beta_between);
    Y(i, dim) = Y(dim, i);
  }

  // Within-cluster edge probability (diagonal element for new cluster)
  Y(dim, dim) = rbeta(rng, beta_alpha, beta_beta);

  return Y;
}



// Function: log_likelihood_mfm_sbm
//
// Computes the log-likelihood contribution for a single node under the
// Mixture of Finite Mixtures Stochastic Block Model (MFM-SBM). Evaluates
// the probability of observed edges between the node and all other nodes
// given their cluster assignments and cluster connection probabilities.
//
// Inputs:
//  - cluster_assign: Vector of cluster assignments for all nodes.
//  - cluster_probs: Matrix of edge probabilities between clusters.
//  - indicator: Upper-triangular matrix of edge indicators (1 = edge present).
//  - node: Index of the node whose contribution is computed.
//  - no_variables: Total number of nodes in the network.
//
// Returns:
//  - Log-likelihood contribution for the specified node.
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
          MY_LOG(cluster_probs(cluster_assign(j), cluster_assign(node)));
        output += (1 - indicator(j, node)) *
          MY_LOG(1 - cluster_probs(cluster_assign(j), cluster_assign(node)));
      } else {
        output += indicator(node, j) *
          MY_LOG(cluster_probs(cluster_assign(node), cluster_assign(j)));
        output += (1 - indicator(node, j)) *
          MY_LOG(1 - cluster_probs(cluster_assign(node), cluster_assign(j)));
      }
    }
  }

  return output;
}

// Function: log_marginal_mfm_sbm
//
// Computes the log-marginal likelihood contribution for a single node under
// the MFM-SBM after integrating out cluster connection probabilities. Used
// when proposing a new cluster for the node, so all interactions are
// between-cluster. Uses Beta-Bernoulli conjugacy.
//
// Inputs:
//  - cluster_assign: Vector of cluster assignments for all nodes.
//  - indicator: Upper-triangular matrix of edge indicators (1 = edge present).
//  - node: Index of the node whose contribution is computed.
//  - no_variables: Total number of nodes in the network.
//  - beta_bernoulli_alpha_between: Alpha hyperparameter for between-cluster edges.
//  - beta_bernoulli_beta_between: Beta hyperparameter for between-cluster edges.
//
// Returns:
//  - Log-marginal likelihood contribution for the specified node.
double log_marginal_mfm_sbm(arma::uvec cluster_assign,
                            arma::umat indicator,
                            arma::uword node,
                            arma::uword no_variables,
                            double beta_bernoulli_alpha_between,
                            double beta_bernoulli_beta_between) {

  arma::uvec indices = arma::regspace<arma::uvec>(0, no_variables-1); // vector of variables indices [0, 1, ..., no_variables-1]
  arma::uvec select_variables = indices(arma::find(indices != node)); // vector of variables indices excluding 'node'
  arma::uvec cluster_assign_wo_node = cluster_assign(select_variables); // vector of cluster labels for all variables but excluding 'node'
  arma::uvec indicator_node = indicator.col(node); // column of indicator matrix corresponding to 'node'
  arma::vec gamma_node = arma::conv_to<arma::vec>::from(indicator_node(select_variables)); // selecting only indicators between 'node' and the remaining variables (thus excluding indicator of node with itself -- that is indicator[node,node])
  arma::uvec table_cluster = table_cpp(cluster_assign_wo_node); // frequency table of clusters excluding node

  double output = 0;
  for(arma::uword i = 0; i < table_cluster.n_elem; i++){
    if(table_cluster(i) > 0){ // if the cluster is empty -- table_cluster(i) = 0 == then it is the previous cluster of 'node' where 'node' was the only member - a singleton, thus skip)
      arma::uvec which_variables_cluster_i = arma::find(cluster_assign_wo_node == i); // which variables belong to cluster i
      int sumG = arma::accu(gamma_node(which_variables_cluster_i)); // sum the indicator variables between node and those variables
      int sumN = static_cast<int>(which_variables_cluster_i.n_elem); // take the size of the group as maximum number of relations

      // Compute log B(alpha + G, beta + N - G) / B(alpha, beta) using the identity
      // Gamma(z+n)/Gamma(z) = prod_{m=0}^{n-1} (z+m), which holds since sumG and sumN are integers.
      // This avoids lbeta/lgamma calls and is faster for small cluster sizes.
      for(int m = 0; m < sumG; m++){
        output += MY_LOG(beta_bernoulli_alpha_between + m);
      }
      for(int m = 0; m < (sumN - sumG); m++){
        output += MY_LOG(beta_bernoulli_beta_between + m);
      }
      for(int m = 0; m < sumN; m++){
        output -= MY_LOG(beta_bernoulli_alpha_between + beta_bernoulli_beta_between + m);
      }
    }
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

// Defined below; the uncorrected sweep samples from max-shifted log-weights.
static arma::uword sample_cluster_log(const arma::vec& log_weights,
                                      SafeRNG& rng);

// ----------------------------------------------------------------------------|
// Sample the block allocations for the MFM - SBM
// Modified to support separate within/between cluster hyperparameters
// ----------------------------------------------------------------------------|
arma::uvec block_allocations_mfm_sbm(arma::uvec cluster_assign,
                                     arma::uword no_variables,
                                     arma::vec log_Vn,
                                     arma::mat block_probs,
                                     arma::umat indicator,
                                     double dirichlet_alpha,
                                     double beta_bernoulli_alpha,
                                     double beta_bernoulli_beta,
                                     double beta_bernoulli_alpha_between,
                                     double beta_bernoulli_beta_between,
                                     SafeRNG& rng) {

  arma::uword old;
  arma::uword cluster;
  arma::uword no_clusters;
  double loglike;
  double logmarg;

  // Generate a randomized order using Rcpp's sample function
  arma::uvec indices = arma_randperm(rng, no_variables);

  for (arma::uword idx = 0; idx < no_variables; idx++) {
    arma::uword node = indices(idx);
    old = cluster_assign(node);

    arma::uvec cluster_size = table_cpp(cluster_assign);
    no_clusters = cluster_size.n_elem;

    if (cluster_size(old) == 1) {
      // Singleton cluster.

      // Cluster sizes without node
      arma::uvec cluster_size_node = cluster_size;

      // Compute log-weights for the sampling process (max-shifted at draw).
      arma::vec log_weights(no_clusters + 1);
      for (arma::uword c = 0; c <= no_clusters; c++) {
        arma::uvec cluster_assign_tmp = cluster_assign;
        cluster_assign_tmp(node) = c;

        if (c < no_clusters) {
          if(c != old){
            loglike = log_likelihood_mfm_sbm(cluster_assign_tmp,
                                             block_probs,
                                             indicator,
                                             node,
                                             no_variables);

            log_weights(c) =
              MY_LOG(dirichlet_alpha + static_cast<double>(cluster_size_node(c))) +
              loglike;
          }
          else{ // if old group, the weight is zero (log-weight -inf)
            log_weights(c) = -arma::datum::inf;
          }

        } else {
          logmarg = log_marginal_mfm_sbm(cluster_assign_tmp,
                                         indicator,
                                         node,
                                         no_variables,
                                         beta_bernoulli_alpha_between,
                                         beta_bernoulli_beta_between);

          log_weights(c) = MY_LOG(dirichlet_alpha) + logmarg +
            (log_Vn(no_clusters - 1) - log_Vn(no_clusters - 2));
        }
      }

      //Choose the cluster number for node
      cluster = sample_cluster_log(log_weights, rng);

      //if the sampled cluster is the new added cluster or the old one
      if (cluster == no_clusters) {
        cluster_assign(node) = old; // new cluster takes the place of the older singleton but doesn't update probabilities, they are kept the same as the old ones
      } else { // otherwise remove old (singleton) empty cluster, redefine cluster_assign and block_probs
        cluster_assign(node) = cluster;
        for (arma::uword i = 0; i < no_variables; i++) {
          if (cluster_assign(i) > old) {
            cluster_assign(i) -= 1;
          }
        }
        // removing row and col index 'old' from block_probs
        block_probs.shed_row(old);
        block_probs.shed_col(old);
      }
    } else {
      // Cluster sizes without node
      arma::uvec cluster_size_node = cluster_size;
      cluster_size_node(old) -= 1;

      // Compute log-weights for the sampling process (max-shifted at draw).
      arma::vec log_weights(no_clusters + 1);
      for (arma::uword c = 0; c <= no_clusters; c++) {
        arma::uvec cluster_assign_tmp = cluster_assign;
        cluster_assign_tmp(node) = c;
        if (c < no_clusters) {
          loglike = log_likelihood_mfm_sbm(cluster_assign_tmp,
                                           block_probs,
                                           indicator,
                                           node,
                                           no_variables);

          log_weights(c) =
            MY_LOG(dirichlet_alpha + static_cast<double>(cluster_size_node(c))) +
            loglike;
        } else {
          logmarg = log_marginal_mfm_sbm(cluster_assign_tmp,
                                         indicator,
                                         node,
                                         no_variables,
                                         beta_bernoulli_alpha_between,
                                         beta_bernoulli_beta_between);

          log_weights(c) = MY_LOG(dirichlet_alpha) + logmarg +
            (log_Vn(no_clusters) - log_Vn(no_clusters-1));
        }
      }


      //Choose the cluster number for node
      cluster = sample_cluster_log(log_weights, rng);

      cluster_assign(node) = cluster;

      if (cluster == no_clusters) {
        block_probs = add_row_col_block_prob_matrix(block_probs,
                                                    beta_bernoulli_alpha,
                                                    beta_bernoulli_beta,
                                                    rng,
                                                    beta_bernoulli_alpha_between,
                                                    beta_bernoulli_beta_between);
      }
    }
  }
  return cluster_assign;

}

// ----------------------------------------------------------------------------|
// Sample the block parameters for the MFM - SBM
// Modified to support separate within/between cluster hyperparameters
// ----------------------------------------------------------------------------|
arma::mat block_probs_mfm_sbm(arma::uvec cluster_assign,
                              arma::umat indicator,
                              arma::uword no_variables,
                              double beta_bernoulli_alpha,
                              double beta_bernoulli_beta,
                              double beta_bernoulli_alpha_between,
                              double beta_bernoulli_beta_between,
                              SafeRNG& rng) {

  arma::uvec cluster_size = table_cpp(cluster_assign);
  arma::uword no_clusters = cluster_size.n_elem;

  arma::mat block_probs(no_clusters, no_clusters);

  double sumG;
  double size;

  for(arma::uword r = 0; r < no_clusters; r++) {
    for(arma::uword s = r; s < no_clusters; s++) {
      sumG = 0;

      if(r == s) {
        // Within-cluster: always use main parameters
        update_sumG(sumG, cluster_assign, indicator, r, r, no_variables);
        size = static_cast<double>(cluster_size(r)) * (static_cast<double>(cluster_size(r)) - 1) / 2;
        block_probs(r, s) = rbeta(rng,
                    sumG + beta_bernoulli_alpha,
                    size - sumG + beta_bernoulli_beta);
      } else {
        // Between-cluster: use between parameters
        update_sumG(sumG, cluster_assign, indicator, r, s, no_variables);
        update_sumG(sumG, cluster_assign, indicator, s, r, no_variables);
        size = static_cast<double>(cluster_size(s)) * static_cast<double>(cluster_size(r));

        block_probs(r, s) = rbeta(rng, sumG + beta_bernoulli_alpha_between, size - sumG + beta_bernoulli_beta_between);
      }
      block_probs(s, r) = block_probs(r, s);
    }
  }

  return block_probs;
}

// ----------------------------------------------------------------------------|
// Normalizing-constant corrections for the MFM-SBM under the determinant-
// tilted GGM prior. The per-graph normalizer Z(Gamma) tilts the graph law;
// the block-model hyperparameter updates read the correction locally off the
// slope curve f'(local density) at each edge's min-endpoint expected degree
// density, and the new-cluster collapsed marginal integrates the per-pair
// curve f(theta).
// ----------------------------------------------------------------------------|

// ----------------------------------------------------------------------------|
// Self-consistent per-edge slopes c_e given labels and block probabilities:
// fixed point of e = th * exp(c) / (1 - th + th * exp(c)) with
// c = f'(min-endpoint expected degree density). Only continuous-continuous
// pairs carry the tilt: other pairs keep slope zero and stay out of the
// expected-degree sums, which run over the continuous subgraph.
// ----------------------------------------------------------------------------|
arma::mat compute_ce_sbm(const arma::uvec& cluster_assign,
                         const arma::mat& block_probs,
                         arma::uword no_variables,
                         const SBMCorrection& correction) {
  arma::uword q = no_variables;
  arma::mat ce(q, q, arma::fill::zeros);
  arma::uword q_cc = correction.num_continuous(q);
  if(q_cc < 2) return ce;
  arma::mat TH(q, q, arma::fill::zeros), ep(q, q, arma::fill::zeros);
  for(arma::uword i = 0; i < q - 1; i++) {
    for(arma::uword j = i + 1; j < q; j++) {
      double th = block_probs(cluster_assign(i), cluster_assign(j));
      TH(i, j) = th; TH(j, i) = th;
      ep(i, j) = th; ep(j, i) = th;
    }
  }
  for(int it = 0; it < 40; it++) {
    arma::vec deg(q, arma::fill::zeros);
    for(arma::uword i = 0; i < q; i++) {
      if(!correction.node_continuous(i)) continue;
      double s = 0;
      for(arma::uword j = 0; j < q; j++) {
        if(j == i || !correction.node_continuous(j)) continue;
        s += ep(i, j);
      }
      deg(i) = s / static_cast<double>(q_cc - 1);
    }
    double max_diff = 0;
    for(arma::uword i = 0; i < q - 1; i++) {
      if(!correction.node_continuous(i)) continue;
      for(arma::uword j = i + 1; j < q; j++) {
        if(!correction.node_continuous(j)) continue;
        double dmin = std::min(deg(i), deg(j));
        double c = correction.fprime_at(dmin);
        ce(i, j) = c; ce(j, i) = c;
        double th = TH(i, j);
        double ec = MY_EXP(c);
        double e = th * ec / (1.0 - th + th * ec);
        max_diff = std::max(max_diff, std::abs(e - ep(i, j)));
        ep(i, j) = e; ep(j, i) = e;
      }
    }
    if(max_diff < 1e-9) break;
  }
  return ce;
}

// ----------------------------------------------------------------------------|
// Baseline expected degree densities under (labels, block probabilities) via
// the self-consistent local-density prediction, over the continuous
// subgraph. Discrete nodes keep density zero; their entries are never read
// by the tilt terms.
// ----------------------------------------------------------------------------|
arma::vec degrees_ld_sbm(const arma::uvec& cluster_assign,
                         const arma::mat& block_probs,
                         arma::uword no_variables,
                         const SBMCorrection& correction) {
  arma::uword q = no_variables;
  arma::vec deg(q, arma::fill::zeros);
  arma::uword q_cc = correction.num_continuous(q);
  if(q_cc < 2) return deg;
  arma::mat ce = compute_ce_sbm(cluster_assign, block_probs, q, correction);
  for(arma::uword i = 0; i < q; i++) {
    if(!correction.node_continuous(i)) continue;
    double s = 0;
    for(arma::uword j = 0; j < q; j++) {
      if(j == i || !correction.node_continuous(j)) continue;
      double th = block_probs(cluster_assign(i), cluster_assign(j));
      double ec = MY_EXP(ce(i, j));
      s += th * ec / (1.0 - th + th * ec);
    }
    deg(i) = s / static_cast<double>(q_cc - 1);
  }
  return deg;
}

// ----------------------------------------------------------------------------|
// Mini thermodynamic integration: logC(z with node in cand) - logC(z with
// node in cur), integrating the predicted per-edge score along the morph of
// the node's incident block-pair probabilities. The node's expected degree
// is updated self-consistently along the morph; the other endpoints' stay
// at deg_base.
// ----------------------------------------------------------------------------|
double miniti_node_sbm(arma::uword node,
                       const arma::uvec& cluster_assign,
                       const arma::mat& block_probs,
                       arma::uword no_variables,
                       arma::uword cur,
                       arma::uword cand,
                       const arma::vec& deg_base,
                       const SBMCorrection& correction) {
  const int T = 4;
  arma::uword q = no_variables;
  arma::uword q_cc = correction.num_continuous(q);
  if(q_cc < 2 || !correction.node_continuous(node)) return 0.0;
  std::vector<double> th_old(q), th_new(q);
  std::vector<arma::uword> js;
  js.reserve(q - 1);
  for(arma::uword j = 0; j < q; j++) {
    if(j == node || !correction.node_continuous(j)) continue;
    th_old[j] = block_probs(cur, cluster_assign(j));
    th_new[j] = block_probs(cand, cluster_assign(j));
    js.push_back(j);
  }
  double dlogC = 0.0, prev = 0.0;
  for(int ti = 0; ti <= T; ti++) {
    double t = static_cast<double>(ti) / T;
    double degk = 0.0;
    for(arma::uword j : js) degk += (1 - t) * th_old[j] + t * th_new[j];
    degk /= static_cast<double>(q_cc - 1);
    for(int it = 0; it < 3; it++) {
      double s = 0;
      for(arma::uword j : js) {
        double th = (1 - t) * th_old[j] + t * th_new[j];
        double dmin = std::min(degk, deg_base(j));
        double ec = MY_EXP(correction.fprime_at(dmin));
        s += th * ec / (1.0 - th + th * ec);
      }
      double ndeg = s / static_cast<double>(q_cc - 1);
      if(std::abs(ndeg - degk) < 1e-9) { degk = ndeg; break; }
      degk = ndeg;
    }
    double sc = 0;
    for(arma::uword j : js) {
      double th = (1 - t) * th_old[j] + t * th_new[j];
      double dth = th_new[j] - th_old[j];
      double dmin = std::min(degk, deg_base(j));
      double ec = MY_EXP(correction.fprime_at(dmin));
      sc += dth * (ec - 1.0) / (1.0 - th + th * ec);
    }
    if(ti > 0) dlogC += 0.5 * (prev + sc) * (1.0 / T);
    prev = sc;
  }
  return dlogC;
}

// ----------------------------------------------------------------------------|
// Mini thermodynamic integration for the node's edge REMOVAL: logC(z with the
// node's incident pair probabilities morphed to zero) - logC(z). Anchors the
// new-cluster weight to the same relative-to-current baseline as the
// existing-cluster mini-TI: the corrected collapsed marginal carries the new
// block's absolute tilt factors, so without this anchor the grow option is
// over-credited by the node's current-state tilt contribution.
// ----------------------------------------------------------------------------|
double miniti_removal_sbm(arma::uword node,
                          const arma::uvec& cluster_assign,
                          const arma::mat& block_probs,
                          arma::uword no_variables,
                          arma::uword cur,
                          const arma::vec& deg_base,
                          const SBMCorrection& correction) {
  const int T = 4;
  arma::uword q = no_variables;
  arma::uword q_cc = correction.num_continuous(q);
  if(q_cc < 2 || !correction.node_continuous(node)) return 0.0;
  std::vector<double> th_old(q);
  std::vector<arma::uword> js;
  js.reserve(q - 1);
  for(arma::uword j = 0; j < q; j++) {
    if(j == node || !correction.node_continuous(j)) continue;
    th_old[j] = block_probs(cur, cluster_assign(j));
    js.push_back(j);
  }
  double dlogC = 0.0, prev = 0.0;
  for(int ti = 0; ti <= T; ti++) {
    double t = static_cast<double>(ti) / T;
    double degk = 0.0;
    for(arma::uword j : js) degk += (1 - t) * th_old[j];
    degk /= static_cast<double>(q_cc - 1);
    for(int it = 0; it < 3; it++) {
      double s = 0;
      for(arma::uword j : js) {
        double th = (1 - t) * th_old[j];
        double dmin = std::min(degk, deg_base(j));
        double ec = MY_EXP(correction.fprime_at(dmin));
        s += th * ec / (1.0 - th + th * ec);
      }
      double ndeg = s / static_cast<double>(q_cc - 1);
      if(std::abs(ndeg - degk) < 1e-9) { degk = ndeg; break; }
      degk = ndeg;
    }
    double sc = 0;
    for(arma::uword j : js) {
      double th = (1 - t) * th_old[j];
      double dth = 0.0 - th_old[j];
      double dmin = std::min(degk, deg_base(j));
      double ec = MY_EXP(correction.fprime_at(dmin));
      sc += dth * (ec - 1.0) / (1.0 - th + th * ec);
    }
    if(ti > 0) dlogC += 0.5 * (prev + sc) * (1.0 / T);
    prev = sc;
  }
  return dlogC;
}

// ----------------------------------------------------------------------------|
// Corrected collapsed log-marginal for assigning a node to a NEW cluster: for
// each existing cluster r the fresh block-pair shares one probability drawn
// from the Beta prior, so the tilt enters as the homogeneous whole-pair
// factor exp(-n_cc * f(theta)) inside the Beta-Bernoulli integral, computed
// by quadrature on the per-pair f curve. Only the node's continuous-
// continuous pairs are tilted: n_cc counts the cluster's continuous members,
// and zero for a discrete node. Cluster terms without tilted pairs use the
// exact conjugate Beta-Bernoulli marginal.
// ----------------------------------------------------------------------------|
double corrected_log_marginal_mfm_sbm(const arma::uvec& cluster_assign,
                                      const arma::umat& indicator,
                                      arma::uword node,
                                      arma::uword no_variables,
                                      double beta_bernoulli_alpha_between,
                                      double beta_bernoulli_beta_between,
                                      const SBMCorrection& correction) {
  const arma::vec& thq = correction.quad_theta();
  const arma::vec& fq = correction.quad_f();
  arma::uword Nq = thq.n_elem;
  double dth = thq(1) - thq(0);
  double lbab = R::lbeta(beta_bernoulli_alpha_between, beta_bernoulli_beta_between);
  arma::uword no_clusters = arma::max(cluster_assign) + 1;
  bool node_cc = correction.node_continuous(node);

  double out = 0;
  std::vector<double> lp(Nq);
  // log(th)/log1p(-th) are cluster-invariant; precompute the grids once.
  std::vector<double> log_thq(Nq), log1p_neg_thq(Nq);
  for(arma::uword k = 0; k < Nq; k++) {
    log_thq[k] = MY_LOG(thq(k));
    log1p_neg_thq[k] = MY_LOG1P(-thq(k));
  }
  for(arma::uword r = 0; r < no_clusters; r++) {
    int nr = 0, mr = 0, nr_cc = 0;
    for(arma::uword j = 0; j < no_variables; j++) {
      if(j == node || cluster_assign(j) != r) continue;
      nr++;
      if(correction.node_continuous(j)) nr_cc++;
      arma::uword a = std::min(node, j), b = std::max(node, j);
      if(indicator(a, b) == 1) mr++;
    }
    if(nr == 0) continue;
    int tilt_n = node_cc ? nr_cc : 0;
    if(tilt_n == 0) {
      out += R::lbeta(beta_bernoulli_alpha_between + mr,
                      beta_bernoulli_beta_between + nr - mr) - lbab;
      continue;
    }
    double mx = -std::numeric_limits<double>::infinity();
    for(arma::uword k = 0; k < Nq; k++) {
      lp[k] = (beta_bernoulli_alpha_between - 1.0 + mr) * log_thq[k] +
        (beta_bernoulli_beta_between - 1.0 + nr - mr) * log1p_neg_thq[k] -
        static_cast<double>(tilt_n) * fq(k);
      if(lp[k] > mx) mx = lp[k];
    }
    double s = 0;
    for(arma::uword k = 0; k < Nq; k++) s += MY_EXP(lp[k] - mx);
    out += (MY_LOG(s) + mx + MY_LOG(dth)) - lbab;
  }
  return out;
}

// ----------------------------------------------------------------------------|
// Sample a cluster from log-weights (max-shifted).
// ----------------------------------------------------------------------------|
static arma::uword sample_cluster_log(const arma::vec& log_weights,
                                      SafeRNG& rng) {
  double mx = log_weights.max();
  arma::vec w = ARMA_MY_EXP(log_weights - mx);
  double u = runif(rng) * arma::accu(w);
  double cum = 0;
  for(arma::uword c = 0; c < w.n_elem; c++) {
    cum += w(c);
    if(u <= cum) return c;
  }
  return w.n_elem - 1;
}

// ----------------------------------------------------------------------------|
// Corrected block-allocation sweep: the collapsed Gibbs weights gain
// -(logC(z with node in c) - logC(z with node in old)) via the mini
// thermodynamic integration, and the new-cluster weight uses the corrected
// collapsed marginal. Expected degrees are estimated once per sweep. Grows
// and shrinks block_probs in place so the subsequent corrected block-
// probability draw centers on the carried values.
// ----------------------------------------------------------------------------|
arma::uvec block_allocations_mfm_sbm_corrected(arma::uvec cluster_assign,
                                               arma::uword no_variables,
                                               const arma::vec& log_Vn,
                                               arma::mat& block_probs,
                                               const arma::umat& indicator,
                                               double dirichlet_alpha,
                                               double beta_bernoulli_alpha,
                                               double beta_bernoulli_beta,
                                               double beta_bernoulli_alpha_between,
                                               double beta_bernoulli_beta_between,
                                               const SBMCorrection& correction,
                                               SafeRNG& rng) {
  if(correction.is_continuous().n_elem > 0 &&
     correction.is_continuous().n_elem != no_variables) {
    // std::runtime_error rather than Rcpp::stop: this runs on worker threads,
    // where constructing an Rcpp exception is not safe.
    throw std::runtime_error("SBM correction: is_continuous mask length does "
                             "not match the number of variables.");
  }
  arma::uvec indices = arma_randperm(rng, no_variables);
  double dir_alpha = static_cast<double>(dirichlet_alpha);

  // One expected-degree baseline per sweep, held fixed across candidates.
  arma::vec deg_base = degrees_ld_sbm(
    cluster_assign, block_probs, no_variables, correction);

  for(arma::uword idx = 0; idx < no_variables; idx++) {
    arma::uword node = indices(idx);
    arma::uword old = cluster_assign(node);

    arma::uvec cluster_size = table_cpp(cluster_assign);
    arma::uword no_clusters = cluster_size.n_elem;
    bool singleton = (cluster_size(old) == 1);

    arma::vec log_weights(no_clusters + 1);
    for(arma::uword c = 0; c < no_clusters; c++) {
      if(singleton && c == old) {
        log_weights(c) = -std::numeric_limits<double>::infinity();
        continue;
      }
      arma::uvec cluster_assign_tmp = cluster_assign;
      cluster_assign_tmp(node) = c;
      double loglike = log_likelihood_mfm_sbm(
        cluster_assign_tmp, block_probs, indicator, node, no_variables);
      double corr = (c == old) ? 0.0 : miniti_node_sbm(
        node, cluster_assign, block_probs, no_variables, old, c,
        deg_base, correction);
      double size_excl = static_cast<double>(cluster_size(c)) -
        ((c == old) ? 1.0 : 0.0);
      log_weights(c) = MY_LOG(dir_alpha + size_excl) + loglike - corr;
    }

    double logmarg = corrected_log_marginal_mfm_sbm(
      cluster_assign, indicator, node, no_variables,
      beta_bernoulli_alpha_between, beta_bernoulli_beta_between, correction);
    double removal = miniti_removal_sbm(
      node, cluster_assign, block_probs, no_variables, old,
      deg_base, correction);
    double vn_ratio = singleton
      ? (log_Vn(no_clusters - 1) - log_Vn(no_clusters - 2))
      : (log_Vn(no_clusters) - log_Vn(no_clusters - 1));
    log_weights(no_clusters) =
      MY_LOG(dir_alpha) + logmarg + vn_ratio - removal;

    arma::uword cluster = sample_cluster_log(log_weights, rng);

    if(singleton) {
      if(cluster == no_clusters) {
        // Keep the node in its own (old) cluster; probabilities unchanged.
        cluster_assign(node) = old;
      } else {
        cluster_assign(node) = cluster;
        for(arma::uword i = 0; i < no_variables; i++) {
          if(cluster_assign(i) > old) cluster_assign(i) -= 1;
        }
        block_probs.shed_row(old);
        block_probs.shed_col(old);
      }
    } else {
      cluster_assign(node) = cluster;
      if(cluster == no_clusters) {
        block_probs = add_row_col_block_prob_matrix(
          block_probs, beta_bernoulli_alpha, beta_bernoulli_beta, rng,
          beta_bernoulli_alpha_between, beta_bernoulli_beta_between);
      }
    }
  }
  return cluster_assign;
}

// ----------------------------------------------------------------------------|
// Corrected block-probability draw for pair (r, s): the conjugate density is
// multiplied by exp(-sum over the pair's edges of log(1 - theta +
// theta * exp(c_e))) with per-edge slopes from the self-consistent expected
// densities, and sampled by inverse CDF on a fine grid centered at the
// current value.
// ----------------------------------------------------------------------------|
static double draw_theta_local_density(int m, int n_pairs, double a, double b,
                                       double theta_current,
                                       const std::vector<double>& ce_edges,
                                       SafeRNG& rng) {
  const int N = 201;
  double mu = (a + m) / (a + b + n_pairs);
  double sd = std::sqrt(mu * (1.0 - mu) / (a + b + n_pairs + 1.0));
  double half_width = std::min(0.49, std::max(1e-4, 10.0 * sd));
  double lo = std::max(1e-7, theta_current - half_width);
  double hi = std::min(1.0 - 1e-7, theta_current + half_width);
  double grid[N], weight[N];
  double step = (hi - lo) / (N - 1);
  // exp(c_e) is grid-invariant; hoist it out of the N-point loop.
  std::vector<double> exp_ce(ce_edges.size());
  for(size_t e = 0; e < ce_edges.size(); e++) {
    exp_ce[e] = MY_EXP(ce_edges[e]);
  }
  double mx = -std::numeric_limits<double>::infinity();
  for(int k = 0; k < N; k++) {
    double th = lo + step * k;
    grid[k] = th;
    double lc = 0;
    for(size_t e = 0; e < exp_ce.size(); e++) {
      lc += MY_LOG(1.0 - th + th * exp_ce[e]);
    }
    weight[k] = (a + m - 1.0) * MY_LOG(th) +
      (b + n_pairs - m - 1.0) * MY_LOG1P(-th) - lc;
    if(weight[k] > mx) mx = weight[k];
  }
  double s = 0;
  for(int k = 0; k < N; k++) { weight[k] = MY_EXP(weight[k] - mx); s += weight[k]; }
  double u = runif(rng) * s;
  double cum = 0;
  int k = 0;
  for(; k < N - 1; k++) { cum += weight[k]; if(cum >= u) break; }
  return grid[k];
}

arma::mat block_probs_mfm_sbm_corrected(const arma::uvec& cluster_assign,
                                        const arma::mat& block_probs_current,
                                        const arma::umat& indicator,
                                        arma::uword no_variables,
                                        double beta_bernoulli_alpha,
                                        double beta_bernoulli_beta,
                                        double beta_bernoulli_alpha_between,
                                        double beta_bernoulli_beta_between,
                                        const SBMCorrection& correction,
                                        SafeRNG& rng) {
  if(correction.is_continuous().n_elem > 0 &&
     correction.is_continuous().n_elem != no_variables) {
    // std::runtime_error rather than Rcpp::stop: this runs on worker threads,
    // where constructing an Rcpp exception is not safe.
    throw std::runtime_error("SBM correction: is_continuous mask length does "
                             "not match the number of variables.");
  }
  arma::uvec cluster_size = table_cpp(cluster_assign);
  arma::uword no_clusters = cluster_size.n_elem;
  arma::mat block_probs(no_clusters, no_clusters);

  arma::mat ce = compute_ce_sbm(
    cluster_assign, block_probs_current, no_variables, correction);

  for(arma::uword r = 0; r < no_clusters; r++) {
    for(arma::uword s = r; s < no_clusters; s++) {
      bool within = (r == s);
      double sumG = 0;
      if(within) {
        update_sumG(sumG, cluster_assign, indicator, r, r, no_variables);
      } else {
        update_sumG(sumG, cluster_assign, indicator, r, s, no_variables);
        update_sumG(sumG, cluster_assign, indicator, s, r, no_variables);
      }
      int n_pairs = within
        ? static_cast<int>(cluster_size(r) * (cluster_size(r) - 1) / 2)
        : static_cast<int>(cluster_size(r) * cluster_size(s));
      std::vector<double> ce_edges;
      ce_edges.reserve(n_pairs);
      for(arma::uword i = 0; i < no_variables - 1; i++) {
        if(!correction.node_continuous(i)) continue;
        for(arma::uword j = i + 1; j < no_variables; j++) {
          if(!correction.node_continuous(j)) continue;
          arma::uword zi = cluster_assign(i), zj = cluster_assign(j);
          bool in_pair = within ? (zi == r && zj == r)
            : ((zi == r && zj == s) || (zi == s && zj == r));
          if(in_pair) ce_edges.push_back(ce(i, j));
        }
      }
      double a = within ? beta_bernoulli_alpha : beta_bernoulli_alpha_between;
      double b = within ? beta_bernoulli_beta : beta_bernoulli_beta_between;
      if(ce_edges.empty()) {
        // No tilted pairs in this block pair: the conjugate draw is exact.
        block_probs(r, s) = rbeta(
          rng, a + sumG, b + static_cast<double>(n_pairs) - sumG);
      } else {
        block_probs(r, s) = draw_theta_local_density(
          static_cast<int>(sumG), n_pairs, a, b,
          block_probs_current(r, s), ce_edges, rng);
      }
      block_probs(s, r) = block_probs(r, s);
    }
  }
  return block_probs;
}
