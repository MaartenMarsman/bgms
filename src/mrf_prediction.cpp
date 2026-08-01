// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils/variable_helpers.h"
using namespace Rcpp;


// ============================================================================
//   GGM Conditional Prediction
// ============================================================================

// Compute conditional Gaussian parameters for a GGM.
//
// For a GGM with precision matrix Omega, the conditional distribution of
// X_j given X_{-j} = x_{-j} is:
//
//   X_j | X_{-j} ~ N( -omega_{jj}^{-1} sum_{k != j} omega_{jk} x_k,
//                       omega_{jj}^{-1} )
//
// Takes an n x p matrix of observed continuous data, 0-based indices of the
// variables to predict, and the p x p precision matrix Omega. Returns a list
// of n x 2 matrices (one per predicted variable), where column 0 is the
// conditional mean and column 1 the conditional SD.
// [[Rcpp::export]]
Rcpp::List compute_conditional_ggm(
    const arma::mat& observations,
    const arma::ivec& predict_vars,
    const arma::mat& precision
) {
  int n = observations.n_rows;
  int num_predict_vars = predict_vars.n_elem;

  Rcpp::List result(num_predict_vars);

  for (int pv = 0; pv < num_predict_vars; pv++) {
    int j = predict_vars[pv];

    double omega_jj = precision(j, j);
    double cond_var = 1.0 / omega_jj;
    double cond_sd  = std::sqrt(cond_var);

    // Compute conditional means for all observations at once:
    // cond_mean_i = -(1/omega_jj) * sum_{k != j} omega_{jk} * x_{ik}
    //             = -(1/omega_jj) * (observations * omega_{.,j} - x_{ij} * omega_{jj})
    arma::vec omega_col_j = precision.col(j);
    arma::vec linear_pred = observations * omega_col_j;   // n x 1
    // Subtract the self-contribution x_{ij} * omega_{jj}
    linear_pred -= observations.col(j) * omega_jj;
    arma::vec cond_means = -cond_var * linear_pred;

    // Build n x 2 output matrix: [mean, sd]
    Rcpp::NumericMatrix out(n, 2);
    for (int i = 0; i < n; i++) {
      out(i, 0) = cond_means[i];
      out(i, 1) = cond_sd;
    }

    result[pv] = out;
  }

  return result;
}


// ============================================================================
//   OMRF Conditional Prediction
// ============================================================================

// Compute conditional probabilities P(X_j = c | X_{-j}) for ordinal MRF
// Uses numerically stable vectorized computation from variable_helpers.h

// [[Rcpp::export]]
Rcpp::List compute_conditional_probs(
    arma::imat observations,          // n x p matrix of observed data
    arma::ivec predict_vars,          // which variables to predict (0-based)
    arma::mat pairwise,           // p x p pairwise matrix
    arma::mat main,             // p x max_cat main effect matrix
    arma::ivec num_categories,         // number of categories per variable
    Rcpp::StringVector variable_type, // "ordinal" or "blume-capel" per variable
    arma::ivec baseline_category      // baseline for blume-capel variables
) {
  int no_persons = observations.n_rows;
  int num_variables = observations.n_cols;
  int num_predict_vars = predict_vars.n_elem;

  // Output is a list of probability matrices, one per predict variable
  Rcpp::List result(num_predict_vars);

  for(int pv = 0; pv < num_predict_vars; pv++) {
    int variable = predict_vars[pv];
    int n_cats = num_categories[variable] + 1;  // Include category 0

    // Compute rest scores for all persons at once (vectorized)
    arma::vec rest_scores(no_persons, arma::fill::zeros);

    for(int vertex = 0; vertex < num_variables; vertex++) {
      if(vertex == variable) continue;  // Skip the variable we're predicting

      arma::vec obs_col = arma::conv_to<arma::vec>::from(observations.col(vertex));

      if(std::string(variable_type[vertex]) != "blume-capel") {
        rest_scores += 2.0 * obs_col * pairwise(vertex, variable);
      } else {
        int ref = baseline_category[vertex];
        rest_scores += 2.0 * (obs_col - double(ref)) * pairwise(vertex, variable);
      }
    }

    // Use numerically stable probability computation
    arma::mat probs;

    if(std::string(variable_type[variable]) == "blume-capel") {
      int ref = baseline_category[variable];
      double lin_eff = main(variable, 0);
      double quad_eff = main(variable, 1);
      arma::vec bound;  // Will be computed inside

      probs = compute_probs_blume_capel(
        rest_scores,
        lin_eff,
        quad_eff,
        ref,
        num_categories[variable],
        bound
      );
    } else {
      // Regular ordinal variable
      // Extract main effect parameters for this variable
      arma::vec main_param = main.row(variable).head(num_categories[variable]).t();

      // Compute bounds for numerical stability: max exponent per person
      arma::vec bound(no_persons, arma::fill::zeros);
      for(int c = 0; c < num_categories[variable]; c++) {
        arma::vec exps = main_param[c] + (c + 1) * rest_scores;
        bound = arma::max(bound, exps);
      }

      probs = compute_probs_ordinal(
        main_param,
        rest_scores,
        bound,
        num_categories[variable]
      );
    }

    // Convert arma::mat to Rcpp::NumericMatrix
    Rcpp::NumericMatrix prob_mat(no_persons, n_cats);
    for(int i = 0; i < no_persons; i++) {
      for(int c = 0; c < n_cats; c++) {
        prob_mat(i, c) = probs(i, c);
      }
    }

    result[pv] = prob_mat;
  }

  return result;
}


// ============================================================================
//   Mixed MRF Conditional Prediction
// ============================================================================

// Compute conditional distributions for a mixed MRF.
//
// For discrete variables: P(x_s = c | x_{-s}, y) using the conditional OMRF.
// For continuous variables: E(y_j | y_{-j}, x) and SD(y_j | y_{-j}, x)
// using the conditional GGM.
//
// Takes discrete data (n x p integer), continuous data (n x q numeric), and
// 0-based predict_vars indices into the combined (p+q) variable list
// (0..p-1 discrete, p..p+q-1 continuous), plus the interaction blocks
// (pairwise_disc: p x p with zero diagonal, pairwise_cross: p x q,
// pairwise_cont: q x q negative-definite), mux (p x max_cats threshold /
// Blume-Capel parameters), muy (q-vector of means), and per-variable
// metadata (num_categories, "ordinal" or "blume-capel" type, baseline
// categories). Returns a list of prediction matrices, one per predicted
// variable: n x (num_cats+1) probabilities for discrete, n x 2 (mean, sd)
// for continuous.
// [[Rcpp::export]]
Rcpp::List compute_conditional_mixed(
    const arma::imat& x_observations,
    const arma::mat& y_observations,
    const arma::ivec& predict_vars,
    const arma::mat& pairwise_disc,
    const arma::mat& pairwise_cross,
    const arma::mat& pairwise_cont,
    const arma::mat& mux,
    const arma::vec& muy,
    const arma::ivec& num_categories,
    const Rcpp::StringVector& variable_type,
    const arma::ivec& baseline_category
) {
  int n = x_observations.n_rows;
  int p = x_observations.n_cols;
  int q = y_observations.n_cols;
  int num_predict_vars = predict_vars.n_elem;

  // Convert discrete to double (centered for rest-score computation)
  arma::mat x_dbl = arma::conv_to<arma::mat>::from(x_observations);

  Rcpp::List result(num_predict_vars);

  for (int pv = 0; pv < num_predict_vars; pv++) {
    int var_idx = predict_vars[pv];

    if (var_idx < p) {
      // --- Discrete variable: P(x_s = c | x_{-s}, y) ---
      int s = var_idx;
      int Cs = num_categories[s];

      // Rest score from discrete neighbours (centered by baseline)
      arma::vec rest_discrete(n, arma::fill::zeros);
      for (int k = 0; k < p; k++) {
        if (k == s) continue;
        arma::vec obs_k = x_dbl.col(k);
        double ref_k = static_cast<double>(baseline_category[k]);
        rest_discrete += 2.0 * (obs_k - ref_k) * pairwise_disc(k, s);
      }

      // Rest score from continuous (factor of 2)
      arma::vec rest_continuous(n, arma::fill::zeros);
      for (int j = 0; j < q; j++) {
        rest_continuous += 2.0 * pairwise_cross(s, j) * y_observations.col(j);
      }

      arma::vec rest_scores = rest_discrete + rest_continuous;

      arma::mat probs;
      if (std::string(variable_type[s]) == "blume-capel") {
        int ref = baseline_category[s];
        double lin_eff = mux(s, 0);
        double quad_eff = mux(s, 1);
        arma::vec bound;
        probs = compute_probs_blume_capel(
          rest_scores, lin_eff, quad_eff, ref, Cs, bound
        );
      } else {
        arma::vec main_param = mux.row(s).head(Cs).t();
        arma::vec bound(n, arma::fill::zeros);
        for (int c = 0; c < Cs; c++) {
          arma::vec exps = main_param[c] + (c + 1) * rest_scores;
          bound = arma::max(bound, exps);
        }
        probs = compute_probs_ordinal(main_param, rest_scores, bound, Cs);
      }

      Rcpp::NumericMatrix prob_mat(n, Cs + 1);
      for (int i = 0; i < n; i++) {
        for (int c = 0; c <= Cs; c++) {
          prob_mat(i, c) = probs(i, c);
        }
      }
      result[pv] = prob_mat;

    } else {
      // --- Continuous variable: y_j | y_{-j}, x ---
      int j = var_idx - p;

      // Continuous interaction matrix is negative-definite: precision = -2 * pairwise_cont
      double precision_jj = -2.0 * pairwise_cont(j, j);
      double cond_var = 1.0 / precision_jj;
      double cond_sd = std::sqrt(cond_var);

      // Contribution from other continuous variables:
      // -sum_{k != j} precision[j,k] * (y_k - muy_k)
      arma::vec lp_continuous(n, arma::fill::zeros);
      for (int k = 0; k < q; k++) {
        if (k == j) continue;
        lp_continuous += 2.0 * pairwise_cont(j, k) * (y_observations.col(k) - muy(k));
      }

      // Contribution from discrete variables (factor of 2):
      // sum_s 2 * cross_interaction(s, j) * x_s_centered
      arma::vec lp_discrete(n, arma::fill::zeros);
      for (int s = 0; s < p; s++) {
        double ref_s = static_cast<double>(baseline_category[s]);
        lp_discrete += 2.0 * pairwise_cross(s, j) * (x_dbl.col(s) - ref_s);
      }

      arma::vec cond_means = muy(j) + cond_var * (lp_continuous + lp_discrete);

      Rcpp::NumericMatrix out(n, 2);
      for (int i = 0; i < n; i++) {
        out(i, 0) = cond_means[i];
        out(i, 1) = cond_sd;
      }
      result[pv] = out;
    }
  }

  return result;
}