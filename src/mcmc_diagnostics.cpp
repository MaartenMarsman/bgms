// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>


// ============================================================================
//   ESS via AR spectral density (matching coda::effectiveSize)
// ============================================================================
//
// Algorithm per chain per parameter:
//   1. Compute biased autocovariance c[0..L], L = floor(10*log10(n))
//   2. Levinson-Durbin recursion for AR orders 1..L
//   3. AIC order selection: argmin_k { n*log(v_k) + 2*k }
//   4. Spectral density at 0: spec0 = v_best / (1 - sum(a_best))^2
//   5. ESS = n * var_unbiased / spec0
//
// Multi-chain ESS = sum of per-chain ESS (same as coda for mcmc.list).
//
// Origin: Plummer et al. (2006), coda package, spectrum0.ar + effectiveSize.
// ============================================================================


// --- Single-column ESS (called from worker) ----------------------------------

static double compute_column_ess(const double* x, int n, int max_order) {
  // Need at least 2 observations to estimate autocorrelation
  if(n <= 1) return NA_REAL;

  // Step 1: mean
  double mean = 0.0;
  for(int i = 0; i < n; i++) mean += x[i];
  mean /= n;

  // Guard against non-finite input (NaN / Inf)
  if(!std::isfinite(mean)) return NA_REAL;

  // Step 2: biased autocovariance c[0..max_order]
  std::vector<double> acov(max_order + 1, 0.0);
  for(int lag = 0; lag <= max_order; lag++) {
    double s = 0.0;
    for(int i = 0; i < n - lag; i++) {
      s += (x[i] - mean) * (x[i + lag] - mean);
    }
    acov[lag] = s / n;
  }

  // Constant chain: no autocorrelation structure to estimate, so it carries no
  // effective draws. Zero rather than NA, so that pooling over chains keeps the
  // ESS of the chains that do vary (a near-saturated edge indicator is often
  // constant in some chains and not in others).
  if(acov[0] < 1e-15) return 0.0;

  // Step 3: Levinson-Durbin + AIC selection
  double best_aic = n * std::log(acov[0]); // AIC for order 0
  double best_var = acov[0];
  double best_sum_ar = 0.0;
  int best_order = 0;

  std::vector<double> a_prev, a_curr;
  double v = acov[0];

  for(int k = 1; k <= max_order; k++) {
    // Reflection coefficient
    double num = acov[k];
    for(int i = 0; i < k - 1; i++) {
      num -= a_prev[i] * acov[k - 1 - i];
    }
    double lambda = num / v;

    // Update AR coefficients
    a_curr.resize(k);
    a_curr[k - 1] = lambda;
    for(int i = 0; i < k - 1; i++) {
      a_curr[i] = a_prev[i] - lambda * a_prev[k - 2 - i];
    }

    // Update prediction variance
    v *= (1.0 - lambda * lambda);
    if(v <= 0.0) break; // Numerical guard

    // AIC for order k
    double aic_k = n * std::log(v) + 2.0 * k;
    if(aic_k < best_aic) {
      best_aic = aic_k;
      best_var = v;
      best_order = k;
      best_sum_ar = 0.0;
      for(int i = 0; i < k; i++) best_sum_ar += a_curr[i];
    }

    a_prev.swap(a_curr);
  }

  // Step 4: spectral density at frequency 0
  // Guard: df correction requires n - order - 1 > 0; if the selected order
  // is too high (can happen when n is tiny), fall back to order 0.
  if(n - best_order - 1 <= 0) {
    best_order = 0;
    best_var = acov[0];
    best_sum_ar = 0.0;
  }
  // Apply df correction: var.pred * n / (n - order - 1), matching R's ar.yw
  double df_corrected_var = best_var * n / (n - best_order - 1.0);
  double denom = 1.0 - best_sum_ar;
  double spec0 = df_corrected_var / (denom * denom);

  // Step 5: ESS = n * var_unbiased / spec0
  double var_unbiased = acov[0] * n / (n - 1.0);
  return (spec0 > 0.0 && std::isfinite(spec0)) ? n * var_unbiased / spec0 : NA_REAL;
}


// --- RcppParallel worker for multi-parameter ESS -----------------------------

struct ESSWorker : public RcppParallel::Worker {
  // Input: 3D array stored as [niter x nchains x nparam] in column-major
  const double* data;
  const int niter;
  const int nchains;
  const int nparam;
  const int max_order;

  // Output
  RcppParallel::RVector<double> ess;

  ESSWorker(const double* data, int niter, int nchains, int nparam,
            int max_order, Rcpp::NumericVector ess)
    : data(data), niter(niter), nchains(nchains), nparam(nparam),
      max_order(max_order), ess(ess) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++) {
      double total_ess = 0.0;
      for(int c = 0; c < nchains; c++) {
        // Column-major 3D: element [i, c, j] is at data[i + c*niter + j*niter*nchains]
        const double* col = data + c * niter + j * niter * nchains;
        total_ess += compute_column_ess(col, niter, max_order);
      }
      // A zero total means every chain was constant; a non-finite one means a
      // chain carried non-finite draws. Neither yields an effective sample size.
      ess[j] = (total_ess > 0.0) ? total_ess : NA_REAL;
    }
  }
};


// ============================================================================
//   Gelman-Rubin Rhat (classic split-Rhat, Gelman et al. 2013 / Stan)
// ============================================================================
//
// Rhat = sqrt(var_plus / W), with var_plus = (n-1)/n * W + B/n, where W is the
// mean within-sub-chain variance and B the between-sub-chain variance. The
// R-level split_chains() halves each chain first, so this is the split-Rhat of
// Gelman et al. (2013) / Vehtari et al. (2021).
//
// The Brooks-Gelman degrees-of-freedom adjustment used by coda::gelman.diag
// (Rhat = sqrt(df_adj * R2), df_adj = (df_V + 3)/(df_V + 1)) is intentionally
// NOT applied. On a nearly-saturated binary indicator chain -- one sub-chain
// carrying a brief excursion while the rest are constant -- df_V collapses to
// exactly 2, so df_adj -> 5/3 and Rhat -> sqrt(5/3) ~ 1.291 regardless of the
// draws, flagging the most decisive edges as unconverged.
// ============================================================================

struct RhatWorker : public RcppParallel::Worker {
  const double* data;
  const int niter;
  const int nchains;
  const int nparam;

  RcppParallel::RVector<double> rhat;

  RhatWorker(const double* data, int niter, int nchains, int nparam,
             Rcpp::NumericVector rhat)
    : data(data), niter(niter), nchains(nchains), nparam(nparam),
      rhat(rhat) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n = niter;
    int m = nchains;

    for(std::size_t j = begin; j < end; j++) {
      // Compute per-chain means and variances
      double grand_mean = 0.0;
      std::vector<double> chain_mean(m), chain_var(m);

      for(int c = 0; c < m; c++) {
        const double* col = data + c * n + j * n * m;
        double s = 0.0;
        for(int i = 0; i < n; i++) s += col[i];
        chain_mean[c] = s / n;
        grand_mean += chain_mean[c];

        double s2 = 0.0;
        for(int i = 0; i < n; i++) {
          double d = col[i] - chain_mean[c];
          s2 += d * d;
        }
        chain_var[c] = s2 / (n - 1.0);
      }
      grand_mean /= m;

      // Guard: if any chain summary is non-finite (NaN/Inf in input),
      // skip this parameter.
      if(!std::isfinite(grand_mean)) {
        rhat[j] = NA_REAL;
        continue;
      }

      // W = mean of within-chain variances
      double W = 0.0;
      for(int c = 0; c < m; c++) W += chain_var[c];
      W /= m;

      // B = n * variance of chain means (unbiased)
      double B = 0.0;
      for(int c = 0; c < m; c++) {
        double d = chain_mean[c] - grand_mean;
        B += d * d;
      }
      B *= n / (m - 1.0);

      // --- Classic split-Rhat (Gelman et al. 2013 / Stan) ---
      // var_plus = (n-1)/n * W + B/n;  Rhat = sqrt(var_plus / W).
      if(W > 0.0) {
        double var_plus = (n - 1.0) / n * W + B / n;
        rhat[j] = std::sqrt(var_plus / W);
      } else {
        // W == 0: every sub-chain is individually constant.
        //   B == 0 -> all sub-chains hold the SAME constant; Rhat is undefined
        //            (0/0), report NA.
        //   B  > 0 -> sub-chains are constant at DIFFERENT values (chains stuck
        //            in different states); the between/within ratio diverges, so
        //            report +Inf as a loud non-convergence alarm rather than the
        //            silent NA that would hide genuinely stuck chains.
        rhat[j] = (B > 0.0) ? R_PosInf : NA_REAL;
      }
    }
  }
};


// ============================================================================
//   Rcpp exports
// ============================================================================

// Compute ESS for a 3D array [niter x nchains x nparam].
// Multi-chain ESS = sum of per-chain ESS.
// [[Rcpp::export(.compute_ess_cpp)]]
Rcpp::NumericVector compute_ess_cpp(Rcpp::NumericVector array3d) {
  Rcpp::IntegerVector dims = array3d.attr("dim");
  int niter   = dims[0];
  int nchains = dims[1];
  int nparam  = dims[2];

  if(niter <= 1) {
    return Rcpp::NumericVector(nparam, NA_REAL);
  }

  int max_order = std::min(niter - 1,
                           (int)std::floor(10.0 * std::log10((double)niter)));

  Rcpp::NumericVector ess(nparam);
  ESSWorker worker(array3d.begin(), niter, nchains, nparam, max_order, ess);
  RcppParallel::parallelFor(0, nparam, worker);
  return ess;
}


// Compute Rhat for a 3D array [niter x nchains x nparam].
// Returns NA for single-chain input.
// [[Rcpp::export(.compute_rhat_cpp)]]
Rcpp::NumericVector compute_rhat_cpp(Rcpp::NumericVector array3d) {
  Rcpp::IntegerVector dims = array3d.attr("dim");
  int niter   = dims[0];
  int nchains = dims[1];
  int nparam  = dims[2];

  if(nchains < 2 || niter <= 1) {
    return Rcpp::NumericVector(nparam, NA_REAL);
  }

  Rcpp::NumericVector rhat(nparam);
  RhatWorker worker(array3d.begin(), niter, nchains, nparam, rhat);
  RcppParallel::parallelFor(0, nparam, worker);
  return rhat;
}


// ============================================================================
//   Binary indicator transition ESS
// ============================================================================
//
// For binary (0/1) MCMC draws, ESS is computed from transition counts
// rather than AR spectral density. For each parameter, pool draws across
// chains into a single vector and count consecutive-pair transitions:
//   n00, n01, n10, n11
// Then:
//   a = n01 / (n00 + n01)
//   b = n10 / (n10 + n11)
//   tau_int = (2 - a - b) / (a + b)
//   n_eff = n_total / tau_int
//
// Returns an [nparam x 8] matrix with columns:
//   mean, sd, mcse, n00, n01, n10, n11, n_eff_mixt
//
// Origin: two-state Markov chain integrated autocorrelation time.
// ============================================================================

struct IndicatorESSWorker : public RcppParallel::Worker {
  const double* data;
  const int niter;
  const int nchains;
  const int nparam;

  // Output: nparam x 8 matrix in column-major order
  RcppParallel::RVector<double> out;

  IndicatorESSWorker(const double* data, int niter, int nchains, int nparam,
                     Rcpp::NumericMatrix out)
    : data(data), niter(niter), nchains(nchains), nparam(nparam),
      out(out) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_total = niter * nchains;

    for(std::size_t j = begin; j < end; j++) {
      // Pool draws across chains: column-major 3D layout
      // element [i, c, j] at data[i + c*niter + j*niter*nchains]
      // Pool order: chain 0 iter 0..niter-1, chain 1 iter 0..niter-1, ...

      // Single pass: accumulate sum and transition counts
      double sum_x = 0.0;
      int c00 = 0, c01 = 0, c10 = 0, c11 = 0;
      bool has_nonfinite = false;

      // First element (no transition from previous)
      const double* base = data + j * niter * nchains;

      // Guard: if the first value is non-finite, skip computation
      if(!std::isfinite(base[0])) {
        has_nonfinite = true;
      }

      int prev = 0;
      if(!has_nonfinite) {
        prev = (int)base[0];
        sum_x += prev;

        for(int c = 0; c < nchains && !has_nonfinite; c++) {
          const double* col = base + c * niter;
          int start = (c == 0) ? 1 : 0;
          for(int i = start; i < niter; i++) {
            double raw = col[i];
            if(!std::isfinite(raw)) {
              has_nonfinite = true;
              break;
            }
            int curr = (int)raw;
            sum_x += curr;
            // Transition from prev to curr
            if(prev == 0) {
              if(curr == 0) c00++; else c01++;
            } else {
              if(curr == 0) c10++; else c11++;
            }
            prev = curr;
          }
        }
      }

      if(has_nonfinite) {
        for(int k = 0; k < 8; k++) out[j + k * nparam] = NA_REAL;
        continue;
      }

      double p_hat = sum_x / n_total;
      double sd = std::sqrt(p_hat * (1.0 - p_hat));

      double n_eff_mixt, mcse;
      if(c01 + c10 == 0) {
        n_eff_mixt = NA_REAL;
        mcse = NA_REAL;
      } else {
        double a = (double)c01 / (c00 + c01);
        double b = (double)c10 / (c10 + c11);
        double tau_int = (2.0 - a - b) / (a + b);
        n_eff_mixt = n_total / tau_int;
        mcse = (n_eff_mixt > 0.0) ? sd / std::sqrt(n_eff_mixt) : NA_REAL;
      }

      // Store in column-major matrix: out[j + col * nparam]
      out[j + 0 * nparam] = p_hat;
      out[j + 1 * nparam] = sd;
      out[j + 2 * nparam] = mcse;
      out[j + 3 * nparam] = (double)c00;
      out[j + 4 * nparam] = (double)c01;
      out[j + 5 * nparam] = (double)c10;
      out[j + 6 * nparam] = (double)c11;
      out[j + 7 * nparam] = n_eff_mixt;
    }
  }
};


// Compute indicator transition ESS for a 3D array [niter x nchains x nparam].
// Returns an [nparam x 8] matrix with columns:
//   mean, sd, mcse, n00, n01, n10, n11, n_eff_mixt
// [[Rcpp::export(.compute_indicator_ess_cpp)]]
Rcpp::NumericMatrix compute_indicator_ess_cpp(Rcpp::NumericVector array3d) {
  Rcpp::IntegerVector dims = array3d.attr("dim");
  int niter   = dims[0];
  int nchains = dims[1];
  int nparam  = dims[2];

  Rcpp::NumericMatrix out(nparam, 8);
  Rcpp::colnames(out) = Rcpp::CharacterVector::create(
    "mean", "sd", "mcse", "n00", "n01", "n10", "n11", "n_eff_mixt"
  );

  if(niter <= 1) {
    // With 0 or 1 iterations we can still compute mean/sd/counts but
    // transition-based mcse and n_eff_mixt are undefined.
    if(niter == 0) {
      std::fill(out.begin(), out.end(), NA_REAL);
      return out;
    }
    // niter == 1: compute mean, sd, zero transitions; mcse/n_eff_mixt = NA
    for(int j = 0; j < nparam; j++) {
      double sum_x = 0.0;
      int n_total = nchains; // niter * nchains with niter == 1
      for(int c = 0; c < nchains; c++) {
        sum_x += array3d[c * niter + j * niter * nchains];
      }
      double p_hat = sum_x / n_total;
      out[j + 0 * nparam] = p_hat;
      out[j + 1 * nparam] = std::sqrt(p_hat * (1.0 - p_hat));
      out[j + 2 * nparam] = NA_REAL; // mcse
      out[j + 3 * nparam] = 0.0;     // n00
      out[j + 4 * nparam] = 0.0;     // n01
      out[j + 5 * nparam] = 0.0;     // n10
      out[j + 6 * nparam] = 0.0;     // n11
      out[j + 7 * nparam] = NA_REAL; // n_eff_mixt
    }
    return out;
  }

  IndicatorESSWorker worker(array3d.begin(), niter, nchains, nparam, out);
  RcppParallel::parallelFor(0, nparam, worker);
  return out;
}
