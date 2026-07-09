/**
 * @file custom_arma_explog.h
 * @brief Element-wise exp/log for Armadillo matrices using OpenLibM.
 *
 * Drop-in replacements for `arma::exp` and `arma::log` that call the
 * portable OpenLibM implementations from `custom_explog.h` instead of
 * the platform's libm. Used on Windows where MSVC's `std::exp` /
 * `std::log` are slower. The `ARMA_MY_EXP` /
 * `ARMA_MY_LOG` macros in `explog_macros.h` resolve to these functions
 * on Windows and to `arma::exp` / `arma::log` elsewhere.
 *
 * @see custom_explog.h    Scalar OpenLibM exp/log
 * @see explog_macros.h    Platform-conditional macro definitions
 */
#ifndef BGMS_CUSTOM_ARMA_EXPLOG_H
#define BGMS_CUSTOM_ARMA_EXPLOG_H

#include "RcppArmadillo.h"

double __ieee754_exp(double x); // forward declaration
double __ieee754_log(double x); // forward declaration

/**
 * Element-wise exponential of an Armadillo matrix expression.
 *
 * Evaluates the input expression, then applies `__ieee754_exp` to
 * every element.
 *
 * @tparam T1  Armadillo expression type
 * @param  X   Matrix expression
 * @return Matrix of the same dimensions with e^x_ij
 */
template<typename T1>
arma::Mat<double> custom_arma_exp(const arma::Base<double, T1>& X)
{
  arma::Mat<double> Xin = X.get_ref();
  arma::Mat<double> out(Xin.n_rows, Xin.n_cols, arma::fill::none);

  const double* in_mem  = Xin.memptr();
  double* out_mem       = out.memptr();
  const arma::uword N   = Xin.n_elem;

  for (arma::uword i = 0; i < N; ++i)
    out_mem[i] = __ieee754_exp(in_mem[i]);

  return out;
}

/**
 * Element-wise natural logarithm of an Armadillo matrix expression.
 *
 * Evaluates the input expression, then applies `__ieee754_log` to
 * every element.
 *
 * @tparam T1  Armadillo expression type
 * @param  X   Matrix expression (all elements must be positive)
 * @return Matrix of the same dimensions with ln(x_ij)
 */
template<typename T1>
arma::Mat<double> custom_arma_log(const arma::Base<double, T1>& X)
{
    arma::Mat<double> Xin = X.get_ref();
    arma::Mat<double> out(Xin.n_rows, Xin.n_cols, arma::fill::none);

    double* out_mem       = out.memptr();
    const double* in_mem  = Xin.memptr();
    const arma::uword N   = Xin.n_elem;

    for (arma::uword i = 0; i < N; ++i)
        out_mem[i] = __ieee754_log(in_mem[i]);

    return out;
}


#endif
