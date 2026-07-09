#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include "math/explog_macros.h"

/**
 * Numerically stable log(exp(a) + exp(b)).
 *
 * Guards -inf inputs (a contributor that is -inf cannot change the sum) and
 * uses the standard m + log1p(exp(-|a - b|)) form to avoid overflow in the
 * exponent and catastrophic cancellation in the log.
 */
inline double log_sum_exp(double a, double b) {
  constexpr double neg_inf = -std::numeric_limits<double>::infinity();
  if (a == neg_inf) return b;
  if (b == neg_inf) return a;
  const double m = std::max(a, b);
  return m + MY_LOG1P(MY_EXP(-std::abs(a - b)));
}
