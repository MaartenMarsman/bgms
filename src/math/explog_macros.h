/**
 * @file explog_macros.h
 * @brief Platform-conditional exp/log macro selection.
 *
 * Defines `MY_EXP`, `MY_LOG` (scalar) and `ARMA_MY_EXP`, `ARMA_MY_LOG`
 * (element-wise matrix) macros that resolve to either the portable
 * OpenLibM implementations or the standard-library versions.
 *
 * Default behaviour:
 *   - **Windows:** Uses OpenLibM (`__ieee754_exp`, `custom_arma_exp`, etc.)
 *     because MSVC's `std::exp` / `std::log` are significantly slower.
 *   - **macOS / Linux:** Uses `std::exp` / `arma::exp` etc.
 *
 * Override at build time:
 * @code
 * Sys.setenv("PKG_CPPFLAGS" = "-DCUSTOM_EXP_LOG=1")  # force OpenLibM
 * Sys.setenv("PKG_CPPFLAGS" = "-DCUSTOM_EXP_LOG=0")  # force std
 * @endcode
 *
 * @see custom_explog.h       Scalar OpenLibM exp/log
 * @see custom_arma_explog.h  Armadillo element-wise wrappers
 */
#ifndef _EXPLOG_MACROS_H_
#define	_EXPLOG_MACROS_H_

// An explicit -DCUSTOM_EXP_LOG value wins on every platform (=1 OpenLibM,
// =0 std); with no value, Windows defaults to OpenLibM and others to std.
#if (defined(CUSTOM_EXP_LOG) && CUSTOM_EXP_LOG != 0) \
  || (!defined(CUSTOM_EXP_LOG) && defined(_WIN32))

#define USE_CUSTOM_LOG 1

#else

#define USE_CUSTOM_LOG 0

#endif

#if USE_CUSTOM_LOG

#include "math/custom_explog.h"
#include "math/custom_arma_explog.h"

#define MY_EXP __ieee754_exp
#define MY_LOG __ieee754_log

#define ARMA_MY_EXP custom_arma_exp
#define ARMA_MY_LOG custom_arma_log

// log1p has no OpenLibM port shipped here yet; keep std::log1p universally.

#define MY_LOG1P std::log1p


#else

#define MY_EXP std::exp
#define MY_LOG std::log
#define MY_LOG1P std::log1p
#define ARMA_MY_EXP arma::exp
#define ARMA_MY_LOG arma::log

#endif

#endif