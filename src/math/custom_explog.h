/**
 * @file custom_explog.h
 * @brief Portable scalar exp() and log() from OpenLibM.
 *
 * The MSVC runtime's `std::exp` and `std::log` are slower
 * than the glibc or macOS libm implementations. This caused bgms to run
 * much slower on Windows. To eliminate the platform gap, the package
 * ships stripped-down versions of the IEEE 754 exp and log from
 * OpenLibM (originally Sun Microsystems / FreeBSD libm). These are
 * used on Windows by default; see `explog_macros.h` for the
 * platform-conditional switching logic.
 *
 * Source: https://github.com/JuliaMath/openlibm  (e_exp.c, e_log.c)
 *
 * @see custom_arma_explog.h  Armadillo element-wise wrappers
 * @see explog_macros.h       Platform-conditional MY_EXP / MY_LOG macros
 */
#ifndef _CUSTOM_EXPLOG_H_
#define	_CUSTOM_EXPLOG_H_

/**
 * IEEE 754 double-precision exponential (OpenLibM).
 * @param x  Exponent
 * @return e^x
 */
double __ieee754_exp(double x);

/**
 * IEEE 754 double-precision natural logarithm (OpenLibM).
 * @param x  Argument (must be positive)
 * @return ln(x)
 */
double __ieee754_log(double x);

#endif
