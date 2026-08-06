#include "math/cholupdate.h"

extern "C" {

// from mgcv: https://github.com/cran/mgcv/blob/1b6a4c8374612da27e36420b4459e93acb183f2d/src/mat.c#L1876-L1883
static inline double hypote(double x, double y) {
/* stable computation of sqrt(x^2 + y^2) */
  double t;
  x = fabs(x);y=fabs(y);
  if (y>x) { t = x;x = y; y = t;}
  if (x==0) return(y); else t = y/x;
  return(x*sqrt(1+t*t));
} /* hypote */

// from mgcv: https://github.com/cran/mgcv/blob/1b6a4c8374612da27e36420b4459e93acb183f2d/src/mat.c#L1956
void chol_up(double *R,double *u, int *n,int *up,double *eps) {
/* Rank 1 update of a cholesky factor. Works as follows:

   [up=1] R'R + uu' = [u,R'][u,R']' =  [u,R']Q'Q[u,R']', and then uses Givens rotations to
   construct Q such that Q[u,R']' = [0,R1']'. Hence R1'R1 = R'R + uu'. The construction
   operates from first column to last.

   [up=0] uses an almost identical sequence, but employs hyperbolic rotations
   in place of Givens. See Golub and van Loan (2013, 4e 6.5.4)

   Givens rotations are of form [c,-s] where c = cos(theta), s = sin(theta).
                                [s,c]

   Assumes R upper triangular, and that it is OK to use first two columns
   below diagonal as temporary strorage for Givens rotations (the storage is
   needed to ensure algorithm is column oriented).

   For downdate returns a negative value in R[1] (R[1,0]) if not +ve definite.
*/
  double c0,s0,*c,*s,z,*x,z0,*c1;
  int j,j1,n1;
  n1 = *n - 1;
  if (*up) for (j1=-1,j=0;j<*n;j++,u++,j1++) { /* loop over columns of R */
    z = *u; /* initial element of u */
    x = R + *n * j; /* current column */
    c = R + 2;s = R + *n + 2; /* Storage for first n-2 Givens rotations */
    for (c1=c+j1;c<c1;c++,s++,x++) { /* apply previous Givens */
      z0 = z;
      z = *c * z - *s * *x;
      *x = *s * z0 + *c * *x;
    }
    if (j) {
      /* apply last computed Givens */
      z0 = z;
      z = c0 * z - s0 * *x;
      *x = s0 * z0 + c0 * *x;
      x++;
      if (j<n1) {*c = c0; *s = s0;} /* store if needed for further columns */
    }

    /* now construct the next rotation u[j] <-> R[j,j] */
    z0 = hypote(z,*x); /* sqrt(z^2+R[j,j]^2) */
    c0 = *x/z0; s0 = z/z0;  /* need to zero z */
    /* now apply this rotation and this column is finished (so no need to update z) */
    *x = s0 * z + c0 * *x;
  } else  for (j1=-1,j=0;j<*n;j++,u++,j1++) { /* loop over columns of R for down-dating */
    z = *u; /* initial element of u */
    x = R + *n * j; /* current column */
    c = R + 2;s = R + *n + 2; /* Storage for first n-2 hyperbolic rotations */
    for (c1=c+j1;c<c1;c++,s++,x++) { /* apply previous hyperbolic */
      z0 = z;
      z =  *c * z - *s * *x;
      *x = -*s * z0 + *c * *x;
    }
    if (j) {
      /* apply last computed hyperbolic */
      z0 = z;
      z = c0 * z - s0 * *x;
      *x = -s0 * z0 + c0 * *x;
      x++;
      if (j<n1) {*c = c0; *s = s0;} /* store if needed for further columns */
    }

    /* now construct the next hyperbolic rotation u[j] <-> R[j,j] */
    z0 = z / *x; /* sqrt(z^2+R[j,j]^2) */
    if (fabs(z0)>=1) { /* downdate not +ve def */
      //Rprintf("j = %d  d = %g ",j,z0);
      if (*n>1) R[1] = -2.0;
      return; /* signals error */
    }
    if (z0 > 1 - *eps) z0 = 1 - *eps;
    c0 = 1/sqrt(1-z0*z0);s0 = c0 * z0;
    /* now apply this rotation and this column is finished (so no need to update z) */
    *x = -s0 * z + c0 * *x;
  }

 /* now zero c and s storage */
  c = R + 2;s = R + *n + 2;
  for (x = c + *n - 2;c<x;c++,s++) *c = *s = 0.0;
} /* chol_up */
}

// for internal use
void cholesky_update(arma::mat& R, arma::vec& u, double eps) {
    int n = R.n_cols;
    int up = 1;
    chol_up(R.memptr(), u.memptr(), &n, &up, &eps);
}

bool cholesky_downdate(arma::mat& R, arma::vec& u, double eps) {
    int n = R.n_cols;
    if (n == 1) {
        // chol_up cannot signal failure for n = 1; handle directly.
        double d = R(0, 0) * R(0, 0) - u(0) * u(0);
        if (d <= 0.0) return false;
        R(0, 0) = std::sqrt(d);
        return true;
    }
    int up = 0;
    chol_up(R.memptr(), u.memptr(), &n, &up, &eps);
    // chol_up signals a non-positive-definite downdate by writing -2.0 into
    // the (1,0) sub-diagonal slot and returning early; on success that slot
    // is always 0.
    return R.memptr()[1] != -2.0;
}
