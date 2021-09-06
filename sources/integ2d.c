#include "relat.h"

double integ2d( double **A, double *r2drdi, double *weights)
{
  register int i,j;
  double *integrand, temp;

  integrand = vector_alloc(nrmax);

  temp = 0.0;
  LLoop(i) {
    RLoop(j) {
      integrand[j] = r2drdi[j] * A[i][j];
    }
    temp += weights[i] * simpson(integrand, nrmax, 1.);
  }

  vector_free(integrand);
  return 2.0*M_PI*temp;
}

