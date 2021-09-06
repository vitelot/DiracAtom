#include "relat.h"

double iexc_nc( Grid2d Rho, Grid2d Mx, Grid2d Mz, 
		double *r2drdi,
	        double *weights)
{
  int i, n;
  double a,b,c,s,r, nup, ndn, enxc, etot, *integrand;

  integrand = vector_alloc(nrmax);

  for( etot=i=0; i<nlmax; i++)
  {
    for( n=0; n<nrmax; n++)
    {
      r = Rho[i][n];
      a = .5*( r + Mz[i][n]);
      b = .5*( r - Mz[i][n]);
      c = .5*Mx[i][n];

      s = sqrt( (a-b)*(a-b)+4.0*c*c);
      nup = .5* ( r + s);
      ndn = .5* ( r - s);
      
      enxc = epsxc( nup, ndn, 1) *2.0; /* Ryd! */

      integrand[n] = r * enxc * r2drdi[n];	/* LDA */
    }
    etot += weights[i] * (simpson( integrand, nrmax, 1.));
  } 

  vector_free(integrand);

return etot * 2.*M_PI; /* 2Pi is for phi integration */
}
	        
double ivexc_nc( Grid2d Rho, Grid2d Mx, Grid2d Mz, 
		 Grid2d   V, Grid2d Bx, Grid2d Bz,
		 double *r2drdi,
		 double *weights)
{
  int i,n;
  double etot, *integrand;

  integrand = vector_alloc(nrmax);

  for( etot=i=0; i<nlmax; i++) {
    for( n=0; n<nrmax; n++) {
      integrand[n] =
	(Rho[i][n]*V[i][n] + Bx[i][n]*Mx[i][n] + Bz[i][n]*Mz[i][n])*r2drdi[n];
    }
    etot += weights[i] * simpson( integrand, nrmax, 1.);
  } 

  vector_free(integrand);

return -etot * 2.*M_PI; /* 2Pi is for phi integration */
}
