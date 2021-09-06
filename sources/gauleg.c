/* GAULEG: Given the lower and upper limits of integration x1, x2, and
           given n, returns arrays x and w of lenght n containing the
	   abscissas and weights of the Gauss-Legendre N-point quadrature 
	   formula.
     cfr Num recipes pg 125   
*/

#include "my_math.h"

void gauleg( double x1, double x2, double *x, double *w, int n)
{
  const double eps=3.0e-14;
  double xm, xl, z, p1, p2, p3, pp, z1;
  int i, j, m;

    m = (n+1)/2;	/* The roots are symm in the interval, so we find 
			   only half of them. */
    xm = .5*(x2+x1);
    xl = .5*(x2-x1);
    
    for( i=0; i<m; i++) {
      z = cos( M_PI * (i+.75)/(n+.5)); /* Starting with this approx to the
					  i-th root we refine with Newton */
      do {
	p1 = 1.;
	p2 = .0;
	for( j=1; j<=n; j++) {
	  p3 = p2;
	  p2 = p1;
	  p1 = ((2.*j-1.)*z*p2-(j-1.)*p3)/j;
	}
	/* p1 is the desired Legendre Polinomial of degree n */

	pp = n*(z*p1-p2)/(z*z-1.);	/* pp is its derivative */
	z1 = z;
	z = z1 - p1/pp;			/* Newton's method */
      } while( fabs(z-z1) > eps);

      x[i] = xm - z*xl;		/* Scale x to the desired interval */
      x[n-1-i] = xm + z*xl;     /* Its symm counterpart */

      w[n-1-i] = w[i] = 2./((1.-z*z)*pp*pp);
    }
}
