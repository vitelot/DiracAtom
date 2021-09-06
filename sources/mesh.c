#include <math.h>

void rmesh( double *r, double *drdi, int imax,
	    double rmin, double rmax,
	    int way)
{
  register i;
  double expdx, dx;
  
  switch (way) {
  case 0:
/*  Radial linear mesh */
    /* rmin ignored */
    dx = rmax/imax;
    for(i=0; i<imax; i++) {
      r[i] = dx*(i+1);
      drdi[i] = dx;
    }
    break;    
  case 1:
/*  Radial quadratic mesh */
    /* rmin ignored */
    dx = rmax/((double)imax*imax);
    for(i=0; i<imax; i++) {
      r[i] = dx*(double)(i+1)*(i+1);
      drdi[i] = 2.0*dx*(i+1);
    }
    break;
  default:
/*  Radial log-mesh */
    r[0] = rmin;
    dx = log(rmax/rmin)/(imax-1);
    drdi[0] = rmin * dx;
    expdx = exp(dx);
    for( i=1; i<imax; i++) {
      r[i] = expdx * r[i-1];
      drdi[i] = r[i]*dx;
    }
    break;
  }
}
