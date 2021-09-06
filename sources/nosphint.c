#include "relat.h"

void nsphint(double *vl, double *vlp, double *rhol,
	     double *rad, double *drdi, int l)
{
register int i;
double r, temp, *a, *b, rl, fact;
/* FILE *fp; */
/* char filename[128]; */

  a = vector_alloc(nrmax);
  b = vector_alloc(nrmax);

  /* REM: dr = dr/di * di with di=1			*/
  for( i=0; i<nrmax; i++) {
    r = rad[i];
    b[i] = rhol[i] * pow(r, 1.-l) * drdi[i];
    a[i] = rhol[i] * pow(r, 2.+l) * drdi[i];
  }

  bodeint( &a[0], nrmax, 1.0, 'o', &a[0]);
  bodeint( &b[0], nrmax, 1.0, 'i', &b[0]);

  temp = 1. / (2*l+1);

  for( i=0; i<nrmax; i++)
  {
    r = rad[i];
    rl = 1.;
    if( l)
      rl = pow( r, (double) l); 
    vl[i]  = temp * ( 1./(rl*r)*a[i] + rl*b[i] );
    vlp[i] = temp * ( -(l+1.)/(rl*r*r)*a[i] + l*rl*b[i]/r );
  }

  for( fact=1., i=l; i>1; i--)
    fact *= i;

  printf("\t[L=%1d]-pole = %.12lg (a.u.)\n", l, temp*a[nrmax-1]*fact);

/*   sprintf( filename, "AB%02d", l); */
/*   fp = fopen( filename, "w"); */
/*   for( i=0; i<nrmax; i++) */
/*     { */
/*       fprintf( fp, "%lg %lg %lg\n", rad[i], a[i], b[i]); */
/*     } */
/*   fclose(fp); */
/*   printf("\nFile %s saved\n", filename); */
  
  vector_free(a);
  vector_free(b);
}
