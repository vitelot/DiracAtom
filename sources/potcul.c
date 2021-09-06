/*
  POTCUL computes of culomb part of potential function, multiplied
  by r, its value for r=0, e-n potential energy,
  e-e potential energy and integral of charge density.
  internal calculated in rydberg units   
  called by 'potat'.
*/

#include "relat.h"

void potcul_ryd( double *rad, double *drdi, double z,
                 double *rho, double *v0, double *ur,
                 double *zb, double *ev, double *q)
/**********************************************************************
 INPUT:
  *rad		Radial mesh
  *drdi		dr(i)/di
  z		Atomic number (nuclear charge)
  *rho		Charge density (not multiplied by r*r)
 OUTPUT:
  v0		Coulomb part of potential evaluated at R=0
  *ur		R* Coulomb part of potential, R* (-z/R) not included
  zb		Potential energy e-nucl.
  ev		Potential energy e-e
  q		Total Charge
 *********************************************************************/
{
  register i;
  double *a, *b, r;

  a = vector_alloc(nrmax);
  b = vector_alloc(nrmax);

  /* REM: dr = dr/di * di with di=1			*/
  for( i=0; i<nrmax; i++) {
    r = rad[i];
    b[i] = r*rho[i] * drdi[i];
    a[i] = r*b[i];
  }

  *q  = bodeint( &a[0], nrmax, 1.0, 'o', &a[0]);
  *v0 = bodeint( &b[0], nrmax, 1.0, 'i', &b[0]) * 2.0;

  *zb = -z * *v0;

  for( i=0; i<nrmax; i++) {
    ur[i] = 2.0 * (a[i]+ rad[i]*b[i]);
    b[i] = .25*ur[i] *rad[i]*rho[i] * drdi[i]; /* its integral gives e-e en */
  }

  *ev = 2.0 * simpson(b, nrmax, 1.0);

  vector_free(a);
  vector_free(b);

}
