#include "relat.h"

/* POTAT computes the potential function from the charge density.
 * called by main.
 * used subroutine: 'potcul', 'simpson'.
 */     

void potat( double *rad, double *drdi, double zn, 
            double *rho, double *rhosp,
            double *vt, double *bt, double *vnuc,
	    double *Vxc,
	    double v0[2], double *q,
	    double *ev, double *evxc, double *eexc
	  )
/**********************************************************************
 INPUT:
  opt.nspin	nspin FALSE for full coupling. Para case otherwise.
  *vnuc		nuclear potential (usually -2Z/r)
  *rad		Radial mesh
  zn		Atomic number (nuclear charge)
  *rho		|Psi|^2
  *rhosp	Spin density
 OUTPUT:
  *Vxc		Spherical Exchange potential
  *vt		Spherical Total Potential (Coul+XC+Nucl) 
  *bt		Spherical Exchange field
  q		Total Charge
  ev		e-e Energy
  evxc		Exchange-Correlation Potential integral		
  eexc		   ,,       ,,       Energy
 *********************************************************************/
{
register i;
Dbp ur[2], bgx[2], y, rhor2, vc, aa, bb;
double temp, v00, ev4, r, rhoup, rhodn;
int is, ismax;

#include "alloc_potat.h"

  if( opt.nspin) {
    for( i=0; i<nrmax; i++) {
      r = rad[i];
      rhor2[i] = r*r*rho[i];
    }
  } else {
    for( i=0; i<nrmax; i++) {
      r = rad[i];
      rhor2[i] = r*r*rho[i];
    }
  }

  potcul_ryd( rad, drdi, zn,
	      rho, &v00, vc,
	      &ev4, ev, q );
  /* if the nucleus is charged then add its potential */
  if(zn>1e-5)
    for( i=0; i<nrmax; i++)
      vc[i] += rad[i]*vnuc[i];

  if ( opt.nspin) {
    for( i=0; i<nrmax; i++) {
      r = rad[i];
      temp = rho[i]/M_4PI;

      bgx[0][i] =    vxc( temp, 0., 0) *2.0; /* *2.0 Ryd! */
      y[i]      =  epsxc( temp, 0., 0) *2.0; /* *2.0 Ryd! */

      aa[i] = rhor2[i] * (-bgx[0][i]) * drdi[i];
      bb[i] = rhor2[i] * y[i] * drdi[i];
    }
    ismax = 0;
  } else {
    for( i=0; i<nrmax; i++) {
      r = rad[i];

      rhoup = .5*(rho[i] + rhosp[i])/M_4PI;
      rhodn = .5*(rho[i] - rhosp[i])/M_4PI;

      bgx[0][i] =    vxc( rhoup, rhodn,  1) *2.0; /* Ryd! */
      bgx[1][i] =    vxc( rhoup, rhodn, -1) *2.0;
      y[i]      =  epsxc( rhoup, rhodn,  1) *2.0;

      aa[i] = ( rhoup*(-bgx[0][i]) + rhodn*(-bgx[1][i]) ) *
	M_4PI * r*r * drdi[i];
      bb[i] = rhor2[i] * y[i] * drdi[i];
    }
    ismax = 1;
  }

  *evxc = simpson(aa, nrmax, 1.0);
  *eexc = simpson(bb, nrmax, 1.0);
  
  for( is=0; is <= ismax; is++) {
    v0[is] = v00 + bgx[is][0];
    for( i=0; i<nrmax; i++) {
      r = rad[i];
      ur[is][i] = vc[i]/r + bgx[is][i];
    }
  }
  
   if ( opt.nspin)
   {
     for( i=0; i<nrmax; i++)
     {
       Vxc[i] = bgx[0][i];
       vt[i] = ur[0][i];
       bt[i] = 0.0;
     }
   } else {
     for( i=0; i<nrmax; i++)
     {
       r = rad[i];
       Vxc[i] = .5*( bgx[0][i]+bgx[1][i] );
       vt[i] = .5*( ur[0][i]+ur[1][i] );
       bt[i] = .5*( ur[0][i]-ur[1][i] );
     }
   }

#include "free_potat.h"
   
}
