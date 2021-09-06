#include "relat.h"

/* DENTOPOT: given a 2*2 spin density matrix, returns a 2*2
	     potential after applying the exchange operator
INPUT:
  light		speed of light
  *Rho		Spatial density
  *Mx		X component of magnetization
  *My		Y component of magnetization
OUTPUT:
  *V		XC non-spherical potential
  *Bx		x component of non-sph XC-field 
  *By		y component of non-sph XC-field 
*/

void dentopot_nc( Grid2d Rho, 
		  Grid2d Mx, Grid2d Mz,
		  Grid2d V, Grid2d Bx, Grid2d Bz)
{
  double U[2][2], c00,c01,c10,c11;
  double a,b,c, rhoud[2], vud[2];
  register i, n;

   for( i=0; i<nlmax; i++) {
     for( n=0; n<nrmax; n++) {
       if(opt.nspin) { /* non spin-polarized */
	 a = .5 * Rho[i][n];
	 b = .5 * Rho[i][n];
	 c = 0.0;
       } else {		/* spin polarized */
	 a = .5 * ( Rho[i][n] + Mz[i][n]);
	 b = .5 * ( Rho[i][n] - Mz[i][n]);
	 c = .5 * Mx[i][n];
       }
       if( a != .0 || b != .0)
	 if( c*c/(a*a+b*b) < 1e-12) c = .0;

/* Diagonalize the density matrix: rho_up=rhoud[0], U unitary */
       diagsm2d( a,b,c, rhoud, U);

/* Apply the exchange operator: */
       vud[0] = vxc( rhoud[0], rhoud[1],  1) *2.0; /* Ryd! */
       vud[1] = vxc( rhoud[0], rhoud[1], -1) *2.0;

/* The new non collinear potential is U*vud*U` */
       c00 = U[0][0];
       c01 = U[0][1];
       c10 = U[1][0];
       c11 = U[1][1];

       a = c00*c00*vud[0] + c01*c01*vud[1];
       b = c10*c10*vud[0] + c11*c11*vud[1];
       V[i][n]  = .5 * (a+b);
       Bz[i][n] = .5 * (a-b);

       if( c!=0.0)
	 Bx[i][n] = vud[0] * c10 * c00 + vud[1] * c11 * c01;
       else
	 Bx[i][n] = 0.0;
     }
   }
}
	      
