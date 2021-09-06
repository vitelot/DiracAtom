#include "relat.h"

static void potcul_nsph_l( Grid2d Rho,
			   Grid2d Ex, Grid2d Ez,
			   double *cos_theta, double *weights,
			   double *r, double *drdi, int maxlexp);

void cleanb( Grid2d Bx, Grid2d Bz,
	     double *cos_theta, double *weights,
	     Grid2d sources,
	     int maxlexp, double *ra, double *drdi)
{
Grid2d r2Br, Ex, Ez;

double **sBt; /* sBt[NRMAX][NLMAX] reverted indexes */
register int nr, nl;

double r, dc, sBt_dc, r2Br_dr, dr, r_dr, c, c_dc, s, *sin_theta;

#include "alloc_cleanb.h"

  for( nl=0; nl<nlmax; nl++) 
  {
    c = cos_theta[nl];
    s = sin_theta[nl] = sqrt( 1-cos_theta[nl]*cos_theta[nl]);
    for( nr=0; nr<nrmax; nr++)
    {
      r2Br[nl][nr] = ra[nr]*ra[nr] * (s*Bx[nl][nr]+c*Bz[nl][nr]);
      sBt[nr][nl] = s*( c*Bx[nl][nr]-s*Bz[nl][nr]);
    }
  }

  dr = (ra[1]-ra[0])/10.;
  dc = 0.0001;
  for( nl=0; nl<nlmax; nl++) 
  {
    c = cos_theta[nl];
    s = sin_theta[nl];
    if( c<0) c_dc = c+dc; else c_dc = c-dc;
    for( nr=0; nr<nrmax-1; nr++)
    {      
      r = ra[nr];
      r_dr = r+dr;

      r2Br_dr = ylag( r_dr, ra, r2Br[nl], nr+1, 6, nrmax-1);
      sBt_dc = ylag( c_dc, cos_theta, sBt[nr], (c<0)?(nl+1):nl, 4, nlmax-1);

      sources[nl][nr] =
	(r2Br_dr-r2Br[nl][nr])/(r*r*dr) -
	((c<0)?(1):(-1))*(sBt_dc-sBt[nr][nl])/(r*dc);

    }   /* End for nl */
  }     /* End for nr */

  potcul_nsph_l( sources, Ex, Ez, cos_theta, weights, ra, drdi, maxlexp);

  for( nl=0; nl<nlmax; nl++) 
    for( nr=0; nr<nrmax-1; nr++)
    {
      Bx[nl][nr] -= Ex[nl][nr];
      Bz[nl][nr] -= Ez[nl][nr];
    }

#include "free_cleanb.h"
}

static void potcul_nsph_l( Grid2d Rho,
			   Grid2d Ex, Grid2d Ez,
			   double *cos_theta, double *weights,
			   double *r, double *drdi, int maxlexp)
{
  register int i,l,n;
  int il;
  double temp, **rhotr, **vctr, **vctrp, **lp, **lp1, rn,
         costheta, sintheta;
/*   char filename[128]; */
/*   FILE *fp; */

  rhotr = matrix_alloc( 1+maxlexp, nrmax);
  vctr  = matrix_alloc( 1+maxlexp, nrmax);
  vctrp = matrix_alloc( 1+maxlexp, nrmax);
  lp    = matrix_alloc( 1+maxlexp, nlmax);
  lp1   = matrix_alloc( 1+maxlexp, nlmax);

  /* find legendre components of Rho storing them in rhotr */
  for( l=1; l<=maxlexp; l+=2) /* only odd components are in Rho */
  {
    il = l;
    for( i=0; i<nlmax; i++)
    {
       lp[il][i] = LegPol( l, 0, cos_theta[i]);
      lp1[il][i] = LegPol( l, 1, cos_theta[i]);
    }

    for( n=0; n<nrmax; n++)
    {
      for( temp=i=0; i<nlmax; i++)
	temp += weights[i] * lp[il][i] * Rho[i][n];
      rhotr[il][n] = (l+.5)*temp;
    }
  }

  for( i=0; i<nlmax; i++)
    for( n=0; n<nrmax; n++)
      Ex[i][n] = Ez[i][n] = .0;
  for( l=1; l<=maxlexp; l+=2) /* only odd components are in Vc */
  {
    il = l;
    nsphint( vctr[il], vctrp[il], rhotr[il], r, drdi, l);
    for( i=0; i<nlmax; i++)
    {
      costheta = cos_theta[i];
      sintheta = sqrt(1-costheta*costheta);
      for( n=0; n<nrmax; n++)
      {
	rn = r[n];
	Ex[i][n] -= sintheta*(vctrp[il][n] * lp[il][i]+
			      vctr[il][n] * lp1[il][i]*costheta/(rn*sintheta));
	Ez[i][n] -= costheta * vctrp[il][n] * lp[il][i] - 
	            vctr[il][n] * lp1[il][i]*sintheta/rn;
	
      }
      
    }
  }

/*   for( l=1; l<=maxlexp; l+=2) { */
/*     sprintf( filename, "SOURCESRHOTR%02d", l);  */
/*     fp = fopen( filename, "w");  */
/*     for( n=0; n<nrmax; n++) { */
/*       fprintf( fp, "%lg %lg %lg %lg\n", */
/* 	       r[n], M_4PI/(2.*l+1.)*rhotr[l][n], vctr[l][n], vctrp[l][n]);  */
/*     }  */
/*     fclose(fp);  */
/*   }   */

  matrix_free( rhotr,  1+maxlexp, nrmax);
  matrix_free(  vctr,  1+maxlexp, nrmax);
  matrix_free( vctrp,  1+maxlexp, nrmax);
  matrix_free(    lp,  1+maxlexp, nlmax);
  matrix_free(   lp1,  1+maxlexp, nlmax);

}

