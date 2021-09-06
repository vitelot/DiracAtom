#include "relat.h"

#define  A(l,mj)  sqrt( (.5+(l)+(mj))/(2.0*(l)+1.0) )
#define  B(l,mj)  sqrt( (.5+(l)-(mj))/(2.0*(l)+1.0) )

struct State { int occ;
		      int l;
		      double mj;
		      int nsol;  };

static void getpsi( double *gck[2], double *fck[2],
		    struct State *, double,
		    double *psi[4]);

void deltaw_nc( int oncalc[TOTST],
		double *cos_th, double *weights,
		double *r2drdi,
		double *vt, double *Vxc, double *bt,
		Grid2d Vc_nsph,  Grid2d Vxc_nsph,
		Grid2d Bx,  Grid2d Bz,
		double *GCK[TOTST][2], double *FCK[TOTST][2],
		double DW_r[TOTST][TOTST]
		)
/**********************************************************************
	 DELTAW: Calculates DWjk simmetric matrix elements 

 INPUT:
  *oncalc    Calculated orbital base: 0=not calculated, 1,2=nsol
  *cos_th    Cosinus table for Gauss-Leg integration
  *weights   Weights table for Gauss-Leg integration
  *r2drdi    (=rc*rc*drdi)
  *Vxc	     Collinear spherical exc-corr potential
  *vt	     Collinear spherical coulomb potential (nuclear included)
  *bt	     Collinear spherical exc field
  **Vc_nsph  Non spherical coulomb potential (nuclear included)
  **Vxc_nsph Non spherical XC potential      
  **Bx		,,
  **Bz		,,
  ***GCK     Large components wf for each eigenfunct. (Base)
  ***FCK     Small components wf for each eigenfunct. (Base)

 OUTPUT:
  **DW_r     Hermitian non collinear correction matrix. Real part.

 *********************************************************************/

{
register int ic, id, l, i,n, s, ilshell;
int nsol, nlshell;

double mj, temp1, temp2, dw11, dw22;

double c11, c22, c21, c21s, dw21, tinteg;
Dbp psij[4], psik[4], rinteg_r;

struct State state[TOTST];

#include "alloc_deltaw.h"

      nlshell = TOTSH;  
      ic = 0;   /* Counter of processed polarized states */
  
      for( ilshell=0; ilshell<nlshell; ilshell++) {
        l   = lqntab[ ord[ilshell] ];

        for( mj=-.5-l; mj<=.5+l; mj += 1.0) {
          nsol = 2;
          if( v_abs(mj) > l)
            nsol = 1;

	  for( s=0; s<nsol; s++) {
	    state[ic].occ = oncalc[ic];
	    state[ic].l = l;
	    state[ic].mj = mj;
	    state[ic].nsol = nsol;
	    
/* 	    printf(" [%03d] on:%d l:%d mj:%+3.1lf nsol:%d\n", */
/* 		   ic, oncalc[ic], l, mj, nsol);	      */

	    ic++;
	  }			/* End 'for s'		*/
	}			/* End 'for mj'		*/
      }				/* End 'for ilshell'	*/

      for( ic=0; ic<TOTST; ic++) {    /*  Rows   */
	if( state[ic].occ) {
	  for( id=0; id<=ic; id++) { /* Columns */
	    if( state[id].occ) {
	     if( (state[ic].mj == state[id].mj) 
		  && !(abs(state[id].l-state[ic].l)&01) ) { /* if Dl even.. */ 
/* Selection rules for mj and l */
/*	      printf("[%02d] % 12.8lf [%02d] % 12.8lf\n",
		     ic, state[ic].mj, id, state[id].mj);
*/
	      for( tinteg=i=0; i<nlmax; i++) {
		getpsi( GCK[ic], FCK[ic], &state[ic], cos_th[i], psij);
		getpsi( GCK[id], FCK[id], &state[id], cos_th[i], psik);
		for( n=0; n<nrmax; n++) {
		  if( opt.nocoll && opt.nosph) {
		    temp1 = Vxc_nsph[i][n]+Vc_nsph[i][n]-vt[n];
		    temp2 = Bz[i][n]-bt[n];
		  }
		  else if( opt.nocoll && !opt.nosph) {
		    temp1 = Vxc_nsph[i][n] - Vxc[n];
		    temp2 = Bz[i][n]-bt[n];
		  }
		  else if( !opt.nocoll && opt.nosph) {
		    temp1 = Vxc_nsph[i][n]+Vc_nsph[i][n]-vt[n];
		    temp2 = Bz[i][n]-bt[n];

/* The following is good if you are not interested in the XC part of the potential */
/* 		    temp1 = Vc_nsph[i][n] - vt[n] + Vxc[n]; */
/* 		    temp2 = 0.0; */
		  }
		  else temp1 = temp2 = 0.0;
		    
		  dw11 = temp1 + temp2; 
		  dw22 = temp1 - temp2;
		  dw21 = (opt.nocoll)? Bx[i][n]: 0.0;

		  c11 = psij[0][n] * psik[0][n] +
			psij[3][n] * psik[3][n];
		  c22 = psij[1][n] * psik[1][n] +
		        psij[2][n] * psik[2][n];
      
		  c21 = psij[1][n] * psik[0][n] -
			psij[3][n] * psik[2][n];
	          c21s= psij[0][n] * psik[1][n] -
			psij[2][n] * psik[3][n];
	  
 	          c11 = dw11 * c11;
 	          c22 = dw22 * c22;
 	          c21 = dw21 * c21;
 	          c21s= dw21 * c21s;
	
		  rinteg_r[n] = (c11 + c22 + c21 + c21s)*r2drdi[n];
		}
		tinteg += weights[i] * simpson( rinteg_r, nrmax, 1.);
	      }
	      DW_r[ic][id] = tinteg*2.0*M_PI; /* 2pi for phi integ */
	      if( ic != id) {
		DW_r[id][ic] = DW_r[ic][id];
	      }
	     }		/* End 'if m==m`'   */ 
	     else {
	       DW_r[id][ic] = DW_r[ic][id] = 0.0;
	     }
	    }		/* End 'if state_id */
	  }		/* End 'for id'     */
	}		/* End 'if state_ic */
      }			/* End 'for ic'     */

/*       perr("# VVVVV\n"); */
/*       for( n=0; n<nrmax; n++) { */
/* 	for( tinteg=i=0; i<nlmax; i++) { */
/* 	  tinteg += weights[i] * (Vxc_nsph[i][n]+Vc_nsph[i][n]-vt[n]); */
/* 	} */
/* 	perr("%d %lg VVVVV\n", n, tinteg); */
/*       } */

#include "free_deltaw.h"

}

static void getpsi( double *gck[2], double *fck[2],
		    struct State *state, double cos_theta, 
		    double *psi[4])
{
  register int n;
  int nsol, l;
  double mj;
  double c1, c2, c3, c4, c5, c6;

    nsol = state->nsol;
    l = state->l;
    mj = state->mj;

    if( nsol==1)
    {
      if( mj>0)
      {
	c1 = Y_nc( l  , l  , cos_theta);
	c2 = Y_nc( l+1, l  , cos_theta);
	c3 = Y_nc( l+1, l+1, cos_theta);
	
	for( n=0; n<nrmax; n++)
	{
	  psi[0][n] =  gck[0][n]	     * c1;
	  psi[1][n] =  0.0;  
	  psi[2][n] = -fck[0][n]*B( l+1, mj) * c2;
	  psi[3][n] =  fck[0][n]*A( l+1, mj) * c3;
	}
      }
      else
      {
	c1 = Y_nc( l  , -l  , cos_theta);
	c2 = Y_nc( l+1, -l-1, cos_theta);
	c3 = Y_nc( l+1, -l  , cos_theta);
	
	for( n=0; n<nrmax; n++)
	{
	  psi[0][n] =  0.0;
	  psi[1][n] =  gck[0][n]             * c1;
	  psi[2][n] = -fck[0][n]*B( l+1, mj) * c2;
	  psi[3][n] =  fck[0][n]*A( l+1, mj) * c3;
	}
      }
    }
    else
    {
      c1 = Y_nc( l  , (int)(mj-.5), cos_theta);
      c2 = Y_nc( l  , (int)(mj+.5), cos_theta);
      if( mj == .5-l)
	c3 = 0.0;
      else
	c3 = Y_nc( l-1, (int)(mj-.5), cos_theta);
      c4 = Y_nc( l+1, (int)(mj-.5), cos_theta);
      if( mj == l-.5)
	c5 = 0.0;
      else
	c5 = Y_nc( l-1, (int)(mj+.5), cos_theta);
      c6 = Y_nc( l+1, (int)(mj+.5), cos_theta);
      
      for( n=0; n<nrmax; n++)
      {
	psi[0][n] = (-gck[1][n]*B(l, mj)+
		      gck[0][n]*A(l, mj) ) * c1;
	psi[1][n] = ( gck[1][n]*A(l, mj)+
		      gck[0][n]*B(l, mj) ) * c2;
	psi[2][n] =   fck[1][n]*A(l-1, mj) * c3 
		     -fck[0][n]*B(l+1, mj) * c4;
	psi[3][n] =   fck[1][n]*B(l-1, mj) * c5 +
		      fck[0][n]*A(l+1, mj) * c6;
      }
    }
}
