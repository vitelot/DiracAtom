#include "relat.h"

#define NPEMAX 20        /* Max. degree of expansion                  */
#define ITMAX  200       /* Max. iterations for Corrector to converge */
#define TOL    (1.0e-9)  /* Convergence tolerance for Corrector       */

#define cgo (*_cgo)  /* cgo is output, so pointer must be passed */

/* Following ones are global vars used in this subroutine */
#define C	(opt.c)
#define Nucl	(opt.nucl) 

void coredir( double z, double En, int l, double mj, char way,
	      Dbp vv, double *bb, double *rc, 
	      double *drdic, double *dovrc, int nmatch, int nzero,
	      double *gc[2][2], double *fc[2][2],
	      double *dp[2][2], double *dq[2][2],
	      double *wp[2][2], double *wq[2][2], 
	      double Pow[2][2], double Qow[2][2], 
	      double Piw[2][2], double Qiw[2][2],
	      double cgd[2], double cgmd[2], double *_cgo, double enstep)

/**********************************************************************
 INPUT:
  c        Light speed
  En       Extimated Energy-eigenvalue
  l        Large component angular momentum
  mj       Total angular momentum z projection
  way      'i' (inward) || 'o' (outward) integration
  *vv      Potential energy in each grid point
  *bb      Magnetic field in each grid point
  *rc      Grid point distance from origin
  *drdic   Deriv. of r with respect to index. As di=1, this is dr.
  *dovrc   drdic/rc
  nmatch   Meeting point between inward and outward solutions
  nzero    Last point from which we start inward integration
  nrc      Dimension of rc, bb, vv vectors
 
 OUTPUT:
  ***gc    G solution (large part) in each point
  ***fc    idem for F (small part)
  ***dp    P Differential in each point (derivative*delta(r))
  ***dq    idem for Q
  ***wp    r*G
  ***wq    r*c*F
  **Pow    Value of wp outward solution in nmatch
  **Qow    idem for wq 
  **Piw    Value of wp inward solution in nmatch
  **Qiw    idem for wq
  *cgd     G-coefficients of B field: <large|sz|large>
  *cgmd    G-coefficients of B field: <small|sz|small>
  cgo      G-coefficient of B field: <large(k1)|sz|large(k2)>. THE COUPLER!
 *********************************************************************/

{
  int i, j, jcorr, m, n;
  int kap1, kap2, nsol;

  double emvqq, emvpp, bqq, bpp, rr, diff;
  double csqr;
  double cg1, cg2, cg4, cg5, cg8;
  double bc[NPEMAX+1], vc[NPEMAX+1], tz;
  double gam[2], kap[2];
  double pnew[2][2], qnew[2][2], qold[2][2], pold[2][2];

  static double pc[2][2][NPEMAX+1], qc[2][2][NPEMAX+1];

# define h24 0.04166666666666666667

      csqr = C*C;

/* Expansion coefficients for the potential and B-field */

      if( Nucl) {
	tz = 2.*z;
	vc[0] = vv[0];
      } else {
	tz = v_nint(-vv[0]*rc[0]);
	vc[0] = vv[0] + tz/rc[0];
      }
      bc[0] = bb[0];

/* Calculate G-coefficients of B-field.
   See "Cortona. Phys. Rew. A 31-5, 2842, (1985)" */

      kap1 = -l-1;
      kap2 = l;

      cg1 = -mj/(kap1+.5);
      cg5 = -mj/(.5-kap1);
      cgd[0] = cg1;
      cgmd[0] = cg5;
      kap[0] = (double) kap1;

      if( Nucl) gam[0] = fabs(kap[0]);
      else gam[0] = sqrt(kap[0]*kap[0] - tz*tz/csqr);
  
      if (fabs(mj) > l)
      { 
	nsol = 1;
	cg2 = cg4 = cg8 = cgd[1] = cgo = cgmd[1] = gam[1] = kap[1] = 0.;
      }
      else
      {
	cg2 = -sqrt(1. - cg1*cg1);
	cg4 = - mj/(kap2+.5);
	cg8 = - mj/(.5-kap2);
	nsol = 2;
	cgd[1] = cg4;
	cgo = cg2;
	cgmd[1] = cg8;
	kap[1] = (double) kap2;
	if( Nucl) gam[1] = fabs(kap[1]);
	else gam[1] = sqrt(kap[1]*kap[1] - tz*tz/csqr);
      }

      if (way == 'i' || way == 'I') goto inward;

/**********************************************************************
 **********************************************************************
                     OUTWARD INTEGRATION
 *********************************************************************/
 {
   const int mps=20;
   double aa12, aa21, aa11, aa22, det, bb1, bb2, rpwgpm, gpm;
   double alpha, beta, w1,w3,w4;
   
/* Determine mps higher expansion coefficients for the wave funct. */

      aa12 = -tz/csqr;
      aa21 = tz;
      emvqq = ( En-vc[0]+csqr)/csqr;
      emvpp = -En +vc[0];
      bqq = bc[0]/csqr;

      if( !Nucl) {
	for (j=0; j<nsol; j++)
	  {
	    i = 1-j;
  	    pc[j][j][0] = sqrt(fabs(kap[j])-gam[j]);
	    qc[j][j][0] = (kap[j]+gam[j])*csqr*pc[j][j][0]/tz;
	    pc[i][j][0] = qc[i][j][0] = 0.;
	  }

	for (j=0; j<nsol; j++)
	  for (m=1; m<=mps; m++)
	    for (i=0; i<nsol; i++)
	      {
		bb1 = (emvqq + bqq*cgmd[i])*qc[i][j][m-1];
		bb2 = (emvpp + bc[0]*cgd[i])*pc[i][j][m-1] + 
	          bc[0]* cgo *pc[1-i][j][m-1];
		aa11 = gam[i]+m+kap[i];
		aa22 = gam[i]+m-kap[i];
		det = aa11*aa22-aa12*aa21;
		pc[i][j][m] = (bb1*aa22-aa12*bb2)/det;
		qc[i][j][m] = (bb2*aa11-aa21*bb1)/det;
	      }
      } else { /* finite size nucleus */
	for (j=0; j<nsol; j++) {
	  i = 1-j;
	  /* arbitrary starting values */
	  if( kap[j]>0) {
	    alpha = 0.0;
	    beta = 1.0;
	  } else {
	    alpha = 1.0;
	    beta = 0.0;
	  }
	  pc[j][j][0] = alpha;
	  pc[i][j][0] = 0.0;
	  qc[j][j][0] = beta;
	  qc[i][j][0] = 0.0;
	}
	w4 = bc[0]*cgo;
	for (m=1; m<=mps; m++)
	  for (j=0; j<nsol; j++)
	    for (i=0; i<nsol; i++) {
	      w1 = emvqq + bqq*cgmd[i];
	      w3 = emvpp + bc[0]*cgd[i];
	      aa11 = gam[j]+kap[i]+m;
	      aa22 = gam[j]-kap[i]+m;
	      /* this happens when gam[1]=k[1]=l and k[0]=-l-1 */
	      if(aa11==0.0) pc[i][j][m] = 0.0;
	      else	    pc[i][j][m] = w1*qc[i][j][m-1]/aa11;
	      qc[i][j][m] = (w3*pc[i][j][m-1] + w4*pc[1-i][j][m-1])/aa22;
	  }
      }

/* Calculates wave functions for the first 4 points 
   using expansion coefficients previously calculated */

      for( n=0; n<4; n++)
      {
	rr = rc[n];
	for( j=0; j<nsol; j++)
	{
	  rpwgpm = pow(rr, gam[j]);
	  for( i=0; i<nsol; i++)
	  {
	    wp[i][j][n] = pc[i][j][0] * rpwgpm; 
	    wq[i][j][n] = qc[i][j][0] * rpwgpm;
	    dp[i][j][n] = wp[i][j][n] * gam[j] * dovrc[n]; 
	    dq[i][j][n] = wq[i][j][n] * gam[j] * dovrc[n]; 
	  }
	  for( m=1; m<=mps; m++)
	  {
	    rpwgpm *= rr;
	    gpm = gam[j] + m;
	    for( i=0; i<nsol; i++)
	    {
	      wp[i][j][n] += pc[i][j][m]*rpwgpm;
	      wq[i][j][n] += qc[i][j][m]*rpwgpm;
	      dp[i][j][n] += pc[i][j][m]*rpwgpm*gpm*dovrc[n];
	      dq[i][j][n] += qc[i][j][m]*rpwgpm*gpm*dovrc[n];
	    }
	  }
	}
      }
	  
/* Calculate all next points by predictor-corrector method:
   (Adams-Moulton-Bashforth). See "Num. Recipes" pg. 571    */

      for( n=4; n<=nmatch; n++)
      {

/* Evaluate predictor */

	for( j=0; j<nsol; j++)
	  for( i=0; i<nsol; i++)
	  {
	    pnew[i][j] = wp[i][j][n-1] + 
	                 h24*(55.0*dp[i][j][n-1]-59.0*dp[i][j][n-2]+
			      37.0*dp[i][j][n-3]- 9.0*dp[i][j][n-4]);
	    qnew[i][j] = wq[i][j][n-1] + 
	                 h24*(55.0*dq[i][j][n-1]-59.0*dq[i][j][n-2]+
			      37.0*dq[i][j][n-3]- 9.0*dq[i][j][n-4]);
	  }
	emvqq =  (En - vv[n] + csqr)*drdic[n]/csqr;
	emvpp = -(En - vv[n]       )*drdic[n];
	bpp = bb[n]*drdic[n];
	bqq = bpp/csqr;

/* Evaluate corrector */

	for(jcorr=0; jcorr<ITMAX; jcorr++)
	{
	  for( j=0; j<nsol; j++)
	    for( i=0; i<nsol; i++)
	    {
	      pold[i][j] = pnew[i][j];
	      qold[i][j] = qnew[i][j];
	      dp[i][j][n] = -kap[i]*pnew[i][j]*dovrc[n]+
		            (emvqq+bqq*cgmd[i])*qnew[i][j];
	      dq[i][j][n] = kap[i]*qnew[i][j]*dovrc[n]+
		            (emvpp+bpp*cgd[i])*pnew[i][j]+
			    bpp* cgo *pnew[1-i][j];

	      pnew[i][j] = wp[i][j][n-1]+
		           h24*(9.0*dp[i][j][n]  +19.0*dp[i][j][n-1]-
				5.0*dp[i][j][n-2]+     dp[i][j][n-3]);
	      qnew[i][j] = wq[i][j][n-1]+
		           h24*(9.0*dq[i][j][n]  +19.0*dq[i][j][n-1]-
				5.0*dq[i][j][n-2]+     dq[i][j][n-3]);
	    }

	  for( j=0; j<nsol; j++)
	    for( i=0; i<nsol; i++)
	    {
	      diff = fabs(pold[i][j]-pnew[i][j]);
	      if( diff > (TOL*fabs(pnew[i][j])) ) goto not_yet_converged;
	      diff = fabs(qold[i][j]-qnew[i][j]);
	      if( diff > (TOL*fabs(qnew[i][j])) ) goto not_yet_converged;
	    }

	  break;

	not_yet_converged:
	  ;
	}

	if (jcorr == ITMAX)
	{
	  printf("\nOUTWARD CORRECTOR DID NOT CONVERGE in %d ITERATIONS (E=%lg L=%d)\n",
		 ITMAX, En, l);
	  if(!opt.wall) exit(1);
	}    
	      
	for( j=0; j<nsol; j++)
	  for( i=0; i<nsol; i++)
	  {
	    wp[i][j][n] = pnew[i][j];
	    wq[i][j][n] = qnew[i][j];
	    dp[i][j][n] = -kap[i]*pnew[i][j]*dovrc[n]+
		            (emvqq+bqq*cgmd[i])*qnew[i][j];
	    dq[i][j][n] = kap[i]*qnew[i][j]*dovrc[n]+
		            (emvpp+bpp*cgd[i])*pnew[i][j]+
			    bpp* cgo *pnew[1-i][j];
	  }
      }

/* Let's transform to the proper wave functions */

      for(n=0; n<=nmatch; n++)
	for(j=0; j<nsol; j++)
	  for(i=0; i<nsol; i++)
	  {
	    gc[i][j][n] = wp[i][j][n]/rc[n];
	    fc[i][j][n] = wq[i][j][n]/(rc[n]*C);
	  }

      for(j=0; j<nsol; j++)
	for(i=0; i<nsol; i++)
	{
	  Pow[i][j] = wp[i][j][nmatch];
	  Qow[i][j] = wq[i][j][nmatch];
	}

  return;
 }

inward:

/**********************************************************************
 **********************************************************************
                     INWARD INTEGRATION
 *********************************************************************/

 {
   double dmue, bova, enold;
   
      enold = En;   /* 220598 */
      En -= enstep;      /* 220598 */    
      dmue = sqrt(-En-En*En/csqr);
      bova = -dmue/(1.0+En/csqr);
      for( n=nzero-4; n<nzero; n++)
      {
	rr = rc[n];

	for( j=0; j<nsol; j++)
	{
	  i = 1-j;

	  wp[j][j][n] = exp(-dmue*rr);
	  dp[j][j][n] = -dmue*drdic[n]*wp[j][j][n];
	  wq[j][j][n] = bova*wp[j][j][n];
	  dq[j][j][n] = bova*dp[j][j][n];
 	  wp[i][j][n] = 0.0; 
 	  wq[i][j][n] = 0.0;
 	  dp[i][j][n] = 0.0; 
 	  dq[i][j][n] = 0.0;
	}
      }

      En = enold;   /* 220598 */

/* Calculate all next points by predictor-corrector method:
   (Adams-Moulton-Bashforth). See "Num. Recipes" pg. 571    */

      for( n=nzero-5; n>=nmatch; n--) {

/* Evaluate predictor */

	for( j=0; j<nsol; j++)
	  for( i=0; i<nsol; i++)
	  {
	    pnew[i][j] = wp[i][j][n+1] - 
	                 h24*(55.0*dp[i][j][n+1]-59.0*dp[i][j][n+2]+
			      37.0*dp[i][j][n+3]- 9.0*dp[i][j][n+4]);
	    qnew[i][j] = wq[i][j][n+1] - 
	                 h24*(55.0*dq[i][j][n+1]-59.0*dq[i][j][n+2]+
			      37.0*dq[i][j][n+3]- 9.0*dq[i][j][n+4]);
	  }

	emvqq =  (En - vv[n] + csqr)*drdic[n]/csqr;
	emvpp = -(En - vv[n]       )*drdic[n];
	bpp = bb[n]*drdic[n];
	bqq = bpp/csqr;

/* Evaluate corrector */

	for(jcorr=0; jcorr<ITMAX; jcorr++) {
	  for( j=0; j<nsol; j++)
	    for( i=0; i<nsol; i++) {
	      pold[i][j] = pnew[i][j];
	      qold[i][j] = qnew[i][j];
	      dp[i][j][n] = -kap[i]*pnew[i][j]*dovrc[n]+
		            (emvqq+bqq*cgmd[i])*qnew[i][j];
	      dq[i][j][n] = kap[i]*qnew[i][j]*dovrc[n]+
		            (emvpp+bpp*cgd[i])*pnew[i][j]+
			    bpp* cgo *pnew[1-i][j];

	      pnew[i][j] = wp[i][j][n+1]-
		           h24*(9.0*dp[i][j][n]  +19.0*dp[i][j][n+1]-
				5.0*dp[i][j][n+2]+     dp[i][j][n+3]);
	      qnew[i][j] = wq[i][j][n+1]-
		           h24*(9.0*dq[i][j][n]  +19.0*dq[i][j][n+1]-
				5.0*dq[i][j][n+2]+     dq[i][j][n+3]);
	    }

	  for( j=0; j<nsol; j++)
	    for( i=0; i<nsol; i++) {
	      diff = fabs(pold[i][j]-pnew[i][j]);
	      if( diff > (TOL*fabs(pnew[i][j])) ) goto not_yet_conv_inw;
	      diff = fabs(qold[i][j]-qnew[i][j]);
	      if( diff > (TOL*fabs(qnew[i][j])) ) goto not_yet_conv_inw;
	    }

	  break;

	not_yet_conv_inw:
	  ;
	}

	if (jcorr == ITMAX) {
	  printf("\nINWARD CORRECTOR DID NOT CONVERGE in %d ITERATIONS\n",
		 ITMAX);
	  exit(1);
	}    
	      
	for( j=0; j<nsol; j++)
	  for( i=0; i<nsol; i++) {
	    wp[i][j][n] = pnew[i][j];
	    wq[i][j][n] = qnew[i][j];
	    dp[i][j][n] = -kap[i]*pnew[i][j]*dovrc[n]+
		            (emvqq+bqq*cgmd[i])*qnew[i][j];
	    dq[i][j][n] = kap[i]*qnew[i][j]*dovrc[n]+
		            (emvpp+bpp*cgd[i])*pnew[i][j]+
			    bpp* cgo *pnew[1-i][j];
	  }
      }

/* Let's transform to the proper wave functions */

      for(n=nmatch; n<nzero; n++)
	for(j=0; j<nsol; j++)
	  for(i=0; i<nsol; i++)
	  {
	    gc[i][j][n] = wp[i][j][n]/rc[n];
	    fc[i][j][n] = wq[i][j][n]/(rc[n]*C);
	  }

      for(j=0; j<nsol; j++)
	for(i=0; i<nsol; i++)
	{
	  Piw[i][j] = wp[i][j][nmatch];
	  Qiw[i][j] = wq[i][j][nmatch];
	}

    return;
 }	
	    
	
}

#undef NPEMAX 
#undef ITMAX
#undef TOL
#undef cgo
#undef h24
#undef C
#undef Nucl
