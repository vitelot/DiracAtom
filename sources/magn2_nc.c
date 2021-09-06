#include "relat.h"

#define  A(l,mj)  sqrt( (.5+(l)+(mj))/(2.0*(l)+1.0) )
#define  B(l,mj)  sqrt( (.5+(l)-(mj))/(2.0*(l)+1.0) )

void nii( double *psi[4], double *Rho,
		 double *Mx, double *Mz);

void nij( double *psii[4], double *psij[4],
		 double *rho, double *mx, 
		 double *mz);

void getpsi( double *gck[2], double *fck[2],
		    int ic, double,
		    double *psi[4]);

int idxmax(double vr[TOTST][TOTST], int col, int dim);

void pure_state(double v[TOTST], int ic);
void init_density( double **Rho, double **Mx, double **Mz);
void total_density( 
	     double **Rho, double **Mx, double **Mz,
	     double **Rhoic, double **Mxic, double **Mzic);

void orbital_density( 
	     double **Rho,
	     double **Mx,
	     double **Mz,
	     double *GCK[TOTST][2], double *FCK[TOTST][2],
	     double v[TOTST],
	     double *cos_theta);

void magn2_nc( int nel,
	       int onp[TOTST], int oncalc[TOTST],
	       double vr[TOTST][TOTST], 
	       double *cos_theta, 
	       double *GCK[TOTST][2], double *FCK[TOTST][2], 
	       double **Rho, double **Mx, double **Mz )
{
register int i, m, ic;
double v[TOTST];
double **Rhoic, **Mxic, **Mzic; /* density matrix of orbital ic */

  Mxic  = matrix_alloc(nlmax, nrmax);
  Mzic  = matrix_alloc(nlmax, nrmax);
  Rhoic = matrix_alloc(nlmax, nrmax);

  init_density( Rho,Mx,Mz );

  /* occupied states not in the basis */
/*   for(ic=0;ic<TOTST;ic++) */
/*     if(onp[ic] && !oncalc[ic]) { */

/*       pure_state( v, ic); */
/*       orbital_density( Rhoic,Mxic,Mzic, GCK,FCK,v, cos_theta); */
/*       total_density( Rho,Mx,Mz, Rhoic,Mxic,Mzic); */
/*       nel--; */
/*     } */

  for( m=0; m<nel; m++) {
    ic = idxmax( vr, m, TOTST);
    /* skip unoccupied levels */
    if( (fl.fixoccnum && !onp[ic]) )
      /* if the state corresponding to the column m is not occupied */
      { nel++; continue; }
    
    for( i=0; i<TOTST; i++)
      v[i] = vr[i][m];
    
    orbital_density( Rhoic,Mxic,Mzic, GCK,FCK,v, cos_theta);
    total_density( Rho,Mx,Mz, Rhoic,Mxic,Mzic);

  }

  matrix_free( Rhoic, nlmax, nrmax);
  matrix_free(  Mxic, nlmax, nrmax);
  matrix_free(  Mzic, nlmax, nrmax);
}

void nii( double *psi[4],
	  double *Rho, double *Mx, double *Mz)
{
  register int n;
  double c0,c1,c2,c3;

  for( n=0; n<nrmax; n++) {
    c0 = psi[0][n];
    c1 = psi[1][n];
    c2 = psi[2][n];
    c3 = psi[3][n];
		
    Rho[n] = c0*c0 + c1*c1 + c2*c2 + c3*c3;
      
    Mx[n] = 2.* (c1*c0 - c3*c2);
    Mz[n] = c0*c0 - c1*c1 - c2*c2 + c3*c3;
  }
  
}

void nij( double *psii[4], double *psij[4],
	  double *rho, double *mx, double *mz)
{
  register int i, n;
  double psis[4], psj[4], ct0, ct1, ct2, ct3;

  for( n=0; n<nrmax; n++) {
    for( i=0; i<4; i++) {
      psis[i] = psii[i][n];
      psj[i]  = psij[i][n];
    }

    ct0 = ( psis[0] * psj[0]);
    ct1 = ( psis[1] * psj[1]);
    ct2 = ( psis[2] * psj[2]);
    ct3 = ( psis[3] * psj[3]);

    rho[n] = ct0 + ct1 + ct2 + ct3;
    mz[n]  = ct0 - ct1 - ct2 + ct3;
      
    ct0 = ( psis[0] * psj[1]);
    ct1 = ( psis[1] * psj[0]);
    ct2 = ( psis[2] * psj[3]);
    ct3 = ( psis[3] * psj[2]);
      
    mx[n]  =  ct0 + ct1 - ct2 - ct3;
  }
}

void getpsi( double *gck[2], double *fck[2],
	     int ic, double cos_theta, 
	     double *psi[4])
{
  register int n;
  int nsol, l, imj;
  double mj;
  double c1, c2, c3, c4, c5, c6;

  l = st[ic].l;
  imj = st[ic].imj;
  if( imj == -l || imj == l+1)  nsol = 1;
  else				nsol = 2;
  mj = imj-0.5;

  if( nsol==1) {
    if( imj>0) {
      c1 = Y_nc( l  , l  , cos_theta);
      c2 = Y_nc( l+1, l  , cos_theta);
      c3 = Y_nc( l+1, l+1, cos_theta);

      for( n=0; n<nrmax; n++) {
	psi[0][n] =  gck[0][n]		   * c1;
	psi[1][n] =  0.0;  
	psi[2][n] = -fck[0][n]*B( l+1, mj) * c2;
	psi[3][n] =  fck[0][n]*A( l+1, mj) * c3;
	}
    } else {
      c1 = Y_nc( l  , -l  , cos_theta);
      c2 = Y_nc( l+1, -l-1, cos_theta);
      c3 = Y_nc( l+1, -l  , cos_theta);
	
      for( n=0; n<nrmax; n++) {
	psi[0][n] =  0.0;
	psi[1][n] =  gck[0][n]             * c1;
	psi[2][n] = -fck[0][n]*B( l+1, mj) * c2;
	psi[3][n] =  fck[0][n]*A( l+1, mj) * c3;
      }
    }
  } else {
    c1 = Y_nc( l  , imj-1, cos_theta);
    c2 = Y_nc( l  ,   imj, cos_theta);
    if( imj == 1-l)
      c3 = 0.0;
    else
      c3 = Y_nc( l-1, imj-1, cos_theta);
    c4 = Y_nc( l+1,   imj-1, cos_theta);
    if( imj == l)
      c5 = 0.0;
    else
      c5 = Y_nc( l-1, imj, cos_theta);
    c6 = Y_nc( l+1,   imj, cos_theta);
      
    for( n=0; n<nrmax; n++) {
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
#undef A
#undef B

int idxmax(double vr[TOTST][TOTST], int col, int dim)
/* return the row where lies the maximum value in the column col of matr vr */
{
  register int i;
  int imax;
  double tmp, max=-1.0;
  
  for(i=0;i<dim;i++) {
    tmp = fabs(vr[i][col]);
    if( tmp>max) {
      max = tmp;
      imax = i;
    }
  }
  if(max<1e-6) v_error("ERROR: max is zero in IDXMAX\n");
  return imax;
}

void orbital_density( 
	     double **Rho,
	     double **Mx,
	     double **Mz,
	     double *GCK[TOTST][2], double *FCK[TOTST][2],
	     double v[TOTST],
	     double *cos_theta)
{
register int i,j, n,l;
int psiigot;
double  cim2, ctemp, cims, cjm, cimscjm;
Dbp  rhoii, mxii, mzii, rhoij, mxij, mzij, psii[4], psij[4];
  
  /* allocations */
  rhoii = vector_alloc(nrmax);
  mxii = vector_alloc(nrmax);
  mzii = vector_alloc(nrmax);
  rhoij = vector_alloc(nrmax);
  mxij = vector_alloc(nrmax);
  mzij = vector_alloc(nrmax);

  for(i=0;i<4;i++) {
    psii[i] = vector_alloc(nrmax);
    psij[i] = vector_alloc(nrmax);
  }
  /***************/

  init_density( Rho,Mx,Mz );

  for( i=0; i<TOTST; i++) {
    psiigot = FALSE;
    cims = v[i];
    cim2 =  cims*cims;
    for( l=0; l<nlmax; l++) {
    /* There are a lot of zeros in the coefficients */
      if( cim2 > 1e-16) { 
	getpsi( GCK[i], FCK[i], i, cos_theta[l], psii);
	psiigot = TRUE;
	nii( psii, rhoii, mxii, mzii);
	for( n=0; n<nrmax; n++) {
	  Rho[l][n] += cim2 * rhoii[n];
	   Mz[l][n] += cim2 *  mzii[n];
	}
	if( opt.nocoll)
	    for( n=0; n<nrmax; n++)
	      Mx[l][n] += cim2 * mxii[n];
      }
      for( j=0; j<i; j++) {
	cjm  = v[j];
	cimscjm = cims * cjm;
      /* There are a lot of zeros in the coefficients */
	if( fabs(cimscjm) > 1e-16) {
	  if( !psiigot) {
	    getpsi( GCK[i], FCK[i], i, cos_theta[l], psii);
	    psiigot = TRUE;
	  }
	
	  getpsi( GCK[j], FCK[j], j, cos_theta[l], psij);
	  nij( psii, psij, rhoij, mxij, mzij);
	  for( n=0; n<nrmax; n++) {
	    ctemp = cimscjm * rhoij[n];
	    Rho[l][n] += 2.*ctemp;
		
	    ctemp = cimscjm * mzij[n];
	    Mz[l][n] += 2.*ctemp;
	  }
	if( opt.nocoll)	      
	    for( n=0; n<nrmax; n++)
	      Mx[l][n] += 2.*cimscjm * mxij[n];
	}
      }
    }
  }

  /* deallocations */
  vector_free(rhoii);
  vector_free(mxii);
  vector_free(mzii);
  vector_free(rhoij);
  vector_free(mxij);
  vector_free(mzij);

  for(i=0;i<4;i++) {
    vector_free(psii[i]);
    vector_free(psij[i]);
  }
  /*****************/
}

void init_density( double **Rho, double **Mx, double **Mz)
{
  register int n,l;

  for( l=0; l<nlmax; l++)
    for( n=0; n<nrmax; n++)
      Rho[l][n] = Mx[l][n] = Mz[l][n] = 0.0;
}

void total_density( 
	     double **Rho, double **Mx, double **Mz,
	     double **Rhoic, double **Mxic, double **Mzic)
{
  register int n,l;

  for( l=0; l<nlmax; l++)
    for( n=0; n<nrmax; n++) {
      Rho[l][n] += Rhoic[l][n];
       Mx[l][n] +=  Mxic[l][n];
       Mz[l][n] +=  Mzic[l][n];
    }
}

void pure_state(double v[TOTST], int ic)
{
  register int i;
  for( i=0;i<TOTST;i++) {
    v[i] = 0.0;
    if(i==ic) v[i]=1.0;
  }
}

