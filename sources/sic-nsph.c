#include "relat.h"

/******************* EXTERN *********************************/
extern void diagrsg( 
	      double As[TOTST][TOTST], double Ns[TOTST][TOTST],
	      int dim,
	      double lambda[TOTST], double vr[TOTST][TOTST]);

extern void getpsi( double *gck[2], double *fck[2],
	       int ic, double cos_theta, 
	       double *psi[4]);
extern void nii( double *psi[4],
	    double *Rho, double *Mx, double *Mz);
extern void nij( double *psii[4], double *psij[4],
	    double *rho, double *mx, double *mz);
extern int idxmax(double vr[TOTST][TOTST], int col, int dim);

extern void orbital_density( 
	     double **Rho,
	     double **Mx,
	     double **Mz,
	     double *GCK[TOTST][2], double *FCK[TOTST][2],
	     double v[TOTST],
	     double *cos_theta);

extern void init_density( double **Rho, double **Mx, double **Mz);
extern void total_density( 
	     double **Rho, double **Mx, double **Mz,
	     double **Rhoic, double **Mxic, double **Mzic);
extern void pure_state(double v[TOTST], int ic);
/************************************************************/

/******************* STATIC *********************************/
static void select_v( double *L, double v[TOTST],
		      int ic, int idx[TOTST],
		      double lambda[TOTST], double vr[TOTST][TOTST], int dim);

static void deltawic( double DW[TOTST][TOTST],
		      double *enee_nsph, double *teevxc,
		      double **Rho, double **Mx, double **Mz,
		      int oncalc[TOTST],
		      double *GCK[TOTST][2], double *FCK[TOTST][2],
		      double *vsic, double *bsic, double *Vxcsic,
		      double *vnuc,
		      double *cos_theta, double *weights,
		      double *r, double *drdi, double *r2drdi);

static void gendiag_matrix( 
		    double A[TOTST][TOTST],
		    double N[TOTST][TOTST],
		    int ic,
		    double *L,
		    double v[TOTST],
		    int oncalc[TOTST]
		    );
static void join_matrix( 
		 double out[TOTST][TOTST],
		 double ectab[TOTST],
		 double ovlp[TOTST][TOTST],
		 double dvsic[TOTST][TOTST],
		 double DW[TOTST][TOTST],
		 double DWic[TOTST][TOTST],
		 int oncalc[TOTST]
		 );
static void delta_sicpot_matrix(
		    double dvsic[TOTST][TOTST],
		    int i,
		    int oncalc[TOTST],
		    double *r2drdi,
		    double *vsic[TOTST], double *bsic[TOTST],
		    double *GCK[TOTST][2], double *FCK[TOTST][2]
		    );
static void overlap_matrix(
		    double ovlp[TOTST][TOTST],
		    int oncalc[TOTST],
		    double *r2drdi,
		    double *GCK[TOTST][2], double *FCK[TOTST][2]
		    );
/************************************************************/

void sic_nsph(
	      /* input */
	      int onp[TOTST],
	      int oncalc[TOTST], double ectab[TOTST],
	      double *r, double *drdi, double *r2drdi,
	      double *vsic[TOTST], double *bsic[TOTST], double *Vxcsic[TOTST],
	      double *vnuc,
	      double *GCK[TOTST][2], double *FCK[TOTST][2],
	      double DW[TOTST][TOTST],
	      double *cos_theta, double *weights,
	      double **Rho, double **Mx, double **Mz,
	      /* output */
	      double vr[TOTST][TOTST],
	      double *tesic, double *teps )
{
  register int ic, jc;
  static int count=0; /* counts the number of execution of this sub */
  static double ovlp[TOTST][TOTST]; /* overlap matrix not unity if SIC is applied */
  double dvsic[TOTST][TOTST]; /* overlap potential matrix not symmetric */
  double Aikl[TOTST][TOTST]; /* final symmetric matrix to diag */
  double DWic[TOTST][TOTST]; /* non spherical self-interaction matrix for state ic */
  double v[TOTST][TOTST];	    /* current list of eigenvectors */
  static double vold[TOTST][TOTST]; /* old list of eigenvectors */
  double lambda[TOTST]; /* eigenvalues */
  double enee_nsphic; /* hartree self-energy for eigenvector ic */
  double eexcsic; /* XC energy of ic */
  double eevxcsic; /* Int of (Vxc_nsph+Bz)*Rho + Mx*Bx for eigenvector ic */
  double seevxcsic; /* sum of integrals of (Vxc_nsph+Bz)*Rhoic + Mxic*Bx */
  double seexcsic; /* sum of XC energies of ic */
  double seneesic; /* sum of Hartree energies of ic */
  double **Rhoic, **Mxic, **Mzic; /* density matrix of orbital ic */

  Mxic  = matrix_alloc(nlmax, nrmax);
  Mzic  = matrix_alloc(nlmax, nrmax);
  Rhoic = matrix_alloc(nlmax, nrmax);

  *teps = seevxcsic = seexcsic = seneesic = 0.0;

  if(!count) { /* it is true only at the first instance */
    overlap_matrix( ovlp, oncalc, r2drdi, GCK, FCK);
    for(ic=0;ic<TOTST;ic++)
      if(oncalc[ic])
	pure_state( vold[ic], ic);

    /* the following may be too memory expensive */
/* 	delta_sicpot_matrix( dvsic[ic], ic, oncalc, r2drdi, vsic, bsic, GCK, FCK); */
  }

  init_density( Rho,Mx,Mz );

  /* occupied states not in the basis */
  for(ic=0;ic<TOTST;ic++)
    if(onp[ic] && !oncalc[ic]) {

      pure_state( v[ic], ic);
      orbital_density( Rhoic,Mxic,Mzic, GCK,FCK,v[ic], cos_theta);
      total_density( Rho,Mx,Mz, Rhoic,Mxic,Mzic);

    } else if( !onp[ic])
      pure_state( v[ic], ic);


  for(ic=0;ic<TOTST;ic++) {
    if(onp[ic] && oncalc[ic]) {

      printf("-------------------------------------------"
	     "-------------------------------------------\n");
      printf("ORBITAL #%03d STATE:[ %s K=%+d Jz=%+d/2 ]\n",
	     ic, st[ic].name, st[ic].k, 2*st[ic].imj-1);

      delta_sicpot_matrix( dvsic, ic, oncalc, r2drdi, vsic, bsic, GCK, FCK);
      orbital_density( Rhoic,Mxic,Mzic, GCK,FCK,vold[ic], cos_theta);

      printf("Averages:\t<Rho>=%.8lf\t<Mx>=%.6lf\t<Mz>=%.6lf\n", 
	     integ2d(Rhoic,r2drdi,weights),
	     integ2d( Mxic,r2drdi,weights),
	     integ2d( Mzic,r2drdi,weights));

      deltawic( DWic,
		&enee_nsphic, &eevxcsic,
		Rhoic, Mxic, Mzic,
		oncalc,
		GCK, FCK,
		vsic[ic], bsic[ic], Vxcsic[ic],
		vnuc,
		cos_theta, weights,
		r, drdi, r2drdi);

      join_matrix( Aikl, ectab, ovlp, dvsic, DW, DWic, oncalc);
      gendiag_matrix( Aikl, ovlp, ic, &lambda[ic], v[ic], oncalc);

      orbital_density( Rhoic,Mxic,Mzic, GCK,FCK,v[ic], cos_theta);
      total_density( Rho,Mx,Mz, Rhoic,Mxic,Mzic);

      vector_cpy( vold[ic], v[ic], TOTST);

      /* sum of the eigenvalues */
      *teps += lambda[ic];
      /* sum of the XC self-energies */
      eexcsic = iexc_nc( Rhoic,Mxic,Mzic, r2drdi, weights);
      seexcsic += eexcsic;
      /* sum of hartree self-energies */
      seneesic += enee_nsphic;
      /* -sum of integrals of (Vxc_nsph+Bz)*Rhoic + Mxic*Bx */
      seevxcsic += eevxcsic;

      printf("\nOrbital SIC contributions:\n");
      printf("   %-12s\t %-12s\t %-12s\t %-12s\t %-18s\n",
	     "Lambda", "Exc", "Evxc", "Ee", "Tot (Exc+Evxc+Ee)");
      printf("   %-12.4lf\t %-12.4lf\t %-12.4lf\t %-12.4lf\t %-18.4lf\n", 
	     lambda[ic], eexcsic, eevxcsic, -enee_nsphic, eexcsic+eevxcsic-enee_nsphic);

      printf("-------------------------------------------"
	     "-------------------------------------------\n");
    }
  }
  count++;

  /* eigenvectors ordered per columns */
  for(ic=0;ic<TOTST;ic++)
    for(jc=0;jc<TOTST;jc++)
      vr[jc][ic] = v[ic][jc];

  *tesic = seexcsic - seneesic + seevxcsic;


  matrix_free( Rhoic, nlmax, nrmax);
  matrix_free(  Mxic, nlmax, nrmax);
  matrix_free(  Mzic, nlmax, nrmax);
}

static void deltawic( double DW[TOTST][TOTST],
		      double *enee_nsph, double *teevxc,
		      double **Rho, double **Mx, double **Mz,
		      int oncalc[TOTST],
		      double *GCK[TOTST][2], double *FCK[TOTST][2],
		      double *vsic, double *bsic, double *Vxcsic,
		      double *vnuc,
		      double *cos_theta, double *weights,
		      double *r, double *drdi, double *r2drdi)
{
  register int i, n;
  double **Vxc_nsph, **Vc_nsph, **Bx, **Bz;
  double *vt, *bt, *Vxcs, *zero;

  zero = vector_alloc(nrmax);
  vt   = vector_alloc(nrmax);
  bt   = vector_alloc(nrmax);
  Vxcs = vector_alloc(nrmax);
  Vxc_nsph = matrix_alloc(nlmax, nrmax);
  Vc_nsph  = matrix_alloc(nlmax, nrmax);
  Bx       = matrix_alloc(nlmax, nrmax);
  Bz       = matrix_alloc(nlmax, nrmax);

  /* PZ81 SIC, i.e. orbitals are fully polarized */
/*   for( n=0; n<nrmax; n++) { */
/*     vt[n]   = vsic[n]; */
/*     bt[n]   = bsic[n]; */
/*     Vxcs[n] = Vxcsic[n]; */
/*     zero[n] = 0.0; */
/*   } */
  for( n=0; n<nrmax; n++) {
    vt[n] = vsic[n]+bsic[n];
    Vxcs[n] = Vxcsic[n]+bsic[n];
    zero[n] = bt[n] = 0.0;
  }

  *enee_nsph = potcul_nsph( Rho, Vc_nsph, cos_theta, weights,
			    r, zero, drdi, opt.maxlexp);

  /* PZ81 SIC, i.e. orbitals are fully polarized */
  dentopot_nc( Rho,Mx,Rho, Vxc_nsph,Bx,Bz);
/*   dentopot_nc( Rho,Mx,Mz, Vxc_nsph,Bx,Bz); */
  for(i=0;i<nlmax;i++) {
    for( n=0; n<nrmax; n++) {
      Vxc_nsph[i][n] += Bz[i][n];
      Bz[i][n] = 0.0;
    }
  }

/*   LLoop(i) { */
/*     RLoop(n) */
/*       perr("%lg %lg XXXXX\n",r[n],Vc_nsph[i][n]-vsic[n]); */
/*     perr(" \n"); */
/*   }  */

  deltaw_nc( oncalc,
	     cos_theta, weights,
	     r2drdi,
	     vt, Vxcs, bt,
	     Vc_nsph,
	     Vxc_nsph, Bx, Bz,
	     GCK, FCK,
	     DW);

  /* PZ81 SIC, i.e. orbitals are fully polarized */
  *teevxc = ivexc_nc( Rho, Mx, Rho, Vxc_nsph, Bx, Bz, r2drdi, weights); 
					/* -Integ( n*vxc) */
/*   *teevxc = ivexc_nc( Rho, Mx, Mz, Vxc_nsph, Bx, Bz, r2drdi, weights);  */

  vector_free(zero);
  vector_free(vt);
  vector_free(bt);
  vector_free(Vxcs);
  matrix_free( Vxc_nsph, nlmax, nrmax);
  matrix_free( Vc_nsph,  nlmax, nrmax);
  matrix_free( Bx,       nlmax, nrmax);
  matrix_free( Bz,       nlmax, nrmax);

}

static void gendiag_matrix( 
		    double A[TOTST][TOTST],
		    double N[TOTST][TOTST],
		    int ic,
		    double *L,
		    double v[TOTST],
		    int oncalc[TOTST]
		    )
/* Generalized diagonalization: Av=LNv
   returns the eigenval and eigenvect (L,v) corresponding to the state ic
   (if the mix is too high it not possible to decide which is which).
   The matrix containing the eigenvectors needs not be unitary any more
*/
{
  register int i,j,dim;
  int idx[TOTST]; /* stores the index to the states */
  double As[TOTST][TOTST], Ns[TOTST][TOTST]; /* smaller matrices to diag */
  double lambda[TOTST], vr[TOTST][TOTST];

  /* build the index table */
  for(dim=i=0;i<TOTST;i++) {
    if( oncalc[i] &&
	(st[i].imj == st[ic].imj) &&
	( !(abs(st[i].l-st[ic].l)&01) ) ) { /* in A there can be non-sph terms */ 
      idx[dim] = i;
      ++dim;
    }}
  /* dim has the dimension of As,Ns now */
  
  /* build the smaller symmetric matrices */
  for(i=0;i<dim;i++) {
    for(j=0;  j<=i;j++) {
      As[i][j] = As[j][i] = A[ idx[i] ][ idx[j] ];
      Ns[i][j] = Ns[j][i] = N[ idx[i] ][ idx[j] ];
    }}
      
  { 
    int jj,ii;
    perr("\n#%03d dim:%d\n", ic, dim);
    for(i=0;i<dim;i++) {
      for(j=0;j<dim;j++) {
	ii=idx[i];
	jj=idx[j];
	perr("[%03d: %s%+d%+d/2] [%03d: %s%+d%+d/2]\t",
	     ii, st[ii].name, st[ii].k, 2*st[ii].imj-1,
	     jj, st[jj].name, st[jj].k, 2*st[jj].imj-1); 
	perr("%+16.8lg %+16.8lg\n", As[i][j], Ns[i][j]);
      }
    }
  }
  /* generalized diagonalization Av=(lambda)Nv */
  diagrsg( As, Ns, dim, lambda, vr);

  select_v( L, v, ic, idx, lambda, vr, dim);


  for(i=0;i<dim;i++) {
    perr("Eigenvalue:%+16.8lg Selected:%+16.8lg\nEigenvector:\n",
	 lambda[i], *L);
    for(j=0;j<dim;j++) {
      perr("%+16.8lg\n", vr[j][i]);
    }
  }

}

static void select_v( double *L, double v[TOTST],
		      int ic, int idx[TOTST],
		      double lambda[TOTST], double vr[TOTST][TOTST], int dim)
/* From the list of eigenvectors in vr, select the one corresponding
   to the state ic. Put it in v and the eigenval in L */
{
  register int i;
  int imax, iccol;

  for(i=0;i<dim;i++) {
    imax = idxmax( vr, i, dim);
    if( idx[imax] == ic ) break;
  }
  /* i is now the column corresp to ic */
  iccol = i;
  *L = lambda[iccol];
    
  for(i=0;i<TOTST;i++)
    v[i] = 0.0;
  for(i=0;i<dim;i++)
    v[ idx[i] ] = vr[i][iccol];
}

static void join_matrix( 
		  double out[TOTST][TOTST],
		  double ectab[TOTST],
		  double ovlp[TOTST][TOTST],
		  double dvsic[TOTST][TOTST],
		  double DW[TOTST][TOTST],
		  double DWic[TOTST][TOTST],
		  int oncalc[TOTST]
		  )
{
  register int ic,jc;

  for( ic=0; ic<TOTST; ic++) {
    if(oncalc[ic]) {
      for( jc=0; jc<=ic; jc++) {
	if(oncalc[jc]) {
	  /* this must be symmetric */
	  out[jc][ic] = out[ic][jc] = 
	    ectab[jc]*ovlp[ic][jc] + 
	    dvsic[ic][jc]          + 
	    DW[ic][jc]             -
	    DWic[ic][jc];

	  perr("[%03d: %s%+d%+d/2] [%03d: %s%+d%+d/2]\t",
	       ic, st[ic].name, st[ic].k, 2*st[ic].imj-1,
	       jc, st[jc].name, st[jc].k, 2*st[jc].imj-1); 
	  perr("%+16.8lg = %+16.8lg + %+16.8lg + %+16.8lg - %+16.8lg VVVV\n",
	       out[ic][jc],
	       ectab[jc]*ovlp[ic][jc],
	       dvsic[ic][jc],
	       DW[ic][jc],
	       DWic[ic][jc] );

	}
      }
    }
  }

}

static void delta_sicpot_matrix(
		    double dvsic[TOTST][TOTST],
		    int i,
		    int oncalc[TOTST],
		    double *r2drdi,
		    double *vsic[TOTST], double *bsic[TOTST],
		    double *GCK[TOTST][2], double *FCK[TOTST][2]
		    )
/* We are interested only in the state ic
   thus we do not calculate the whole matrix
   but only those elements which are coupled with ic.
   For a spherical potential the states coupled are
   those with the same spinangular functions.
   The matrix dvsic can be calculated only once.
*/
{
  register int ic,jc,n;
  double *integrand;

  integrand = vector_alloc(nrmax);

  if( opt.usesic==1) {
    /*    ic = i; we are interested in the orbital i only */
    for( ic=0; ic<TOTST; ic++) {
      if(oncalc[ic]) {
	for( jc=0; jc<TOTST; jc++) {
	  if(oncalc[jc]) {
	    if(st[ic].k==st[jc].k && st[ic].imj==st[jc].imj) {
	      /* if the states have the same angular part ... */
	      if( st[ic].imj == -st[ic].l || st[ic].imj == 1+st[ic].l) {
		/* nsol = 1 */
		for(n=0;n<nrmax;n++) {
		  integrand[n] = ( GCK[ic][0][n] * GCK[jc][0][n] +
				   FCK[ic][0][n] * FCK[jc][0][n] ) * 
		    (vsic[jc][n]+bsic[jc][n] - vsic[i][n]-bsic[i][n]) *
		    r2drdi[n];
		}
	      } else { 
		/* nsol = 2; */
		for(n=0;n<nrmax;n++) {
		  integrand[n] = ( GCK[ic][0][n] * GCK[jc][0][n] +
				   GCK[ic][1][n] * GCK[jc][1][n] +
				   FCK[ic][0][n] * FCK[jc][0][n] +
				   FCK[ic][1][n] * FCK[jc][1][n] ) *
		    (vsic[jc][n]+bsic[jc][n] - vsic[i][n]-bsic[i][n]) *
		    r2drdi[n];		  
		}
	      }
	      /* no more symmetric */
	      dvsic[ic][jc] = simpson( integrand, nrmax, 1.0);
	    } else {
	      /* the angular part is orthogonal because vsic,bsic are spherical*/
	      dvsic[ic][jc] = 0.0;
	    }
/* 	    perr("[%03d: %s%+d%+d/2] [%03d: %s%+d%+d/2]\t", */
/* 		   ic, st[ic].name, st[ic].k, 2*st[ic].imj-1,  */
/* 		   jc, st[jc].name, st[jc].k, 2*st[jc].imj-1); */
/* 	    perr("%+16.8lg VVVVV\n", dvsic[ic][jc]); */
	  }
	}
      }
    }
  } else {
    printf("Non spherical SIC is possible only with the PZ81 SIC, for the moment\n");
    printf("opt.usesic: %d\n",opt.usesic);
  }

  vector_free(integrand);

}

static void overlap_matrix(
		    double ovlp[TOTST][TOTST],
		    int oncalc[TOTST],
		    double *r2drdi,
		    double *GCK[TOTST][2], double *FCK[TOTST][2]
		    )
{
  register int ic,jc,n;
  double *integrand;

  integrand = vector_alloc(nrmax);

  for( ic=0; ic<TOTST; ic++) {
    if(oncalc[ic]) {
      for( jc=0; jc<ic; jc++) {
	if(oncalc[jc]) {
	  if(st[ic].k==st[jc].k && st[ic].imj==st[jc].imj) {
	    /* if the states have the same angular part ... */
	    if( st[ic].imj == -st[ic].l || st[ic].imj == 1+st[ic].l) {
	      /* nsol = 1 */
	      for(n=0;n<nrmax;n++) {
		integrand[n] = ( GCK[ic][0][n] * GCK[jc][0][n] +
				 FCK[ic][0][n] * FCK[jc][0][n] ) * r2drdi[n];
	      }
	    } else { 
	      /* nsol = 2; */
	      for(n=0;n<nrmax;n++) {
		integrand[n] = ( GCK[ic][0][n] * GCK[jc][0][n] +
				 GCK[ic][1][n] * GCK[jc][1][n] +
				 FCK[ic][0][n] * FCK[jc][0][n] +
				 FCK[ic][1][n] * FCK[jc][1][n] ) * r2drdi[n];
	      }
	    }
	    ovlp[ic][jc] = ovlp[jc][ic] = simpson( integrand, nrmax, 1.0);
	  } else {
	    /* the angular part is orthogonal */
	    ovlp[ic][jc] = ovlp[jc][ic] = 0.0;
	  }
/* 	  printf("[%03d: %s%+d%+d/2] [%03d: %s%+d%+d/2]\t", */
/* 		ic, st[ic].name, st[ic].k, 2*st[ic].imj-1, */
/* 		jc, st[jc].name, st[jc].k, 2*st[jc].imj-1);  */
/* 	  printf("%+16.8lg VVVVV\n", ovlp[ic][jc]); */
	}
      }
      ovlp[ic][ic] = 1.0;
    }
  }

  vector_free(integrand);
}
