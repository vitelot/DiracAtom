#include "relat.h"
#include "core.h"

void core( 
  double *vv_nosic, double *bb_nosic, double z,
  double *rc, double *r2drdic, double *drdic,
  double onp[TOTST], int totiter,
  double *rhochr, double *rhospn, double *grad,
  double *orbdens[TOTST], double *orbspin[TOTST], double *orbgrad[TOTST],
  double *rhoconvrg,
  double ectab[TOTST], double befftab[TOTST], double SOtab[TOTST], 
  FILE *enrg,
  double *vsic[TOTST], double *bsic[TOTST],
  double *GCK[TOTST][2], double *FCK[TOTST][2]
  )
     /**********************************************************************
 It extracts also each eigenfunction and puts it in:
   GCK[][0][] = Large negativ k component (if with coupling) or gck no cpld.
   GCK[][1][] = Large positiv k component (if with coupling) or nothing
   FCK[][0][] = Small negativ k component (if with coupling) or fck no cpld.
   FCK[][1][] = Small positiv k component (if with coupling) or nothing

 INPUT:
  c          Light speed: c>>137 for non relativistic limit
  *vv_nosic  Potential without sic corrections
  *bb        Magnetic field
  z          Atomic number
  *rc	     Radial mesh
  *r2drdic   (=rc*rc*drdic)
  *drdic     Derivative of r respect to index (actually 'dr')     
  *onr       Non relativistic occupation numbers
  *onp       Polarized relativistic occupation numbers
  *enrg_file Name of Energy lvls file
  **vsic     Orbital dependent sic correction

 OUTPUT:
 *rhochr    Charge density
 *rhospn    Sz spin density
  bcor       Hyperfine b summed over all orbitals
  bcors      Idem but only 's' orbitals
  *ectab     List of all eigenvalues
  *befftab   List of the mean magn. field for each orbital
  *********************************************************************/

{ /* sub core start */
int i,  s, ic, ish, n, l, nqn, noBst[2];
int kap, nsol, nzero, nmatch, node, imj;
double bsum, ec, temp, sz, cgd[2], cgmd[2], cgo, var[4];
double spnchar, ecnob[2], deB1[2], ecold;
Dbp rhochrcur, rhospncur, integrand, gradcur, rho_old, b00, vv, bb, dovrc;
/* static int optmatchingvalue; */
double spinorbit;
  /*  
 ***gck,fck Actual normalized solutions of Dirac eq. (k-summed)
 ***gc,fc   Solutions to Dirac eq. projected on different k's (large,small)
 ***dp,dq   Differential of P,Q
 ***wp,wq   r*G, r*c*F = P, Q 
 */
Dbp gck[2][2], fck[2][2],
  gc[2][2],  fc[2][2],
  dp[2][2],  dq[2][2],
  wp[2][2],  wq[2][2],
  dpk[2][2], dqk[2][2];

#include "alloc_core.h"

  bsum = 0.0;
    
  switch( opt.matching) {
  case 1:
    Newton = newton; break;	/* Vitus version */
  case 2:
    Newton = newton2; break;	/* E.Engel version */
  default:
    Newton = newton_old;	/* H.Ebert version */
  }
	    
	      
  for( n=0; n<nrmax; n++) {
    vv[n] = vv_nosic[n];
    bb[n] = bb_nosic[n];
    dovrc[n] = drdic[n]/rc[n];
    rho_old[n] = rhochr[n];
    rhochr[n] = rhospn[n] = grad[n] = 0.0;
    b00[n] = 1e-7;
    bsum += fabs(bb[n]);
  }

  
  if( bsum < 1.0e-8) {
    if( !opt.nspin) {
      printf("No XC-field: Paramagnetic case forced (NSPIN=1)\n");
      opt.nspin=1;
    }
  }

  nqn = l = -1;
  for( ic=0; ic<TOTST; ic++) { /* 'for ic' */
    
    if( onp[ic]) {     /* 'if onp' :Jump if not occupied */

      if( nqn != st[ic].n || l != st[ic].l) {
	ish = 0;     /* ish: Counter of processed shell states */
      }

      nqn = st[ic].n;
      l   = st[ic].l;
      kap = st[ic].k;

      imj = st[ic].imj;
      if( imj == -l || imj == l+1) nsol = 1;
      else			   nsol = 2;

      if( kap < 0) s = 0;
      else	   s = 1;

      if( opt.usesic)
	sicpot( vv_nosic, bb_nosic, vsic[ic], bsic[ic], vv, bb);

      ec = ectab[ic];
      /* precalculation with Bxc=0 */
      if( l!=0 && ish==0 && opt.nspin==0) { /* 'if' */

	for(i=0;i<TOTST;i++) /* find the state in the shell with k<0 */
	  if( st[i].n==nqn && st[i].l==l && st[i].k<0) break;
	noBst[1] = i;
	for(i=0;i<TOTST;i++) /* find the state in the shell with k>0 */
	  if( st[i].n==nqn && st[i].l==l && st[i].k>0) break;
	noBst[0] = i;

	/* calculate eigenval for these two states */
	for(i=0;i<2;i++) { /* 'for i=1,2 ' */
	  ecold = ec;
	  do {
	    /* scan the energy to find the correct number of nodes */
/*  	    printf("State: [N:%d L:%d K:%+d Mj:%+d/2] ", */
/*  		   nqn, l, kap, 2*imj-1); */
/*  	    printf("ec1: %lf\n",ec); vvvvv */
	    ec = setnodes( z, ec, noBst[i], vv, b00, rc,
			   drdic, dovrc, &nmatch, &nzero,
			   gc, fc,
			   dp, dq,
			   wp, wq);
	    /* refine the energy matching inw and outw solutions */
	    ec = (*Newton)( z, ec, noBst[i], vv, b00, rc, 
			    drdic, dovrc, &nmatch, &nzero,
			    gc, fc,
			    dp, dq,
			    wp, wq, 
			    cgd, cgmd, &cgo, var);
	    
	    node = nodes( gc[s][s], 1, nmatch);
	    
	    if( node==nqn-l-1) break;  /* good solution */
	    if( node>nqn-l-1) ec = (ecold = ecorrect(ecold,2));
	    else              ec = (ecold = ecorrect(ecold,3));

	    /* perr("corevnew: NoB: Bad number of nodes %d of %d\n",
		   node,nqn-l-1);  vvvvv */
	    
	  } while(1);
	  ectab[noBst[i]] = ecnob[i] = ec; /* i==0 : no B energy with k>0 */

	  normalize(rc, drdic, r2drdic, 
		    nmatch, nzero, var, noBst[i],
		    gc, fc,
		    dp, dq,
		    gck, fck,
		    dpk, dqk,
		    cgd, cgmd,
		    rhochrcur, rhospncur, gradcur,
		    integrand);
	  /* calculate the en shift due to B. 1st ord. perturb.  */
	  for( n=0; n<nrmax; n++)
	    integrand[n] = rhochrcur[n] * r2drdic[n] * bb[n];
	  deB1[i] = simpson( integrand, nrmax, 1.);

	  if( fl.printenrg) {
	    fprintf( enrg, "Iter:%03d +++\t", totiter);
	    fprintf( enrg, "E[#%03d %s K:%+d  J: %d/2] = %14.8lf\t",
		     noBst[i], st[noBst[i]].name, st[noBst[i]].k,
		     (2*abs(st[noBst[i]].k)-1), ecnob[i]);
	    fprintf( enrg, "Beff: ++++++++++++++");
	    if(i==0) fprintf( enrg, "\n");
	  }
	}
	spinorbit = fabs(ecnob[1]-ecnob[0]);

	if( fl.printenrg)
	  fprintf( enrg, "\tS.O.: %13.8lf\n", spinorbit);
      }

      ish++;
      
      if( opt.nspin==0)
	SOtab[ic] = spinorbit;

      ec = ectab[ic];
      if( nsol==2 && opt.nspin==0) 	  /* <sz> = -2mj/(2k+1) */
	ec = deB1[(kap<0)?1:0] * ((kap>0)?-1.0:1.0) + ecnob[(kap<0)?1:0];
      /* it is better to overshoot, instead of the following */
/* 	ec = deB1[(kap<0)?1:0] *  (-2.0*st[ic].imj + 1.0)/(2.0*kap+1) */
/* 	  + ecnob[(kap<0)?1:0]; */

      ecold = ec;
      do {
/*  	printf("State: [N:%d L:%d K:%+d Mj:%+d/2] ", */
/*  	       nqn, l, kap, 2*imj-1); */
/*  	printf("ec2: %lf\n",ec); vvvvv */
	/* scan the energy to find the correct number of nodes */
	ec = setnodes( z, ec, ic, vv, bb, rc,
		       drdic, dovrc, &nmatch, &nzero,
		       gc, fc,
		       dp, dq,
		       wp, wq);
	/* refine the energy matching inw and outw solutions */
	ec = (*Newton)( z, ec, ic, vv, bb, rc, 
			drdic, dovrc, &nmatch, &nzero,
			gc, fc,
			dp, dq,
			wp, wq, 
			cgd, cgmd, &cgo, var);

	node = nodes( gc[s][s], 1, nmatch);

	if( node==nqn-l-1) break;  /* good solution */
	if( node>nqn-l-1) ec = (ecold = ecorrect(ecold,2));
	else              ec = (ecold = ecorrect(ecold,3));
	
	printf("corevnew: Bad number of nodes (%d instead of %d)\n",
	       node, nqn-l-1); /* vvvvv */

      } while(1);

      ec = var[0];
      /* normalize the wave functs */
      normalize(rc, drdic, r2drdic, 
	        nmatch, nzero, var, ic,
	        gc, fc,
	        dp, dq,
	        gck, fck,
	        dpk, dqk,
		cgd, cgmd,
		rhochrcur, rhospncur, gradcur,
		integrand);

      if( nsol == 2) {
	for( n=0; n<nrmax; n++) {
	  temp = 2.*cgo*gck[0][s][n]*gck[1][s][n];
	  rhospncur[n] += temp;
	}
      }

      /* Calculate total charge && spin density && gradient */
      for( n=0; n<nrmax; n++) {
	rhochr[n] += onp[ic]*rhochrcur[n];
	rhospn[n] += onp[ic]*rhospncur[n];
	grad[n] += onp[ic]*gradcur[n];
      }

      /* distance beetween the old and new density */
      for( *rhoconvrg=n=0; n<nrmax; n++) {
	temp  = fabs(rhochr[n]-rho_old[n]);
	if( temp > *rhoconvrg) *rhoconvrg=temp;
      }

      if( nsol==1) {
	for( n=0; n<nrmax; n++) {
	  GCK[ic][0][n] = gck[0][0][n];
	  FCK[ic][0][n] = fck[0][0][n];
	}
      }
      else {
	for( n=0; n<nrmax; n++) {
	  GCK[ic][0][n] = gck[0][s][n];
	  GCK[ic][1][n] = gck[1][s][n];
	  FCK[ic][0][n] = fck[0][s][n];
	  FCK[ic][1][n] = fck[1][s][n];
	}
      }

/* VVVVV */
/* 	for( n=0; n<nrmax; n++)  */
/* 	  printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t VVV\n", */
/* 		 rc[n], */
/* 		 rhochrcur[n], */
/* 		 GCK[ic][0][n], */
/* 		 GCK[ic][1][n],	      */
/* 		 FCK[ic][0][n], */
/* 		 FCK[ic][1][n]); */

      /* 	      if( usesic) */
      for( n=0; n<nrmax; n++) {
	orbdens[ic][n] = onp[ic]*rhochrcur[n];
	orbspin[ic][n] = onp[ic]*rhospncur[n];
	orbgrad[ic][n] = onp[ic]*gradcur[n];
      }

      ectab[ic] = ec;

      /*	Calculate the mean of B for each orbital */
      for( n=0; n<nrmax; n++)
	integrand[n] = rhochrcur[n] * r2drdic[n] * bb[n]; 
      befftab[ic] = simpson( integrand, nrmax, 1.); 

      /*	Calculate the spin character of each orbital */
      for( n=0; n<nrmax; n++)
	integrand[n] = rhospncur[n] * r2drdic[n]; 
      spnchar = onp[ic]*simpson( integrand, nrmax, 1.); 
      if( fl.printenrg) {
	fprintf( enrg, "Iter:%03d\t", totiter);
	fprintf( enrg, "E[#%03d %s K:%+d Mj:%+d/2] = %14.8lf\t",
		 ic, st[ic].name, kap, 2*imj-1, ec);
	fprintf( enrg, "Beff: %14.8lf\tBSz: %14.8lf\n",
		 befftab[ic], spnchar);
      }
      /* check if the solution has not collapsed to a wrong spin */
      if( opt.nspin==0 && nsol==2 && fabs(befftab[ic])>SOtab[ic])
	if( (kap>0 && spnchar*befftab[ic]>0) ||
	    (kap<0 && spnchar*befftab[ic]<0)   ) {
	  puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
	  puts("WARN: the following state might be with wrong spin");
	  printf("State #%d: [N:%d L:%d K:%+d Mj:%+d/2]\n",
		 ic, nqn, l, kap, (int)(2.0*imj-1));
	  puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
	}

    }
  }

  /*			Calculate total spin character		*/
  if( fl.printenrg) {
    for( n=0; n<nrmax; n++)
      integrand[n] = rhospn[n] * r2drdic[n];
    sz = simpson( integrand, nrmax, 1.);

    fprintf( enrg, "<Beta Sigma_z>: %lf\n", sz);
    fprintf( enrg, "----------------------------------------\n");
  }

#include "free_core.h"

}
