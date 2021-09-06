#include "relat.h"

double Hec( double z, int n, int k)
{
  int n1;
  int c=opt.c;
  double temp;

  n1 = n - v_abs(k);
  temp = n1*.5*c+sqrt(k*k*c*c*.25-z*z);
  return .5*c*c* (-1.+1./sqrt(1.+z*z/(temp*temp)));
}

void sicpot(double *vv_nosic, double *bb_nosic,
	    double *vsic, double *bsic, /* these are orbital dependent */ 
	    double *vv, double *bb)
{
  register int n;
  
  if( opt.usesic==1) {
    /************************************************/
    /*          Perdew-Zunger-81 SIC		    */

    for( n=0; n<nrmax; n++) {
      vv[n] = vv_nosic[n] - vsic[n] - bsic[n];   
      bb[n] = bb_nosic[n];   
    } 
    /************************************************/
  } else if( opt.usesic==2) {
    for( n=0; n<nrmax; n++) {
      vv[n] = vv_nosic[n] - vsic[n];
      bb[n] = bb_nosic[n] - bsic[n];
    }
  }

  /* the following hybrid version is commented out */

  /*       		if( ic<48) */
  /*     		  for( n=0; n<nrmax; n++)    */
  /*     		  {    */
  /*  		    vv[n] = vv_nosic[n] - vsic[ic][n]; */
  /* 		    bb[n] = bb_nosic[n] - bsic[ic][n];  */
  /*  		  } */
  /*       		else */
  /*   		{ */
  /*       		  for( n=0; n<nrmax; n++)     */
  /*       		  {     */
  /*    		    vv[n] = vv_nosic[n] - vsic[ic][n] - bsic[ic][n];    */
  /*       		    bb[n] = bb_nosic[n];    */
  /*       		  } */
  /*   		} */
}

int nodes( double *f, int start, int end)
{
  register int n;
  int node=0;
  double f0,f1;

  for( n=start; n<end; n++) {
    f0 = f[n-1];
    f1 = f[ n ];
    if( (f0<0.0 && f1>0.0) || (f0>0.0 && f1<0.0) ) ++node;
    if( f1==0.0) { ++node; ++n; }
  }

  return node;
}    

double elimit( double *v, double *r, int l)
{
  int n;
  double val, elim, lll;

  lll = l*(l+1);
  elim = v[0] + lll/(r[0]*r[0]);
  for( n=1; n<nrmax; n++)
    if( (val=v[n]+ lll/(r[n]*r[n])) < elim)
      elim = val;

  return elim;
}

#define UNEND    600.0    /* Infinite for nzero finding			*/
int n0val( double *v, double *r, double e)
{
  register int n;

  for( n=0; n<nrmax-1; n++)
    if( (v[n]-e)*r[n]*r[n] > UNEND)
      break;

  return n ;
}
#undef UNEND

int turning_pt( double *v, double *r, int l, double e, int nzero)
{
  register int n;
  double lll;
  
  lll=l*(l+1);
  for( n=nzero; n>0; n--)
    if( (v[n] + lll/(r[n]*r[n]) - e) <0)
      break;

  return n;  
}

double ecorrect( double ec, int c)
{
  if(!opt.wall) {
    switch(c) {
    case 1: return ec*0.9;
    case 2: return ec*1.2;
    case 3: return ec*0.8;
    default: return ec*.9;
    }
  } else {    /* there can be positive energies */
    switch(c) {
    case 1: if(ec<-.1) return ec*0.9; else return ec+0.1;
    case 2: if(ec<-.1) return ec*1.2; else return ec-0.15;
    case 3: if(ec<-.1) return ec*0.8; else return ec+0.15;
    default: return ec*.9;
    }
  }
}

void outward(double z, double ec, int ic,
		 double *vv, double *bb, double *rc, 
		 double *drdic, double *dovrc, int *nmatch, int *nzero,
		 double *gc[2][2], double *fc[2][2],
		 double *dp[2][2], double *dq[2][2],
		 double *wp[2][2], double *wq[2][2] )
{
  double mj;
  int l,nqn;
  double cgd[2], cgmd[2], cgo, 
    Piw[2][2], Qiw[2][2],
    Pow[2][2], Qow[2][2];

  l   = st[ic].l;
  nqn = st[ic].n;
  mj  = (double)st[ic].imj - 0.5;


  /************************************************************
   *                        Find  nzero			*
   ************************************************************/
  *nzero = n0val( vv, rc, ec);
  if( !opt.wall && *nzero >= nrmax-1) {
    puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
    printf("Warn: E=%lg Nqn=%d L=%d K%d Mj=%+d/2: nzero set to %d\n",
	   ec, nqn, l, st[ic].k, (int)(mj*2.0), *nzero);
    puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
  }
  
  /************************************************************
   *                        Find  nmatch			*
   ************************************************************/
  *nmatch = turning_pt( vv, rc, l, ec, *nzero);

  if( !opt.wall && (*nmatch==0 || *nmatch>=*nzero-8)) {
    v_error("ERROR: St:%d E=%lg N=%d L=%d K=%d Mj=%+d/2 n0=%d: Nmatch not found\n",
	    ic, ec, nqn, l, st[ic].k, 2*st[ic].imj-1, *nzero);
  }
  
  if( opt.wall && (*nmatch==0 || *nmatch>=*nzero-8)) *nmatch = *nzero;
  
  coredir(z, ec, l, mj, 'o',	/* Outward integration */
	  vv, bb, rc,
	  drdic, dovrc, *nmatch, *nzero,
	  gc, fc, dp, dq,
	  wp , wq,
	  Pow, Qow,
	  Piw, Qiw,
	  cgd, cgmd, &cgo, 0.0);
  
  return;
}

/* looks for a reasonable starting energy, which gives */
/* the correct number of nodes and derivatives at the turning point */
#define INFINITELOOPERROR 1000
double setnodes( double z, double ec, int ic,
		 double *vv, double *bb, double *rc, 
		 double *drdic, double *dovrc, int *nmatch, int *nzero,
		 double *gc[2][2], double *fc[2][2],
		 double *dp[2][2], double *dq[2][2],
		 double *wp[2][2], double *wq[2][2] )
{
  double rat, elim;
  double ecwide, ecleft, ecright; /* for bisection */
  int l, nqn, node, s;
  int i=0;
  
  if( st[ic].k < 0) s = 0;
  else		    s = 1;
  l   = st[ic].l;
  nqn = st[ic].n;
  
  /************************************************************
   *                        Find  potential minimum		*
   ************************************************************/
  elim = elimit( vv, rc, l);
  if(ec<elim) ec = elim*0.7;
  
  ecwide = fabs(ec*.3);
  ecleft = ec-ecwide;
  ecright= ec+ecwide;

  if(ecleft<elim) ecleft = elim*0.8;
  
  do { /* infinite loop */
    i++;
/*     perr("#%d elim=%lg ecleft=%lg ecright=%lg\n",i,elim,ecleft,ecright); */
    outward( z, ec, ic,
	     vv, bb, rc, 
	     drdic, dovrc, nmatch, nzero,
	     gc, fc,
	     dp, dq,
	     wp, wq);

    /* count nodes */
    node = nodes( gc[s][s], 1, *nmatch);
    /* slope at turning point */
    rat = gc[s][s][*nmatch]/gc[s][s][*nmatch-1];
	
    if( node==nqn-l-1 && rat>0.0 && rat<1.0) break;
    if(  node>nqn-l-1) ecright=ec;
    if( (node<nqn-l-1) || (node==nqn-l-1 && (rat<0.0 || rat>1.0)) ) ecleft =ec;

    if(fabs(ecright-ecleft)<1e-8) { /* wrong initial limits */
/*       perr("Wrong init: eleft=eright=%lg\n",ecleft); */
      if(  node>nqn-l-1) {
	if(opt.wall) ecleft=ec-ecwide;
	else         ecleft *= 1.1; /* should remain negative */
	if(ecleft<elim) ecleft = elim*0.8;
      }
      if( (node<nqn-l-1) || (node==nqn-l-1 && (rat<0.0 || rat>1.0)) ) {
	if(opt.wall) ecright =ec+ecwide;
	else         ecright *= .9; /* should remain negative */
      }
    }

/*     if( node==nqn-l-1) {			 correct # of nodes */
/*       if( rat<0.0 || rat>1.0) ec = (ecold += ecstepbadnodes); ecorrect(ec,1); but not the slope */
/*       else break;				  loop exit */
/*     } else if( node>nqn-l-1)  ec = (ecold -= ecstepbadnodes*2.3); ecorrect(ec,2); */
/*     else		      ec = (ecold += ecstepbadnodes*2.1); ecorrect(ec,3); */

/*     perr("setnodes: Bad number of nodes. %s ec=%lg node=%d rat=%lg\n", */
/*  	   st[ic].name, ec,node,rat); vvvvv */

/*     for(i=0;i< nrmax;i++) */
/*       printf("%lg\t%lg\t%lg   VVVVV\n",rc[i],gc[s][s][i],vv[i]); */

    ec = .5*(ecleft+ecright); /* bisection */

    if( i>INFINITELOOPERROR) {
      printf("WARNING: Loop interrupted in SETNODE\n");
      break;
    }

  } while(1);
  
  /*  printf("Bisection steps in SETNODE: %d\n", i); */
  
  return ec;
}
#undef INFINITELOOPERROR

#define DVSTEP   0.01	  /* Step for Newton-Raphson			*/
#define TOLVAR   1e-8	  /* Convergence criterion for Newton-Raphson	*/
#define ITERMAX  50       /* Max iterations to converge Newton-Raphson	*/
#define TRYMIX   0.01
/**************************************************************/
double newton_old(double z, double ec, int ic,
	      double *vv, double *bb, double *rc, 
	      double *drdic, double *dovrc, int *_nmatch, int *_nzero,
	      double *gc[2][2], double *fc[2][2],
	      double *dp[2][2], double *dq[2][2],
	      double *wp[2][2], double *wq[2][2], 
	      double cgd[2], double cgmd[2], double *cgo, double var[4])
{
#define nmatch (*_nmatch)
#define nzero  (*_nzero)
  int setbreak;
  int l, kap, nsol, s,t, iter, nqn, imj, nvar, j,n, iv,jv,ie;
  double err[4], errnew[4], varnew[4], dv[4], dedv[4][4], dvde[4][4];
  double Pow[2][2], Qow[2][2], Piw[2][2],  Qiw[2][2], now[2], niw[2];
  double mj, rr, ratt;
  int node;
  double ecold, ecstepbadnodes;
    
  ecold=ec;
  ecstepbadnodes = .01*fabs(ecold);

  l   = st[ic].l;
  imj = st[ic].imj;
  if( imj == -l || imj == l+1)  nsol = 1;
  else				nsol = 2;
  nvar = 2*nsol;

  kap = st[ic].k;
  if( kap < 0) s=0;
  else         s=1;
  t = 1-s;

  nqn = st[ic].n;
  mj  = (double)imj - 0.5;

  coredir(z, ec, l, mj, 'o',	/* Outward integration */
	  vv, bb, rc,
	  drdic, dovrc, nmatch, nzero,
	  gc, fc, dp, dq,
	  wp , wq,
	  Pow, Qow,
	  Piw, Qiw,
	  cgd, cgmd, cgo, 0.0);

  var[0] = ec;

  if(nmatch != nzero) {		/* no wall */

    coredir(z, ec, l, mj, 'i',     /* Inward integration */
	    vv, bb, rc,
	    drdic, dovrc, nmatch, nzero,
	    gc, fc, dp, dq,
	    wp , wq,
	    Pow, Qow,
	    Piw, Qiw,
	    cgd, cgmd, cgo, 0.0);

    /************************************************************
                       Start values for
                       matching parameters
    ************************************************************/

    var[1] = Pow[s][s]/Piw[s][s];
    if( nsol==2) {
      for( j=0; j<nsol; j++)
	now[j] = niw[j] = .0;
      for( n=0; n<nmatch; n++) {
	rr = rc[n]*rc[n]*rc[n]*drdic[n];	    /* it was r*r*r */
	for( j=0; j<nsol; j++) 
	  now[j] += gc[j][j][n]*gc[j][j][n]*rr;
      }

      for( n=nmatch; n<nzero; n++) {
	rr = rc[n]*rc[n]*rc[n]*drdic[n];	    /* it was r*r*r */
	for( j=0; j<nsol; j++)
	  niw[j] += gc[j][j][n]*gc[j][j][n]*rr;
      }

      ratt = Pow[t][t]/Piw[t][t];
      var[2] = TRYMIX*(now[s]+niw[s]*var[1])/(now[t]+niw[t]*ratt);
      var[3] = ratt * var[2]/var[1];
    } else
      for( iv=2; iv<4; iv++)
	err[iv]=errnew[iv]=var[iv]=varnew[iv]=dv[iv]=.0;

    coreerr( err, var, s, nsol, Pow,Qow, Piw,Qiw);

    for( iv=0; iv<nvar; iv++) 
      dv[iv] = var[iv]; 

    iter = 0;
    while( ++iter < ITERMAX) {
      for( iv=1; iv<nvar; iv++) {
	for( jv=0; jv<nvar; jv++)
	  varnew[jv] = var[jv];
      
	varnew[iv] = var[iv] + dv[iv]*DVSTEP;
      
	if( v_abs(dv[iv]/var[iv]) < TOLVAR)
	  varnew[iv] = var[iv]*(1.0+v_sign(DVSTEP*TOLVAR, dv[iv]));
      
	coreerr( errnew, varnew, s, nsol, Pow,Qow, Piw,Qiw);
      
	for( ie=0; ie<nvar; ie++) {
	  if( errnew[ie]-err[ie] == 0.0) {
	    dedv[ie][iv] = .0;
	    if( ie==iv && opt.nspin==1)
	      dedv[ie][iv] = 1.;
	  } else
	    dedv[ie][iv] = (errnew[ie]-err[ie])/(varnew[iv]-var[iv]);
	}
      }
					
      for( jv=0; jv<nvar; jv++)
	varnew[jv] = var[jv];
					
      varnew[0] = var[0] + dv[0]*DVSTEP;
      if( v_abs(dv[0]/var[0]) < TOLVAR)
	varnew[0] = var[0]*(1.0+v_sign(DVSTEP*TOLVAR, dv[0]));

      coredir(z, varnew[0], l, mj, 'o',   /* Outward integration */
	      vv, bb, rc,
	      drdic, dovrc, nmatch, nzero,
	      gc, fc, dp, dq,
	      wp , wq,
	      Pow, Qow,
	      Piw, Qiw,
	      cgd, cgmd, cgo, 0.0);
      coredir(z, varnew[0], l, mj, 'i',   /* Inward integration */
	      vv, bb, rc,
	      drdic, dovrc, nmatch, nzero,
	      gc, fc, dp, dq,
	      wp , wq,
	      Pow, Qow,
	      Piw, Qiw,
	      cgd, cgmd, cgo, 0.0);
    
      coreerr( errnew, varnew, s, nsol, Pow,Qow, Piw,Qiw);
    
      for( ie=0; ie<nvar; ie++)
	dedv[ie][0] = (errnew[ie]-err[ie])/(varnew[0]-var[0]);
    
      rinvgj( dvde, dedv, nvar);
    
      for( iv=0; iv<nvar; iv++) {
	dv[iv] = .0;
	for( ie=0; ie<nvar; ie++)
	  dv[iv] += dvde[iv][ie]*err[ie];
	var[iv] -= dv[iv];
      }

      if( var[0]>0.0 ) {
	puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
	printf("Warning: E=%16.8lf > 0 !!!\n", var[0]);
	printf("State: [N:%d L:%d K:%+d Mj:%+d/2]\n",
	       nqn, l, kap, (int)(2.0*mj));
	puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
 	var[0] = (ecold*=1.1);
/* 	var[0] = -1.0; */
      }

/*        printf("State: [N:%d L:%d K:%+d Mj:%+d/2] ", */
/*  	     nqn, l, kap, (int)(2.0*mj)); */
/*        printf("var0: %lf\n", var[0]); */

      coredir(z, var[0], l, mj, 'o',   /* Outward integration */
	      vv, bb, rc,
	      drdic, dovrc, nmatch, nzero,
	      gc, fc, dp, dq,
	      wp , wq,
	      Pow, Qow,
	      Piw, Qiw,
	      cgd, cgmd, cgo, 0.0);
      coredir(z, var[0], l, mj, 'i',   /* Inward integration */
	      vv, bb, rc,
	      drdic, dovrc, nmatch, nzero,
	      gc, fc, dp, dq,
	      wp , wq,
	      Pow, Qow,
	      Piw, Qiw,
	      cgd, cgmd, cgo, 0.0);
    
      coreerr( err, var, s, nsol, Pow,Qow, Piw,Qiw);
    
      ec = var[0];

      setbreak = TRUE;
      for( iv=0; iv<nvar; iv++) {
	if( fabs(dv[iv]/var[iv])>TOLVAR && fabs(var[iv])>1e-30)
	  setbreak=FALSE;
      }
		
      if( nodes( gc[s][s], 1, nmatch) != nqn-l-1) return ec; /* bad solution */
      if( setbreak) break;
    }			 /* End while */
  } else { /* possible wall */

    var[1]=var[3]=var[2]=0;
    dv[0] = var[0];

    iter = 0;
    while( ++iter < ITERMAX) { /* while */
      
      /* only the large components are forced to zero at nmatch=nzero */
      if( nsol==2) {
	var[2] = -Pow[t][s]/Pow[t][t];
	err[0] = Pow[s][s] + var[2]*Pow[s][t];
      } else { err[0] = Pow[s][s]; }


      if( fabs(dv[0]/var[0])<TOLVAR ) break;

      varnew[0] = var[0] + dv[0]*DVSTEP;
      if( fabs(dv[0]/var[0]) < TOLVAR)
	varnew[0] = var[0]*(1.0+v_sign(DVSTEP*TOLVAR, dv[0]));

      coredir(z, varnew[0], l, mj, 'o',   /* Outward integration */
	      vv, bb, rc,
	      drdic, dovrc, nmatch, nzero,
	      gc, fc, dp, dq,
	      wp , wq,
	      Pow, Qow,
	      Piw, Qiw,
	      cgd, cgmd, cgo, 0.0);

      if( nsol==2) {
	var[2] = -Pow[t][s]/Pow[t][t];
	errnew[0] = Pow[s][s] + var[2]*Pow[s][t];
      } else { errnew[0] = Pow[s][s]; }

      dv[0] = err[0]*(varnew[0]-var[0])/(errnew[0]-err[0]);

/*       printf("ITER:%d\tEC=%lg\tDV=%lg\tERR=%lg\tV2=%lg\n", vvvvv*/
/* 	     iter,var[0], dv[0], err[0], var[2]); */


      var[0] -= dv[0];

      coredir(z, var[0], l, mj, 'o',   /* Outward integration */
	      vv, bb, rc,
	      drdic, dovrc, nmatch, nzero,
	      gc, fc, dp, dq,
	      wp , wq,
	      Pow, Qow,
	      Piw, Qiw,
	      cgd, cgmd, cgo, 0.0);
      
      ec = var[0];
      node = nodes( gc[s][s], 1, nmatch);
      if( node < nqn-l-1) var[0] = (ecold += ecstepbadnodes);
      if( node > nqn-l-1) var[0] = (ecold -= ecstepbadnodes*.57);
    }
  }

  if( iter>=ITERMAX) {
    puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
    printf("Newton_old: Not converged after ITERMAX=%d steps:\n", ITERMAX);
    printf("     E=%16.8lg var[1]=%16.8lg var[2]=%16.8lg var[3]=%16.8lg\n",
	   ec, var[1], var[2], var[3]);
    printf(" dv[0]=%16.8lg  dv[1]=%16.8lg  dv[2]=%16.8lg  dv[3]=%16.8lg\n",
	   dv[0], dv[1], dv[2], dv[3]);
    printf("err[0]=%16.8lg err[1]=%16.8lg err[2]=%16.8lg err[3]=%16.8lg\n",
	   err[0], err[1], err[2], err[3]);
    puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
  }
  return ec;
}
#undef nmatch
#undef nzero

#undef DVSTEP
#undef TOLVAR
#undef ITERMAX
#undef TRYMIX

#define DVSTEP   0.0001	  /* Step for Newton-Raphson			*/
#define TOLVAR   1e-8	  /* Convergence criterion for Newton-Raphson	*/
#define ITERMAX  20       /* Max iterations to converge Newton-Raphson	*/
double newton(double z, double ec, int ic,
	      double *vv, double *bb, double *rc, 
	      double *drdic, double *dovrc, int *_nmatch, int *_nzero,
	      double *gc[2][2], double *fc[2][2],
	      double *dp[2][2], double *dq[2][2],
	      double *wp[2][2], double *wq[2][2], 
	      double cgd[2], double cgmd[2], double *cgo, double var[4])
{
#define nmatch (*_nmatch)
#define nzero  (*_nzero)
  extern double matcherr( double var[4], int s, int nsol,
			  double Pow[2][2], double Qow[2][2], 
			  double Piw[2][2], double Qiw[2][2]);

  int l, nsol, s, iter, nqn, imj;
  double mj, err, errnew, varnew[4], derrde, dv;
  double Pow[2][2], Qow[2][2], Piw[2][2],  Qiw[2][2];
  
  l   = st[ic].l;
  imj = st[ic].imj;
  if( imj == -l || imj == l+1)  nsol = 1;
  else				nsol = 2;

  if( st[ic].k < 0) s=0;
  else              s=1;

  nqn = st[ic].n;
  mj  = (double)imj - 0.5;

  coredir(z, ec, l, mj, 'o',	/* Outward integration */
	  vv, bb, rc,
	  drdic, dovrc, nmatch, nzero,
	  gc, fc, dp, dq,
	  wp , wq,
	  Pow, Qow,
	  Piw, Qiw,
	  cgd, cgmd, cgo, 0.0);
  coredir(z, ec, l, mj, 'i',     /* Inward integration */
	  vv, bb, rc,
	  drdic, dovrc, nmatch, nzero,
	  gc, fc, dp, dq,
	  wp , wq,
	  Pow, Qow,
	  Piw, Qiw,
	  cgd, cgmd, cgo, 0.0);
  
  err = matcherr( var, s, nsol, Pow,Qow, Piw,Qiw);

  iter = 0; var[0] = ec;
  while( ++iter < ITERMAX) { /* 'while iter' */
					
    varnew[0] = var[0] + DVSTEP;

    coredir(z, varnew[0], l, mj, 'o',   /* Outward integration */
	    vv, bb, rc,
	    drdic, dovrc, nmatch, nzero,
	    gc, fc, dp, dq,
	    wp , wq,
	    Pow, Qow,
	    Piw, Qiw,
	    cgd, cgmd, cgo, 0.0);
    coredir(z, varnew[0], l, mj, 'i',   /* Inward integration */
	    vv, bb, rc,
	    drdic, dovrc, nmatch, nzero,
	    gc, fc, dp, dq,
	    wp , wq,
	    Pow, Qow,
	    Piw, Qiw,
	    cgd, cgmd, cgo, 0.0);
    
    errnew = matcherr( varnew, s, nsol, Pow,Qow, Piw,Qiw);

    derrde = (errnew-err)/DVSTEP;

    dv = -err/derrde;
    if( fabs(dv)>fabs(var[0]/2.0)) dv /= 10.0;
    var[0] += dv;

    if( var[0] >= -1e-10) {
      puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
      printf("Warning: E=%16.8lf >= 0 !!!\n", var[0]);
      printf("State: [N:%d L:%d K:%+d Mj:%+d/2]\n",
	     nqn, l, st[ic].k, (int)(2.0*mj));
      puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
      var[0] = -1.0;
    }

    if( fabs(dv/var[0]) < TOLVAR) break;

    coredir(z, var[0], l, mj, 'o',   /* Outward integration */
	    vv, bb, rc,
	    drdic, dovrc, nmatch, nzero,
	    gc, fc, dp, dq,
	    wp , wq,
	    Pow, Qow,
	    Piw, Qiw,
	    cgd, cgmd, cgo, 0.0);

    if( nodes( gc[s][s], 1, nmatch) != nqn-l-1) return ec; /* bad solution */

    coredir(z, var[0], l, mj, 'i',   /* Inward integration */
	    vv, bb, rc,
	    drdic, dovrc, nmatch, nzero,
	    gc, fc, dp, dq,
	    wp , wq,
	    Pow, Qow,
	    Piw, Qiw,
	    cgd, cgmd, cgo, 0.0);
		
    err = matcherr( var, s, nsol, Pow,Qow, Piw,Qiw);

    ec = var[0];
  }
  if( iter>=ITERMAX) {
    puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
    printf("Newton: Not converged after ITERMAX=%d steps:\n",ITERMAX);
    printf("E=%16.8lf var[1]=%16.8lf var[2]=%16.8lf var[3]=%16.8lf\n",
	   ec, var[1], var[2], var[3]);
    printf(" dv=%16.8lf\t", dv);
    printf("err=%16.8lf\n", err);
    puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
  }
  return ec;
}
#undef nmatch
#undef nzero

#undef DVSTEP
#undef TOLVAR
#undef ITERMAX

#define DVSTEP   0.0001	  /* Step for Newton-Raphson			*/
#define TOLVAR   1e-12	  /* Convergence criterion for Newton-Raphson	*/
#define ITERMAX  20       /* Max iterations to converge Newton-Raphson	*/
double newton2(double z, double ec, int ic,
	      double *vv, double *bb, double *rc, 
	      double *drdic, double *dovrc, int *_nmatch, int *_nzero,
	      double *gc[2][2], double *fc[2][2],
	      double *dp[2][2], double *dq[2][2],
	      double *wp[2][2], double *wq[2][2], 
	      double cgd[2], double cgmd[2], double *cgo, double var[4])
{
#define nmatch (*_nmatch)
#define nzero  (*_nzero)

extern double matcherr( double var[4], int s, int nsol,
			double Pow[2][2], double Qow[2][2], 
			double Piw[2][2], double Qiw[2][2]);
void glue(int s, int nsol, int, int,
	  double var[4],
	  double *gc[2][2],  double *fc[2][2],
	  double *dp[2][2],  double *dq[2][2] );

int l, nsol, s,t, iter, nqn, imj, k,j,n, node;
double mj, err, ecold;
double Pow[2][2], Qow[2][2], Piw[2][2], Qiw[2][2];
double norm2, gck,fck, de, *integrand;
/* FILE *tmpf; */
/* char fname[128]; */

  integrand= vector_alloc(nrmax);
  
  l   = st[ic].l;
  imj = st[ic].imj;
  if( imj == -l || imj == l+1)  nsol = 1;
  else				nsol = 2;

  if( st[ic].k < 0) s=0;
  else              s=1;
  t = 1-s;

  nqn = st[ic].n;
  mj  = (double)imj - 0.5;

/*   printf("nsol=%d\n",nsol);  vvvvv */
  iter = 0; ecold=var[0] = ec;
  while( ++iter < ITERMAX) { /* 'while iter' */

/*     nzero =  n0val( vv, rc, ec); */
/*     nmatch = turning_pt( vv, rc, l, ec, nzero); */

    coredir(z, ec, l, mj, 'o',	/* Outward integration */
	    vv, bb, rc,
	    drdic, dovrc, nmatch, nzero,
	    gc, fc, dp, dq,
	    wp , wq,
	    Pow, Qow,
	    Piw, Qiw,
	    cgd, cgmd, cgo, 0.0);
    coredir(z, ec, l, mj, 'i',     /* Inward integration */
	    vv, bb, rc,
	    drdic, dovrc, nmatch, nzero,
	    gc, fc, dp, dq,
	    wp , wq,
	    Pow, Qow,
	    Piw, Qiw,
	    cgd, cgmd, cgo, 0.0);

    err = matcherr( var, s, nsol, Pow,Qow, Piw,Qiw);

/*     if(nsol==2) { */
/*     printf("e1=%lg e2=%lg e3=%lg e4=%lg\n", */
/* 	   Pow[s][s]+var[3]*Pow[s][t]-var[1]*Piw[s][s]-var[2]*Piw[s][t], */
/* 	   Qow[s][s]+var[3]*Qow[s][t]-var[1]*Qiw[s][s]-var[2]*Qiw[s][t], */
/* 	   Pow[t][s]+var[3]*Pow[t][t]-var[1]*Piw[t][s]-var[2]*Piw[t][t], */
/* 	   Qow[t][s]+var[3]*Qow[t][t]-var[1]*Qiw[t][s]-var[2]*Qiw[t][t]); */
/*     printf("%lg %lg %lg %lg\n%lg %lg %lg %lg\n%lg %lg %lg %lg\n%lg %lg %lg %lg\n\n", */
/* 	   Pow[s][s], Pow[s][t], Piw[s][s], Piw[s][t], */
/* 	   Pow[t][s], Pow[t][t], Piw[t][s], Piw[t][t], */
/* 	   Qow[s][s], Qow[s][t], Qiw[s][s], Qiw[s][t], */
/* 	   Qow[t][s], Qow[t][t], Qiw[t][s], Qiw[t][t]); */
/*     } */

    glue( s, nsol, nmatch, nzero, var, gc, fc, dp, dq);
    
/*     { */
/*       sprintf( fname, "N%dL%dK%+dIMJ%+d", nqn, l, st[ic].k,st[ic].imj); */
/*       tmpf = fopen(fname,"w"); */
/*       fprintf(tmpf,"# EC=%lg  Rmatch=%lg  R0=%lg\n", */
/* 	      ec, rc[nmatch], rc[nzero]); */
/*       fprintf(tmpf,"# v0=%lg  v1=%lg  v2=%lg  v3=%lg  ERR=%lg\n", */
/* 	      var[0], var[1], var[2], var[3], err); */
/*       if(nsol==2) */
/* 	for( n=0; n<nzero; n++) */
/* 	  fprintf(tmpf,"%-12lg %16lg %16lg %16lg %16lg %16lg %16lg %16lg %16lg\n", */
/* 		  rc[n], */
/* 		  gc[s][s][n], gc[t][s][n], */
/* 		  gc[s][t][n], gc[t][t][n], */
/* 		  fc[s][s][n], fc[t][s][n], */
/* 		  fc[s][t][n], fc[t][t][n]); */
/*       else */
/* 	for( n=0; n<nzero; n++) */
/* 	  fprintf(tmpf,"%-12lg %16lg %16lg\n", */
/* 		  rc[n], gc[s][s][n], fc[s][s][n]); */
/*       fclose(tmpf); */
/*     } */

    /* Calculate  norm^2      */
    norm2 = .0;
    for( k=0; k<nsol; k++) {
      for( n=0; n<nzero; n++) {
	gck = fck = 0.0;
	for( j=0; j<nsol; j++) {
	  gck += gc[k][j][n];
	  fck += fc[k][j][n];
	}
	integrand[n] = drdic[n] * rc[n] * rc[n] * (gck*gck + fck*fck);
      }
      norm2 += simpson(integrand, nzero, 1.0);
    }
    /* the following is E-<Psi|H|Psi>/<Psi|Psi>
       In deriving it, note that the derivative of a Heavyside
       function is a delta function */
    if( nsol==2)
      de = -(Qow[s][s] + Qow[s][t]*var[3]) *
	( (Pow[s][s]+var[3]*Pow[s][t]) - (var[1]*Piw[s][s]+var[2]*Piw[s][t]) )/norm2;
    else
      de = -Qow[s][s] * 
	(  Pow[s][s] - var[1]*Piw[s][s] )/norm2;

    /* vvvvv */
/*     if(nsol==2) */
/*       printf("ec=%lg\t de=%16lg norm2=%16lg Nm=%4d N0=%4d\n", */
/* 	     ec,de,norm2,nmatch,nzero); */

    ec += de;

    var[0] = ec;

    if( var[0] >= -1e-10) {
      puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
      printf("Warning: E=%16.8lf >= 0 !!!\n", var[0]);
      printf("State: [N:%d L:%d K:%+d Mj:%+d/2]\n",
	     nqn, l, st[ic].k, (int)(2.0*mj));
      puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
      ec = var[0] = (ecold*=2.0);
      continue;
    }
    
    node = nodes( gc[s][s], 1, nmatch);
    if( node != nqn-l-1) {
      printf("newton2: Bad number of nodes\n"); /* vvvvv */
      if( node > nqn-l-1) 
	ec = var[0] = (ecold *= 1.05);
      else
	ec = var[0] = (ecold *= 0.95);

      continue;
    }

    if( fabs(de/var[0]) < TOLVAR) break;
    if( fabs(de/ecold) > 0.7) {
      ecold *= 0.9;
      var[0] = ec = ecold;
      iter=0;
    /* vvvvv */
/*       printf("INVERSION\n"); */
      continue;
    }
  }
  if( iter>=ITERMAX) {
    puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
    printf("Newton2: Not converged after ITERMAX=%d steps:\n",ITERMAX);
    printf("E=%16.8lf var[1]=%16.8lf var[2]=%16.8lf var[3]=%16.8lf\n",
	   ec, var[1], var[2], var[3]);
    printf(" de=%16.8lf\t", de);
    printf("err=%16.8lf\n", err);
    puts("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
  }

  /* vvvvv */
/*   puts(""); */

  /* the following calls reload the gc and fc without var multiplication */
  coredir(z, ec, l, mj, 'o',	/* Outward integration */
	  vv, bb, rc,
	  drdic, dovrc, nmatch, nzero,
	  gc, fc, dp, dq,
	  wp , wq,
	  Pow, Qow,
	  Piw, Qiw,
	  cgd, cgmd, cgo, 0.0);
  coredir(z, ec, l, mj, 'i',     /* Inward integration */
	  vv, bb, rc,
	  drdic, dovrc, nmatch, nzero,
	  gc, fc, dp, dq,
	  wp , wq,
	  Pow, Qow,
	  Piw, Qiw,
	  cgd, cgmd, cgo, 0.0);

  vector_free(integrand);
  return ec;
}
#undef nmatch
#undef nzero

#undef DVSTEP
#undef TOLVAR
#undef ITERMAX

void glue_old(int s, int nsol, int nmatch, int nzero,
	      double var[4],
	      double *gc[2][2],  double *fc[2][2],
	      double *dp[2][2],  double *dq[2][2] )
{
  int t, n, k, j;

  t = 1-s;

  /************************************************************
   *           Glue wavefunctions according              *
   *                to matching conditions                    *
   ************************************************************/
  
  
  /*  g[kj]: k gives the angular behaviour, 
      j the two indipendent solutions.
      g[ss] can be summed only with g[st] which belongs to
      the other solution. */
  
  /*              Inward solution                             */
  for( n=nmatch; n<nzero; n++)
    for( j=0; j<nsol; j++)
      for( k=0; k<nsol; k++) {
	gc[k][j][n] *= var[1];
	fc[k][j][n] *= var[1];
	dp[k][j][n] *= var[1];
	dq[k][j][n] *= var[1];
      }
  if( nsol==2) {
    /*                                      Outward solution     */
    for( n=0; n<nmatch; n++)
      for( k=0; k<nsol; k++) {
	gc[k][t][n] *= var[2];
	fc[k][t][n] *= var[2];
	dp[k][t][n] *= var[2];
	dq[k][t][n] *= var[2];
      }
    /*                                       Inward solution     */
    for( n=nmatch; n<nzero; n++)
      for( k=0; k<nsol; k++) {
	gc[k][t][n] *= var[3];
	fc[k][t][n] *= var[3];
	dp[k][t][n] *= var[3];
	dq[k][t][n] *= var[3];
      }
  }
}

void glue(int s, int nsol, int nmatch, int nzero,
	  double var[4],
	  double *gc[2][2],  double *fc[2][2],
	  double *dp[2][2],  double *dq[2][2] )
{
  int t, n, k;

  t = 1-s;

  /************************************************************
   *           Glue wavefunctions according              *
   *                to matching conditions                    *
   ************************************************************/
  
  /*  g[kj]: k gives the angular behaviour, 
      j the two indipendent solutions.
      g[ss] can be summed only with g[st] which belongs to
      the other solution. */
  
  /*					Inward solution    */
  if( nsol==1) {
    for( n=nmatch; n<nzero; n++) {
      gc[s][s][n] *= var[1];
      fc[s][s][n] *= var[1];
      dp[s][s][n] *= var[1];
      dq[s][s][n] *= var[1];
    }
  } else {
    /*                                      Outward solution     */
    for( n=0; n<nmatch; n++)
      for( k=0; k<nsol; k++) {
	gc[k][t][n] *= var[3];
	fc[k][t][n] *= var[3];
	dp[k][t][n] *= var[3];
	dq[k][t][n] *= var[3];
      }
    /*                                       Inward solution     */
    for( n=nmatch; n<nzero; n++)
      for( k=0; k<nsol; k++) {
	gc[k][t][n] *= var[2];
	fc[k][t][n] *= var[2];
	dp[k][t][n] *= var[2];
	dq[k][t][n] *= var[2];
	
	gc[k][s][n] *= var[1];
	fc[k][s][n] *= var[1];
	dp[k][s][n] *= var[1];
	dq[k][s][n] *= var[1];
      }
  }
}

void normalize(double *rc, double *drdic, double *r2drdic, 
	       int nmatch, int nzero, double var[4], int ic,
	       double *gc[2][2],  double *fc[2][2],
	       double *dp[2][2],  double *dq[2][2],
	       double *gck[2][2], double *fck[2][2],
	       double *dpk[2][2], double *dqk[2][2],
	       double cgd[2], double cgmd[2],
	       double *rhochrcur, double *rhospncur, double *gradcur,
	       double *integrand)
{
void (*Glue)(int s, int nsol, int nmatch, int nzero,
	     double var[4],
	     double *gc[2][2],  double *fc[2][2],
	     double *dp[2][2],  double *dq[2][2] );

  int n, s, k, j, l, nsol;
  double norm, gck2, fck2, temp;

  l   = st[ic].l;

  if( st[ic].k < 0)  s=0;
  else               s=1;
  
  if( st[ic].imj == -l || st[ic].imj == l+1) nsol = 1;
  else					     nsol = 2;

  if( opt.matching) Glue = glue;
  else		    Glue = glue_old;

  (*Glue)( s, nsol, nmatch, nzero, var, gc, fc, dp, dq );

/*       sprintf( fname, "N%dL%dK%+d", nqn, l, kap[s]); */
/*       tmpf = fopen(fname,"w"); */
/*       for( n=0; n<nzero; n++) */
/* 	fprintf(tmpf,"%-12lg %16lg %16lg %16lg %16lg %16lg %16lg %16lg %16lg\n", */
/* 		rc[n], */
/* 		gc[s][s][n], gc[t][s][n], */
/* 		gc[s][t][n], gc[t][t][n], */
/* 		fc[s][s][n], fc[t][s][n], */
/* 		fc[s][t][n], fc[t][t][n]); */

/*       fclose(tmpf); */
  
  /*                              Sum for each K               */
  for( n=0; n<nzero; n++)
    for( k=0; k<nsol; k++) {
      gck[k][s][n]=fck[k][s][n]=.0;
      dpk[k][s][n]=dqk[k][s][n]=.0;
      for( j=0; j<nsol; j++) {
	gck[k][s][n] += gc[k][j][n];
	fck[k][s][n] += fc[k][j][n];
	dpk[k][s][n] += dp[k][j][n];
	dqk[k][s][n] += dq[k][j][n];
      }
    }

  /* Calculate  norm  and normalize to 1       */
  norm = .0;
  for( k=0; k<nsol; k++) {
    for( n=0; n<nzero; n++) {
      integrand[n] = r2drdic[n] *
	(gck[k][s][n]*gck[k][s][n]+fck[k][s][n]*fck[k][s][n]);
    }
    norm += simpson(integrand, nzero, 1.0);
  }
  norm = 1./sqrt(norm);
  for( n=0; n<nzero; n++)
    for( k=0; k<nsol; k++) {
      gck[k][s][n] *= norm;
      fck[k][s][n] *= norm;
      dpk[k][s][n] *= norm/drdic[n]; /* dP/dr */
      dqk[k][s][n] *= norm/drdic[n]; /* dQ/dr */
    }
  
  for( n=nzero; n<nrmax; n++)
    for( k=0; k<nsol; k++)
      gck[k][s][n] = fck[k][s][n] = 0.;

  /*	    Calculate charge && spin density && gradient	*/
  for( n=0; n<nrmax; n++) {
    rhochrcur[n] = .0;
    rhospncur[n] = .0;
    gradcur[n] = .0;
    for( k=0; k<nsol; k++) {
      gck2 = gck[k][s][n]*gck[k][s][n];
      fck2 = fck[k][s][n]*fck[k][s][n];
      temp = gck2 + fck2;
      rhochrcur[n] += temp;
      temp = gck2*cgd[k] - fck2*cgmd[k];
      rhospncur[n] += temp;
      gradcur[n] += (gck[k][s][n]*dpk[k][s][n]+
		     fck[k][s][n]*dqk[k][s][n]/opt.c-
		     gck2-fck2)/rc[n]*2.0;
    }
  }
}
