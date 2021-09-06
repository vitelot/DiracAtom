#include "my_math.h"
#include "relat.h"

static void pwgcorr(double, double, double, double,
		    double, double, double, double,
		    double *, double *);
static double pzcpot(double, double, double, double,
		     double, double, double, double);
static double pwcpot(double, double, double, double,
		     double *, double *);
static double hlcpot( double, double, double, double);
static double pzce(double, double, double, double,
                   double, double, double, double);
static double pwce(double, double, double);
static double hlce(double, double, double);
static double vxrel(double);
static double exrel(double);

/************************************************************************ 
 *     vxc (rho  (r), rho    (r))                                       *
 *             up        down                                           *
 *                                                                      *
 *     exchange and correlation potential of a homogeneous electron     *
 *     liquid at density rho                                            *
 *     if no spin-polarization: nspol=0                                 *
 *     if spin-polarization:                                            *
 *        nspol=-1: vxc=vxcminus                                        *
 *        else:     vxc=vxcplus                                         *
 *                                                                      *
 *     nvers=2 and 3 questionable (parameters!)                         *
 *     if option 'vxonly  ':   exchange only (vxc=vx) used              *
 ************************************************************************/

#define CUBSQR2	  1.2599210499	/* 2^(1/3) */
#define TWO4_3	  0.5198420998	/* 2^(4/3)-2 */
#define THRD	  0.3333333333  /* 1/3 */   
#define THRD4	  1.3333333333  /* 4/3 */
#define XFACTOR   0.4581652933  /* 3/4 * (9/(4*Pi^2))^(1/3) */

/*     correlation potential versions:                  */

/*     Perdew-Zunger (Ceperley-Alder) 1981 (nvers=1):   */

#define gamp  -0.1423
#define  b1p   1.0529
#define  b2p   0.3334
#define  aap   0.0311
#define  bbp  -0.048
#define  ccp   0.002
#define  ddp  -0.0116
#define gamf  -0.0843
#define  b1f   1.3981
#define  b2f   0.2611
#define  aaf   0.01555
#define  bbf  -0.0269
#define  ccf   0.0007
#define  ddf  -0.0048

/*     Gunnarsson-Lundqvist 1976 (nvers=2):     */

#define ap2  0.0333
#define bp2  11.4
#define af2  0.0203
#define bf2  15.9
#define vap2  0.0333
#define vaf2  0.0203

/*     v.Barth-Hedin 1972 (nvers=3):            */

#define ap3  0.0252
#define bp3  30.0
#define af3  0.0127
#define bf3  75.0
#define vap3  0.0252
#define vaf3  0.0127

/*     Hedin-Lundqvist 1971 (nvers=4, nspol=0): */

#define ap4  0.0225
#define bp4  21.0
#define vap4 0.0246

/*	Perdew-Wang 1992 (nvers=5) */
#define FZZ	1.709921

double vxc(double rhoup, double rhodown, int nspol)
{

  double zeta, vxcf, vxcp, f, epscf, epscp, epsxf, f1, 
	 f2, epsxp, af, bf, fd, ap, bp, rs, epsxcf, fd1, fd2,
         epsxcp, vaf, vcf, vap, rho, vcp, vxf, vxp, ret_val;

    if( opt.nvers==-1) return 0.0; /* only Hartree */

    if (rhoup < 0.) rhoup = 0.;
    if (rhodown < 0.) rhodown = 0.;

    rho = rhoup + rhodown;
    if (rho <= 1e-16) return 0.;
    
    rs = pow( .75/(M_PI*rho), THRD);
    vxp = -THRD4 * XFACTOR / rs;

    if( opt.xchng==1 ) vxp *= vxrel(rho);
    else if ( opt.xchng==-1 ) vxp=0.0;   /* No exchange */

    if( opt.nvers==0 && nspol==0) /* Exchange only, non s.p. */
      return vxp;

    epsxp = .75 * vxp;

    if (nspol != 0) {
      zeta = (rhoup-rhodown) / (rhoup+rhodown);

      f1 = pow( 1.+zeta, THRD4);
      f2 = pow( 1.-zeta, THRD4);
      f = (f1+f2-2.) / TWO4_3;
      
      fd1 = pow( 1.+zeta, THRD);
      fd2 = pow( 1.-zeta, THRD);
      fd = (fd1-fd2) /.75 /TWO4_3;
      vxf = CUBSQR2 * vxp;
      epsxf = CUBSQR2 * epsxp;

      if( opt.nvers==0) {		/* Exchange only */
	if (nspol == -1) {
	  return (vxp+(vxf-vxp)*f+(epsxf-epsxp)*(-1-zeta)*fd);
	} else {
	  return (vxp+(vxf-vxp)*f+(epsxf-epsxp)*( 1-zeta)*fd);
	}
      }
    }

    switch(opt.nvers) {
    case 1:
      vcp = pzcpot(gamp, b1p, b2p, aap, bbp, ccp, ddp, rs);
      epscp = pzce(gamp, b1p, b2p, aap, bbp, ccp, ddp, rs);
      break;
    case 2:
      vap = vap2;
      ap = ap2;
      bp = bp2;
      vaf = vaf2;
      af = af2;
      bf = bf2;
      vcp = hlcpot(vap, bp, rs, vxp);
      epscp = hlce(ap, bp, rs);
      break;
    case 3:
      vap = vap3;
      ap = ap3 *1.09; /* corrected value. prefactor is (2Pi/3)^(2/3) 2/3 */
      bp = bp3;
      vaf = vaf3;
      af = af3 *1.09; /* see above */
      bf = bf3;
      vcp = hlcpot(vap, bp, rs, vxp);
      epscp = hlce(ap, bp, rs);
      break;
    case 4:
      vap = vap4;
      ap = ap4;
      bp = bp4;
      vcp = hlcpot(vap, bp, rs, vxp);
      epscp = hlce(ap, bp, rs);
      break;
    case 5:
      if( nspol==0) {
	pwcpot( rs, 0.0, 0.0, 0.0, &vcp, &vcf);
	return vcp+vxp;
      } else {
	pwcpot( rs, zeta, f, fd, &vcp, &vcf);
	if(nspol==-1) return vcf + vxp + (vxf-vxp)*f + (epsxf-epsxp)*(-1-zeta)*fd;
	else          return vcp + vxp + (vxf-vxp)*f + (epsxf-epsxp)*(+1-zeta)*fd;
      }
      break;
    default:
      v_error("Correlation functional number %d not existent. Aborting\n",
	      opt.nvers);
    }

    vxcp = vxp + vcp;
    epsxcp = epsxp + epscp;
    if (nspol != 0) {
      switch(opt.nvers) {
      case 1:
        vcf = pzcpot(gamf, b1f, b2f, aaf, bbf, ccf, ddf, rs);
        epscf = pzce(gamf, b1f, b2f, aaf, bbf, ccf, ddf, rs);
	break;
      case 2:
      case 3:
      case 4:
        vcf = hlcpot(vaf, bf, rs, vxf);
        epscf = hlce(af, bf, rs);
	break;
      }
      vxcf = vxf + vcf;
      epsxcf = epsxf + epscf;
      
      if (nspol == -1) /* vxcminus */
	ret_val = vxcp + (vxcf-vxcp)*f + (epsxcf-epsxcp)*(-1-zeta)*fd; 
      else	       /* vxcplus  */
	ret_val = vxcp + (vxcf-vxcp)*f + (epsxcf-epsxcp)*(+1-zeta)*fd;
	
    }
    else
      ret_val = vxcp;
  
  return ret_val;

} /* vxc */

/************************************************************************ 
 *                                                                      *
 *     epsxc (rho  (r), rho    (r))                                     *
 *               up        down                                         *
 *                                                                      *
 *     exchange and correlation energy (in Hartree) per electron        *
 *     of a homogeneous electron liquid at density rho                  *
 *                                                                      *
 *     if no spin-polarization: nspol=0                                 *
 *                                                                      *
 ************************************************************************/

double epsxc(double rhoup, double rhodown, int nspol)
{
  double zeta, f, epscf, epscp, epsxf, f1, f2, epsxp, af, bf, ap,
         bp, rs, epsxcf, epsxcp, epscprs, rho, ret_val;

    if( opt.nvers==-1) return 0.0; /* only Hartree */

    f = 0.;

    if (rhoup < 0.) rhoup = 0.;
    if (rhodown < 0.) rhodown = 0.;
    
    rho = rhoup + rhodown;
    if (rho <= 1e-16) return .0;
    
    rs = pow( .75/(M_PI*rho), THRD);
    epsxp = -XFACTOR/rs;

    if( opt.xchng==1 ) epsxp *= exrel(rho);
    else if ( opt.xchng==-1 ) epsxp=0.0;   /* No exchange */

    if( opt.nvers==0 && nspol==0) { /* X only, non spin-polarized */
      return epsxp;
    }

    switch (opt.nvers) {
    case 0:
      break;
    case 1:
      epscp = pzce( gamp, b1p, b2p, aap, bbp, ccp, ddp, rs);
      break;
    case 2:
      ap = ap2;
      bp = bp2;
      af = af2;
      bf = bf2;
      epscp = hlce(ap, bp, rs);
      break;
    case 3:
      ap = ap3 *1.09; /* corrected value. prefactor is (2Pi/3)^(2/3) 2/3 */
      bp = bp3;
      af = af3 *1.09; /* see above */
      bf = bf3;
      epscp = hlce(ap, bp, rs);
      break;
    case 4:
      ap = ap4;
      bp = bp4;
      epscp = hlce(ap, bp, rs);
      break;
    case 5:
      pwgcorr( 0.031091, 0.21370,  7.5957, 3.5876, 1.6382, 0.49294, 1.0,
	       rs, &epscp, &epscprs);
    break;
    default:
      v_error("Correlation functional number %d not existent. Aborting\n",
	      opt.nvers);
    }

    epsxcf = 0.;
    if (nspol != 0) {
      zeta = (rhoup-rhodown) / (rhoup+rhodown);
      
      f1 = pow( 1.+zeta, THRD4);
      f2 = pow( 1.-zeta, THRD4);
      f = ( f1+f2-2.) / TWO4_3;

      epsxf = CUBSQR2 * epsxp;

      switch( opt.nvers) {
      case 0: /* X only */
	return (epsxp+(epsxf-epsxp)*f);
      case 1:
        epscf = pzce( gamf, b1f, b2f, aaf, bbf, ccf, ddf, rs);
	break;
      case 2:
      case 3:
      case 4:
        epscf = hlce( af, bf, rs);
	break;
      case 5:
	return pwce( rs, zeta, f) + epsxp + f*(epsxf-epsxp);
      }
      epsxcf = epsxf + epscf;
    }
    epsxcp = epsxp + epscp;

    ret_val = epsxcp + (epsxcf-epsxcp)*f;

  return  ret_val;

} /* epsxc */

static double pzcpot(double gam, double b1, double b2, double aa,
                     double bb, double cc, double dd, double rs)
{
  double h, p1, p2, ret_val;

    if (rs >= 1.)
    {
      p1 = 1. + b1*7.*sqrt(rs)/6. + b2*4.*rs/3.;
      p2 = 1. + b1*sqrt(rs) + b2*rs;
      ret_val = pzce(gam, b1, b2, aa, bb, cc, dd, rs)*p1/p2;
    }
    else
    {
      h = log(rs);
      ret_val =  aa*h + (bb-aa/3.) +cc*2.*rs*h/3. + (dd*2.-cc)*rs/3.;
    }

  return ret_val;
} /* pzcpot */

static double pwcpot(double rs, double zeta, double f, double fz,
		     double *vcup, double *vcdn)
{
  /* input:	rs, zeta, f, f'(zeta)
     output:	vcup, vcdown, spin stiffness */

  double eu, ep, eurs, eprs, alfm, alfrsm, z4, ec, ecrs, eczet, comm;

  pwgcorr( 0.031091, 0.21370,  7.5957, 3.5876, 1.6382, 0.49294, 1.0,
	   rs, &eu, &eurs);
  pwgcorr( 0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517, 1.0,
	   rs, &ep, &eprs);
  pwgcorr( 0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671, 1.0,
	   rs, &alfm, &alfrsm);

  z4 = zeta*zeta*zeta*zeta;

  ec = eu*(1.0-f*z4) + ep*f*z4 - alfm*f*(1-z4)/FZZ;
  ecrs = eurs*(1.0-f*z4) + eprs*f*z4 - alfrsm*f*(1.0-z4)/FZZ;
  eczet = 4.0*zeta*zeta*zeta*f*(ep-eu+alfm/FZZ) + 
          fz*(z4*ep-z4*eu-(1.0-z4)*alfm/FZZ);
  comm = ec - rs*ecrs/3.0 - zeta*eczet;
  *vcup = comm + eczet;
  *vcdn = comm - eczet;
  return -alfm;
}

static double hlcpot( double a, double b, double rs, double vx)
{
  double alph, x;

    x = rs/b;
    alph = a*rs* log( 1./x + 1.) + 2./3.;
    return vx*( alph*3./2. - 1.);
} /* hlcpot */



static double pzce(double gam, double b1, double b2, double aa,
                   double bb, double cc, double dd, double rs)
{
  double ch, ret_val;

    if (rs >= 1.)
    {
      ret_val = gam / ( b1*sqrt(rs) + 1. + b2*rs);
    }
    else
    {
      ch = log(rs);
      ret_val = aa*ch + bb + cc*rs*ch + dd*rs;
    }

  return ret_val;
} /* pzce */

static double pwce(double rs, double zeta, double f)
{
  double eu, ep, eurs, eprs, alfm, alfrsm, z4;

  pwgcorr( 0.031091, 0.21370,  7.5957, 3.5876, 1.6382, 0.49294, 1.0,
	   rs, &eu, &eurs);
  pwgcorr( 0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517, 1.0,
	   rs, &ep, &eprs);
  pwgcorr( 0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671, 1.0,
	   rs, &alfm, &alfrsm);

  z4 = zeta*zeta*zeta*zeta;
  
  return ( eu*(1.0-f*z4) + ep*f*z4 - alfm*f*(1.0-z4)/FZZ);

}

static void pwgcorr(double a, double a1, double b1, double b2,
		    double b3, double b4, double p, double rs,
		    double *gg, double *ggrs)
{
  /* PRB 45, 13244 (1992) */

  double p1,q0,q1,q2,q3,rs12,rs32,rsp;

  p1 = p + 1.0;
  q0 = -2.0*a*(1.0+a1*rs);
  rs12 = sqrt(rs);
  rs32 = rs12*rs12*rs12;
  rsp = pow(rs, p);
  q1 = 2.0*a*( b1*rs12 + b2*rs + b3*rs32 + b4*rs*rsp);
  q2 = log(1.0+1.0/q1);
  *gg = q0*q2;
  q3 = a*( b1/rs12 + 2.0*b2 + 3.0*b3*rs12 + 2.0*b4*p1*rsp);
  *ggrs = -2.0*a*a1*q2 - q0*q3/(q1*q1+q1); /* d(gg)/d(rs) */
  return;
}

static double hlce(double a, double b, double rs)
{
  double g, x, g1, g2;

    x = rs/b;

    g1 = 1. + x*x*x;
    g2 = log( 1./x + 1.);

    g = g1*g2 - x*x + x/2. - THRD;
    return -a*g;
} /* hlce */

#define CUB3PI2	  6.1873354526	/* 2*(3*Pi^2)^(1/3) ; 2* is for Ryd. units */

static double vxrel(double rho)
{
  double beta, eta;

    beta = CUB3PI2 * pow(rho, THRD)/ opt.c;
    eta = sqrt(1+beta*beta);
    return .5*(3*log(beta+eta)/(beta*eta) -1.0);
}

static double exrel(double rho)
{
  double eta, beta, temp;

    beta = CUB3PI2 * pow(rho, THRD)/ opt.c;
    eta = sqrt(1+beta*beta);
    temp = (beta*eta - log(beta+eta))/(beta*beta);
    return (1.0 - 1.5*temp*temp);
}
#undef CUBSQR2
#undef TWO4_3
#undef CUB3PI2
#undef THRD
#undef THRD4
#undef XFACTOR

#undef gamp
#undef  b1p
#undef  b2p
#undef  aap
#undef  bbp
#undef  ccp
#undef  ddp
#undef gamf
#undef  b1f
#undef  b2f
#undef  aaf
#undef  bbf
#undef  ccf
#undef  ddf

#undef ap2
#undef bp2
#undef af2
#undef bf2
#undef vap2
#undef vaf2

#undef ap3
#undef bp3
#undef af3
#undef bf3
#undef vap3
#undef vaf3

#undef ap4
#undef bp4
#undef vap4

#undef FZZ
