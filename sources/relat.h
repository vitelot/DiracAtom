#ifndef _RELATH
#define _RELATH

#define CURRENT_VERSION 3.15

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_math.h"
#include "tensor.h"

/* typedef double Grid2d[NLMAX][NRMAX]; */
typedef double *Dbp;
typedef double **Grid2d;
typedef int bool;
 
#define  TOTST	  280	      /* Total number of KS states		  */
#define  TOTSH	  28	      /* Total number of non-relativistic shells  */

#define RLoop(x) for((x)=0;(x)<nrmax;(x)++)
#define LLoop(x) for((x)=0;(x)<nlmax;(x)++)


struct Options {
  bool nspin;
  bool nocoll;
  bool nosph;
  bool force_sph;
  bool clean_field;
  bool nucl;
  int usesic;
  int nsphsic;
  int mesh;
  int nvers;
  int xchng;
  int wall;
  double bconst;
  double rmin, rmax;
  double p;
  double accrcy;
  char basis[128];
  int basis_recalc;
  int maxlexp;
  int matching;
  double c; /* = 2.0*137.03599; */
  int nrrectpt;	/* Number of points to calc in rect */
  int maxiter;
};

struct Flags {
  bool print_base;
  int  write_avect;
  int  write_deltaw;
  bool nonintoccnum;
  bool onlyocc;
  bool fixoccnum;	/* If TRUE use the occupied lvls not the lwst in enrgy */
  bool nsphorbsic;	/* Display SIC for each Nsph-orbital */
  bool compton_profile; /* Fourier transforms of the magnetization */
  bool printrho;
  bool printmagn;
  bool printpot;
  bool printenrg;
  bool r_powers;
  int iprint;		/* Level of infos printed  0->3	*/
};

struct State_qn {
  int n;	/* principal quantum number (non-relat analogous) */
  int l;	/* large component angular momentum */
  int k;	/* spin-orbit quantum number */
  int imj;	/* Jz=imj-0.5 */
  /* J=|k|-0.5 ; nsol=1 if k<0, nsol=2 if k>0 */
  char name[3]; /* name of the orbital belonging to */
};
  
/**************************************************************************/
/*			Here some extern declarations			  */

extern int nrmax, nlmax;
extern FILE *fperr;

extern struct Options opt;
extern struct Flags fl;
extern struct State_qn st[];
extern const int nqntab[];
extern const int lqntab[];
extern const int    ord[];
extern const char *orbname[];

extern void gauleg( double, double, double *, double *, int);

extern double integ2d( double **A, double *r2drdi, double *weights);

extern int prepat( int *, int *, char name[128],
		   int *, char atom_name[128]);

extern double LegPol( int, int, double);

extern double Y_nc( int, int, double);

extern double epsxc( double, double, int);

extern double   vxc( double, double, int);

extern void potcul_ryd( double *rad, double *drdi, double,
		       double *rho, double *, double *ur,
		       double *, double *, double *);

extern double simpson( double *, int, double);
extern double bodeint( double *, int, double, char, double *);

extern void diagsm2d( double, double, double,
		     double lambda[2], double U[2][2]);

extern void core( 
  double *vv_nosic, double *bb_nosic, double z,
  double *rc, double *r2drdic, double *drdic,
  double onp[TOTST], int totiter,
  double *rhochr, double *rhospn, double *grad,
  double *orbdens[TOTST], double *orbspin[TOTST], double *orbgrad[TOTST],
  double *rhoconvrg,
  double ectab[TOTST], double befftab[TOTST], double SOtab[TOTST], 
  FILE *enrg,
  double *vsic[TOTST], double *bsic[TOTST],
  double *GCK[TOTST][2], double *FCK[TOTST][2] );

extern void coredir( double, double, int, double, char,
	      double *vv, double *bb, double *rc, 
	      double *drdic, double *dovrc, int, int,
	      double *gc[2][2], double *fc[2][2],
	      double *dp[2][2], double *dq[2][2],
	      double *wp[2][2], double *wq[2][2], 
	      double Pow[2][2], double Qow[2][2], 
	      double Piw[2][2], double Qiw[2][2],
	      double cgd[2], double cgmd[2], double *, double);

extern void coreerr( double err[4], double var[4], int, int,
		     double Pow[2][2], double Qow[2][2], 
		     double Piw[2][2], double Qiw[2][2]);

extern double ferr( double var[4], int s, int nsol,
		    double Pow[2][2], double Qow[2][2], 
		    double Piw[2][2], double Qiw[2][2]);

extern void rinvgj( double a[4][4], double b[4][4], int);

extern void basis_xtrct(
   double *vv_nosic, double *bb_nosic, double z, int zat,
   double *rc, double *r2drdic, double *drdic,
   int onp[TOTST], int oncalc[TOTST],
   double *orbdens[TOTST], double *orbspin[TOTST], double *orbgrad[TOTST],
   double ectab[TOTST], double befftab[TOTST], double SOtab[TOTST], 
   FILE *enrg,
   double *vsic[TOTST], double *bsic[TOTST],
   double *GCK[TOTST][2], double *FCK[TOTST][2]
   );

extern void potat( double *rad, double *drdi, double,
		   double *rho, double *rhosp,
		   double *vt, double *bt, double *vnuc,
		   double *Vxc,
		   double v0[2], double *,
		   double *, double *, double *);

extern void dentopot_nc( Grid2d Rho, 
			 Grid2d Mx, Grid2d Mz,
			 Grid2d V,
			 Grid2d Bx, Grid2d Bz );

extern void deltaw_nc( int oncalc[TOTST],
		   double *cos_th, double *weights,
		   double *r2drdi,
		   double *vt, double *Vxc, double *bt,
		   Grid2d Vc_nsph, Grid2d V, 
		   Grid2d Bx, Grid2d Bz,
		   double *GCK[TOTST][2], double *FCK[TOTST][2],
		   double DW_r[TOTST][TOTST]);

extern void diagrs( double hr[TOTST][TOTST], int,
		   double lambda[TOTST], double vr[TOTST][TOTST]);

extern void magn2_nc( int nel,
		      int onp[TOTST], int oncalc[TOTST],
		      double vr[TOTST][TOTST], 
		      double *cos_theta, 
		      double *GCK[TOTST][2], double *FCK[TOTST][2], 
		      double **Rho, double **Mx, double **Mz );

extern double iexc_nc( Grid2d Rho, Grid2d Mx,
		       Grid2d Mz, double *r2drdi,
		       double *weights);

extern double ivexc_nc( Grid2d Rho, Grid2d Mx,
			Grid2d Mz, Grid2d V,
			Grid2d Bx, Grid2d Bz,
			double *r2drdi, double *weights);

extern double potcul_nsph( Grid2d Rho, Grid2d Vc,
		  double *cos_theta, double *weights,
		  double *rad, double *vnuc, double *drdi, int);

extern void cleanb( Grid2d Bx, Grid2d Bz,
		    double *cos_theta, double *weights,
		    Grid2d sources,
		    int maxlexp, double *ra, double *drdi);

extern void plotfield( Grid2d Bx, Grid2d Bz,
		       double *r, double *cos_theta,
		       double rmean, int NRRECTPT, char *name);

extern void plotscalar( Grid2d V,
			double *r, double *cos_theta,
			double rmean, int NRRECTPT, char *name);

extern double ylag(double, double *, double *, int, int, int);

extern void rmesh( double *r, double *drdi, int imax,
		   double rmin, double rmax,
		   int way);

extern void nsphint(double *vl, double *vlp, double *rhol,
		    double *rad, double *drdi, int l);

extern char l_to_c(int);

extern void perr(const char *fmt, ...);

extern void operr(char *errfilename);

extern void clerr(void);

extern void angular_avrg( double **, double *, int, int, double *);

#endif /* _RELATH */
