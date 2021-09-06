#include "relat.h"
#include "globalvar.h"

#define  prbool(x) ((x)? "TRUE": "FALSE")
#define  COEFFCUTOFF 1e-16

/******************* EXTERN *********************************/
extern int idxmax(double vr[TOTST][TOTST], int col, int dim);
extern double Hec( double z, int n, int k);
extern void input_param( char *_file);
extern void printinfo(void);
extern double niocc( int onr[TOTSH], int onp[TOTST], double onpre[TOTST], 
		     char occ_file[128] );
extern int  scan_cmd_line_params( int, char *argv[], char *, int *, char *);
extern void print_rmeans( Dbp r, Dbp r2drdi,
			  int *onp, double *onp_real,
			  Dbp orbdens[TOTST], double *);
extern void print_base( char name[128], int *onp, int *printme, double *r,
		 Dbp orbdens[TOTST], double *onp_real,
		 Dbp rhochr, Dbp vt, Dbp vsic[TOTST],
		 Dbp orbgrad[TOTST]);
extern void print_avect(char filename[128], double Coeffcutoff, int row,
			int nel, int DWidx[TOTST], int onp[TOTST],
			double vr[TOTST][TOTST], double lambda[TOTST]);
extern void print_deltaw(char filename[128], int row, int onp[TOTST],
			 int DWidx[TOTST], double DW_r[TOTST][TOTST]);
extern void prepare_qn_tables(void);

extern void diag_matrix( 
		 double A[TOTST][TOTST],
		 double L[TOTST],
		 double vr[TOTST][TOTST],
		 int oncalc[TOTST]);

extern void sic_nsph(
	      int onp[TOTST],
	      int oncalc[TOTST], double ectab[TOTST],
	      double *r, double *drdi, double *r2drdi,
	      double *vsic[TOTST], double *bsic[TOTST], double *Vxcsic[TOTST],
	      double *vnuc,
	      double *GCK[TOTST][2], double *FCK[TOTST][2],
	      double DW[TOTST][TOTST],
	      double *cos_theta, double *weights,
	      double **Rho, double **Mx, double **Mz,
	      double vr[TOTST][TOTST],
	      double *tesic, double *teps );
/************************************************************/

/******************* STATIC *********************************/

static int maxindex( double vr[TOTST][TOTST], int onp[TOTST], 
		     int DWidx[TOTST], int m, int row);
/************************************************************/
