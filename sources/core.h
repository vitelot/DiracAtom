#define IPRINT   fl.iprint

extern int nodes( double *, int, int);

/* extern int turning_pt( double *v, double *r, int l, double e, int nzero); */

/* extern double ecorrect( double ec, int c); */

extern double Hec( double z, int n, int k);

extern void sicpot(
  double *, double *, double *, double *,
  double *, double *);

extern int n0val( double *v, double *r, double e);

extern double ecorrect( double ec, int c);

/* extern double elimit( double *v, double *r, int l); */

extern void outward(
  double z, double ec, int ic,
  double *vv, double *bb, double *rc, 
  double *drdic, double *dovrc, int *nmatch, int *nzero,
  double *gc[2][2], double *fc[2][2],
  double *dp[2][2], double *dq[2][2],
  double *wp[2][2], double *wq[2][2] );

extern double setnodes( 
  double z, double ec, int ic,
  double *vv, double *bb, double *rc, 
  double *drdic, double *dovrc, int *nmatch, int *nzero,
  double *gc[2][2], double *fc[2][2],
  double *dp[2][2], double *dq[2][2],
  double *wp[2][2], double *wq[2][2] );

double (*Newton)(
  double z, double ec, int ic,
  double *vv, double *bb, double *rc, 
  double *drdic, double *dovrc, int *nmatch, int *nzero,
  double *gc[2][2], double *fc[2][2],
  double *dp[2][2], double *dq[2][2],
  double *wp[2][2], double *wq[2][2], 
  double cgd[2], double cgmd[2], double *cgo, double var[4]);

extern double newton_old(
  double z, double ec, int ic,
  double *vv, double *bb, double *rc, 
  double *drdic, double *dovrc, int *nmatch, int *nzero,
  double *gc[2][2], double *fc[2][2],
  double *dp[2][2], double *dq[2][2],
  double *wp[2][2], double *wq[2][2], 
  double cgd[2], double cgmd[2], double *cgo, double var[4]);

extern double newton2(
  double z, double ec, int ic,
  double *vv, double *bb, double *rc, 
  double *drdic, double *dovrc, int *nmatch, int *nzero,
  double *gc[2][2], double *fc[2][2],
  double *dp[2][2], double *dq[2][2],
  double *wp[2][2], double *wq[2][2], 
  double cgd[2], double cgmd[2], double *cgo, double var[4]);

extern double newton(
  double z, double ec, int ic,
  double *vv, double *bb, double *rc, 
  double *drdic, double *dovrc, int *nmatch, int *nzero,
  double *gc[2][2], double *fc[2][2],
  double *dp[2][2], double *dq[2][2],
  double *wp[2][2], double *wq[2][2], 
  double cgd[2], double cgmd[2], double *cgo, double var[4]);

extern void normalize(
  double *rc, double *drdic, double *r2drdic, 
  int nmatch, int nzero, double var[4], int ic,
  double *gc[2][2],  double *fc[2][2],
  double *dp[2][2],  double *dq[2][2],
  double *gck[2][2], double *fck[2][2],
  double *dpk[2][2], double *dqk[2][2],
  double cgd[2], double cgmd[2],
  double *rhochrcur, double *rhospncur, double *gradcur,
  double *integrand);

/* void glue(int s, int nsol, int nmatch, int nzero, */
/* 	  double var[4], */
/* 	  double *gc[2][2],  double *fc[2][2], */
/* 	  double *dp[2][2],  double *dq[2][2] ); */

/* void glue_old(int s, int nsol, int nmatch, int nzero, */
/* 	      double var[4], */
/* 	      double *gc[2][2],  double *fc[2][2], */
/* 	      double *dp[2][2],  double *dq[2][2] ); */
