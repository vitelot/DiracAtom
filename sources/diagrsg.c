#include "relat.h"
#include "my_f2c.h"

#ifdef LINUX
#  define RSG rsg_
#else
#  define RSG rsg
#endif

#define _TRUE 1

void diagrsg( double As[TOTST][TOTST], double Ns[TOTST][TOTST],
	      int dim,
	      double lambda[TOTST], double vr[TOTST][TOTST])
{
  extern void RSG( int *, int *, double *, double *, double *, int *,
		   double *, double *, double *, int *); 

  double *tmp1, *tmp2;
  int alsovect=_TRUE, errcode;

  double **ar, **br, **vvr, *la, *ar_f, *br_f, *vr_f;
  register int i, j;

  ar  = matrix_alloc( dim, dim);
  br  = matrix_alloc( dim, dim);
  vvr = matrix_alloc( dim, dim);
  la    = vector_alloc( dim);
  tmp1  = vector_alloc( dim);
  tmp2  = vector_alloc( dim);

  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++) {
      ar[i][j] = As[i][j];
      br[i][j] = Ns[i][j];
    }

  ar_f = matrix_to_fortran( ar, dim, dim);
  br_f = matrix_to_fortran( br, dim, dim);
  vr_f = matrix_to_fortran( vvr, dim, dim);

  RSG( &dim, &dim, ar_f, br_f, la, &alsovect,
       vr_f, tmp1, tmp2, &errcode); 

  if( errcode)
    printf("WARNING: error code = %d returned by RSG\n",errcode);

  matrix_to_c( vvr, vr_f, dim, dim, _TRUE);
  free((void *) ar_f);
  free((void *) br_f);

  for (i=0; i<dim; i++) {
    lambda[i] = la[i];
    for (j=0; j<dim; j++) {
      vr[i][j] = vvr[i][j];
    }
  }

  matrix_free(ar,dim,dim);
  matrix_free(br,dim,dim);
  matrix_free(vvr,dim,dim);
  vector_free(la);
  vector_free(tmp1);
  vector_free(tmp2);
}

#undef _TRUE 
#undef RSG
