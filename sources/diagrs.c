#include "relat.h"
#include "my_f2c.h"

#ifdef LINUX
#  define RS rs_
#else
#  define RS rs
#endif

#define _TRUE 1

void diagrs( double hr[TOTST][TOTST], int row,
	     double lambda[TOTST], double vr[TOTST][TOTST])
{
  extern void RS( int *, int *, double *, double *, int *,
		  double *, double *, double *, int *); 

  double *tmp1, *tmp2;
  int alsovect=_TRUE, errcode;

  double **ar, **vvr, *la, *ar_f, *vr_f;
  int i, j;

  ar  = matrix_alloc( row, row);
  vvr = matrix_alloc( row, row);
  la    = vector_alloc( row);
  tmp1  = vector_alloc( row);
  tmp2  = vector_alloc( row);

  for (i=0; i<row; i++)
    for (j=0; j<row; j++)
    {
      ar[i][j] = hr[i][j];
    }

  ar_f = matrix_to_fortran( ar, row, row);
  vr_f = matrix_to_fortran( vvr, row, row);

  RS( &row, &row, ar_f, la, &alsovect,
      vr_f, tmp1, tmp2, &errcode); 

  if( errcode)
    printf("WARNING: error code = %d returned by RS\n",errcode);

  matrix_to_c( vvr, vr_f, row, row, _TRUE); 
  free((void *) ar_f);

  for (i=0; i<row; i++)
  {
    lambda[i] = la[i];
    for (j=0; j<row; j++)
    {
      vr[i][j] = vvr[i][j];
    }
  }

  matrix_free(ar,row,row);
  matrix_free(vvr,row,row);
  vector_free(la);
  vector_free(tmp1);
  vector_free(tmp2);
}

#undef _TRUE 
#undef RS
