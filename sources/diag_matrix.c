#include "relat.h"

static void convert( double  v[TOTST][TOTST], 
		     double vr[TOTST][TOTST],
		     int idx[TOTST], int dim)
/* the indexed eigen vector list v is converted
   into a direct non indexed one vr */
{
  register int i,j;

  for(i=0;i<TOTST;i++)
    for(j=0;j<TOTST;j++)
      vr[i][j] = 0.0;

  for(j=0;j<TOTST;j++)
    vr[j][j] = 1.0;

  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      vr[ idx[i] ][ idx[j] ] = v[i][j];

}

void diag_matrix( 
		 double A[TOTST][TOTST],
		 double L[TOTST],
		 double vr[TOTST][TOTST],
		 int oncalc[TOTST]
		 )
/* Symmetric diagonalization: Av=Lv
   returns the eigenvals and eigenvects (L,v)
*/
{
  register int i;
  int dim;
  int idx[TOTST];
  double v[TOTST][TOTST];

  for( dim=i=0; i<TOTST; i++)
    if( oncalc[i]) {
      idx[dim++] = i;
    }
  /* dim has the dimension of A now */

  diagrs( A, dim, L, v);

  convert( v, vr, idx, dim);
}
