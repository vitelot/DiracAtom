/* MY_F2C.H: Header file containing routines useful for
 *	     c & fortran interfacing
 *   V.D.P. Servedio
 */

#ifndef _V_MY_F2C
#define _V_MY_F2C

#include <stdlib.h>
#include "tensor.h"

/*****************************************************************
 * MATRIX_TO_FORTRAN:                                            *
 *     This function creates the transpose matrix to be          *
 *     passed to fortran subroutines.                            *
 *     Actually the returned matrix is a vector with (col * row) *
 *     elements, so DON'T use it as a matrix in C codes!         *
 *     The 'a' pointer MUST be created with malloc.              *
 *     Otherwise the number of columns must be specified         *
 *     (for ex. 'double (*a)[10]').                              *
 *****************************************************************/
static double *matrix_to_fortran(double **a, int row, int col)
{
  double *f;
  int i, j;
      
      f = (double *) malloc( (size_t) (col * row * sizeof(double)));
      if (!f) v_error("Allocation ERROR in function MATRIX_TO_FORTRAN\n");

      if( row==col) { /* in case of square matrices, do transpose */
	for (i=0; i<row; i++)
	  for (j=0; j<col; j++)
	    f[i+j*row] = a[i][j];
      } else {
	for (i=0; i<row; i++)
	  for (j=0; j<col; j++)
	    f[j+i*col] = a[i][j]; /* j is the fastest index: A[i][j]=F(j,i) */
      }
  return f;
}

/*****************************************************************
 * MATRIX_TO_C:                                                  *
 *      Converts back fortran matrixes, previously created       *
 *      by 'matrix_to_fortran' function, into C matrixes doing   *
 *      simple transposition.                                    *
 *      If the 'doIfree' flag is set, memory is freed from       *
 *      'f_matrix'.                                              *
 *****************************************************************/
static void matrix_to_c(double **c_matrix, double *f_matrix, 
		 int row, int col, int doIfree)
{
  int i, j;

      if( row==col) { /* in case of square matrices, do transpose */
	for (i=0; i<row; i++)
	  for (j=0; j<col; j++)
	    c_matrix[j][i] = f_matrix[j+i*col];
      } else {
	for (i=0; i<row; i++)
	  for (j=0; j<col; j++)
	    c_matrix[i][j] = f_matrix[j+i*col];
      }

      if (doIfree)
	free( (void *) f_matrix);	  
}

#endif /*  _V_MY_F2C  */
