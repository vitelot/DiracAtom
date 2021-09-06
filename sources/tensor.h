/* TENSOR.H: utilities per tensori, matrici e vettori   */
/* HEADER file				                */
/* V.D.P. Servedio */

#ifndef _V_TENSOR
#define _V_TENSOR

#include <stdlib.h>

extern void v_error(const char *fmt, ...);

extern double *vector_alloc(int dimx);

extern void matrix_transpose( double **A, int dim);

extern void vector_cpy(double *dest, double *source, int dim);

extern void vector_normalize( double *v, int dim);

#define vector_free(x)	free( (double *)(x))

extern double **matrix_alloc(int row, int col);

extern void matrix_free(double **x, int row, int col);

extern void matrix_cpy( double **dest, double **source, int dmx, int dmy);

extern double **matrix_load(const char *f_name, int *dimx, int *dimy);

extern double ***tensor_alloc(int dimx, int dimy, int dimz);

extern void tensor_free(double ***x, int dimx, int dimy, int dimz);

extern double ***tensor_load(const char *f_name, int *dimx, int *dimy, int *dimz);

extern void tensor_cpy( double ***dest, double ***source, int dmx, int dmy, int dmz);

extern void tensor_print( double ***pot, int dmx, int dmy, int dmz, char *f_name);

#endif  /* _V_TENSOR  */
