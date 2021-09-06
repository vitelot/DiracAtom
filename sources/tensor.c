/* TENSOR.C: utilities per tensori, matrici e vettori   */
/* V.D.P. Servedio */

#ifndef _V_TENSORC
#define _V_TENSORC

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "my_math.h"

void v_error(const char *fmt, ...)
{
  va_list ap;
	va_start(ap, fmt);
	vfprintf( stderr, fmt, ap);
	va_end(ap);
	exit(1);
}

double *vector_alloc(int dimx)
{
  double *t;

	t = (double *) malloc( (size_t) dimx * sizeof(double) );
	if (!t) v_error("Allocation ERROR in  funct VECTOR_ALLOC\n\n");

  return t;
}

void vector_cpy(double *dest, double *source, int dim)
{
  register int i;
  for(i=0;i<dim;i++) dest[i]=source[i];
}

void vector_normalize( double *v, int dim)
{
  register int i;
  double norm=0;

  for(i=0;i<dim;i++)
    norm += v[i]*v[i];

  norm = sqrt(norm);

  for(i=0;i<dim;i++)
    v[i] /= norm;
}

double **matrix_alloc(int dimx, int dimy)
{
  double **t;
  register int i;

	t = (double **) malloc( (size_t) dimx * sizeof(double *) );
	if (!t) v_error("Allocation ERROR 1 in  funct MATRIX_ALLOC\n\n");
	for (i=0; i<dimx; i++)
	{
	  t[i] = (double *) malloc( (size_t) dimy * sizeof(double) );
	  if (!t[i]) v_error("Allocation ERROR 2 in funct MATRIX_ALLOC\n\n");
	}

  return t;
}

void matrix_free(double **t, int dimx, int dimy)
{
  register int i;

	for (i=0; i<dimx; i++)
	  free( (double *) t[i]);

	free( (double **) t);
}

void matrix_transpose( double **A, int dim)
{
  register int i, j;
  double temp;

	for (i=0; i<dim; i++)
	  for (j=0; j<i; j++) {
	    temp = A[i][j];
	    A[i][j] = A[j][i];
	    A[j][i] = temp;
	  }
}

void matrix_cpy( double **dest, double **source, int dmx, int dmy)
{
  register int i, j;
	for (i=0; i<dmx; i++)
	  for (j=0; j<dmy; j++)
	    dest[i][j] = source[i][j];
}

double **matrix_load(const char *f_name, int *dimx, int *dimy)
{
  FILE *in;
  double **t;
  register int i, j; 
  long count=0;
  char fake[128];

	in = fopen( f_name, "r");
	if (!in) v_error( "Errore apertura file %s in function %s()\n\n",
			  f_name, "matrix_load");

/* Ignore the first line assuming it is a comment */
	fgets( fake, 127, in);

	i = fscanf( in, "%d %d\n", dimx, dimy);
	if (i!=2)
	  v_error("Error: matrix dimensions not found in file %s\n\n",
			  f_name);

	t = matrix_alloc( *dimx, *dimy );

	while( !feof(in) )
	{
	  for (i=0; i<*dimx; i++)
	  {
	    for (j=0; j<*dimy; j++)
	    {
	      if ( fscanf(in, "%lf\n", &t[i][j]) == 1 ) count++;
	    }
	    fscanf(in, "\n");
	  }
	}
	if (count != *dimx * *dimy )
	  v_error("Error: %ld datas loaded do not corresp. to all elem. in file %s\n\n", count, f_name);

	fclose(in);
  return t;
}

double ***tensor_alloc(int dimx, int dimy, int dimz)
{
  double ***t;
  register int i,j;

	t = (double ***) malloc( (size_t) dimx * sizeof(double **) );
	if (!t) v_error("Errore: malloc failure 1 \n\n");
	for (i=0; i<dimx; i++)
	{
	  t[i] = (double **) malloc( (size_t) dimy * sizeof(double *) );
	  if (!t[i]) v_error("Errore: malloc failure 2 \n\n");
	}
	for (i=0; i<dimx; i++)
	  for (j=0; j<dimy; j++)
	  {
	    t[i][j] = (double *) malloc( (size_t) dimz * sizeof(double) );
	    if (!t[i][j]) v_error("Errore: malloc failure 3 \n\n");
	  }
  return t;
}

void tensor_free(double ***t, int dimx, int dimy, int dimz)
{
  register int i,j;

	for (i=dimx-1; i>=0; i--)
	{
	  for (j=dimy-1; j>=0; j--)
	    free((double *) t[i][j]);

	  free((double **)t[i]);
	}

	free((double ***)t);
}

double ***tensor_load(const char *f_name, int *dimx, int *dimy, int *dimz)
{
  FILE *in;
  double ***t;
  register int i, j, k; 
  long count=0;

	in = fopen( f_name, "r");
	if (!in) v_error( "Errore apertura file %s in function %s()\n\n",
			  f_name, "tensor_load");

/* Explorer per motivi che ignoriamo preferisce le dimensioni specificate come 
 *	    Z Y X.
 */
	i = fscanf( in, "%d %d %d\n", dimz, dimy, dimx);
	if (i!=3) v_error("Errore: dimensioni tensore non trovate in file %s\n\n",
			  f_name);

	t = tensor_alloc( *dimx, *dimy, *dimz );

	while( !feof(in) )
	{
	  for (i=0; i<*dimx; i++)
	  {
	    for (j=0; j<*dimy; j++)
	    {
	      for (k=0; k<*dimz; k++)
	      {
		if ( fscanf(in, "%lf\n", &t[i][j][k]) == 1 ) count++;
	      }
	      fscanf(in, "\n");
	    }
	    fscanf(in, "\n");
	  }
	}
	if (count != *dimx * *dimy * *dimz)
	  v_error("Error: %ld dati caricati non corrisp. a tutti gli elem. in file %s\n\n", count, f_name);

	fclose(in);
  return t;
}

void tensor_cpy( double ***dest, double ***source, int dmx, int dmy, int dmz)
{
  register int i, j, k;
	for (i=0; i<dmx; i++)
	  for (j=0; j<dmy; j++)
	    for (k=0; k<dmz; k++)
	      dest[i][j][k] = source[i][j][k];
}


void tensor_print( double ***pot, int dmx, int dmy, int dmz, char *f_name)
{
  register int i, j, k;
  FILE *f_pot;

	f_pot = fopen( f_name, "w");
	if (!f_pot) v_error("Error opening file %s in funct tensor_print\n\n", f_name);
	  fprintf(f_pot, "%d %d %d\n", dmz, dmy, dmx);
	  for (i=0; i<dmx; i++)
	  {
	    for (j=0; j<dmy; j++)
	    {
	      for (k=0; k<dmz; k++)
	      {
		fprintf( f_pot, "%lf\n", pot[i][j][k]);
	      }
	      fprintf(f_pot, "\n");
	    }
	    fprintf(f_pot, "\n");
	  }
	fclose(f_pot);
}

#endif  /* _V_TENSORC  */
