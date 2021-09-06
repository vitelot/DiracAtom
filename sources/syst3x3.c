#include <stdio.h>
#include <stdlib.h>

static double det3x3( double a[3][3]);

void syst3x3(double a[3][3], double b[3])
     /* Solves AX=B; returns B=X at exit. */
{
  int i,j;
  double x,y,z, det;
  double p[3][3];

  for( i=0;i<3;i++)
    for( j=0;j<3;j++)
      p[i][j] = a[i][j]; /* copy */

  det = det3x3( p);

  if( det==0.0) {
    fprintf(stderr, "DET=0 in syst3x3.c; Exit.\n");
    exit(1);
  }

  for( i=0;i<3;i++) p[i][0] = b[i];
  x = det3x3( p)/det;

  for( i=0;i<3;i++) { p[i][1] = b[i]; p[i][0]=a[i][0]; }
  y = det3x3( p)/det;

  for( i=0;i<3;i++) { p[i][2] = b[i]; p[i][1]=a[i][1]; }
  z = det3x3( p)/det;

  b[0]=x; b[1]=y; b[2]=z;
  return;
}

static double det3x3( double a[3][3])
{
  return 
    ( a[0][0]*a[1][1]*a[2][2] +
      a[0][1]*a[1][2]*a[2][0] +
      a[0][2]*a[1][0]*a[2][1] -
      a[0][0]*a[1][2]*a[2][1] -
      a[0][2]*a[1][1]*a[2][0] -
      a[0][1]*a[1][0]*a[2][2] );
}
