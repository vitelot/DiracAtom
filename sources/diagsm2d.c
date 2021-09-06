/* DIAGHM2D: finds eigenvalues and the unitary matrix tha diagonalize
             a 2*2 hermitian matrix
*/

#include <math.h>

void diagsm2d( double a, double b, double c,
	       double lambda[2],
	       double U[2][2])
/********************************************************
       || a   c ||
INPUT: ||       ||
       || c   b ||

OUTPUT: 
	lambda[0], lambda[1], eigenvalues
	lambda[0] >= lambda[1] if a>=b
	lambda[0]  < lambda[1] if a <b

	U, orthogonal matrix with corresponding eigenvectors on columns:
	   first column corresponds to lambda[0].
*********************************************************/
{
  double temp;
  int i;

  if( c==.0 )  /* if already diagonal, do nothing */
  {
    lambda[0] = a;
    lambda[1] = b;

    for( i=0; i<2; i++)
    {
      U[0][i] = (i==0)?1.:0.;
      U[1][i] = (i==1)?1.:0.;
    }
  }
  else
  {
    temp = sqrt( (a-b)*(a-b) + 4*c*c );
    if( a>=b )
    {
      lambda[0] = .5*( a + b + temp);
      lambda[1] = .5*( a + b - temp);
    }
    else
    {
      lambda[0] = .5*( a + b - temp);
      lambda[1] = .5*( a + b + temp);
    }
    for( i=0; i<2; i++)
    {
      temp = sqrt( 1. + ((lambda[i]-a)*(lambda[i]-a))/(c*c) );
      
      U[0][i] = 1./temp;
      U[1][i] = (lambda[i]-a)/(c * temp);
    }
  }
  
}

