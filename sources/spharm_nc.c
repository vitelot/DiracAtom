#include "my_math.h"
#include <stdio.h>
#include <stdlib.h>

double LegPol( int l, int m, double x)
{
  double fact, pll, pmm, pmmp1, somx2;
  int i, ll;

  if( m<0 || m>l || fabs(x)>1.)
  { puts( "Bad args in funct LegPol"); exit(1); }

  pmm = 1.0;		/* Compute Pmm */
  if( m > 0)
  {
    somx2 = sqrt((1.-x)*(1.+x));
    fact = 1.0;
    for( i=1; i<=m; i++)
    {
      pmm *= -fact *somx2;
      fact += 2.0;
    }
  }
  if( l == m)
    return pmm;
  else			/* Compute Pm m+1 */
  {
    pmmp1 = x*(2*m+1)*pmm;
    if( l == (m+1))
      return pmmp1;
    else		/* Compute Pml  l>m+1 */
    {
      for( ll=m+2; ll<=l; ll++)
      {
	pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	pmm = pmmp1;
	pmmp1 = pll;
      }
      return pll;
    }
  }
}

double Y_nc( int l, int m, double cos_theta)
/* Gives the usual spherical harm. function but without exp(im*phi) */
{
  double x, pfact, absY;
  int i, sgn, mp;

  mp = m;
  if( m<0) mp = -m;
  if( mp > l) 
  { puts( "(|m|>l) in funct Y\n"); exit(1); }
   
  x = cos_theta;

  pfact=1.0; 
  if( mp>0)
    for( i=1; i<=2*mp; i++)
      pfact *= (double) (l-mp+i); 

  sgn = ((m<0 && (mp&01))? -1: +1);
  absY = sgn*sqrt( (2.*l+1.)/(4.0*M_PI*pfact))*LegPol( l, mp, x);

  return absY;
}
