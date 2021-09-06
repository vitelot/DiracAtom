#define labs(x) ( ((x)>0)? (x): (-(x)) )

/********************************************************************
 *                                                                  *
 * Lagrangian interpolation                                         *
 * xi is interpolated entry into x-array                            *
 * n is the order of lagrangian interpolation                       *
 * y is array from which ylag is obtained by interpolation          *
 * ind is the min-i for x(i).gt.xi                                  *
 * if ind=-1, x-array will be searched                              *
 * imax is max index of x-and y-arrays                              *
 *                                                                  *
 ********************************************************************/

double ylag(double xi, double *x, double *y, int ind, int n, int imax)
{
  register i, j; 
  int inu, inl;
  long jl, ju, jm;
  double p, d, s, xd, tmp;

    if( n>imax) n=imax;
    if( ind<0) {
      /* apply bisection to find the nearest point */
      jl = -1;
      ju = imax;
      while( ju-jl > 1) {
	jm = (ju+jl) >> 1;
	if( xi >= x[jm]) jl = jm;
	else ju = jm;
      }
      ind = (int) ju;
    }

    inl = ind - (n+1)/2;
    if( inl<0) inl = 0;

    inu = inl+n;
    
    if( inu>imax) {
      inl = imax-n;
      inu = imax;
    }

    s = 0.0;
    p = 1.0;
    for( j=inl; j<inu; j++)
    {
      /* if xi== x[j] then you get 0/0 otherwise */
      if( labs(tmp=xi-x[j])<1e-12) return y[j];
      p *= tmp;
      d = 1.0;
      for( i=inl; i<inu; i++)
      {
	xd = x[j];
	if( i==j) xd = xi;

	d *= xd - x[i];
      }
      s += y[j]/d;
    }

  return s*p;
}
