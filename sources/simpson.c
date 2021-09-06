double simpson( double *f, int nmax, double h)
{
  register int n;
  int simp, m;
  double sum;

  if( nmax%2 ) { /* if nmax is odd */
    m = nmax-1;
    sum = f[0];
    simp = -1;
    for ( n=1; n<m; n++) {
      simp = -simp;
      sum += (3+simp) * f[n];
    }
    sum += f[n];
    sum *= h/3.0;
  } else {

    double l0,l1,l2;

    l2 = f[nmax-1];
    l1 = f[nmax-2];
    l0 = f[nmax-3];

    sum = f[0];
    simp = -1;
    m = nmax-2;
    for ( n=1; n<m; n++) {
      simp = -simp;
      sum += (3+simp) * f[n];
    }
    sum += f[n] - l0/4 + 2*l1 + l2*1.25;
    sum *= h/3.0;
  }

 return sum;  
}  
