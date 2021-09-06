static double vdet( double a[4], double b[4], double c[4]);

double ferr( double var[4], int s, int nsol,
             double Pow[2][2], double Qow[2][2], 
             double Piw[2][2], double Qiw[2][2])
{
  int t;
  double abd, a[4], b[4], c[4], d[4];

  if( nsol==1)
    {
      var[1] = Qow[s][s]/Qiw[s][s];
      var[2] = var[3] =0.0;
      return (Pow[s][s]-var[1]*Piw[s][s]);
    }

  t = 1-s;
  a[1] = Qow[s][s];
  a[2] = Pow[t][s];
  a[3] = Qow[t][s];
  b[1] = Qow[s][t];
  b[2] = Pow[t][t];
  b[3] = Qow[t][t];
  c[1] = Qiw[s][s];
  c[2] = Piw[t][s];
  c[3] = Qiw[t][s];
  d[1] = Qiw[s][t];
  d[2] = Piw[t][t];
  d[3] = Qiw[t][t];

  abd = vdet( a,b,d);
  var[3] = - vdet(a,b,c)/abd;
  var[1] = - abd/vdet(b,c,d);
  var[2] = (var[1]*(c[1]+var[3]*d[1])-a[1])/b[1];
  return (Pow[s][s]+var[2]*Pow[s][t]-var[1]*(Piw[s][s]+var[3]*Piw[s][t]));
}

static double vdet( double a[4], double b[4], double c[4])
{
  return ( a[3]*b[2]*c[1] - a[2]*b[3]*c[1] - a[3]*b[1]*c[2]
         + a[1]*b[3]*c[2] + a[2]*b[1]*c[3] - a[1]*b[2]*c[3] );
}
