double matcherr( double var[4], int s, int nsol,
		 double Pow[2][2], double Qow[2][2], 
		 double Piw[2][2], double Qiw[2][2])
{
  extern void syst3x3(double a[3][3], double b[3]);

  int t;
  double  a[3][3], b[3];

  
  if( nsol==1) {
    var[1] = Qow[s][s]/Qiw[s][s];
    var[2] = var[3] =0.0;
    return (Pow[s][s]-var[1]*Piw[s][s]);
  }

  t = 1-s;

  a[0][0]= Qiw[s][s]; a[0][1]= Qiw[s][t]; a[0][2]= -Qow[s][t];
  a[1][0]= Piw[t][s]; a[1][1]= Piw[t][t]; a[1][2]= -Pow[t][t];
  a[2][0]= Qiw[t][s]; a[2][1]= Qiw[t][t]; a[2][2]= -Qow[t][t];

  b[0]= Qow[s][s];
  b[1]= Pow[t][s];
  b[2]= Qow[t][s];

  syst3x3(a,b);

  var[1] = b[0];
  var[2] = b[1];
  var[3] = b[2];

  return
    ( Pow[s][s] + var[3]*Pow[s][t] - var[1]*Piw[s][s] - var[2]*Piw[s][t] );
}
