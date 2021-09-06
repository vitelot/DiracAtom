/* INPUT:
     V[nlmax][nrmax]	Vector to be averaged
     weights[nlmax]	Gauss-Legendre weights
   OUTPUT:
     Out[nrmax]
*/

void angular_avrg( double **V, double *weights, int nlmax, int nrmax,
		   double *Out)
{
  register int n,i;
  double tmp;

  for( n=0; n<nrmax; n++) {
    for( tmp=i=0; i<nlmax; i++) {
      tmp += weights[i] * V[i][n];
    }
    Out[n] = tmp/2.0; /* Integral[d(cos t)] = 2 {norm} */
  }
  return;
}
