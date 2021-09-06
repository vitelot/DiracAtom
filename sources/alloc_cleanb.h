/* static void alloc_cleanb(Grid2d r2Br, Grid2d Ex, Grid2d Ez, */
/* 		  double **sBt, double *sin_theta ) */
{
  r2Br = matrix_alloc(nlmax, nrmax);
  Ex   = matrix_alloc(nlmax, nrmax);
  Ez   = matrix_alloc(nlmax, nrmax);
  sBt  = matrix_alloc(nrmax, nlmax);

  sin_theta = vector_alloc(nlmax);
}
