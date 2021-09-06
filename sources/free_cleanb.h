/* static void free_cleanb(Grid2d r2Br, Grid2d Ex, Grid2d Ez, */
/* 		  double **sBt, double *sin_theta ) */
{
  matrix_free(r2Br, nlmax, nrmax);
  matrix_free(Ex, nlmax, nrmax);
  matrix_free(Ez, nlmax, nrmax);
  matrix_free(sBt, nrmax, nlmax);

  vector_free(sin_theta);
}
