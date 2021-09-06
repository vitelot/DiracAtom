#include "relat.h"

void plotfield( Grid2d Bx, Grid2d Bz,
		double *r, double *cos_theta,
		double rmean, int nrrectpt, char *name)
{
register i;
double step, x, z, Bxint, Bzint, pr, pcos_t;
double *Bxang, *Bzang;
FILE *fp;

  Bxang = vector_alloc(nlmax);
  Bzang = vector_alloc(nlmax);

  fp = fopen( name, "w");

  step = 2*rmean/nrrectpt;
  for( z = rmean; z>= -rmean; z-=step)
  {
    for( x=step; x<2*rmean; x+=step) 
    {
      pr = sqrt(x*x + z*z);
      pcos_t = z/pr;
      for( i=0; i<nlmax; i++)
      {
	Bxang[i] = ylag( pr, r, Bx[i], -1, 6, nrmax-1);
	Bzang[i] = ylag( pr, r, Bz[i], -1, 6, nrmax-1);
      }
      Bxint = ylag( pcos_t, cos_theta, Bxang, -1, 4, nlmax-1);
      Bzint = ylag( pcos_t, cos_theta, Bzang, -1, 4, nlmax-1);
      fprintf( fp, " %lf %lf", Bxint, Bzint);
    }
    fprintf( fp, "\n");
  }
  fclose( fp);

  vector_free(Bxang);
  vector_free(Bzang);
}

void plotscalar( Grid2d V,
		 double *r, double *cos_theta,
		 double rmean, int nrrectpt, char *name)
{
register i;
double step, x, z, Vint, pr, pcos_t;
double *Vang;
FILE *fp;

  Vang = vector_alloc(nlmax);

  fp = fopen( name, "w");

  step = 2*rmean/nrrectpt;
  for( z = rmean; z>= -rmean; z-=step)
  {
    for( x=step; x<2*rmean; x+=step) 
    {
      pr = sqrt(x*x + z*z);
      pcos_t = z/pr;
      for( i=0; i<nlmax; i++)
	Vang[i] = ylag( pr, r, V[i], -1, 6, nrmax-1);
      Vint = ylag( pcos_t, cos_theta, Vang, -1, 4, nlmax-1);
      fprintf( fp, " %lf", Vint);
    }
    fprintf( fp, "\n");
  }
  fclose( fp);

  vector_free(Vang);
}
