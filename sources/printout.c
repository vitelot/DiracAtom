#include "relat.h"

void print_rmeans( double *r, double *r2drdi,
		   int *onp, double *onp_real, double *orbdens[TOTST],
		   double *rmean)
{
/* calculate the averages of powers of r for all orbitals */
/* gives the value of <r> anyway			  */
FILE *fp;
register int i,j;
double temp, r2mean, r4mean, r6mean, r_3mean;
double *integrand;

  integrand = vector_alloc(nrmax);

  if( fl.r_powers) {
    fp = fopen("R_Powers","w");
    fprintf(fp,"Power of r averages:\n");
    fprintf(fp,"         STATE              <r^-3>        <r>         <r^2>        <r^4>        <r^6>\n");
    for( *rmean=j=0; j<TOTST; j++)
      if( onp[j])
	{
	  /* average of r */
	  for( i=0; i<nrmax; i++)
	    integrand[i] = r[i] * orbdens[j][i] * r2drdi[i];
	  temp = simpson( integrand, nrmax, 1.)/onp_real[j];
	  if( temp > *rmean) *rmean=temp;
	  /* average of r^-3 */
	  for( i=0; i<nrmax; i++)
	    integrand[i] = orbdens[j][i] * r2drdi[i]/(r[i]*r[i]*r[i]);
	  r_3mean= simpson( integrand, nrmax, 1.)/onp_real[j];
	  /* average of r^2 */
	  for( i=0; i<nrmax; i++)
	    integrand[i] = orbdens[j][i] * r2drdi[i] *r[i]*r[i];
	  r2mean= simpson( integrand, nrmax, 1.)/onp_real[j];
	    /* average of r^4 */
	  for( i=0; i<nrmax; i++)
	    integrand[i] = orbdens[j][i] * r2drdi[i] *r[i]*r[i]*r[i]*r[i];
	  r4mean= simpson( integrand, nrmax, 1.)/onp_real[j];
	  /* average of r^6 */
	  for( i=0; i<nrmax; i++)
	    integrand[i] = orbdens[j][i] * r2drdi[i] *r[i]*r[i]*r[i]*r[i]*r[i]*r[i];
	  r6mean= simpson( integrand, nrmax, 1.)/onp_real[j];
	  
	  fprintf(fp,"[#%03d %d%c K=%+d Jz=%+d/2] ",
		  j, st[j].n, l_to_c(st[j].l), st[j].k, 2*st[j].imj -1);
	  fprintf(fp,"%12.3lE %12.3lE %12.3lE %12.3lE %12.3lE\n",
		  r_3mean, temp, r2mean, r4mean, r6mean);
	}
    fprintf(fp,"\nAtomic radius: %.8lf\n\n", *rmean);
    printf("\nAtomic radius: %.8lf\n\n", *rmean);
    fclose(fp);
  } else {
    /* calculate the atomic radius as the maximum <r> among all orbitals */
    for( *rmean=j=0; j<TOTST; j++)
      if( onp[j]) {
	/* average of r */
	for( i=0; i<nrmax; i++)
	  integrand[i] = r[i] * orbdens[j][i] * r2drdi[i];
	temp = simpson( integrand, nrmax, 1.)/onp_real[j];
	if( temp > *rmean) *rmean=temp;
      }
  }
  
  vector_free(integrand);
}

void print_base( char name[128], int *onp, int *printme, double *r,
		 double *orbdens[TOTST], double *onp_real,
		 double *rhochr, double *vt, double *vsic[TOTST],
		 double *orbgrad[TOTST])
{
register int i,j;
char orb_files[128];
FILE *fp;

     if( fl.print_base==2) {/* for Photoemission */
       for( i=0; i<TOTST; i++) { 
 	if( onp[i] && printme[i]) { 
 	  sprintf(orb_files, "%s.Cooper.%03d", name, i);
	  /* Save: r, |Psi|^2, Sum|Psi|^2, rV */
 	  fp = fopen( orb_files, "w"); 
 	    for( j=0; j<nrmax; j++) {
	      if( opt.usesic)
		fprintf( fp, "% 16.8lg % 16.8lg % 16.8lg % 16.8lg\n",
			 r[j], orbdens[i][j]/onp_real[i],
			 rhochr[j], r[j]*(vt[j]-vsic[i][j]));
	      else
		fprintf( fp, "% 16.8lg % 16.8lg % 16.8lg % 16.8lg\n", 
			 r[j], orbdens[i][j]/onp_real[i],
			 rhochr[j], r[j]*vt[j]); 
	    }
 	  fclose( fp); 
 	  printf("%s KS orbital density saved for Cooper calcs\n", orb_files); 
 	} 
       } 
     } else if( fl.print_base==3) {
       for( i=0; i<TOTST; i++) { 
 	if( onp[i]) { 
 	  sprintf(orb_files, "%s_%d%c_%+d_%+3.1f.%03d",
		  name, st[i].n, l_to_c(st[i].l), st[i].k, st[i].imj-0.5, i);
	  /* Save: r, |Psi|^2, grad(|Psi|^2), rV */
 	  fp = fopen( orb_files, "w"); 
 	    for( j=0; j<nrmax; j++) {
	      if( opt.usesic)
		fprintf( fp, "% 16.8lg % 16.8lg % 16.8lg % 16.8lg\n",
			 r[j], orbdens[i][j]/onp_real[i],
			 orbgrad[i][j], r[j]*(vt[j]-vsic[i][j]));
	      else
		fprintf( fp, "% 16.8lg % 16.8lg % 16.8lg % 16.8lg\n", 
			 r[j], orbdens[i][j]/onp_real[i],
			 orbgrad[i][j], r[j]*vt[j]); 
	    }
 	  fclose( fp); 
 	  printf("%s KS orbital density saved\n", orb_files); 
 	} 
       } 
     }       

}

void print_avect(char filename[128], double COEFFCUTOFF, int row,
		 int nel, int DWidx[TOTST], int onp[TOTST],
		 double vr[TOTST][TOTST], double lambda[TOTST])
{
  FILE *fp;
  double temp,temp2;
  register int i,j,jj;
  
  fp = fopen( filename, "w");
  temp2 = sqrt(COEFFCUTOFF);
  for( i=0; i<row; i++) {
    if(i==nel)
      fprintf(fp, "############### Next ones are not occupied ###############\n\n");
    fprintf(fp, "Eigenvector#%03d \tEigenvalue: %+16.8lg Ryd\n", i,lambda[i]);
    for( temp=j=0; j<row; j++) {
      jj = DWidx[j];
      if( fabs(vr[j][i])>=temp2 )
	{
	  fprintf(fp, "[%03d:%s%s%+d%+d/2]\t",
		  jj, (onp[jj])?"  ":"U ", st[jj].name, st[jj].k, 2*st[jj].imj-1);
	  fprintf(fp,"%+16.8lE ---> %13.8lf %%\n", vr[j][i], 100*vr[j][i]*vr[j][i]);
	    temp += vr[j][i]*vr[j][i];
	}
    }
    fprintf(fp, "  Checksum:   1 =%15.12lf\n",temp); 
    fprintf(fp, "\n"); 
  }
  fclose( fp);
  printf("\n%s saved\n", filename);
}

void print_deltaw(char filename[128], int row, int onp[TOTST],
		 int DWidx[TOTST], double DW_r[TOTST][TOTST])
{
  int ii,jj, i,j;
  double temp;
  FILE *fp;

  fp = fopen( filename, "w");
  for( i=0; i<row; i++) {
    for( j=0; j<row; j++) {
      ii = DWidx[i];
      jj = DWidx[j];
      if( st[ii].imj == st[jj].imj ) {
	fprintf(fp, "[%03d:%s%s%+d%+d/2] [%03d:%s%s%+d%+d/2]\t",
		ii, (onp[ii])?"  ":"U ", st[ii].name, st[ii].k, 2*st[ii].imj-1,
		jj, (onp[jj])?"  ":"U ", st[jj].name, st[jj].k, 2*st[jj].imj-1); 
	temp = DW_r [ii] [jj]; 
	fprintf(fp,"%+16.8lg\n", temp);
      }  
    } 
    fprintf(fp, "\n"); 
  }  
  fclose( fp);
  printf("\n%s saved\n", filename);
}
