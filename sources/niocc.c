#include "relat.h"

/* static struct State { int occ; */
/* 		      int l; */
/* 		      int mj2; */
/* 		      int k; */
/* 		      int ord; */
/* 		      const char *name; }; */

static void create_occ( char *, int onp[TOTST]);

double niocc( int onr[TOTSH], int onp[TOTST], double onpre[TOTST], 
	      char occ_file[128] )
{
  const int nlshell=TOTSH;
  register ic, ilshell, imj, s;
  int l, nsol, kap[2];
  int stord[TOTST];
  double tmp;
  char line[128], *p;
  FILE *fp;

  for( ic=ilshell=0; ilshell<nlshell; ilshell++) {
    l = lqntab[ ord[ilshell] ];
    kap[0] = -l-1;
    kap[1] =   l  ;
    for( imj= -l; imj<=l+1; imj++) {
      nsol = 2;
      if( imj == -l || imj == l+1)
	nsol = 1;
      
      for( s=0; s<nsol; s++) {
/* 	state[ic].l = l; */
/* 	state[ic].mj2 = 2*imj-1; */
/* 	state[ic].k = kap[s]; */
	stord[ic] = ord[ilshell];
/* 	state[ic].name = orbname[ ord[ilshell] ]; */

	onr[ ord[ilshell] ] = 0;

	ic++;
      }			/* End 'for s'		*/
    }			/* End 'for mj'		*/
  }			/* End 'for ilshell'	*/
  
  fp = fopen( occ_file, "r");
  if( !fp) {
    printf("Error opening file %s\n", occ_file);
    printf("Trying to create a default one.\n\n");
    create_occ( occ_file, onp);
  }
  
  while( !feof(fp)) {
    fgets( line, 126, fp);
    if( line[0] == '#') continue; /* lines starting with # are comments */
    /* read up to # or line end */
    for( p=line; *p!='#' && *p!='\0'; p++) ;
    *p = '\0'; /* skip inline comments */

    sscanf(line, "%d %lf", &l, &tmp);
    onpre[ic=l] = tmp;
    onp[ic] = 0;
    if( tmp > 0.0) {
      onp[ic] = 1;
      ++onr[ stord[ic] ];
    }
  }
  
  for( tmp=ic=0; ic<TOTST; ic++)
	 tmp += onpre[ic];

  fclose(fp);

  return tmp; /* total charge */
}

static void create_occ( char *file, int onp[TOTST])
{
  register ic;
  char *p="xx";
  FILE *fel;

  fel = fopen( file, "w");
  if( !fel) {
    printf("Cannot create default %s occ-file. Aborting.\n", file);
    exit(1);
  }
  fprintf( fel, "# Occupation number file:\n"
	   "#   1st col.: state number\n"
	   "#   2nd col.: occupation number\n"
	   "#   3rd col.: description\n"
	   "#\n#\n");
  for( ic=0; ic<TOTST; ic++) {
    if( strcmp( p, st[ic].name)!=0 ) {
      p = st[ic].name;
      fprintf( fel, "#*************** %s **********************\n", p);
    }
    fprintf( fel, "%3d %8.4lf\t\t# %s K=%+d Mj=%+d/2 \n",
	     ic, (double) onp[ic],
	     st[ic].name, st[ic].k, 2*st[ic].imj-1);
  }

  fclose(fel);
  printf("Default %s occ-file created. You can edit it!\n", file);
  printf("Exiting...\n");
  exit(1);
}
