#include "relat.h"

static void print_manual(void);

static void print_short_help(void);

int scan_cmd_line_params( int argc, char *argv[],
			  char *atom_name, int *printme, char *ini_file)
{
  register i;
  int state;

     for( i=0; i<TOTST; i++)
       printme[i] = FALSE;

    while( argc--) { /* scan command line params */
      if( *argv[argc]==':') {
	strcpy( atom_name, argv[argc]+1);
	printf("Atom name %s found in the command line\n"
	       "Applying defaults...\n", atom_name);
      }
      if( *argv[argc]=='-') {
	switch(*(++argv[argc])) {
	case 'f':
	  if( strlen(argv[argc+1]))
	    sprintf( ini_file, "%s", argv[argc+1]);
	  break;
	case 'h':
	  print_short_help();
	  exit(0);
	case 'M':
	  print_manual();
	  exit(0);
	case 'p':
	  sscanf( argv[argc]+1, "%d", &state); /* row is used temporary */
	  if( state<TOTST && state>=0) {
	    printme[state] = TRUE;
	    printf("Orbital #%d will be printed\n", state);
	  }
	  break;
	default:
	  break;
	}
      }
    }
    return 0;
}

static void print_short_help(void)
{
  printf("Short help:\n");
  printf(" flags:\n"
	 "   -f filename : uses filename as ini file.\n"
	 "   -p# : output the orbital density number #.\n"
	 "   -M : Produce the postscript manual.\n"
	 "   -h : print this help and exit.\n"
	 " others:\n"
	 "   :elem : calculate element elem (Example :Na+ )\n");
}

static void print_manual(void)
{
FILE *fp;
char *file="MANUAL.ps";

  fp=fopen(file,"w");
#include "manual.h"
  fclose(fp);
  printf("\nThe file %s has been created.\n\n",file);
}
