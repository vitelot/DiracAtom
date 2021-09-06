#include "relat.h"

#ifndef CURRENT_DATE
#  define CURRENT_DATE "not available"
#endif

static void create_ini( char *);

void input_param( char *ini_file)
{
FILE *fp;
char line[128], *p;
double version=0.0;
int totlin=40;

    fp = fopen( ini_file, "r");
    if( !fp) {
      printf("Error opening file %s\n", ini_file);
      printf("Trying to create a default one.\n\n");
      create_ini( ini_file);
    }

    while( !feof(fp)) {
      fgets( line, 126, fp);
      if( line[0] == '#') continue; /* lines starting with # are comments */
      /* read up to # or line end */
      for( p=line; *p!='#' && *p!='\0'; p++) ;
      *p = '\0'; /* skip inline comments */

      if(sscanf(line, "VERSION %lf", &version)==1) {
	if( version != CURRENT_VERSION) {
	  fprintf( stderr, 
		   "The initialization file belongs to an old release.\n"
		   "Please delete it and run the program again to create\n"
		   "the new version.\n"
		   "Exiting.\n");
	  exit(1);
	}
	continue;
      }
      /*###################################################################*/
      if(sscanf(line, "NRMAX   = %d", &nrmax)==1) {totlin--; continue;}
      if(sscanf(line, "MESH    = %d", &opt.mesh)==1) {totlin--; continue;}
      if(sscanf(line, "NVERS   = %d", &opt.nvers)==1) {totlin--; continue;}
      if(sscanf(line, "XCHNG   = %d", &opt.xchng)==1) {totlin--; continue;}
      if(sscanf(line, "NSPIN   = %d", &opt.nspin)==1) {totlin--; continue;}
      if(sscanf(line, "USESIC  = %d", &opt.usesic)==1) {totlin--; continue;}
      if(sscanf(line, "NUCLEUS = %d", &opt.nucl)==1) {totlin--; continue;}
      if(sscanf(line, "NOSPH   = %d", &opt.nosph)==1) {totlin--; continue;}
      if(sscanf(line, "NOCOLL  = %d", &opt.nocoll)==1) {totlin--; continue;}
      if(sscanf(line, "CLEAN_FIELD  = %d", &opt.clean_field)==1) {totlin--; continue;}
      if(sscanf(line, "NONINTOCCNUM = %d", &fl.nonintoccnum)==1) {totlin--; continue;}
      /*###################################################################*/
      if(sscanf(line, "C = %lf", &opt.c)==1) {totlin--; continue;}
      if(sscanf(line, "BCONST = %lf", &opt.bconst)==1) {totlin--; continue;}
      if(sscanf(line, "RMIN   = %lf", &opt.rmin)==1) {totlin--; continue;}
      if(sscanf(line, "RMAX   = %lf", &opt.rmax)==1) {totlin--; continue;}
      if(sscanf(line, "WALL   = %d", &opt.wall)==1) {totlin--; continue;}
      if(sscanf(line, "MIXING = %lf", &opt.p)==1) {totlin--; continue;}
      if(sscanf(line, "ACCRCY = %lf", &opt.accrcy)==1) {totlin--; continue;}
      if(sscanf(line, "MATCHING = %d", &opt.matching)==1) {totlin--; continue;}
      if(sscanf(line, "MAXITER  = %d", &opt.maxiter)==1) {totlin--; continue;}
      /*###################################################################*/
      if(sscanf(line, "NLMAX   = %d", &nlmax)==1) {totlin--; continue;}
      if(sscanf(line, "BASIS = %s", opt.basis)==1) {totlin--; continue;}
      if(sscanf(line, "MAXLEXP = %d", &opt.maxlexp)==1) {totlin--; continue;}
      if(sscanf(line, "NSPHSIC = %d", &opt.nsphsic)==1) {totlin--; continue;}
      if(sscanf(line, "ONLYOCC = %d", &fl.onlyocc)==1) {totlin--; continue;}
      if(sscanf(line, "FIXOCCNUM = %d", &fl.fixoccnum)==1) {totlin--; continue;}
      if(sscanf(line, "BASIS_RECALC = %d", &opt.basis_recalc)==1) {totlin--; continue;}
      if(sscanf(line, "FORCE_SPH = %d", &opt.force_sph)==1) {totlin--; continue;}
      if(sscanf(line, "NRRECTPT = %d", &opt.nrrectpt)==1) {totlin--; continue;}
      /*###################################################################*/
      if(sscanf(line, "PRINT_BASE = %d", &fl.print_base)==1) {totlin--; continue;}
      if(sscanf(line, "WRITE_AVECT = %d", &fl.write_avect)==1) {totlin--; continue;}
      if(sscanf(line, "WRITE_DELTAW = %d", &fl.write_deltaw)==1) {totlin--; continue;}
      if(sscanf(line, "NSPHORBSIC = %d", &fl.nsphorbsic)==1) {totlin--; continue;}
      if(sscanf(line, "COMPTON_PROFILE = %d", &fl.compton_profile)==1) {totlin--; continue;}
      if(sscanf(line, "PRINTRHO = %d", &fl.printrho)==1) {totlin--; continue;}
      if(sscanf(line, "PRINTMAGN = %d", &fl.printmagn)==1) {totlin--; continue;}
      if(sscanf(line, "PRINTPOT = %d", &fl.printpot)==1) {totlin--; continue;}
      if(sscanf(line, "PRINTENRG = %d", &fl.printenrg)==1) {totlin--; continue;}
      if(sscanf(line, "R_POWERS = %d", &fl.r_powers)==1) {totlin--; continue;}
      if(sscanf(line, "IPRINT = %d", &fl.iprint)==1) {totlin--; continue;}
      /* if you reach here, the line is not recognized */
      fprintf( stderr, "The following line in %s was not recognized:\n%s\n"
	       "Exiting.\n", ini_file, line);
      exit(1);
    }
     fclose(fp);

     if(totlin > 0) {
       printf("There are %d lines less in your %s file. Check for missing ones. Exiting.\n",
	      totlin, ini_file);
       exit(1);
     }
     if(totlin < 0) {
       printf("There are %d lines more in your %s file. Check for duplicates. Exiting.\n",
	      -totlin, ini_file);
       exit(1);
     }

     /* no need of memory if only sph */
     if( !opt.nosph && !opt.nocoll) nlmax = 1;
}

#define  prbool(x) ((x)? "TRUE": "FALSE")

void printinfo(void)
{
  puts("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
  /* 	 printf("RHO_FILE  = \"%s\"\n", rho_file); */
  /* 	 printf("ENRG_FILE = \"%s\"\n", enrg_file); */
  /* 	 printf("MAGN_FILE = \"%s\"\n", magn_file); */
  /* 	 printf("POT_FILE  = \"%s\"\n", pot_file); */
  printf("IIIIIIII Main computational parameters IIIIIIIIIII\n");
  printf("NRMAX  = %d\n", nrmax);
  printf("NLMAX  = %d\n", nlmax);
  printf("IIIIIIII Master flags IIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
  printf("MESH  = %d\n", opt.mesh);
  printf("NVERS = %d\n", opt.nvers);
  printf("XCHNG = %d\n", opt.xchng);
  printf("NSPIN ------------- %s\n", prbool(opt.nspin));
  printf("USESIC ............ %s\n", prbool(opt.usesic));
  printf("NUCLEUS ----------- %s\n", prbool(opt.nucl));
  printf("NOSPH ............. %s\n", prbool(opt.nosph));
  printf("NOCOLL ------------ %s\n", prbool(opt.nocoll));
  printf("CLEAN_FIELD ....... %s\n", prbool(opt.clean_field));
  printf("NONINTOCCNUM ------ %s\n", prbool(fl.nonintoccnum));
  printf("IIIIIIII Basic parameters IIIIIIIIIIIIIIIIIIIIIIII\n");
  printf("C = %lf\n", opt.c);
  printf("BCONST = %lf\n", opt.bconst);
  printf("RMIN   = %lg\n", opt.rmin);
  printf("RMAX   = %lf\n", opt.rmax);
  printf("WALL   = %d\n", opt.wall);
  printf("MIXING = %lf\n", opt.p);
  printf("ACCRCY = %lg\n", opt.accrcy);
  printf("MATCHING .......... %s\n", ((opt.matching)?"NEW":"OLD"));
  printf("MAXITER = %d\n", opt.maxiter);
  printf("IIIIIIII Non-spherical stuff IIIIIIIIIIIIIIIIIIIIII\n");
  printf("BASIS = %s\n", opt.basis);
  printf("MAXLEXP = %d\n", opt.maxlexp);
  printf("NSPHSIC = --------- %s\n", prbool(opt.nsphsic));
  printf("ONLYOCC = ......... %s\n", prbool(fl.onlyocc));
  printf("FIXOCCNUM --------- %s\n", prbool(fl.fixoccnum));
  printf("BASIS_RECALC ...... %s\n", prbool(opt.basis_recalc));
  printf("FORCE_SPH --------- %s\n", prbool(opt.force_sph));
  printf("NRRECTPT = %d\n", opt.nrrectpt);
  printf("IIIIIIII Output control IIIIIIIIIIIIIIIIIIIIIIIIIII\n");
  printf("PRINT_BASE ........ %s\n", prbool(fl.print_base));
  printf("WRITE_AVECT ------- %s\n", prbool(fl.write_avect));
  printf("WRITE_DELTAW ...... %s\n", prbool(fl.write_deltaw));
  printf("NSPHORBSIC -------- %s\n", prbool(fl.nsphorbsic));
  printf("COMPTON_PROFILE ... %s\n", prbool(fl.compton_profile));
  printf("PRINTRHO ---------- %s\n", prbool(fl.printrho));
  printf("PRINTMAGN ......... %s\n", prbool(fl.printmagn));
  printf("PRINTPOT ---------- %s\n", prbool(fl.printpot));
  printf("PRINTENRG ......... %s\n", prbool(fl.printenrg));
  printf("R_POWERS ---------- %s\n", prbool(fl.r_powers));
  printf("IPRINT = %d\n", fl.iprint);
  puts("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");

  switch( opt.mesh) {
  case 0:
    puts("A linear radial mesh will be used."); break;
  case 1:
    puts("A quadratic radial mesh will be used."); break;
  default:
    puts("A logarithmic radial mesh will be used.");
  }

  if(opt.nspin==1)
    puts("Non-spin-polarized calculation.");
  else
    puts("Spin-polarized calculation.");

  switch( opt.nvers) {
  case -1:
    puts("Hartree only, no exchange-correlation functional at all.");
    break;	   
  case 0:
    puts("Exchange only, no correlation functional.");
    break;	   
  case 1: 
    puts("Perdew-Zunger 1981 parametrization of the Correlation functional.");
    break;
  case 2:
    puts("Gunnarson-Lundqvist 1976 parametrization of the Correlation functional.");
    break;
  case 3:
    puts("v.Barth-Hedin 1972 parametrization of the Correlation functional.");
    break;
  case 4:
    puts("Hedin-Lundqvist 1971 parametrization of the Correlation functional.");
    if( opt.nspin==0) {
      puts("Only non-spin-polarized calculation is possible."
	   "Change to NSPIN=1. Exiting.\n");
      exit(1);
    }
    break;
  case 5:
    puts("Perdew-Wang 1992 parametrization of the Correlation functional.");
    break;
  default:
    puts("Unexistent parametrization of the Correlation functional. Exiting.\n");
    exit(1);
  }
  if( opt.xchng==1)
    puts("Relativistic correction to the exchange switched on.");
  if( opt.xchng==-1)
    puts("No exchange energy functional.");
	   
  if( opt.usesic)
    puts("SIC will be applied.");
  if( opt.usesic==2)
    puts("The unconventional Vitus' SIC will be used.");
  
  if( opt.nsphsic && (opt.nosph || opt.nocoll)) {
    if( !opt.usesic) {
      puts("ATTENTION: I'll disable the full non spherical SIC because USESIC is switched off.");
      opt.nsphsic=0;
    }
      puts("Fully non spherical SIC will be applied.");
  }

  if( opt.matching)
    puts("A newer and faster (but not so stable) matching procedure will be applied.");
  else
    puts("The old and slow Ebert's matching procedure will be applied.");

  if( opt.wall)
    puts("A wall will be placed at RMAX");

  if( opt.force_sph) {
    puts("The density will be spherically averaged.");
    puts("   This option has sense if one calculates the spin-polarized case");
    puts("   with a variational method using the non spin polarized states");
    puts("   as basis set.\n");
  }

  if( fl.print_base==2)
    puts("Wavefunctions will be saved in a format suitable for Cooper minima calculations");

  if( fl.write_avect==2)
    puts("Eigenvectors and coupling will be written at any iteration.");

  if( fl.write_deltaw==2)
    puts("Matrix elements will be written at any iteration.");

  if( nlmax<3)
    printf("NLMAX has been forced to %d to save memory\n",nlmax);
  puts("");

  fflush(stdout);
} 

static void create_ini( char *file)
{
  FILE *fel;

    fel = fopen( file, "w");
    if( !fel)
      v_error("Cannot create default %s ini-file. Aborting.\n", file);
    fprintf( fel, 
     "############################################################################\n"
     "### ATOM's default initialization file #####################################\n"
     "VERSION %.2lf       # Last update: %s\n"
     "############################################################################\n",
	     CURRENT_VERSION, CURRENT_DATE);
    fprintf( fel, 
     "# Use the flag -h to have a list of command line options\n"
     "################### Master flags ###########################################\n"
     "NRMAX   =    500			# Number of radial points\n"
     "MESH    =      2			# 0=linear 1=quadratic 2=log\n"
     "NVERS   =      1			# Corr. en. parametrization:\n"
     "#... -1=NoXC 0=Xonly 1=PZ81 2=GL76 3=BH72 4=HL71 5=PW92.\n"
     "XCHNG   =      1			# -1=No Xchng 0=non rel. 1=rel. Xchng \n"
     "NSPIN   =      0			# 0=Spin polarized; 1=Para case\n"
     "USESIC  =      0			# 1=Perdew-Zunger-81 SIC; 2=Vitus SIC\n"
     "NUCLEUS =      0			# 0=point-like nucl. 1=homog. charged sphere\n"
     "NOSPH   =      0			# 1=do non spherical calculations\n"
     "NOCOLL  =      0			# 1=do non collinear calculations\n"
     "CLEAN_FIELD  = 0			# 1=removes sources from xc-field\n"
     "NONINTOCCNUM = 0			# 1=Allow non integer occ numbers\n"
     "################### Basic parameters #######################################\n"
     "C = 137.035999	                # Light speed\n"
     "BCONST = 0.0001			# Initial B to split degeneracies\n"
     "RMIN   = 0.0			# Minimum r in the grid. if 0.0 it is estimated\n"
     "RMAX   = 100.0			# Maximum r in the grid\n"
     "WALL   = 0			# 1=Wall with infinite height at RMAX\n"
     "MIXING = 0.30			# Mix parameter, low for heavy atoms\n"
     "ACCRCY = 1e-5			# Density accuracy for convergence\n"
     "MATCHING = 0			# Matching procedure: 1=newer,faster,unstable\n"
     "MAXITER  = 99			# Maximum number of iterations. Then Exits.\n"
     "################### Non-spherical stuff ####################################\n"
     "NLMAX   = 24			# Number of angular mesh points\n"
     "MAXLEXP = 10			# Max L for non spherical expansion\n"
     "BASIS = 1s2s2p3s3p3d4s4p4d4f5s5p  # Basis set: add or remove what you need \n"
     "NSPHSIC   =    0			# 1=apply full non spherical SIC\n"
     "ONLYOCC   =    0			# 1=only occupied states in the basis\n"
     "FIXOCCNUM =    1                  # 1=use the occupied lvls not the lwst in enrgy\n"
     "BASIS_RECALC = 0                  # 1=recalculate the basis set at each nosph step\n"
     "FORCE_SPH =    0                  # 1=force the density to be spherically averaged\n"
     "NRRECTPT  =    0			# Number of points to calc in rect\n"
     "################### Output control #########################################\n"
     "PRINT_BASE =      0		# Print base orbital densities: 2=for Cooper calc 3=all occupied\n"
     "WRITE_AVECT =     0		# Print the composition of eigenstates\n"
     "#...				  1=at the end 2=at each iteration\n"
     "WRITE_DELTAW =    0		# Write DELTAW matrix (see above)\n"
     "NSPHORBSIC =      0		# Display SIC for each Nsph-orbital\n"
     "COMPTON_PROFILE = 0		# Fourier transforms of the magnetization\n"
     "PRINTRHO =        0		# Write total rho and potentials\n"
     "PRINTMAGN =       0		# Write 2d magnetisation\n"
     "PRINTPOT =        0		# Write pot, bx, bz\n"
     "PRINTENRG =       0		# Print KSD eigenval: 0=last iter, 1=any iter\n"
     "R_POWERS =        0		# Averaged powers of r. 1=out to R_Powers file\n"
     "IPRINT    =       0		# Level of infos printed  0->3\n"
     "######################## END ###############################################\n"
     );
    fclose( fel);
    printf("Default %s ini-file created. You MUST edit it!\n", file);
    printf("Exiting...\n");
    exit(1);
}
