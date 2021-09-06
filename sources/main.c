#include "main.h"

int main(int argc, char *argv[])
{
FILE *fp, *enrg;
register int i,j, counter=0;
int totiter, row; 
int onr[TOTSH], onp[TOTST], oncalc[TOTST], printme[TOTST], DWidx[TOTST], colocc[TOTST];
double onp_real[TOTST];
int nat, zat, nel, base_el;		/* Atomic number, Valency and # of electrons*/
int nsphstate, nsphstart, nsphend;

double *vt_p, *bt_p, *Vxc_p;
Dbp vt, bt, r, drdi, r2drdi, vnuc, rhochr, rhospn;
Dbp grad, vto, integrand, bto, Vxc;
Dbp vt0, bt0, Vxc0, Vc0;
Dbp GCK[TOTST][2], FCK[TOTST][2];
Dbp orbdens[TOTST], orbspin[TOTST], orbgrad[TOTST];
Dbp Vxcsic[TOTST], vsic[TOTST], bsic[TOTST];
double ectab[TOTST], befftab[TOTST], SOtab[TOTST];

const double Ryd = 13.6056981; /* eV */
double rhoconvrg, enconvrg;
double v0[2], q, ev, evxc, eexc, sum, temp, temp2, x,z, sin_theta;
double tebs, tebsic, teeexc, teevxc, tesic, esic, teee, te, te_old;
double rmean, te0, eta, eta0, enee_nsph, evs, evxcs, eexcs, rnuc, pmix;

double DW_r[TOTST][TOTST], hr[TOTST][TOTST], vr[TOTST][TOTST];
double lambda[TOTST], lambda_old[TOTST];

Dbp cos_theta, weights;
Grid2d Rho, Vxc_nsph, Rho_old, Vc_nsph;
Grid2d Mx, Mz, Bx, Bz, Mx_old, Mz_old, sources;

char ini_file[128]="relat.ini", 
     occ_file[128]="relat.occ",
       name[128],
       filename[128],
       rho_file[128],
       enrg_file[128],
       magn_file[128],
       orb_files[128],
       pot_file[128],
       atom_name[128]="";

  printf("*******************************************\n");
  printf("******* Program ATOM, version %.2lf ********\n",CURRENT_VERSION);
  printf("*******************************************\n\n");

      scan_cmd_line_params( argc, argv,
			    atom_name, printme, ini_file);

      input_param( ini_file);

      zat = prepat(onp, onr, name, &nat, atom_name);
      nel = nat-zat;

      sprintf(filename, "%s.errors", name);
      operr( filename);

      /* estimate rmin */
      if(opt.rmin==0.0) opt.rmin = 0.01/pow(nat,1.33);

      printinfo();
      
#include "alloc_main.h"

      prepare_qn_tables();

/* prepare for non int occupation numbers */
     for( i=0; i<TOTST; i++)
       onp_real[i] = (double) onp[i];

     if( fl.nonintoccnum) {
       temp = niocc( onr, onp, onp_real, occ_file);       
       printf("Total charge check for non integer occ.num.: %lf\n",temp);
     }
/****************************************/

     sprintf(rho_file, "%s.Rho", name);
     sprintf(enrg_file, "%s.en", name);
     sprintf(magn_file, "%s.M", name);

  opt.c *= 2.0; /* Required by Rydberg units */

     if( opt.nspin) opt.bconst=.0;

     for( i=0; i<TOTST; i++)
     {
       DWidx[i] = -1;
       oncalc[i] = SOtab[i] = ectab[i] = .0;
       befftab[i] = opt.bconst;
       for( j=0; j<TOTST; j++)
	 DW_r[i][j] = .0;
     }

/*  Radial mesh */
     rmesh( r, drdi, nrmax, opt.rmin, opt.rmax, opt.mesh);
     for(i=0; i<nrmax; i++) r2drdi[i] = r[i]*r[i]*drdi[i];
    
      if( opt.nucl) {
      /* nuclear radius is calculated with this empirical
	 formula. it gives the compton proton radius if z=1 .
	 2.6 is to fixed such to have Uranium atomic mass.
	 the last factor converts fermi to a.u.               */
	rnuc = 0.96 * pow( 2.6*nat, 0.33333) * 1.8897268e-5;
	printf("Estimated nuclear radius: %lg (a.u.)\n\n", rnuc);
	if( rnuc <= r[3]) {
	  fprintf(stderr,
		  "Nuclear radius falls inside the first three\n"
		  "radial points. Please decrease the value of RMIN\n"
		  "in the initialization file or switch off the NUCLEUS\n"
		  "flag. Exiting.\n");
	  exit(1);
	}
      } else rnuc = 0.0;

/* Potential start */
     for( i=0; i<nrmax; i++)
     {
       vt[i] = -2.0*nat/r[i];
       bt[i] = opt.bconst;
       if(opt.nucl && r[i]<rnuc)
	 vt[i] = nat/rnuc*(r[i]*r[i]/(rnuc*rnuc)-3.0); /* already Ryd */
       vnuc[i] = vt[i];
     }
     if( opt.nspin)
       for( i=0; i<nrmax; i++)
	 bt[i] = .0;

     for( i=0; i<TOTST; i++)
       RLoop(j)
	 vsic[i][j] = bsic[i][j] = Vxcsic[i][j] = 0.0;
     
     enrg = fopen(enrg_file, "w");

     /* initialize energy eigenvals to the rel case */
     for( i=0;i<TOTST;i++)
       ectab[i] = Hec( nat, st[i].n, st[i].k);

    totiter = 0;
    pmix = opt.p;
    te_old=0;
    do { /* selfconsistent cycle */
      for( i=0; i<nrmax; i++) {
        vto[i]= vt[i];
        bto[i]= bt[i];
      }

      core( vt, bt, (double) nat,
	    r, r2drdi, drdi,
	    onp_real, totiter,
	    rhochr, rhospn, grad,
	    orbdens, orbspin, orbgrad,
	    &rhoconvrg,
	    ectab, befftab, SOtab,
	    enrg,
	    vsic, bsic,
	    GCK, FCK
	  );

      fflush(enrg);

      potat( r, drdi, (double) nat,
	     rhochr, rhospn,
	     vt, bt, vnuc, Vxc,
             v0, &q,
	     &ev, &evxc, &eexc
	   );

      for( i=0; i<nrmax; i++)
      {
        vt[i] = vt[i]*pmix + vto[i]*(1.-pmix);
        bt[i] = bt[i]*pmix + bto[i]*(1.-pmix);
      }

      esic = 0.;
      /************************************************/
      /*		  SIC			      */
      /************************************************/
      if( opt.usesic) {
	for( i=0; i<TOTST; i++)	{
	  if( onp[i]) {
   	    if( opt.usesic==1) {
	      /************************************************/
	      /*          Perdew-Zunger-81 SIC		      */
  	      potat( r, drdi, 0.0,
  		     orbdens[i], orbdens[i],  
  		     vsic[i], bsic[i], vnuc, Vxcsic[i],
  		     v0, &q,  
  		     &evs, &evxcs, &eexcs
  		   );  	      
	      /************************************************/
	    } else if( opt.usesic==2) {
  	      potat( r, drdi, 0.0,
  		     orbdens[i], orbspin[i],
  		     vsic[i], bsic[i], vnuc, Vxcsic[i],
  		     v0, &q,
  		     &evs, &evxcs, &eexcs
  		   );
   	    }
	    esic += (eexcs+evxcs-evs)*onp[i]; /* not onp_real */
	  }
	}
      }
      /************************************************/
      /*		END SIC			      */
      /************************************************/

      for( tebs=i=0; i<TOTST; i++)
	if( onp[i]) tebs += onp_real[i]*ectab[i] ;
      teeexc = eexc;
      teevxc = evxc;
      tesic  = esic;
      teee   = - ev;
      te     = te0 = tebs + teeexc + teevxc - tesic + teee;

      enconvrg = te-te_old;

      puts("-------------------------------------------------------------------------------");
      printf("Iter: %03d  Total Enrg (Se+Exc+Evxc+Ee-Esic):%14.8lf Ry = %14.8lf eV\n",
	      totiter, te, te*Ryd);
      printf("%-12s\t %-12s\t %-12s\t %-12s\t %-12s\t \n",
	     "Se", "Exc", "Evxc", "Ee", "Esic (Ryd)");
      printf("%-12.4lf\t %-12.4lf\t %-12.4lf\t %-12.4lf\t %-12.4lf\n", 
	     tebs, eexc, evxc, teee, esic);
      if( opt.nspin==0) {
	for( i=0; i<nrmax; i++)  integrand[i] = rhospn[i]*r2drdi[i];
	printf("\t\t\t\t\tTotal <Beta Sigma_z> = %lf\n",
	       simpson( integrand, nrmax, 1.));
      }
      printf("Convergence -> Density: %8.2le  Total Energy: %8.2le\n",
	     rhoconvrg, enconvrg);
      puts("-------------------------------------------------------------------------------\n\n");

      totiter++;
      if( rhoconvrg<opt.accrcy && fabs(enconvrg)<opt.accrcy && !fl.printenrg) {
	fl.printenrg = TRUE;
	rhoconvrg = 1;
      }

      te_old = te;

      if( totiter>opt.maxiter) {
	printf("Maximum number of iteration reached.\n\n");
	break;
      }

    } while( rhoconvrg > opt.accrcy || fabs(enconvrg) > opt.accrcy); 

/* calculate the averages of powers of r for all orbitals */
    print_rmeans( r, r2drdi, onp, onp_real, orbdens, &rmean);

    if( fl.printrho) {
         fp = fopen( rho_file, "w");    
         fprintf( fp, "# R          Rhocrg       Rhospn       Grad                 R*V(R)       B\n");     
          for( i=0; i<nrmax; i++)    
    	 fprintf( fp, "%12.8lf %12.8lf %12.8lf %12.8lf %20.10lf %12.8lf\n",    
    		 r[i], rhochr[i], rhospn[i], grad[i],
    		 r[i]*vt[i], bt[i]);    
         fclose( fp);    
    }

    fp = fopen("TOTAL_SPH", "a");
      fprintf(fp,"%-4s: NVERS=%d XCHNG=%d NSPIN=%-5s SIC=%-5s \tEn=%10.4lf (eV)\t<r>=%7.4lf\n", 
	      name, opt.nvers, opt.xchng, prbool(opt.nspin), prbool(opt.usesic),
	      te*Ryd, rmean);
    fclose(fp);

    /* print the base functions */
    print_base(name, onp, printme, r, orbdens, onp_real, rhochr, vt, vsic, orbgrad);

    fclose(enrg);

    if( !opt.nosph && !opt.nocoll) {
      printf("Nothing else to do... What about a short chess match?\n");
      clerr();
      return 0;
    }

    sprintf(enrg_file, "%s.basis", name);
    enrg = fopen(enrg_file, "w");

/*  Angular mesh */
    gauleg( -1., 1., cos_theta, weights, nlmax);

    for(i=0;i<2;i++) {
      printf("Calculating the basis set: ITER %d\n",i);
      basis_xtrct( vt, bt, (double) nat, (double) zat,
		   r, r2drdi, drdi,
		   onp, oncalc,
		   orbdens, orbspin, orbgrad,
		   ectab, befftab, SOtab,
		   enrg,
		   vsic, bsic,
		   GCK, FCK
		   );
    }

    puts("\n\n");

/* Let the potential be the one generated by the actual density, 
   not the one of the previous step */
    potat( r, drdi, (double) nat,
	   rhochr, rhospn,
	   vt, bt, vnuc, Vxc,
	   v0, &q,
	   &ev, &evxc, &eexc
	  );

    fflush( enrg);
    fflush( stdout);

    printf("##################### PRECALCULATIONS #####################\n\n");

    if( fl.print_base==1)
      for( i=0; i<TOTST; i++)
      {
	if( oncalc[i] || onp[i]) 
	{
	  sprintf(orb_files, "%s.%03d", rho_file, i+1);
	  fp = fopen( orb_files, "w");
	    RLoop(j)
	      fprintf( fp, "% 12.8lg % 12.8lg\n", r[j], orbdens[i][j]);
	  fclose( fp);
	  printf("%s KS orbital density saved\n", orb_files);
	}
      }

    for( base_el=i=0; i<TOTST; i++)
      if( onp[i] && oncalc[i])  /* if occupied and in the basis */
	++base_el;		/* electrons in the basis */

    for( row=i=0; i<TOTST; i++)
      if( oncalc[i]) {
	lambda[row] = ectab[i];
	DWidx[row++] = i;
      }

/* Order the eigenvalues in ascending order */
/* and prepare starting eigenvectors        */
    for( i=0; i<row; i++) {
      hr[i][i] = lambda[i];
      for( j=0; j<i; j++) {
	hr[i][j] = hr[j][i] = 0.0;
      }
    }

    diag_matrix( hr, lambda, vr, oncalc);

  if( fl.fixoccnum)
      for( i=0; i<TOTST; i++)
	colocc[i] = onp[ idxmax(vr,i,TOTST)  ];
    /* colocc[i] = TRUE if the eigenvector in the i-th column
       corresponds to an occupied state */
/*********************************************/

    if( fl.write_avect==2) {
      sprintf( filename, "%s.AVECT.%03d", name, 0);
      print_avect( filename, COEFFCUTOFF, row, nel, DWidx, onp, vr, lambda);
    }

    magn2_nc( nel, 
	      onp, oncalc,
	      vr, cos_theta,
	      GCK, FCK, 
	      Rho, Mx, Mz);

    if( opt.force_sph) {
      angular_avrg( Rho, weights, nlmax, nrmax, integrand);
      LLoop(i)
	RLoop(j)
	  Rho[i][j] = integrand[j];
      angular_avrg( Mz, weights, nlmax, nrmax, integrand);
      LLoop(i)
	RLoop(j)
	  Mz[i][j] = integrand[j];
    }

    LLoop(i)   
      RLoop(j) {
	Rho_old[i][j] =  Rho[i][j];   
	 Mz_old[i][j] =   Mz[i][j];   
	 Mx_old[i][j] =   Mx[i][j];   
      }
    dentopot_nc( Rho, Mx, Mz, Vxc_nsph, Bx, Bz);

/*       LLoop(i)    */
/*         RLoop(j)    */
/*   	{ */
/*   	  x = r[j]*sqrt(1-cos_theta[i]*cos_theta[i]); */
/*   	  z = r[j]*cos_theta[i]; */
/* 	  temp = r[j] * r[j] * r[j] * r[j] * r[j]; */
/*   	  Bx[i][j] = 3*x*z/temp; */
/*   	  Bz[i][j] = (3*z*z-r[j]*r[j])/temp; */
/*   	} */

/*     plotfield( Bx, Bz, r, cos_theta, rmean, opt.nrrectpt, "Try_Field"); */

    if( opt.clean_field && opt.nocoll)
    {
      puts("\nI WILL REMOVE THE SOURCES FROM THE XC-FIELD\n");
      cleanb( Bx, Bz, cos_theta, weights, 
	      sources,
	      opt.maxlexp, r, drdi);
    }

/*     plotfield( Bx, Bz, r, cos_theta, rmean, opt.nrrectpt, "Try_Field_clean"); */
/*     plotscalar( sources, r, cos_theta, rmean, opt.nrrectpt, "Try_Field_sources"); */

/*     if( opt.clean_field && opt.nocoll) */
/*       cleanb( Bx, Bz, cos_theta, weights,  */
/* 	      sources, */
/* 	      opt.maxlexp, r, dx); */
/*     plotscalar( sources, r, cos_theta, rmean, opt.nrrectpt, "Try_Field_clean_sources"); */

/*      exit(8);  */

    if( opt.nosph) {
      enee_nsph = potcul_nsph( Rho, Vc_nsph, cos_theta, weights, r, vnuc, drdi, opt.maxlexp); 
      printf("El-El Non sph. En (0th step): %.12lf\n", enee_nsph);
    }

    teeexc = iexc_nc( Rho, Mx, Mz, r2drdi, weights); 
					/* Exch-Corr. Energy */

/* The nsph correction is not yet in the sum of the eigenvalues, 
   so one must add it. For the same reason we don't need a new n*Vxc */
    if( opt.nosph) teee = -ev+(enee_nsph-ev);
    te     = tebs + teeexc + teevxc - tesic + teee;

    puts("\n///////////////////////////////////////////////////////////////////////////////");
    printf("NSC Preliminary step:  ");
    if( opt.nosph && opt.nocoll)
      printf("Non-Coll Non-Sph Tot Enrg:");
    else if( !opt.nosph && opt.nocoll)
      printf("Non-Coll Tot Enrg:");
    else if( opt.nosph && !opt.nocoll)
      printf("Non-Sph Tot Enrg:");
    printf("%14.8lf Ry = %14.8lf eV\n", te, te*Ryd);
    printf("%-12s\t %-12s\t %-12s\t %-12s\t %-12s\t \n",
	   "Se", "Exc", "Evxc", "Ee", "Esic (Ryd)");
    printf("%-12.4lf\t %-12.4lf\t %-12.4lf\t %-12.4lf\t %-12.4lf\n", 
	   tebs, teeexc, teevxc, teee, tesic);
    puts("///////////////////////////////////////////////////////////////////////////////");
    fflush( stdout);

    /* eta */
    for( temp2=i=0; i<nlmax; i++) {
      RLoop(j) {
	if( (temp=Bx[i][j]*Bx[i][j] + Bz[i][j]*Bz[i][j]) ==0.0)
	  integrand[j] = 0.0;
	else
	  integrand[j] = Rho[i][j] * r2drdi[j] * fabs(Bx[i][j])/sqrt(temp);
      }
      temp2 += weights[i] * simpson(integrand, nrmax, 1.);
    }
    printf("\n\nc= %16.8lf eta0= %16.8lf Z= %3d\n\n",
	   opt.c/2.0, eta0=2*M_PI*temp2/nel, nat); 
    
    /* Average magnetization */
    temp2 = integ2d( Mx, r2drdi, weights);
    temp  = integ2d( Mz, r2drdi, weights);
    printf("c= %16.8lf Mx0= %16.8lf Mz0= %16.8lf\n\n",
	   opt.c/2.0, temp2, temp); 

  if( fl.write_deltaw >= 1) {
    sprintf(filename,"%s.LAMBDAOLD",name);
    fp = fopen( filename, "w");
    for( i=0; i<row; i++)
      fprintf(fp, "Row:%03d\t% 12.8lg\n", i,lambda[i]); 
    fclose( fp);
  }

  /* Start self consistent cycle for the non spherical, non collinear part */
  do {
    
    counter++;

    puts("\n\n\n");
    printf("###############################################################\n");
    printf("######################## NSC  Iter %03d ########################\n",
	   counter);
    printf("###############################################################\n");

    if( opt.basis_recalc) {
      totiter++;

      angular_avrg( Vxc_nsph, weights, nlmax, nrmax, Vxc0);

      if(opt.nosph) {
	angular_avrg(  Vc_nsph, weights, nlmax, nrmax,  Vc0);
	for( i=0;i<nrmax;i++) vt0[i] = Vc0[i]+Vxc0[i];
      }
      else
	for( i=0;i<nrmax;i++) vt0[i] = vt[i]-Vxc[i]+Vxc0[i];

/*       angular_avrg( Vc_nsph, weights, nlmax, nrmax, bt0); */
/*       for( i=0;i<nrmax;i++) */
/* 	perr("%lg %lg  %lg\n", */
/* 	     r[i], r[i]*(vt[i]-vnuc[i]-Vxc[i]), r[i]*(bt0[i]-vnuc[i])); */
      
      angular_avrg( Bz, weights, nlmax, nrmax, bt0);

/* the next is necessary for the states that are not in the base */
      core( vt0, bt0, (double) nat,
	    r, r2drdi, drdi,
	    onp_real, totiter,
	    rhochr, rhospn, grad,
	    orbdens, orbspin, orbgrad,
	    &temp,
	    ectab, befftab, SOtab,
	    enrg,
	    vsic, bsic,
	    GCK, FCK
	    );

      basis_xtrct( vt0, bt0, (double) nat, (double) zat,
		   r, r2drdi, drdi,
		   onp, oncalc,
		   orbdens, orbspin, orbgrad,
		   ectab, befftab, SOtab,
		   enrg,
		   vsic, bsic,
		   GCK, FCK
		   );

      puts("\n\n");

      vt_p = vt0;
      bt_p = bt0;
      Vxc_p = Vxc0;
    } else {
      vt_p = vt;
      bt_p = bt;
      Vxc_p = Vxc;
    }

    deltaw_nc( oncalc, cos_theta, weights,
	       r2drdi, vt_p, Vxc_p, bt_p,
	       Vc_nsph, Vxc_nsph, Bx, Bz,
	       GCK, FCK,
	       DW_r);
    
    if( fl.write_deltaw==2) {
      sprintf(pot_file, "%s.DELTAW.%03d", name,counter);
      print_deltaw( pot_file, row, onp, DWidx, DW_r);
    }

    if(opt.nsphsic) {
      printf("\n***************** APPLYING FULL SIC ********************\n");
      sic_nsph( onp, oncalc, ectab,
		r, drdi, r2drdi,
		vsic, bsic, Vxcsic,
		vnuc,
		GCK, FCK,
		DW_r,
		cos_theta, weights,
		Rho,Mx,Mz,
		vr,
		&tesic, &tebsic
		);
      printf("********************************************************\n\n");

    } else {

	/****** Build the hermitian matrix to diag. ******/
	for( i=0; i<row; i++)
	  for( j=0; j<row; j++)
	    hr[i][j] = DW_r[ DWidx[i] ][ DWidx[j] ];
	for( i=0; i<row; i++)
	  hr[i][i] += ectab[ DWidx[i] ];

	/*************************************************/

	for( i=0; i<row; i++)
	  lambda_old[i] = lambda[i];

	diag_matrix( hr, lambda, vr, oncalc);
	if( fl.fixoccnum)
	  for( i=0; i<TOTST; i++)
	    colocc[i] = onp[ idxmax(vr, i, TOTST) ];

	if( fl.iprint>0)
	  /* check to have found the eigenvectors, ordered by columns in vr */
	  {
	    register int ii,jj,kk,ll;
	    double aa;
	    puts("\nCheck the diagonalization procedure: the following must be the unity matrix\n");
	    for(ii=0;ii<row;ii++) {
	      for(aa=jj=0;jj<=ii;jj++) {
		for(kk=0;kk<row;kk++)
		  for(ll=0;ll<row;ll++)
		    aa += vr[kk][ii]*hr[kk][ll]*vr[ll][jj];
		printf("%8.5lf ",(ii==jj)?aa/lambda[ii]:aa);
		/* 	    if( aa!=0.0 && ii!=jj) */
		/* 	      v_error("Error: element OtHO[%d][%d] != 0 = %lg\n",ii,jj,aa); */
		/* 	    if( ii==jj && aa/lambda[ii]!=1.0 ) */
		/* 	      v_error("Error: element OtHO[%d][%d] != 1 = %lg\n",ii,jj,aa/lambda[ii]); */
	      }
	      puts("");
	    }
	    printf("Orthogonal matrix vr checked\n\n");
	  }
      }

      if( fl.write_avect==2) {
	sprintf( filename, "%s.AVECT.%03d", name,counter);
	print_avect( filename, COEFFCUTOFF, row, nel, DWidx, onp, vr, lambda);
      }

      magn2_nc( nel,
		onp, oncalc,
		vr, cos_theta,
		GCK, FCK, Rho, Mx, Mz);

    if( opt.force_sph) {
      angular_avrg( Rho, weights, nlmax, nrmax, integrand);
      LLoop(i)
	RLoop(j)
	  Rho[i][j] = integrand[j];
      angular_avrg( Mz, weights, nlmax, nrmax, integrand);
      LLoop(i)
	RLoop(j) {
	  Mz[i][j] = integrand[j];
	  Mx[i][j] = 0.0;
	}
    }

      /* Check the change in density */
      rhoconvrg=0.0;
      LLoop(i)
	RLoop(j) {   
	  if( (temp=fabs( Rho[i][j]-Rho_old[i][j])) > rhoconvrg) rhoconvrg=temp; 
	  if( (temp=fabs(  Mz[i][j]- Mz_old[i][j])) > rhoconvrg) rhoconvrg=temp; 
	   Rho_old[i][j] =  Rho[i][j];   
	  Mz_old[i][j] = Mz[i][j];   
	}
      if( opt.nocoll)
	LLoop(i) {
	  RLoop(j) {
	    if( (temp=fabs( Mx[i][j]-Mx_old[i][j]))> rhoconvrg) rhoconvrg=temp; 
	    Mx_old[i][j] = Mx[i][j];   
	  }
      }

      dentopot_nc( Rho, Mx, Mz, Vxc_nsph, Bx, Bz);

      if( opt.clean_field && opt.nocoll)
	cleanb( Bx, Bz, cos_theta, weights,
		sources,
		opt.maxlexp, r, drdi);

      if( opt.nosph) {
	enee_nsph = potcul_nsph( Rho, Vc_nsph, cos_theta, weights, r, vnuc, drdi, opt.maxlexp);
	printf("El-El Non sph. En: %14.8lf\n", enee_nsph);
      }
      
      /* Calculate non-collinear and/or non-spherical total energy */
      tebs = temp = 0.0;
      for( j=nel, i=0; i<TOTST; i++) {
	if( onp[i] && !oncalc[i]) {    /* if occupied but not in the basis */
	  tebs += onp_real[i]*ectab[i];/* sum eignvals of the spherical part */
	  --j;
	  if(opt.nsphsic) {
	    /*          Perdew-Zunger-81 SIC		      */
	    potat( r, drdi, 0.0,
		   orbdens[i], orbdens[i],  
		   vsic[i], bsic[i], vnuc, Vxcsic[i],
		   v0, &q,  
		   &evs, &evxcs, &eexcs );  	      
	    temp += (eexcs+evxcs-evs)*onp[i]; /* not onp_real */
	  }
	}
      }
      tesic += temp;
      printf("\nSum of eigenvalues of occupied states not in basis: %lg\n", tebs);
      printf("SIC contribution from occupied states not in basis: %lg\n", temp);

      if(opt.nsphsic)
	tebs += tebsic;
      else
	for( i=0; i<j; i++) {		/* sum what remains */
	  if( fl.fixoccnum)
	    if( !colocc[i] ) {
	      j++;
	      printf("Skipped column %dth in Tebs calculation\n", i); 
	      continue;
	    }
	  tebs += lambda[i];		/* Sum of the lowest eigenvalues */
	}

      /* If collinear then automatically no Mx */

      teeexc = iexc_nc( Rho, Mx, Mz, r2drdi, weights); 
      /* Exch-Corr. Energy */

      teevxc = ivexc_nc( Rho, Mx, Mz, Vxc_nsph, Bx, Bz, r2drdi, weights); 
      /* -Integ( n*vxc) */

      if( opt.nosph) teee = -enee_nsph;
      te     = tebs + teeexc + teevxc - tesic + teee;

      enconvrg = te-te_old;

      puts("\n///////////////////////////////////////////////////////////////////////////////");
      printf("NSC Iter: %03d  ", counter);
      if( opt.nosph && opt.nocoll)
	printf("Non-Coll Non-Sph Tot Enrg:");
      else if( !opt.nosph && opt.nocoll)
	printf("Non-Coll Tot Enrg:");
      else if( opt.nosph && !opt.nocoll)
	printf("Non-Sph Tot Enrg:");
      printf("%14.8lf Ry = %14.8lf eV\n", te, te*Ryd);
      printf("%-12s\t %-12s\t %-12s\t %-12s\t %-12s\t \n",
	     "Se", "Exc", "Evxc", "Ee", "Esic (Ryd)");
      printf("%-12.4lf\t %-12.4lf\t %-12.4lf\t %-12.4lf\t %-12.4lf\n", 
	     tebs, teeexc, teevxc, teee, tesic);
      printf("NSC Convergence -> Density: %8.2le  Total Energy: %8.2le\n",
	     rhoconvrg, enconvrg);
      puts("///////////////////////////////////////////////////////////////////////////////");
      fflush( stdout);

      te_old = te;

      if( counter>opt.maxiter) {
	printf("Maximum number of iteration reached.\n\n");
	break;
      }

    } while( rhoconvrg > opt.accrcy || fabs(enconvrg) > opt.accrcy );

    fclose(enrg);

  for( temp2=i=0; i<nlmax; i++)
  {
    RLoop(j)
    {
      if( (temp=Bx[i][j]*Bx[i][j] + Bz[i][j]*Bz[i][j]) ==0.0)
	integrand[j] = 0.0;
      else
	integrand[j] = Rho[i][j] * r2drdi[j] * fabs(Bx[i][j])/sqrt(temp);
    }
    temp2 += weights[i] * simpson(integrand, nrmax, 1.);
  }
  printf("\nC=%.8lf <eta>=%.8lf Z=%d\n\n",
	 opt.c/2.0, eta=2*M_PI*temp2/nel, nat); 
 
  sprintf(filename,"%s.Integ_Rho",name);
  if( fl.printrho) fp = fopen( filename, "w");
  RLoop(j) {
    for( integrand[j]=i=0; i<nlmax; i++)
      integrand[j] += 2.0*M_PI*(weights[i] * Rho[i][j]);
    if( fl.printrho) fprintf( fp, "%lg %lg\n", r[j], integrand[j]);
    integrand[j] *= r2drdi[j];
  }
  if( fl.printrho) fclose( fp); 
  temp2 = simpson(integrand, nrmax, 1.);
  printf("Z=%d    <Rho>=%.12lf ", nat, temp2); 

    /* Average magnetization */
    temp2 = integ2d( Mx, r2drdi, weights);
    temp  = integ2d( Mz, r2drdi, weights);  
    printf("<Mx>=%.8lf <Mz>=%.8lf\n\n", temp2, temp); 

  if( fl.write_deltaw >= 1) {
    sprintf(filename,"%s.LAMBDA",name);
    fp = fopen( filename, "w");
    for( i=0; i<row; i++)
      fprintf(fp, "Row:%03d\t% 12.8lg\n", i,lambda[i]); 
    fclose( fp);
  }

  if( fl.write_avect >= 1) {
    sprintf(filename,"%s.AVECT",name);
    print_avect( filename, COEFFCUTOFF, row, nel, DWidx, onp, vr, lambda);
  }

  if( fl.write_deltaw>=1) {
    sprintf(pot_file, "%s.DELTAW",name);
    print_deltaw( pot_file, row, onp, DWidx, DW_r);
  }

  if( fl.printpot) {
    sprintf(pot_file, "%s.POTXC",name);
    fp = fopen( pot_file, "w");  
    for( j=0; j<nlmax; j++) 
      for( i=0; i<nrmax; i++) 
 	fprintf(fp, "% 12.8lg\n", Vxc_nsph[j][i]);  
    fclose( fp);  
    
    sprintf(pot_file, "%s.BX",name);
    fp = fopen( pot_file, "w"); 
    for( j=0; j<nlmax; j++) 
      for( i=0; i<nrmax; i++) 
 	fprintf(fp, "% 12.8lg\n", Bx[j][i]);  
    fclose( fp); 
    
    sprintf(pot_file, "%s.BZ",name);
    fp = fopen( pot_file, "w"); 
    for( j=0; j<nlmax; j++) 
      for( i=0; i<nrmax; i++) 
 	fprintf(fp, "% 12.8lg\n", Bz[j][i]);  
    fclose( fp); 
  }
  
  if( fl.printmagn) {
    fp = fopen( magn_file, "w");
    for( j=0; j<nlmax; j++) {
      sin_theta = sqrt( 1-cos_theta[j]*cos_theta[j]);
      for( i=0; i<nrmax; i++) {
	temp = r[i];
	x = temp*sin_theta;
	z = temp*cos_theta[j];
	fprintf(fp,
		"% 12.8lg % 12.8lg % 12.8lg % 12.8lg % 12.8lg % 12.8lg\n",
 	        r[i], x,z, Rho[j][i], Mx[j][i], Mz[j][i]); 
      }
    }  
    fclose( fp);
  }

  fp = fopen( "ZCE0Ee0e", "a");
    fprintf(fp, "%s.%1d%1d%1d %3d % 16.8lf % 16.8lf % 16.8lf % 16.8lf % 16.8lf\n",
	        name, opt.nocoll, opt.nosph, opt.usesic, nat, opt.c/2.0, te0, te, eta0, eta); 
  fclose( fp);

  for( temp=i=0; i<nlmax; i++)
  {
    RLoop(j)
    {
      integrand[j] = r2drdi[j] * Rho[i][j] * sqrt(Bz[i][j]*Bz[i][j]+Bx[i][j]*Bx[i][j]);
    }
    temp += weights[i] * simpson(integrand, nrmax, 1.);
  }
  printf("<(|B|)> = %16.8lf \n\n", 2*M_PI*temp/nel); 

  if( fl.compton_profile)
  { 
    printf("Sorry, no compton profile possible, for now.\n");
    fprintf(stderr,"Sorry, no compton profile possible, for now.\n");
/*     fp = fopen("Compton", "w"); */
/*     LLoop(i) */
/*     { */
/*      if( cos_theta[i]>0) */
/*      { */
/*       fprintf( fp, "# costheta=%lf\n", cos_theta[i]); */
/*       for( temp=0; temp<20; temp+=.2) */
/*       { */
/* 	RLoop(j) */
/* 	{ */
/* 	  integrand[j] = drdi[j] * Mz[i][j] * cos(r[j]*temp); */
/* 	} */
/* 	fprintf( fp, "%lf\t%lf\n", temp, simpson(integrand, nrmax, 1.)); */
/*       } */
/*      } */
/*       fprintf( fp, "\n"); */
/*     } */
/*     fclose( fp); */
  }


  if( opt.nrrectpt >0)
  {
/*     dentopot_nc( opt.c, opt.nvers, Rho, Mx, Mz, V, Bx, Bz); */
/*     plotfield( Bx, Bz, r, cos_theta, rmean, opt.nrrectpt, "Rect_B_before"); */
/*     cleanb( Bx, Bz, cos_theta, weights, */
/* 	    sources, */
/* 	    opt.maxlexp, r, dx); */
    plotfield( Bx, Bz, r, cos_theta, rmean, opt.nrrectpt, "Rect_B");
    plotfield( Mx, Mz, r, cos_theta, rmean, opt.nrrectpt, "Rect_magn");
    /* Just clean the field to find the sources !!!!!! */
    /* WARNING the field is being changed */
	cleanb( Bx, Bz, cos_theta, weights,
		sources,
		opt.maxlexp, r, drdi);

    plotscalar( sources, r, cos_theta, rmean, opt.nrrectpt, "Rect_divB");
    plotscalar( Rho, r, cos_theta, rmean, opt.nrrectpt, "Rect_Rho");
  }

    if( opt.nosph && opt.usesic && fl.nsphorbsic)
    {
      v_error("NSPHORBSIC not possible for now\n");

      fp = fopen("Nsph_orb.txt", "w");

      /* Select the orbital to process */
      nsphstart = 0; nsphend = 279;
      for( sum = 0, nsphstate=nsphstart; nsphstate<=nsphend; nsphstate++)
      {
	if( !onp[nsphstate]) continue;
	for( i=0; i<row; i++)
	  if( DWidx[ maxindex(vr, onp, DWidx, i, row) ] == nsphstate)
	    colocc[i] = onp[ DWidx[ maxindex(vr, onp, DWidx, i, row) ] ];
	  else colocc[i]=0;
	/* Calculate the spin-density only for this orbital */
	magn2_nc( 1,
		  onp, oncalc,
		  vr, cos_theta,
		  GCK, FCK, Rho, Mx, Mz);
	  
	/* the Coulomb non-sph. self-interaction */
	printf("----------------------------------------------------\n");
	printf("Selected orbital #%d\n", nsphstate+1);
	enee_nsph = potcul_nsph( Rho, Vc_nsph, cos_theta, weights, r, vnuc, drdi, opt.maxlexp);
	/* and the exch-corr. energy */
	teeexc = iexc_nc( Rho, Mx, Rho, r2drdi, weights); /* Joerg SIC */ 
	fprintf(fp, "----------------------------------------------------\n");
	fprintf(fp, "Selected orbital #%d\n", nsphstate+1);
	fprintf(fp, "Coul. Non sph. self-En for selected orb.: %14.8lf\n", enee_nsph);
	fprintf(fp, "XC Non sph. self-En for selected orb.: %14.8lf\n", teeexc);
	/* And the sph self-en */
	potat( r, drdi, 0.0,
	       orbdens[nsphstate], orbdens[nsphstate],  /* Joerg SIC */
	       vsic[nsphstate], bsic[nsphstate], vnuc, Vxcsic[nsphstate],  
	       v0, &q,  
	       &evs, &evxcs, &eexcs
	     );  	      
	fprintf(fp, "Coul. sph. self-En for selected orb.: %14.8lf\n", evs);
	fprintf(fp, "XC sph. self-En for selected orb.: %14.8lf\n\n", eexcs);
	fprintf(fp, "Delta(nsph-sph) for selected orb.: %14.8lf\n\n",
		enee_nsph-evs+teeexc-eexcs);

	sum += enee_nsph-evs+teeexc-eexcs;
      }
      fprintf(fp, "++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      fprintf(fp, "++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      fprintf(fp, "Sum of Delta(nsph-sph) for all considered orbitals: %14.8lf\n\n",
	      sum);
      fclose( fp);
    }
    
    clerr();

  return 0;
}

static int maxindex( double vr[TOTST][TOTST], int DWidx[TOTST], 
		     int onp[TOTST], int m, int row)
{
  register int i, idx;
  double temp, max;

  for( max=i=0; i<row; i++)
  {
    temp = vr[i][m];
    temp *= temp;
    if( temp > max)
      { max = temp; idx = i; }
  }
  if( max < .7 && onp[DWidx[idx]])
    printf("Indecision in maxindex::main max=%lf, STATE#%d\n", max, idx);

  return idx;
}

#undef prbool
