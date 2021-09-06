/* void alloc1d( Dbp vt, Dbp bt, Dbp r, Dbp drdi, Dbp r2drdi, */
/* 	      Dbp vnuc, Dbp rhochr, Dbp rhospn, Dbp grad, Dbp vto, */
/* 	      Dbp integrand, Dbp bto, Dbp Vxc, Dbp Vxcs, */
/* 	      double *cos_theta, Dbp weights, */
/* 	      Dbp GCK[TOTST][2], Dbp FCK[TOTST][2], */
/* 	      Dbp orbdens[TOTST], Dbp orbspin[TOTST], Dbp orbgrad[TOTST], */
/* 	      Dbp vsic[TOTST], Dbp bsic[TOTST]) */
{
#define  R_Alloc vector_alloc(nrmax)
#define LR_Alloc matrix_alloc(nlmax, nrmax)

register int ia,ja;

  puts("Allocating one dimensional arrays for main"); 
  /* Allocate one dimensional arrays */
  vt		= R_Alloc;
  bt		= R_Alloc;
  r		= R_Alloc;
  drdi		= R_Alloc;
  r2drdi	= R_Alloc;
  vnuc		= R_Alloc;
  rhochr	= R_Alloc;
  rhospn	= R_Alloc;
  grad		= R_Alloc;
  vto		= R_Alloc;
  integrand	= R_Alloc;
  bto		= R_Alloc;
  Vxc		= R_Alloc;
  Vxc0		= R_Alloc;
  Vc0		= R_Alloc;
  vt0		= R_Alloc;
  bt0		= R_Alloc;
  cos_theta	= vector_alloc(nlmax);
  weights	= vector_alloc(nlmax);
  /***********************************/

  for(ia=0; ia<TOTST; ia++)
    for(ja=0; ja<2; ja++) {
      GCK[ia][ja] = R_Alloc;
      FCK[ia][ja] = R_Alloc;
    }

  for(ia=0; ia<TOTST; ia++) {
    orbdens[ia] = R_Alloc;
    orbspin[ia] = R_Alloc;
    orbgrad[ia] = R_Alloc;
    vsic[ia]    = R_Alloc;
    bsic[ia]    = R_Alloc;
    Vxcsic[ia]  = R_Alloc;
  }

}

/* void alloc2d(Grid2d Rho, Grid2d V, Grid2d V_old, */
/* 	     Grid2d Vc_nsph, Grid2d Vc_nsph_old, */
/* 	     Grid2d Mx, Grid2d Mz, */
/* 	     Grid2d Bx, Grid2d Bz, Grid2d Bx_old, Grid2d Bz_old, */
/* 	     Grid2d sources) */
{
  /* Allocate two dimensional arrays */

  puts("Allocating two dimensional arrays for main"); 

  Rho      = LR_Alloc;
  Vxc_nsph = LR_Alloc;
  Vc_nsph  = LR_Alloc;
  Mx       = LR_Alloc;
  Mz       = LR_Alloc;
  Bx       = LR_Alloc;
  Bz       = LR_Alloc;
  sources  = LR_Alloc;
  Rho_old  = LR_Alloc;
  Mx_old   = LR_Alloc;
  Mz_old   = LR_Alloc;

  puts("Allocation of arrays for main finished\n"); 
}

#undef R_Alloc
#undef LR_Alloc
