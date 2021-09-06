#include "relat.h"

void prepare_qn_tables(void)
{
  int ilshell, imj, ic, l, n, kap[2], s, nsol;

  /*   printf("State list:\n"); */
  for( ic=ilshell=0; ilshell<TOTSH; ilshell++) {
    n = nqntab[ ord[ilshell] ];
    l = lqntab[ ord[ilshell] ];

    for( imj= -l; imj<=l+1; imj++) {
      if( imj == -l || imj == l+1) nsol = 1;
      else nsol = 2;
      kap[0] = -l-1;
      kap[1] =   l ;

      for( s=0; s<nsol; s++) {
	st[ic].n  = n;
	st[ic].l  = l;
	st[ic].k  = kap[s];
	st[ic].imj = imj;

	st[ic].name[0] = n + '0';
	st[ic].name[1] = l_to_c(l);
	st[ic].name[3] = '\0';
	/* 	  printf("IC=%03d N=%d L=%d K=%+d Jz=%+d/2\n", */
	/* 		 ic,n,l,kap[s],2*imj -1); */
	ic++;
      }			/* End 'for s'		*/
    }			/* End 'for imj'	*/
  }			/* End 'for ilshell'	*/
  /*   puts(""); */
}
