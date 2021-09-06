/* Integration routine using Bode method */
/* VDPS, March 2002 after magru's suggestion */

#include <stdio.h>
#include <stdlib.h>

double bodeint( double *f, int n, double h, char way, double *intf)
/* 
   f   : vector array to be integrated
   n   : its dimension
   h   : step
   way : 'o' for outward integration, 'i' for inward one
   intf: the primitive of f; 
         if outward, intf[i] = INTEG[f[x], {x,0,i}]
         if  inward, intf[i] = INTEG[f[x], {x,i,n}]
   bodeint: the value of INTEG[f[x], {x,0,n}]
*/
{
  register i;
  int k;
  double tmp;
  double f0,f1,f2,f3,f4, l0,l1,l2,l3,l4;
  double h720, h90, h380, h245;

  h720 =     h/720.0;
  h90  =     h/90.0;
  h380 = 3.0*h/80.0;
  h245 = 2.0*h/45.0;

  k = (n-1)%4;

  if( way == 'o' ) {

    f0 = f[0];
    f1 = f[1];
    f2 = f[2];
    f3 = f[3];
    f4 = f[4];

    l4 = f[n-1];
    l3 = f[n-2];
    l2 = f[n-3];
    l1 = f[n-4];
    l0 = f[n-5];

    intf[0] = 0.0;
    
    for( i=4; i<n-k; i+=4) {
      tmp = intf[i-4];

      intf[i-3] = tmp +
	          (251.0*f0 + 646.0*f1 - 264.0*f2 + 106.0*f3 - 19.0*f4)*h720;
      intf[i-2] = tmp +
	          (29.0*f0 + 124.0*f1 + 24.0*f2 + 4.0*f3 - f4)*h90;
      intf[i-1] = tmp +
	          (9.0*f0 + 34.0*f1 + 24.0*f2 + 14.0*f3 - f4)*h380;
      intf[ i ] = tmp +
	          (7.0*f0 + 32.0*f1 + 12.0*f2 + 32.0*f3 + 7.0*f4)*h245;
      f0 = f4;
      f1 = f[i+1];
      f2 = f[i+2];
      f3 = f[i+3];
      f4 = f[i+4];
    }
    switch(k) { /* note the absence of breaks */
    case 3:
      intf[n-3] = intf[n-5] + (29.0*l0 + 124.0*l1 + 24.0*l2 + 4.0*l3 - l4)*h90;
    case 2:
      intf[n-2] = intf[n-5] + (9.0*l0 + 34.0*l1 + 24.0*l2 + 14.0*l3 - l4)*h380;
    case 1: 
      intf[n-1] = intf[n-5] + (7.0*l0 + 32.0*l1 + 12.0*l2 + 32.0*l3 + 7.0*l4)*h245;
    default: ;
    }

    return intf[n-1];
    
  } else if( way == 'i' ) {
    f0 = f[n-k-1];
    f1 = f[n-k-2];
    f2 = f[n-k-3];
    f3 = f[n-k-4];
    f4 = f[n-k-5];
    
    l0 = f[n-1];
    l1 = f[n-2];
    l2 = f[n-3];
    l3 = f[n-4];
    l4 = f[n-5];

    intf[n-1] = intf[n-k-1] = 0.0;
    
    switch(k) {
    case 3:
      intf[n-4] = intf[n-1] + (9.0*l0 + 34.0*l1 + 24.0*l2 + 14.0*l3 - l4)*h380;
    case 2:
      intf[n-3] = intf[n-1] + (29.0*l0 + 124.0*l1 + 24.0*l2 + 4.0*l3 - l4)*h90;
    case 1: 
      intf[n-2] = intf[n-1] + (251.0*l0 + 646.0*l1 - 264.0*l2 + 106.0*l3 - 19.0*l4)*h720;
    default: ;
    }
    
    for( i=n-k-5; i>=0; i-=4) {
      tmp = intf[i+4];

      intf[ i ] = tmp +
        (7.0*f0 + 32.0*f1 + 12.0*f2 + 32.0*f3 + 7.0*f4)*h245;
      intf[i+1] = tmp +
        (9.0*f0 + 34.0*f1 + 24.0*f2 + 14.0*f3 - f4)*h380;
      intf[i+2] = tmp +
        (29.0*f0 + 124.0*f1 + 24.0*f2 + 4.0*f3 - f4)*h90;
      intf[i+3] = tmp +
        (251.0*f0 + 646.0*f1 - 264.0*f2 + 106.0*f3 - 19.0*f4)*h720;

      f0 = f4;
      f1 = f[i-1];
      f2 = f[i-2];
      f3 = f[i-3];
      f4 = f[i-4];
    }
      
    return intf[0];

  } else {
    fprintf(stderr, "Direction of integration unknown in bodeint\n");
    exit(1);
  }
  return 1e99;
}
