/********************************************************************
 *                                                                  *
 *                      AINV = A^(-1)                               *
 *                                                                  *
 *  Invert A using the GAUSS-JORDAN - algorithm                     *
 *  the 1- matrix is not set up and use is made of its structure    *
 *                                                                  *
 *                    (double) Version                              *
 *                                                                  *
 ********************************************************************
 *		ainv and pa MUST be created with malloc		    *
 ********************************************************************/

void rinvgj( double ainv[4][4], double pa[4][4], int n)
{
  int icol, l, ll;
  double t1, t;
  double a[4][4];

/*						make local copy of a   */
     for( l=0; l<n; l++)
       for( ll=0; ll<n; ll++)
	 a[l][ll] = pa[l][ll];

/*                                                      scan columns   */
     for( icol=0; icol<n; icol++)
     {
/*                                               make A(ICOL,ICOL) = 1 */
       t1 = 1.0/a[icol][icol];
       for( l=icol+1; l<n; l++)
         a[icol][l] *= t1;

       for( l=0; l<icol; l++)
         ainv[icol][l] *= t1;
       ainv[icol][icol] = t1;

/*                                    make A(LL,ICOL) = 0 for LL<>ICOL */
       for( ll=0; ll<n; ll++)
       {
         if( ll != icol)
         {
           t = a[ll][icol];
           for( l=icol+1; l<n; l++)
             a[ll][l] -= a[icol][l]*t;

           for( l=0; l<icol; l++)
             ainv[ll][l] -= ainv[icol][l]*t;
           ainv[ll][icol] = -t1*t;
         }
       }
     }
}
