/* COREERR:
       Calculates the mismatch of the radial wave functions at
       the point nmatch for out- and inward integration
*/

void coreerr( double err[4], double var[4], int s, int nsol,
              double Pow[2][2], double Qow[2][2], 
              double Piw[2][2], double Qiw[2][2])
{
  int t;

      err[0] = Pow[s][s] - Piw[s][s]*var[1];      
      err[1] = Qow[s][s] - Qiw[s][s]*var[1];
      
      if (nsol==1) return;
      
      t = 1-s;
      
      err[0] += Pow[s][t]*var[2] - Piw[s][t]*var[1]*var[3];
      err[1] += Qow[s][t]*var[2] - Qiw[s][t]*var[1]*var[3];
      err[2]  = Pow[t][s]        - Piw[t][s]*var[1] 
               +Pow[t][t]*var[2] - Piw[t][t]*var[1]*var[3];
      err[3]  = Qow[t][s]        - Qiw[t][s]*var[1] 
               +Qow[t][t]*var[2] - Qiw[t][t]*var[1]*var[3];

  return;
}

 

  
