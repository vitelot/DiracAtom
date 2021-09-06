/* static void free_core( Dbp rhochrcur, Dbp rhospncur, Dbp integrand, */
/* 		       Dbp gradcur, Dbp rho_old, Dbp b00, Dbp vv, Dbp bb, */
/* 		       Dbp gck[2][2], Dbp fck[2][2], */
/* 		       Dbp  gc[2][2],  Dbp fc[2][2], */
/* 		       Dbp  dp[2][2],  Dbp dq[2][2], */
/* 		       Dbp  wp[2][2],  Dbp wq[2][2], */
/* 		       Dbp dpk[2][2], Dbp dqk[2][2]) */
{
register int ia,ja;
  free ((Dbp) rhochrcur);
  free ((Dbp) rhospncur);
  free ((Dbp) integrand);
  free ((Dbp) gradcur);
  free ((Dbp) rho_old);
  free ((Dbp) b00);
  free ((Dbp) vv);
  free ((Dbp) bb);
  free ((Dbp) dovrc);

  for(ia=0;ia<2;ia++) {
    for(ja=0;ja<2;ja++) {
      free ((Dbp) gck[ia][ja]);
      free ((Dbp) fck[ia][ja]);
      free ((Dbp)  gc[ia][ja]);
      free ((Dbp)  fc[ia][ja]);
      free ((Dbp)  dp[ia][ja]);
      free ((Dbp)  dq[ia][ja]);
      free ((Dbp)  wp[ia][ja]);
      free ((Dbp)  wq[ia][ja]);
      free ((Dbp) dpk[ia][ja]);
      free ((Dbp) dqk[ia][ja]);
    }
  }
}
