/* static void alloc_core(Dbp rhochrcur, Dbp rhospncur, Dbp integrand, */
/* 		       Dbp gradcur, Dbp rho_old, Dbp b00, Dbp vv, Dbp bb, */
/* 		       Dbp gck[2][2], Dbp fck[2][2], */
/* 		       Dbp  gc[2][2],  Dbp fc[2][2], */
/* 		       Dbp  dp[2][2],  Dbp dq[2][2], */
/* 		       Dbp  wp[2][2],  Dbp wq[2][2], */
/* 		       Dbp dpk[2][2], Dbp dqk[2][2]) */
{
register int ia,ja;

  rhochrcur = vector_alloc(nrmax);
  rhospncur = vector_alloc(nrmax);
  integrand = vector_alloc(nrmax);
  gradcur   = vector_alloc(nrmax);
  rho_old   = vector_alloc(nrmax);
  b00       = vector_alloc(nrmax);
  vv        = vector_alloc(nrmax);
  bb        = vector_alloc(nrmax);
  dovrc     = vector_alloc(nrmax);

  for(ia=0;ia<2;ia++) {
    for(ja=0;ja<2;ja++) {
      gck[ia][ja] = vector_alloc(nrmax);
      fck[ia][ja] = vector_alloc(nrmax);
       gc[ia][ja] = vector_alloc(nrmax);
       fc[ia][ja] = vector_alloc(nrmax);
       dp[ia][ja] = vector_alloc(nrmax);
       dq[ia][ja] = vector_alloc(nrmax);
       wp[ia][ja] = vector_alloc(nrmax);
       wq[ia][ja] = vector_alloc(nrmax);
      dpk[ia][ja] = vector_alloc(nrmax);
      dqk[ia][ja] = vector_alloc(nrmax);
    }
  }
}
