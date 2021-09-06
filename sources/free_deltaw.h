/* static void free_deltaw(Dbp psij[4], Dbp psik[4], Dbp rinteg_r) */
{
  register int i;

  vector_free(rinteg_r);
  for(i=0;i<4;i++) {
    vector_free(psij[i]);
    vector_free(psik[i]);
  }
}

