/* static void alloc_deltaw(Dbp psij[4], Dbp psik[4], Dbp rinteg_r) */
{
  register int i;

  rinteg_r = vector_alloc(nrmax);
  for(i=0;i<4;i++) {
    psij[i] = vector_alloc(nrmax);
    psik[i] = vector_alloc(nrmax);
  }
}
