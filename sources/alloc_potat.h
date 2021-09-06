/* static void alloc_potat(Dbp ur[2], Dbp bgx[2], */
/* 			Dbp y, Dbp rhor2, Dbp vc,  */
/* 			Dbp aa, Dbp bb) */
{
register i;

  y=rhor2=vc=aa=bb=NULL;

  y = vector_alloc(nrmax);
  rhor2 = vector_alloc(nrmax);
  vc = vector_alloc(nrmax);
  aa = vector_alloc(nrmax);
  bb = vector_alloc(nrmax);

  for(i=0;i<2;i++) {
    ur[i] = vector_alloc(nrmax);
    bgx[i] = vector_alloc(nrmax);
  }
}

