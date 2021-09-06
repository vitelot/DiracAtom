/* static void free_potat(Dbp ur[2], Dbp bgx[2], */
/* 			Dbp y, Dbp rhor2, Dbp vc,  */
/* 			Dbp aa, Dbp bb) */
{
register i;

  vector_free(y);
  vector_free(rhor2);
  vector_free(vc);
  vector_free(aa);
  vector_free(bb);

  for(i=0;i<2;i++) {
    vector_free(ur[i]);
    vector_free(bgx[i]);
  }
}
