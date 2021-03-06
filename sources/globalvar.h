/**************************************************************************/
/*			Global variables				  */

const int nqntab[] = { 1, 2,2, 3,3,3, 4,4,4,4, 5,5,5,5,5, 6,6,6,6,6,6, 7,7,7,7,7,7,7 };
const int lqntab[] = { 0, 0,1, 0,1,2, 0,1,2,3, 0,1,2,3,4, 0,1,2,3,4,5, 0,1,2,3,4,5,6 };
const char *orbname[] = 
{"1s","2s","2p","3s","3p","3d","4s","4p","4d","4f","5s","5p","5d","5f",
 "5g","6s","6p","6d","6f","6g","6h","7s","7p","7d","7f","7g","7h","7i"}; 

const int ord[]=
  { 0, 1, 2, 3, 4, 6, 5, 7, 10, 8, 11, 15, 9, 12, 16, 21, 13, 17, 22,
/* 1s,2s,2p,3s,3p,4s,3d,4p, 5s,4d, 5p, 6s,4f, 5d, 6p, 7s, 5f, 6d, 7p, */
   14, 18, 19, 20, 23, 24, 25, 26, 27};
/* 5g, 6f, 6g, 6h, 7d, 7f, 7g, 7h, 7i */

struct Options opt;
struct Flags fl;
struct State_qn st[TOTST];

int nrmax, nlmax;
FILE *fperr;

/*			End global variables				  */
/**************************************************************************/
