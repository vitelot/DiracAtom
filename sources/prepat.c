#include "relat.h"

/* PREPAT:
	  Reads electronic structure from file and converts
	  it into occupation numbers.
			VDPS
*/

#define nat (*_nat)

static void fillin(const int ord[TOTSH], const int ltb[TOTSH], int on[TOTSH], int ntot);
static char Tolower( char);
static char Toupper( char);
static void expandsq( char *, char *, int *);
static void create_el_str(void);

int prepat( int onp[TOTST], int on[TOTSH], char name[128],
	    int *_nat, char atom_name[128])
/**********************************************************************
 OUTPUT:
  *on	   Occupation numbers in the not relativistic case
  *onp	   Occupation numbers in the relativistic case
  nat      Atomic number

 VALUE RETURNED:
  Atomic valency ( 0 for neutral atoms)
 *********************************************************************/
{
  FILE *in;
  char line[128], es[1024]="", lerr[1024], insq[128], *p;
  int nl=0,				/* Lines read */
       i, j, n, l, ll, totocc, zat, occ, ocu, default_nat;
  int ix[] ={ 0, 1, 3, 6, 10, 15, 21};   /* Index to nqntbl */
  const int *ltb=lqntab;

     name[0] = '\0';
     nat = 0;
     for( i=0; i<TOTSH; i++)
       on[i] = 0;
     for( i=0; i<TOTST; i++)
       onp[i] = 0;

     if( atom_name[0] == '\0') { /* if no atom in command line... */
       in = fopen("El_structure", "r");
       if (!in) {
	 printf("Error opening El_structure file in PREPAT\n");
	 printf("Trying to create a default one.\n\n");
	 create_el_str();
       }
     
       i=0;
       while ( !feof(in) ) {
	 if ( fgets( line, 128, in) == NULL ) break;
	 nl++;
	 p = line;
	 if ( *p != '#' ) {        /* If it isn't a comment ... */
	   if ( sscanf( line, "Name %s\n", name) == 1)
	     continue;
	   if ( sscanf( line, "Atomic number %d\n", &nat) == 1)
	     continue;
	   while( *p != ';') {
	     if( *p == ' ' || *p == '\t') p++;  /* Ignore space chars */
	     else if( *p=='\n' || *p=='\0') break;
	     else {
	       lerr[i] = ' ';
	       es[i++] = *p++;
	     }
	   }
	 }
	 if( *p==';') {  /* End of structure reached */
	   es[i++] = *p; 
	   es[i] = lerr[i] = '\0';
	   break;  
	 }
       }
       fclose(in);
     
       if( name[0] == '\0') {
	 printf("Name not read from El_structure file. Assuming \"X\" as default\n\n");
	 strcpy( name, "X");
       }

       if( es[i-1] != ';')
	 v_error("Error: no semicolon found in El_structure\n\n");
     } else {
       sprintf( es, "[%s];", atom_name);
       sprintf( name, "%s", atom_name);
     }

     printf("Electronic Structure: %s\n", es); 

     for( i=0,p=es; *p!=';'; p++)
       if( *p=='[') {
	 ++p;
	 for( ; *p!=']'; p++) {
	   if( *p==';') {
	     lerr[ p-es ] = '^';
	     v_error("%s\nSemicolon inside square brackets in El_structure\n",lerr);
	   }
	   insq[i++] = Tolower(*p);
	 }
	 if( i==0) {
	   v_error("Empty square brackets in El_structure. Exit.\n");
	 }
       }
     if( i) {	/* If square brackets are found ... */
       insq[i]='\0';
       expandsq( insq, es, &default_nat);
     }

     if( nat==0) {
       printf("Atomic number not read in El_structure file, assuming %d as default\n\n", default_nat);
       nat = default_nat;
     }
     printf("Expanded El_structure:\n%s Atomic number:%d Name:%s\n",
	    es, nat, name);

     for( p=es; *p!=';'; p++) {
       if     ( *p == '[' ) {
	 ++p;
	 if      ( *p=='H' || *p=='h') fillin(ord,ltb, on, 2); /* Helium  */
	 else if ( *p=='N' || *p=='n') fillin(ord,ltb, on,10); /* Neon    */
	 else if ( *p=='A' || *p=='a') fillin(ord,ltb, on,18); /* Argon   */
	 else if ( *p=='K' || *p=='k') fillin(ord,ltb, on,36); /* Kripton */
	 else if ( *p=='X' || *p=='x') fillin(ord,ltb, on,54); /* Xenon   */
	 else if ( *p=='R' || *p=='r') fillin(ord,ltb, on,86); /* Radon   */
	 else if ( *p=='1' && p[1]=='1' && p[2]=='8')
	   fillin(ord,ltb, on,118);     /* Element 118   */
	 else {
	   lerr[ p-es ] = '^';  
	   v_error("%s\nSymbol ['%c'?] unknown in El_structure file \n\n",
		   lerr,*p);
	 }
	    
	 while( *++p != ']')   /* Go after brackets */
	   if( *p==';')
	     v_error("Error: no ']' found in El_structure\n\n");
       }
       else if( *p == '(' ) {
	 n = *(++p)-'0';
	 if( n<1 || n>7)
	 {
	   lerr[ p-es ] = '^'; 
	   v_error("%s\nQuantum number N=%d out of range in El_str. file\n\n",
		   lerr, n);
	 }
	 
	 while( *++p != ')') {
	   if( *p==';')
	     v_error("Error: no ')' found in El_structure\n\n");
	   
	   if     ( *p=='s' || *p=='S') l=0;
	   else if( *p=='p' || *p=='P') l=1; 
	   else if( *p=='d' || *p=='D') {
	       l=2;
	       if( n>6) {
		 lerr[ p-es ] = '^'; 
		 v_error("%s\nNo orbital 'd' with n = %d allowed\n\n",
		     lerr, n);
	       }
	   }
	   else if( *p=='f' || *p=='F') {
	       l=3;
	       if( n>5) {
		 lerr[ p-es ] = '^';
		 v_error("%s\nNo orbital 'f' with n = %d allowed\n\n",
		     lerr, n);
	       }
	   }
	   else if( *p=='g' || *p=='G') {
	       l=4;
	       if( n>5) {
		 lerr[ p-es ] = '^';
		 v_error("%s\nNo orbital 'g' with n = %d allowed\n\n",
		     lerr, n);
	       }
	   }
 
	   else 
	   {
	     lerr[ p-es ] = '^'; 
	     v_error("%s\nNo orbital 'spdf' found in El_structure file\n\n",
		     lerr);
	   }
	   if( n <= l)
	   {
	     lerr[ p-es ] = '^'; 
	     v_error("%s\nNo orbital '%c' with n = %d allowed!!!\n\n",
		     lerr, *p, n);
	   }
	     
	   occ = *++p-'0';
	   if( occ<0 || occ>9)
	   {
	     lerr[ p-es ] = '^'; 
	     v_error("%s\nNumber aspected after orbital, '%c' instead\n\n",
		     lerr, *p);
	   }
	   ocu = *(p+1)-'0';
	   if( ocu>=0 && ocu<=9) { occ = 10*occ+ocu; p++;}
	   if( occ>4*l+2)
	   {
	     lerr[ p-es ] = '^'; 
	     v_error("%s\nOccupation number too large %d in El_struct.\n\n",
		     lerr, occ);
	   }

	   on[ ix[n-1]+l ] = occ;
	 }
       }
       else
       {
	 lerr[ p-es ] = '^'; 
	 v_error("%s\nCould not find open bracket in El_struct.\n\n",
		 lerr);
       }
     }
  
     for( i=0,totocc=0; i<TOTSH; i++)
       totocc += on[i];

     printf("Atomic number: %d\tTotal number of electrons: %d\tValence: %+d\n\n",
	     nat, totocc, (zat=nat-totocc));

/*******************************************************************
                 Now let's fill out onp
 *******************************************************************/

     for( n=0, i=0; n<TOTSH && totocc; n++) /* n & i indexes of on & onp */
     {
       l = ltb[ ord[n] ];
       ll = 4*l+2;
       if( (occ = on[ ord[n] ]) == 0)
       {
	 i += ll;
	 continue;
       }

       totocc -= occ;

       for( j=ll-2; j>=0 && occ; j-=2, --occ)
	 onp[i+j] = 1;
       for( j=1; j<ll && occ; j+=2, --occ)
	 onp[i+j] = 1;

       i += ll;
     }
  return zat; 
}
#undef nat

static void fillin(const int ordi[TOTSH], const int ltb[TOTSH], int on[TOTSH], int ntot)
{
 int i, ll;

    for( i=0; (i<TOTSH && ntot); i++)
    {
      ll = 4*ltb[ ordi[i] ] +2;
      on[ ordi[i] ] = ll;
      ntot -= ll;
    }
 return;
} 

static char Tolower( char c)
{
  if( c>='A' && c<='Z') 
    return (c-'A'+'a');
  else return c;
}

static char Toupper( char c)
{
  if( c>='a' && c<='z') 
    return (c-'a'+'A');
  else return c;
}

static void expandsq( char *sq, char *es, int *n)
{
  char *p, bef[128], aft[128], in[128];
  int i;

  /* Store things before [] in bef */
  for( i=0,p=es; *p!='['; p++)
    bef[i++] = *p;
  bef[i] = '\0';
  /* Skip [*sq] */ 
  while( *p++ != ']')
    ;
  /* Store things after [] in aft */
  for( i=0; *p!='\0'; p++)
    aft[i++] = *p;
  aft[i] = '\0';

#define Elif(x,y,z)  else if(strcmp(sq,x)==0){strcpy(in,y);*n=z;}
  /* Start of symbol parsing */
  if( strcmp( sq, "h" ) == 0 ) { strcpy( in, "(1s1)"); *n=1; }
  Elif("he", "[He]",			   2)
  Elif("li", "[He](2s1)",		   3)
  Elif("be", "[He](2s2)",		   4)
  Elif("b" , "[He](2s2p1)",		   5)
  Elif("c" , "[He](2s2p2)",		   6)
  Elif("n" , "[He](2s2p3)",		   7)
  Elif("o" , "[He](2s2p4)",		   8)
  Elif("f" , "[He](2s2p5)",	           9)
  Elif("ne", "[Ne]",		          10)
  Elif("na", "[Ne](3s1)",		  11)
  Elif("mg", "[Ne](3s2)",		  12)
  Elif("al", "[Ne](3s2p1)",	          13)
  Elif("si", "[Ne](3s2p2)",	          14)
  Elif("p" , "[Ne](3s2p3)",	          15)
  Elif("s" , "[Ne](3s2p4)",	          16)
  Elif("cl", "[Ne](3s2p5)",	          17)
  Elif("ar", "[Ar]",		          18)
  Elif("k" , "[Ar](4s1)",		  19)
  Elif("ca", "[Ar](4s2)",		  20)
  Elif("sc", "[Ar](3d1)(4s2)",	          21)
  Elif("ti", "[Ar](3d2)(4s2)",	          22)
  Elif("v" , "[Ar](3d3)(4s2)",	          23)
  Elif("cr", "[Ar](3d5)(4s1)",	          24)
  Elif("mn", "[Ar](3d5)(4s2)",	          25)
  Elif("fe", "[Ar](3d6)(4s2)",	          26)
  Elif("co", "[Ar](3d7)(4s2)",	          27)
  Elif("ni", "[Ar](3d8)(4s2)",	          28)
  Elif("cu", "[Ar](3d10)(4s1)",	          29)
  Elif("zn", "[Ar](3d10)(4s2)",	          30)
  Elif("ga", "[Ar](3d10)(4s2p1)",         31)
  Elif("ge", "[Ar](3d10)(4s2p2)",         32)
  Elif("as", "[Ar](3d10)(4s2p3)",         33)
  Elif("se", "[Ar](3d10)(4s2p4)",         34)
  Elif("br", "[Ar](3d10)(4s2p5)",         35)
  Elif("kr", "[Kr]",		          36)
  Elif("rb", "[Kr](5s1)",		  37)
  Elif("sr", "[Kr](5s2)",		  38)
  Elif("y" , "[Kr](4d1)(5s2)",	          39)
  Elif("zr", "[Kr](4d2)(5s2)",	          40)
  Elif("nb", "[Kr](4d4)(5s1)",	          41)
  Elif("mo", "[Kr](4d5)(5s1)",	          42)
  Elif("tc", "[Kr](4d5)(5s2)",	          43)
  Elif("ru", "[Kr](4d7)(5s1)",	          44)
  Elif("rh", "[Kr](4d8)(5s1)",		  45)
  Elif("pd", "[Kr](4d10)",		  46)
  Elif("ag", "[Kr](4d10)(5s1)",	          47)
  Elif("cd", "[Kr](4d10)(5s2)",           48)
  Elif("in", "[Kr](4d10)(5s2p1)",         49)
  Elif("sn", "[Kr](4d10)(5s2p2)",         50)
  Elif("sb", "[Kr](4d10)(5s2p3)",         51)
  Elif("te", "[Kr](4d10)(5s2p4)",         52)
  Elif("i" , "[Kr](4d10)(5s2p5)",         53)
  Elif("xe", "[Xe]",			  54)
  Elif("cs", "[Xe](6s1)",                 55)
  Elif("ba", "[Xe](6s2)",                 56)
  Elif("la", "[Xe](5d1)(6s2)",            57)
  Elif("ce", "[Xe](4f1)(5d1)(6s2)",       58)
  Elif("pr", "[Xe](4f3)(6s2)",            59)
  Elif("nd", "[Xe](4f4)(6s2)",            60)
  Elif("pm", "[Xe](4f5)(6s2)",            61)
  Elif("sm", "[Xe](4f6)(6s2)",            62)
  Elif("eu", "[Xe](4f7)(6s2)",            63)
  Elif("gd", "[Xe](4f7)(5d1)(6s2)",       64)
  Elif("tb", "[Xe](4f9)(6s2)",            65)
  Elif("dy", "[Xe](4f10)(6s2)",           66)
  Elif("ho", "[Xe](4f11)(6s2)",           67)
  Elif("er", "[Xe](4f12)(6s2)",           68)
  Elif("tm", "[Xe](4f13)(6s2)",           69)
  Elif("yb", "[Xe](4f14)(6s2)",           70)
  Elif("lu", "[Xe](4f14)(5d1)(6s2)",      71)
  Elif("hf", "[Xe](4f14)(5d2)(6s2)",      72)
  Elif("ta", "[Xe](4f14)(5d3)(6s2)",      73)
  Elif("w" , "[Xe](4f14)(5d4)(6s2)",      74)
  Elif("re", "[Xe](4f14)(5d5)(6s2)",      75)
  Elif("os", "[Xe](4f14)(5d6)(6s2)",      76)
  Elif("ir", "[Xe](4f14)(5d7)(6s2)",      77)
  Elif("pt", "[Xe](4f14)(5d9)(6s1)",      78)
  Elif("au", "[Xe](4f14)(5d10)(6s1)",     79)
  Elif("hg", "[Xe](4f14)(5d10)(6s2)",     80)
  Elif("tl", "[Xe](4f14)(5d10)(6s2p1)",   81)
  Elif("pb", "[Xe](4f14)(5d10)(6s2p2)",   82)
  Elif("bi", "[Xe](4f14)(5d10)(6s2p3)",   83)
  Elif("po", "[Xe](4f14)(5d10)(6s2p4)",   84)
  Elif("at", "[Xe](4f14)(5d10)(6s2p5)",   85)
  Elif("rn", "[Rn]",		          86)
  Elif("fr", "[Rn](7s1)",		  87)
  Elif("ra", "[Rn](7s2)",		  88)
  Elif("ac", "[Rn](6d1)(7s2)",		  89)
  Elif("th", "[Rn](6d2)(7s2)",		  90)
  Elif("pa", "[Rn](5f2)(6d1)(7s2)",       91)
  Elif("u", "[Rn](5f3)(6d1)(7s2)",        92)
  Elif("np", "[Rn](5f4)(6d1)(7s2)",       93)
  Elif("pu", "[Rn](5f6)(7s2)",		  94)
  Elif("am", "[Rn](5f7)(7s2)",		  95)
  Elif("cm", "[Rn](5f7)(6d1)(7s2)",       96)
  Elif("bk", "[Rn](5f9)(7s2)",		  97)
  Elif("cf", "[Rn](5f10)(7s2)",		  98)
  Elif("es", "[Rn](5f11)(7s2)",	          99)
  Elif("fm", "[Rn](5f12)(7s2)",	         100)
  Elif("md", "[Rn](5f13)(7s2)",	         101)
  Elif("no", "[Rn](5f14)(7s2)",	         102)
  Elif("lr", "[Rn](5f14)(6d1)(7s2)",     103)
  Elif("rf", "[Rn](5f14)(6d2)(7s2)",     104)
  Elif("db", "[Rn](5f14)(6d3)(7s2)",     105)
  Elif("sg", "[Rn](5f14)(6d4)(7s2)",     106)
  Elif("bh", "[Rn](5f14)(6d5)(7s2)",     107)
  Elif("hs", "[Rn](5f14)(6d6)(7s2)",     108)
  Elif("mt", "[Rn](5f14)(6d7)(7s2)",     109)
  Elif("110","[Rn](5f14)(6d8)(7s2)",     110)
  Elif("111","[Rn](5f14)(6d9)(7s2)",     111)
  Elif("112","[Rn](5f14)(6d10)(7s2)",    112)
  Elif("113","[Rn](5f14)(6d10)(7s2p1)",  113)
  Elif("114","[Rn](5f14)(6d10)(7s2p2)",  114)
  Elif("115","[Rn](5f14)(6d10)(7s2p3)",  115)
  Elif("116","[Rn](5f14)(6d10)(7s2p4)",  116)
  Elif("117","[Rn](5f14)(6d10)(7s2p5)",  117)
  Elif("118","[118]",  118)
/*   Elif("118","[Rn](5f14)(6d10)(7s2p6)",  118) */
     /* and now some +1 ions follow */
  Elif("he+", "(1s1)",			   2)
  Elif("li+", "[He]",			   3)
  Elif("be+", "[He](2s1)",		   4)
  Elif("b+" , "[He](2s2)",		   5)
  Elif("c+" , "[He](2s2p1)",		   6)
  Elif("n+" , "[He](2s2p2)",		   7)
  Elif("o+" , "[He](2s2p3)",		   8)
  Elif("f+" , "[He](2s2p4)",		   9)
  Elif("ne+", "[He](2s2p5)",		  10)
  Elif("na+", "[Ne]",			  11)
  Elif("mg+", "[Ne](3s1)",		  12)
  Elif("al+", "[Ne](3s2)",		  13)
  Elif("si+", "[Ne](3s2p1)",		  14)
  Elif("p+" , "[Ne](3s2p2)",		  15)
  Elif("s+" , "[Ne](3s2p3)",		  16)
  Elif("cl+", "[Ne](3s2p4)",		  17)
  Elif("ar+", "[Ne](3s2p5)",		  18)
  Elif("k+" , "[Ar]",			  19)
  Elif("ca+", "[Ar](4s1)",		  20)
  Elif("sc+", "[Ar](3d1)(4s1)",		  21)
  Elif("ti+", "[Ar](3d2)(4s1)",		  22)
  Elif("v+" , "[Ar](3d3)(4s1)",		  23)
  Elif("cr+", "[Ar](3d5)",		  24)
  Elif("mn+", "[Ar](3d5)(4s1)",		  25)
  Elif("fe+", "[Ar](3d6)(4s1)",		  26)
  Elif("co+", "[Ar](3d7)(4s1)",		  27)
  Elif("ni+", "[Ar](3d8)(4s1)",		  28)
  Elif("cu+", "[Ar](3d10)",		  29)
  Elif("zn+", "[Ar](3d10)(4s1)",	  30)
  Elif("ga+", "[Ar](3d10)(4s2)",	  31)
  Elif("ge+", "[Ar](3d10)(4s2p1)",	  32)
  Elif("as+", "[Ar](3d10)(4s2p2)",	  33)
  Elif("se+", "[Ar](3d10)(4s2p3)",	  34)
  Elif("br+", "[Ar](3d10)(4s2p4)",	  35)
  Elif("kr+", "[Ar](3d10)(4s2p5)",	  36)
  Elif("rb+", "[Kr]",			  37)
  Elif("sr+", "[Kr](5s1)",		  38)
  Elif("y+" , "[Kr](4d1)(5s1)",		  39)
  Elif("zr+", "[Kr](4d2)(5s1)",		  40)
  Elif("nb+", "[Kr](4d4)",		  41)
  Elif("mo+", "[Kr](4d5)",		  42)
  Elif("tc+", "[Kr](4d5)(5s1)",		  43)
  Elif("ru+", "[Kr](4d7)",		  44)
  Elif("rh+", "[Kr](4d8)",		  45)
  Elif("pd+", "[Kr](4d9)",		  46)
  Elif("ag+", "[Kr](4d10)",		  47)
  Elif("cd+", "[Kr](4d10)(5s1)",	  48)
  Elif("in+", "[Kr](4d10)(5s2)",	  49)
  Elif("sn+", "[Kr](4d10)(5s2p1)",	  50)
  Elif("sb+", "[Kr](4d10)(5s2p2)",	  51)
  Elif("te+", "[Kr](4d10)(5s2p3)",	  52)
  Elif("i+" , "[Kr](4d10)(5s2p4)",	  53)
  Elif("xe+", "[Kr](4d10)(5s2p5)",	  54)
  Elif("cs+", "[Xe]",			  55)
  Elif("ba+", "[Xe](6s1)",		  56)
  Elif("la+", "[Xe](4f1)(6s1)",		  57)
  Elif("ce+", "[Xe](4f1)(5d2)",		  58)
  Elif("pr+", "[Xe](4f3)(6s2)",		  59)
  Elif("nd+", "[Xe](4f4)(6s1)",		  60)
  Elif("pm+", "[Xe](4f5)(6s1)",		  61)
  Elif("sm+", "[Xe](4f6)(6s1)",		  62)
  Elif("eu+", "[Xe](4f7)(6s1)",		  63)
  Elif("gd+", "[Xe](4f7)(5d1)(6s1)",	  64)
  Elif("tb+", "[Xe](4f9)(6s1)",		  65)
  Elif("dy+", "[Xe](4f10)(6s1)",	  66)
  Elif("ho+", "[Xe](4f11)(6s1)",	  67)
  Elif("er+", "[Xe](4f12)(6s1)",	  68)
  Elif("tm+", "[Xe](4f13)(6s1)",	  69)
  Elif("yb+", "[Xe](4f14)(6s1)",	  70)
  Elif("lu+", "[Xe](4f14)(6s2)",	  71)
  Elif("hf+", "[Xe](4f14)(5d2)(6s1)",	  72)
  Elif("ta+", "[Xe](4f14)(5d3)(6s1)",	  73)
  Elif("w+" , "[Xe](4f14)(5d4)(6s1)",	  74)
  Elif("re+", "[Xe](4f14)(5d5)(6s1)",	  75)
  Elif("os+", "[Xe](4f14)(5d6)(6s1)",	  76)
  Elif("ir+", "[Xe](4f14)(5d7)(6s1)",	  77)
  Elif("pt+", "[Xe](4f14)(5d9)",	  78)
  Elif("au+", "[Xe](4f14)(5d10)",	  79)
  Elif("hg+", "[Xe](4f14)(5d10)(6s1)",	  80)
  Elif("tl+", "[Xe](4f14)(5d10)(6s2)",	  81)
  Elif("pb+", "[Xe](4f14)(5d10)(6s2p1)",  82)
  Elif("bi+", "[Xe](4f14)(5d10)(6s2p2)",  83)
  Elif("po+", "[Xe](4f14)(5d10)(6s2p3)",  84)
  Elif("at+", "[Xe](4f14)(5d10)(6s2p4)",  85)
  Elif("rn+", "[Xe](4f14)(5d10)(6s2p5)",  86)
  Elif("fr+", "[Rn]",			  87)
  Elif("ra+", "[Rn](7s1)",		  88)
  Elif("ac+", "[Rn](6d1)(7s1)",		  89)
  Elif("th+", "[Rn](6d2)(7s1)",		  90)
  Elif("pa+", "[Rn](5f2)(6d1)(7s1)",      91)
  Elif("u+" , "[Rn](5f3)(6d1)(7s1)",      92)
  Elif("np+", "[Rn](5f4)(6d1)(7s1)",      93)
  Elif("pu+", "[Rn](5f6)(7s1)",		  94)
     /* Some -1 ions follow */  
  Elif("h-" , "[He]",			   1)
  Elif("f-" , "[Ne]",			   9)
  Elif("cl-", "[Ar]",			  17)
     /* Lanthanides in the trivalent state follow */  
  Elif("la_iii", "[Xe](5d1)(6s2)",	  57)
  Elif("ce_iii", "[Xe](4f1)(5d1)(6s2)",	  58)
  Elif("pr_iii", "[Xe](4f2)(5d1)(6s2)",	  59)
  Elif("nd_iii", "[Xe](4f3)(5d1)(6s2)",	  60)
  Elif("pm_iii", "[Xe](4f4)(5d1)(6s2)",	  61)
  Elif("sm_iii", "[Xe](4f5)(5d1)(6s2)",	  62)
  Elif("eu_iii", "[Xe](4f6)(5d1)(6s2)",	  63)
  Elif("gd_iii", "[Xe](4f7)(5d1)(6s2)",	  64)
  Elif("tb_iii", "[Xe](4f8)(5d1)(6s2)",	  65)
  Elif("dy_iii", "[Xe](4f9)(5d1)(6s2)",	  66)
  Elif("ho_iii", "[Xe](4f10)(5d1)(6s2)",  67)
  Elif("er_iii", "[Xe](4f11)(5d1)(6s2)",  68)
  Elif("tm_iii", "[Xe](4f12)(5d1)(6s2)",  69)
  Elif("yb_iii", "[Xe](4f13)(5d1)(6s2)",  70)
  Elif("lu_iii", "[Xe](4f14)(5d1)(6s2)",  71)
     /* Lanthanum in the bivalent state follows */
  Elif("la_ii",  "[Xe](4f1)(6s2)",	  57)
  else
    v_error("Symbol [%s] not recognized in El_structure. Use noble gases in []. Exit.\n", sq);

  sprintf(es, "%s%s%s", bef,in,aft);
}
#undef Elif

static void create_el_str(void)
{
  FILE *fel;

  fel = fopen("El_structure", "w");
  if( !fel)
    v_error("Cannot create default El_structure file. Aborting.\n");
  fprintf( fel, 
	   "# Electronic structure file:\n"
	   "#\n"
	   "# - Comment lines begin with '#'.\n"
	   "# - Empty lines are ignored.\n"
	   "# - Letter case is irrelevant.\n"
	   "# - Electronic structure MUST end with ';'.\n"
	   "# - Electronic structure may be splitted into several lines.\n"
	   "#   any space can be introduced.\n"
	   "# - A line may contain: \"Name nice-name\" to be used creating files.\n"
	   "# - A line may contain: \"Atomic number ##\", \n"
	   "#   if otherwise absent, the default value corresponding\n"
	   "#   to the element in square brackets will be used.\n"
	   "# - Any structure may contain elements into [].\n"
	   "# - Any shell MUST be enclosed in ().\n"
	   "# - Any structure MUST have the principal quantum number after the '('.\n"
	   "# - It is possible to modify noble elements reference. \n"
	   "#  \n"
	   "# EXAMPLES: \n"
	   "#    [Ne](1s1); is Ne- with only one 1s electron. \n"
	   "#    [Xe] (4f7) (5d1) (6s2); is Gd.\n"
	   "#    (1s2)(2s2 p6); is the same as [Ne];\n"
	   "#\n"
	   "\n"
	   "Name Gd+\n"
	   "\n"
	   "Atomic number 64\n"
	   "       \n"
	   "# and this is its electronic structure:\n"
	   "\n"
	   "    [Xe](4f7)(5d1)(6s1);\n"
	   " \n"
	   "# End\n");
  fclose( fel);
  printf("Default El_structure file created. You MUST edit it!\n");
  printf("Exiting...\n");
  exit(1);
}
