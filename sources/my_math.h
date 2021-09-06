#ifndef _MY_MATH
#define _MY_MATH

#include <math.h>

#define v_abs(x)	( ( (x) >= 0 )? (x): -(x) )
#define v_nint(x)   (((x) >=0)? floor((x)+.5): ceil((x)-.5))
#define v_max(x,y)      ( ( (x)>=(y) )? (x): (y) )
#define v_min(x,y)      ( ( (x)<=(y) )? (x): (y) )
#define v_sign(x,y)     ( ( (y) >= 0 )? (x): -(x) )

typedef int Logical;

#define TRUE       1
#define FALSE      0

#define M_PI       3.14159265358979323846
#define M_4PI	   12.5663706143591729539

#endif  /*  _MY_MATH  */
