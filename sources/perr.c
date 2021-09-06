#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

extern FILE *fperr;

/* print on file *fperr opened before as global */
void perr(const char *fmt, ...)
{
va_list ap;
 va_start(ap, fmt);
 vfprintf( fperr, fmt, ap);
 va_end(ap);
}

void operr(char *errfilename)
{
  fperr = fopen( errfilename, "w");
  if(!fperr) {
    fprintf(stderr,"Error opening file %s as error log file \n", errfilename);
    exit(1);
  }
}

void clerr(void)
{ fclose( fperr); }
