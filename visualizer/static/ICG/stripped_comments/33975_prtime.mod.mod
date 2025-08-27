NEURON {
    SUFFIX nothing
}

VERBATIM
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <stdio.h>
ENDVERBATIM


FUNCTION prtime () {
VERBATIM
  double prt;
  static double PRTIME;
  prt = (clock()-PRTIME)/CLOCKS_PER_SEC;
  
  if (prt<0) prt += UINT_MAX/CLOCKS_PER_SEC; 
  PRTIME=clock();
  _lprtime = prt;
ENDVERBATIM
}