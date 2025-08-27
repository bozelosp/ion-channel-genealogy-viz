INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX nothing
}

VERBATIM
#include "misc.h"
ENDVERBATIM

FUNCTION istmpobj () {
VERBATIM
  _listmpobj=hoc_is_tempobj_arg(1);
ENDVERBATIM  
}