INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX nothing
}


PROCEDURE vseed (seed) {
  VERBATIM
  srand48((unsigned)_lseed);
  set_seed(_lseed);
  ENDVERBATIM
}