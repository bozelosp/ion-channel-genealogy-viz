NEURON {
        SUFFIX atten     
        RANGE v_atten
}

ASSIGNED {
        v (millivolt)
        v_atten (millivolt)
}

INITIAL {
        v_atten = 0
}

BREAKPOINT {

VERBATIM
  v_atten = (v + 60.0) / (60.0);
ENDVERBATIM
}