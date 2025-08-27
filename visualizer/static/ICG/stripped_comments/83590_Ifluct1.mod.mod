INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS Ifluct1
    RANGE m, s, tau, x
    RANGE new_seed
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp) 
    (mV) = (millivolt)
}

PARAMETER {
     dt   (ms)
     m   = 0. (nA)      
     s   = 0. (nA)      
     tau = 2. (ms)      
}

ASSIGNED {
    i     (nA)          
    x                   
}

INITIAL {
    x = m               
}


BREAKPOINT {
    SOLVE oup
    if (tau <= 0) {  x = m + s  * normrand(0,1) }  
    i = - x
}


PROCEDURE oup() {       
if (tau > 0) {  x = x + (1. - exp(-dt/tau)) * (m - x) + sqrt(1.-exp(-2.*dt/tau)) * s  * normrand(0,1) }
}


PROCEDURE new_seed(seed) {      
    set_seed(seed)
    VERBATIM
      printf("Setting random generator with seed = %g\n", _lseed);
    ENDVERBATIM
}