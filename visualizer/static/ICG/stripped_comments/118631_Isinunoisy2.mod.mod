INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS Isinunoisy2
    RANGE amp, i, freq, m, s, tau, x, new_seed
    ELECTRODE_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    dt   (ms)
    m   = 0. (nA)       
    s   = 0.5 (nA)       
    tau = 0. (ms)       
    amp = 0. (nA)       
    freq= 0. (Hz)       
    phas= 0. (HZ)       
    fr2 = 0.(Hz)        
}

ASSIGNED { 
    i (nA)              
    x                   
}

INITIAL {
    i = m
    x = 0               
}

BREAKPOINT {  
    SOLVE oup
    
    if (tau <= 0) {  x =  normrand(0,1)  }
     
    
     i = x * (s + amp * sin(0.0062831853071795866 * freq * t)) + m

}


PROCEDURE oup() {       
if (tau > 0) {  x = x + (1. - exp(-dt/tau)) * ( - x) + sqrt(1.-exp(-2.*dt/tau)) * normrand(0,1)}
}

PROCEDURE new_seed(seed) {      
    set_seed(seed)
    VERBATIM
      printf("Setting random generator with seed = %g\n", _lseed);
    ENDVERBATIM
}