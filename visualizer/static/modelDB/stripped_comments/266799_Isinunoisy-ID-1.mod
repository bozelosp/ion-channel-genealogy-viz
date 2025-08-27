INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS Isinunoisy
    RANGE amp, i, freq, m, s, tau, x, new_seed,onset,np
    ELECTRODE_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    dt   (ms)
    m   = 0. (nA)       
    s   = 0. (nA)       
    tau = 2. (ms)       
    amp = 0. (nA)       
    freq= 0. (Hz)       
    phas= 0. (HZ)       
    fr2 = 5.(Hz)       
    onset = 10
    np = 100000000000   
}

ASSIGNED { 
    i (nA)              
    x                   
}

INITIAL {
    i = m
    x = m               
}

BREAKPOINT {  
    SOLVE oup
    
    if (tau <= 0) {  x = m + s  * normrand(0,1) }  
    if ((t<onset)||(t>onset+(1000/freq)*np)) {
        i = x
    } else {
        i = x + amp * sin(0.0062831853071795866 * freq * t)
    }

}


PROCEDURE oup() {       
if (tau > 0) {  x = x + (1. - exp(-dt/tau)) * (m - x) + sqrt(1.-exp(-2.*dt/tau)) * s  * normrand(0,1)}
}

PROCEDURE new_seed(seed) {      
    set_seed(seed)

}