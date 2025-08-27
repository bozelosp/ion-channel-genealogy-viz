INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX gDKkin3
    USEION k READ ek WRITE ik
    RANGE  gk, gbar
    GLOBAL VoltageOffset, gv, tau, tstart
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
}

PARAMETER {
        gv[300]                             (ps/um2)
        tau[300]                            (ms)        
        gbar            = 5                 (pS/um2)    
        v                                   (mV)
        VoltageOffset   = 100.5             (mV)        
        tstart          = 10                (ms)        
        dt                                  (ms)
        measTime                            (ms)        
         celsius                             (degC)
}


ASSIGNED {
        ik      (mA/cm2)
        gk      (pS/um2)
        ek      (mV)
       }



BREAKPOINT {
        SOLVE states
        ik = (1e-4) * gk * (v - ek)
}


PROCEDURE states() {        
                          
        gbar = gv[v+VoltageOffset]
        if(t<tstart) { gk=0 }
        if (t>tstart) {gk = gbar}
        VERBATIM
        return 0;
        ENDVERBATIM
}