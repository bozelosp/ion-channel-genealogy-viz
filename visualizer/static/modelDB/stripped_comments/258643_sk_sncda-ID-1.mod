NEURON {
    SUFFIX sk
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gbar, oinf, n, ik
    GLOBAL km
    EXTERNAL apc_metap, fpc_metap
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)    
    (molar) = (1/liter)
    (mM) = (millimolar)
    (S) = (siemens)
}

PARAMETER {
    gbar = 1e-05    (S/cm2)
    n = 4 
    km = 0.0002 (mM) 
    cai (mM)
    dt  (ms)
    ek  (mV)

    
    

    v   (mV)
}
ASSIGNED {
    ik  (mA/cm2)
    oinf
}

BREAKPOINT {
    oinf = 1/(1 + pow(km/cai,n))
    ik = oinf*gbar*(v-ek)
}