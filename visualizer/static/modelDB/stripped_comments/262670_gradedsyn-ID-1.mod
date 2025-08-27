NEURON {
    POINT_PROCESS GradedSyn
    POINTER vpre
    RANGE e, g, i, gbar, vmid, vslope, tau, sinf
    NONSPECIFIC_CURRENT i
}
UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}

PARAMETER {
    e = -80 (mV)  
    vslope = 5.0  (mV)  
    vmid = -40  (mV)    
    gbar = 0.0  (uS)
    tau = 4.0  (ms)
}

ASSIGNED {
    v (mV)
    vpre  (mV)  
    g  (uS)
    i  (nA)
    sinf
}

STATE {
    s
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gbar * s
    i = g * (v - e)
}

INITIAL {
    s = 0.0
}

DERIVATIVE states {
    sinf = 1 / (1 + exp((vmid - vpre) / vslope))
    s' = (sinf - s) / tau
}