NEURON {
    POINT_PROCESS syn_inh_sat
    RANGE tau_o                           
    RANGE tau_c                          
    RANGE erev                              
    RANGE i                                 
    RANGE syn_step
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (S) = (siemens)
}

PARAMETER {
    tau_o = 1.5 (ms)
    tau_c = 4.0 (ms)
    erev = -75 (mV)
    syn_step = 3.0
    inh_sat = 1.0
}

ASSIGNED {
    v (mV)
    i (nA)
}

STATE {
    o (1/ms)
    c (1/ms)
}

INITIAL {
    o = 0
    c = 0
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = (c - o) * (v-erev)
}

NET_RECEIVE(weight (uS)) {
    
    o = o + syn_step*(1-(c-o)/inh_sat)*weight
    c = c + syn_step*(1-(c-o)/inh_sat)*weight
}

DERIVATIVE states {
    o' = -o/tau_o
    c' = -c/tau_c
}