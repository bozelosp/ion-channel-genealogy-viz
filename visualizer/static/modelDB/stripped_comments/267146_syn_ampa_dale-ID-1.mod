NEURON {
    POINT_PROCESS syn_ampa_dale
    RANGE tau_o                           
    RANGE tau_c                          
    RANGE erev                              
    RANGE syn_step
    RANGE i                                 
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
    erev = 0.0 (mV)
    syn_step = 1.25
}

ASSIGNED {
    v (mV)
    i (nA)
}

STATE {
    o
    c
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
    o = o + syn_step*weight
    c = c + syn_step*weight
}

DERIVATIVE states {
    o' = -o/tau_o
    c' = -c/tau_c
}