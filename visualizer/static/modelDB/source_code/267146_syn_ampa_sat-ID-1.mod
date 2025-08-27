NEURON {
    POINT_PROCESS syn_ampa_sat
    RANGE tau_o, tau_c, erev, syn_step, ampa_sat, alpha, count
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
    tau_o = 0.2 (ms)
    tau_c = 3.0 (ms)
    erev = 0.0 (mV)
    syn_step = 1.25
    ampa_sat = 100 :0.001
    alpha = 0.95
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

NET_RECEIVE(weight (uS), w, count) {
    w=weight*pow(alpha,count)
    :o = o + syn_step*w
    :c = c + syn_step*w
    o = o + syn_step*(1-(c-o)/ampa_sat)*w
    c = c + syn_step*(1-(c-o)/ampa_sat)*w
    :printf("%g \n", c-o)
    count=count+1
}


DERIVATIVE states {
    o' = -o/tau_o
    c' = -c/tau_c
}
