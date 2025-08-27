NEURON {
    POINT_PROCESS syn_inh_std
    RANGE tau_o, tau_c, erev, i
    RANGE c1, c2
    RANGE Gtau_o, Gtau_c, Ginc
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
    Gtau_o = 40.0 (ms)
    Gtau_c = 41.05 (ms)
    Ginc = 3.0
}

ASSIGNED {
    v (mV)
    i (nA)
    factor
    Gfactor
}

STATE {
    o
    c
}

INITIAL {
    LOCAL tp
    o = 0
    c = 0
    tp = (tau_o*tau_c)/(tau_c-tau_o)*log(tau_c/tau_o)
    factor = -exp(-tp/tau_o)+exp(-tp/tau_c)
    factor = 1/factor
    tp = (Gtau_o*Gtau_c)/(Gtau_c-Gtau_o)*log(Gtau_c/Gtau_o)
    Gfactor = -exp(-tp/Gtau_o)+exp(-tp/Gtau_c)
    Gfactor = 1/Gfactor
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = (c - o) * (v-erev)
}

NET_RECEIVE(weight (uS), w, G_o, G_c, t0 (ms)) {
    printf("%g |%g | ", t-t0, G_c-G_o)
    G_o = G_o*exp(-(t-t0)/Gtau_o)
    G_c = G_c*exp(-(t-t0)/Gtau_c)
    G_o = G_o + Ginc*Gfactor
    G_c = G_c + Ginc*Gfactor
    t0=t

    w=weight*(1+G_c-G_o)
    printf("w = %g \n", w)
    
    o = o + w*1.25
    c = c + w*1.25
}


DERIVATIVE states {
    o' = -o/tau_o
    c' = -c/tau_c
}

UNITSON