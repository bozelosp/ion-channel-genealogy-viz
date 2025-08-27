UNITS {
    (molar) = (1/liter)
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (millimolar)
}

NEURON {
    SUFFIX sk
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gbar, ik, I
}

PARAMETER {
    gbar = 0.0 (mho/cm2)
    q = 1
}

ASSIGNED {
    v (mV)
    ik (mA/cm2)
    cai (mM) 
    ek (mV)
    oinf
    otau (ms)
    I (mA/cm2)
}

STATE { o }

BREAKPOINT {
    SOLVE state METHOD cnexp
    I  = gbar*o*(v-ek)
    ik = I 
}

DERIVATIVE state {
    rate(v, cai)
    o' = (oinf-o)/otau*q
}

INITIAL {
    rate(v, cai)
    o = oinf
}

PROCEDURE rate(v (mV), ca (mM)) {
    LOCAL a
    
    a = (ca/0.57e-3)^5.2
    oinf = a/(1+a)
    otau = 4.9
}