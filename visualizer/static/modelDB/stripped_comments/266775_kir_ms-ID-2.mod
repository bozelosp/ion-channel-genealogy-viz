NEURON {
    SUFFIX kir_ms
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik, shift
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 	(S/cm2) 
    shift = 0.0 (mV)
    q = 1 	
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    gk (S/cm2)
    minf
    mtau (ms)
}

STATE { m }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gbar*m
    ik = gk*(v-ek)
}

DERIVATIVE states {
    rates()
    m' = (minf-m)/mtau*q
}

INITIAL {
    rates()
    m = minf
}

PROCEDURE rates() {
    UNITSOFF
    minf = 1/(1+exp((v-(-82)-shift)/13))
    mtau = 1/(exp((v-(-103))/(-14.5))+0.125/(1+exp((v-(-35))/(-19))))
    UNITSON
}