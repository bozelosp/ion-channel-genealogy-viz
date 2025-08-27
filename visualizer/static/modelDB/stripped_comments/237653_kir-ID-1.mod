NEURON {
    THREADSAFE
    SUFFIX kir
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik, base, factor
    POINTER pka
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 (S/cm2) 
    q = 3
    base   = 0.0      
	factor = 0.0      
} 

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    gk (S/cm2)
    minf
    mtau (ms)
    pka (1)
}

STATE { m }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = modulation() * gbar*m
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
    LOCAL alpha, beta, sum
    UNITSOFF
    minf = 1/(1+exp((v-(-102))/13))
    alpha = 0.1*exp((v-(-60))/(-14))
    beta = 0.27/(1+exp((v-(-31))/(-23)))
    sum = alpha+beta
    mtau = 1/sum
    UNITSON
}


FUNCTION modulation() {
    
    
    
    modulation = 1 + factor * (pka - base)
    
}