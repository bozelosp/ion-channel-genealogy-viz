NEURON {
    THREADSAFE
    SUFFIX kaf
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
    q = 2
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
    hinf
    htau (ms)
    pka (1)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = modulation() * gbar*m*m*h
    ik = gk*(v-ek)
}

DERIVATIVE states {
    rates()
    m' = (minf-m)/mtau*q
    h' = (hinf-h)/htau*q
}

INITIAL {
    rates()
    m = minf
    h = hinf
}

PROCEDURE rates() {
    LOCAL alpha, beta, sum
    UNITSOFF
    alpha = 1.5/(1+exp((v-4)/(-17)))
    beta = 0.6/(1+exp((v-10)/9))
    sum = alpha+beta
    minf = alpha/sum
    mtau = 1/sum

    alpha = 0.105/(1+exp((v-(-121))/22))
    beta = 0.065/(1+exp((v-(-55))/(-11)))
    sum = alpha+beta
    hinf = alpha/sum
    htau = 1/sum
    UNITSON
}


FUNCTION modulation() {
    
    
    
    modulation = 1 + factor * (pka - base)
    
}