NEURON {
    THREADSAFE
    SUFFIX kaf
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik
    RANGE damod, maxMod, level, max2, lev2, modShift
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 (S/cm2) 
    q = 2
    damod = 0
    maxMod = 1
    level = 0
    max2 = 1
    lev2 = 0
    modShift = 0
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
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gbar*m*m*h *modulation()
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
    alpha = 1.5/(1+exp((v-4+modShift)/(-17)))
    beta = 0.6/(1+exp((v-10+modShift)/9))
    sum = alpha+beta
    minf = alpha/sum
    mtau = 1/sum
    
    
    alpha = 0.105/(1+exp((v-(-121)+modShift)/22))
    beta = 0.065/(1+exp((v-(-55)+modShift)/(-11)))
    sum = alpha+beta
    hinf = alpha/sum
    htau = 1/sum
    
    UNITSON
}

FUNCTION modulation() {
    
    
    modulation = 1 + damod * ( (maxMod-1)*level + (max2-1)*lev2 ) 
    if (modulation < 0) {
        modulation = 0
    }  
}