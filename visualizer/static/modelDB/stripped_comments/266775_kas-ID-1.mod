NEURON {
    THREADSAFE
    SUFFIX kas
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik
    RANGE damod, maxMod, level, max2, lev2
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 (S/cm2) 
    q = 3
    a = 0.2
    damod = 0
    maxMod = 1
    level = 0
    max2 = 1
    lev2 = 0
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
    gk = gbar*m*m*h*modulation()
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
    alpha = 0.25/(1+exp((v-50)/(-20)))
    beta = 0.05/(1+exp((v-(-90))/35))
    sum = alpha+beta
    minf = alpha/sum
    mtau = 1/sum

    alpha = 0.0025/(1+exp((v-(-95))/16))
    beta = 0.002/(1+exp((v-50)/(-70)))
    sum = alpha+beta
    hinf = a+(alpha/sum)*(1-a)
    htau = 1/sum
    UNITSON
}

FUNCTION modulation() {
    
    
    modulation = 1 + damod * ( (maxMod-1)*level + (max2-1)*lev2 ) 
    if (modulation < 0) {
        modulation = 0
    }  
}