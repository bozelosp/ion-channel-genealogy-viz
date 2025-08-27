NEURON {
    SUFFIX kaf_ms
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik, q
    RANGE damod, maxMod, level
    
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 (S/cm2) 
    q = 1	
    
    
    damod = 0
    maxMod = 1
    level = 0
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
    UNITSOFF
    minf = 1/(1+exp((v-(-10))/(-17.7)))
    mtau = 0.9+1.1/(1+exp((v-(-30))/10))
    hinf = 1/(1+exp((v-(-75.6))/11.8))
    htau = 14
    UNITSON

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}


FUNCTION modulation() {
    
    
    modulation = 1 + damod*(maxMod-1)*level 
}