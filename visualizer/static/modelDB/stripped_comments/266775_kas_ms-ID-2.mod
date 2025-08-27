NEURON {
    SUFFIX kas_ms
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 	(S/cm2) 
    a = 0.8
    
    q = 3	
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
    gk = gbar*m*m*(h*a+1-a)
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
    minf = 1/(1+exp((v-(-27))/(-16)))
    mtau = 3.4+89.2*exp(-((v-(-34.3))/30.1)^2)
    hinf = 1/(1+exp((v-(-33.5))/21.5))
    htau = 548.7*6/(exp((v-(-96))/(-29.01))+exp((v-(-96))/100))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    UNITSON
}