NEURON {
    SUFFIX naf_ms
    USEION na READ ena WRITE ina
    RANGE gbar, gna, ina
    RANGE damod, maxMod, level
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 (S/cm2) 
    
    q = 1.8	
    damod = 0
    maxMod = 1
    level = 0
}

ASSIGNED {
    v (mV)
    ena (mV)
    ina (mA/cm2)
    gna (S/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gbar*m*m*m*h*modulation()
    ina = gna*(v-ena)
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
    
    
    
    
    minf = 1/(1+exp((v-(-25))/(-10)))
    mtau = 0.33+1/(exp((v-(-62))/14)+exp((v-(-60))/(-17)))
    hinf = 1/(1+exp((v-(-62))/6))
    htau = 0.6+1/(exp((v-(-44))/8)+exp((v-(-99))/(-44)))
    UNITSON
}

FUNCTION modulation() {
    
    
    modulation = 1 + damod*(maxMod-1)*level 
}