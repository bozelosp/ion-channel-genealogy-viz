NEURON {
    THREADSAFE
    SUFFIX naf
    USEION na READ ena WRITE ina
    RANGE gbar, gna, ina, mVhalf, hVhalf, mSlope, hSlope, taum, tauh, taun, base, factor
    POINTER pka
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 (S/cm2) 
    q = 1.8
    mVhalf     = -25.0 (mV)
    hVhalf     = -62.0 (mV)
    mSlope     =  -9.2 (mV)
    hSlope     =   6.0 (mV)
    taum       =   0.09 (ms)
    taun       =   0.34 (ms)
    tauh       =   0.34 (ms)
    base   = 0.0      
	factor = 0.0      
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
    pka (1)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = modulation() * gbar*m*m*m*h
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
    minf = 1 / (1 + exp( (v-mVhalf) / mSlope ) )
    hinf = 1 / (1 + exp( (v-hVhalf) / hSlope ) )
    
    mtau = 0.38 + 1/( 0.6*exp((v-(-58.0))/8.0) + 1.8*exp((v-(-58.0))/(-35.0))  )
    
    if (v < - 60) {
        htau = 3.4 + 0.015*v
    }else{
        htau = 0.56 + 1.1/(1+exp((v-(-48))/15.0)) + 1.2/(1+exp((v-(-48))/4.0))
    }
        
    
    
    UNITSON
}

FUNCTION modulation() {
    
    
    
    modulation = 1 + factor * (pka - base)
    
}