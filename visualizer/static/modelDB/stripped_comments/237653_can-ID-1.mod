UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

NEURON {
    THREADSAFE
    SUFFIX can
    USEION ca READ cai, cao WRITE ica VALENCE 2
    RANGE pbar, ica, base, factor
    POINTER pka
}

PARAMETER {
    pbar = 0.0 (cm/s)
    a = 0.21
    
    q = 2	
    base   = 0.0      
	factor = 0.0      
} 

ASSIGNED { 
    v (mV)
    ica (mA/cm2)
    eca (mV)
    celsius (degC)  
    cai (mM)
    cao (mM)
    minf
    mtau (ms)
    hinf
    htau (ms)
    pka (1)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = modulation() * pbar*m*m*(h*a+1-a)*ghk(v, cai, cao)
}

INITIAL {
    rates()
    m = minf
    h = hinf
}

DERIVATIVE states { 
    rates()
    m' = (minf-m)/mtau*q
    h' = (hinf-h)/htau*q
}

PROCEDURE rates() {
    UNITSOFF
    minf = 1/(1+exp((v-(-3))/(-8)))
    
    mtau = (0.06+1/(exp((v-25)/18)+exp((v-(-31))/-44)))*2
    hinf = 1/(1+exp((v-(-74.8))/6.5))
    htau = 70
    UNITSON
}

FUNCTION ghk(v (mV), ci (mM), co (mM)) (.001 coul/cm3) {
    LOCAL z, eci, eco
    z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
    if(z == 0) {
        z = z+1e-6
    }
    eco = co*(z)/(exp(z)-1)
    eci = ci*(-z)/(exp(-z)-1)
    ghk = (1e-3)*2*FARADAY*(eci-eco)
}


FUNCTION modulation() {
    
    
    
    modulation = 1 + factor * (pka - base)
    
}