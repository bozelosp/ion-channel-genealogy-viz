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
    SUFFIX cal12
    USEION cal READ cali, calo WRITE ical VALENCE 2
    RANGE pbar, ical, base, factor
    POINTER pka
}

PARAMETER {
    pbar = 0.0 (cm/s)
    a = 0.17
    
    q = 2	          
    base   = 0.0      
	factor = 0.0      
} 

ASSIGNED { 
    v (mV)
    ical (mA/cm2)
    ecal (mV)
    celsius (degC)
    cali (mM)
    calo (mM)
    minf
    mtau (ms)
    hinf
    htau (ms)
    pka (1)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ical = modulation() * pbar*m*(h*a+1-a)*ghk(v, cali, calo)
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
    minf = 1/(1+exp((v-(-8.9))/(-6.7)))
    
    mtau = 0.06+1/(exp((v-10)/20)+exp((v-(-17))/-48))
    hinf = 1/(1+exp((v-(-13.4))/11.9))
    htau = 44.3
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