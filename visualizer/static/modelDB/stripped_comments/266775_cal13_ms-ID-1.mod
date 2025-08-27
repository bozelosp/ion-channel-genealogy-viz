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
    SUFFIX cal13_ms
    USEION cal READ cali, calo WRITE ical VALENCE 2
    RANGE pbar, ical
    RANGE damod, maxMod, level
}

PARAMETER {
    pbar = 0.0 (cm/s)
    
    q = 2	
    damod = 0
    maxMod = 1
    level = 0
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
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ical = pbar*m*m*h*ghk(v, cali, calo)*modulation()
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
    minf = 1/(1+exp((v-(-33))/(-6.7)))
    mtau = 0.06+1/(exp((v-10)/20)+exp((v-(-17))/-48))
    hinf = 1/(1+exp((v-(-13.4))/11.9))
    htau = 44.3
    UNITSON
}

FUNCTION ghk(v (mV), ci (mM), co (mM)) (.001 coul/cm3) {
    LOCAL z, eci, eco
    z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
    eco = co*efun(z)
    eci = ci*efun(-z)
    ghk = (1e-3)*2*FARADAY*(eci-eco)
}

FUNCTION efun(z) {
    if (fabs(z) < 1e-4) {
        efun = 1-z/2
    }else{
        efun = z/(exp(z)-1)
    }
}


FUNCTION modulation() {
    
    
    modulation = 1 + damod*(maxMod-1)*level 
}