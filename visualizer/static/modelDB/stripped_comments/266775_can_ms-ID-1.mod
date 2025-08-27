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
    SUFFIX can_ms
    USEION ca READ cai, cao WRITE ica VALENCE 2
    RANGE pbar, ica
    RANGE damod, maxMod, level
}

PARAMETER {
    pbar = 0.0 	(cm/s)
    a = 0.21
    
    q = 3	
    damod = 0
    maxMod = 1
    level = 0
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
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = pbar*m*m*(h*a+1-a)*ghk(v, cai, cao)*modulation()
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
    mtau = 0.06+1/(exp((v-25)/18)+exp((v-(-31))/(-44)))
    hinf = 1/(1+exp((v-(-74.8))/6.5))
    htau = 70
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