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
    SUFFIX cav32
    USEION cal READ cali, calo WRITE ical VALENCE 2
    RANGE pbar, ical, a, perm, I
}

PARAMETER {
    pbar = 6.7e-6   (cm/s)
    mvhalf = -61.5  (mV)
    mslope = -8.0   (mV)
    hvhalf = -73.7  (mV)
    hslope = 9.1   (mV) 
    a      = 0.9
}

ASSIGNED { 
    v (mV)
    ical (mA/cm2)
    ecal (mV)
    celsius (degC)
    cali (mM)
    calo (mM)
    minf
    hinf
    mtau  (ms)
    htau  (ms)
    htau2 (ms)
    htot  (ms)
    perm
    I
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    perm = pbar*m*m*m*h
    ical = ghk(v, cali, calo)*perm
    I    = ical
}

INITIAL {
    rates(v)
    m = minf
    h = hinf
}

DERIVATIVE states { 
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htot
}

PROCEDURE rates(v (mV)) {
    minf  = 1/(1+exp((v-mvhalf)/mslope))
    hinf  = 1/(1+exp((v-hvhalf)/hslope))    
    mtau  = 6.0/(1+exp((v+66.0)/15.0  ))+0.6
    htau  = 4.3/(1+exp(0.06*(v)))+8
    htau2 = 95*exp(-(v+58.0)/25.0)+20
    htot  = a*htau + (1-a)*htau2
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