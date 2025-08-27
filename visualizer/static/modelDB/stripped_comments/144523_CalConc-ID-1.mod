NEURON {
    SUFFIX CalConc
    USEION cal READ ical WRITE cali VALENCE 2    
    RANGE cali, kCal, depth
    GLOBAL tauCal
}

UNITS {
    (molar) = (1 / liter)
    (mM) = (millimolar)
    (mV) = (millivolt)
    (mA) = (milliamp)
    PI = (pi) (1)
}

PARAMETER {
    
    

    
    
    
    
    
    kCal = 3.45e-7 (1/coulomb)
    tauCal = 70 (ms)
    caliBase = 50e-6 (mM) 
    depth = 0.2 (micron)
}

ASSIGNED {
    C (kilo / m3 / s)
    D (kilo / m3 / s)
    ical (mA/cm2)
}

STATE {
    cali (mM)
}

INITIAL {
    cali = caliBase
}

BREAKPOINT {
    SOLVE states METHOD cnexp
}

DERIVATIVE states {
    C = (cali - caliBase) / tauCal
	D = - kCal / depth * ical * (1e4)
    cali' = D - C
}