NEURON {
    SUFFIX CaConc
    USEION ca READ ica WRITE cai
    RANGE cai, kCa, depth
    GLOBAL tauCa
}

UNITS {
    (molar) = (1 / liter)
    (mM) = (millimolar)
    (mV) = (millivolt)
    (mA) = (milliamp)
    PI = (pi) (1)
}

PARAMETER {
    
    

    
    
    
    
    
    kCa = 3.45e-7 (1/coulomb)
    tauCa = 70 (ms)
    caiBase = 50e-6 (mM) 
    depth = 0.2 (micron)
}

ASSIGNED {
    C (kilo / m3 / s)
    D (kilo / m3 / s)
    ica (mA/cm2)
}

STATE {
    cai (mM)
}

INITIAL {
    cai = caiBase
}

BREAKPOINT {
    SOLVE states METHOD cnexp
}

DERIVATIVE states {
    C = (cai - caiBase) / tauCa
	D = - kCa / depth * ica * (1e4)
    cai' = D - C
}