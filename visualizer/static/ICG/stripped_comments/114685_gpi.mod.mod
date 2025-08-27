UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX GPi
    RANGE gpas, epas
    NONSPECIFIC_CURRENT ipas
}

PARAMETER {
    v (mV)
    dt (ms)
    gpas  = 6.8027e-5 (mho/cm2) <0,1e9>     
    epas  = -60 (mV)                        
    celsius
}

ASSIGNED { 
    ipas (mA/cm2)
}

BREAKPOINT {
    ipas = gpas*(v - epas)
}