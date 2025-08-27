NEURON {
    SUFFIX ConeIh
    NONSPECIFIC_CURRENT ih
    RANGE ghbar, gh, eh 
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mS) = (millimho)
}

PARAMETER {
    
    ghbar = 3.5 (mS/cm2) <0,1e9>
    eh  = -32.5 (mV)
}

STATE {
    nh
}

ASSIGNED {
    v (mV)
    
    ih (mA/cm2)
    
    infh
    tauh   (ms)
    
    gh (mho/cm2)
}

INITIAL {
    rate(v)
    nh  = infh
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gh  = (0.001)*ghbar*(1-(1+3*nh)*(1-nh)^3)
    ih  = gh*(v - eh)
    
    
}

DERIVATIVE states {
    rate(v)
    nh'  = (infh  - nh )/tauh
}

FUNCTION alphah(v(mV)) (/ms) { 
    alphah = 0.001*18/( exp  (  ( v+88)/12 ) + 1 )
}

FUNCTION betah(v(mV)) (/ms) { 
    betah = 0.001*18/( exp  ( - ( v+18)/19 ) + 1 )
}

PROCEDURE rate(v (mV)) { LOCAL a, b
    a = alphah(v)
    b = betah(v)
    tauh = 1/(a + b)
    infh = a/(a + b)
}