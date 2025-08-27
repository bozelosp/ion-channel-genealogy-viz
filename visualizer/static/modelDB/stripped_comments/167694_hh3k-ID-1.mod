NEURON {
    SUFFIX hh3k
    USEION k READ ek WRITE ik
    RANGE gbar
    GLOBAL vmin, vmax
    GLOBAL taun, pown
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    v (mV)
    gbar = .12 (mho/cm2)
    vmin = -120 (mV)
    vmax = 100 (mV)
    ek (mV)
    taun = 2 (ms)
    pown = 2
}

STATE {    
    n   <1e-1> 
}

ASSIGNED {
    ik (mA/cm2)
}

INITIAL {
    n = ssn(v)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gbar*(v - ek)*n^pown
}

DERIVATIVE states {
    n' = (ssn(v) - n)/taun
}

FUNCTION ssn(v(mV)) {  
    TABLE FROM vmin TO vmax WITH 200
    ssn = 1/(1 + exp((v + 40 (mV))/-3 (mV)))
}