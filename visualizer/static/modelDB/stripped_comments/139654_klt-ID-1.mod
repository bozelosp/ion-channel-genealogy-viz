UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX klt
        USEION k READ ek WRITE ik
        RANGE gbar, g, ik
        GLOBAL winf, wtau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        ek = -90 (mV)
        gbar = 0 (mho/cm2) <0,1e9>
        tfac = 0.5
}

STATE {
        w
}

ASSIGNED {
    ik (mA/cm2) 
    g (mho/cm2)
    winf 
    wtau (ms)
    }

BREAKPOINT {
        SOLVE states METHOD cnexp 

        g = gbar*(w^4)
        ik = g*(v - ek)
}
 
 
INITIAL {
    rates(v)
    w = winf
}

DERIVATIVE states {  
        rates(v)
        w' =  (winf-w)/wtau
}


UNITSOFF

PROCEDURE rates(v) {  
                      

    winf = (1 / (1 + exp(-(v + 48) / 6)))^0.25
    wtau =  tfac*((100 / (6*exp((v+60) / 6) + 16*exp(-(v+60) / 45))) + 1.5)
}

UNITSON