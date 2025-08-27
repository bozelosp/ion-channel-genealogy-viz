INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

NEURON { 
        SUFFIX kv
        USEION k READ ek WRITE ik
        RANGE gbar, ik, g
}

UNITS {
        (S) = (siemens)
        (mV) = (millivolt) 
        (mA) = (milliamp) 
}
 
PARAMETER { 
        gbar = 0.0      (mho/cm2)
}
 
ASSIGNED { 
        ek      (mV)
        v       (mV)
        ik      (mA/cm2)
        g       (S/cm2)
        minf
        mtau    (ms)
}
 
STATE {
    m
}

BREAKPOINT { 
        SOLVE states METHOD cnexp 
        g = gbar * m * m * m * m
        ik = g * ( v - ek ) 
}
 
INITIAL { 
        settables(v)
        m = minf
} 

DERIVATIVE states { 
        settables(v) 
        m' = (minf - m ) / mtau
}











PROCEDURE settables(v (mV)) { 
UNITSOFF
        TABLE minf, mtau FROM -120 TO 40 WITH 641
        minf  = 1.0 / (1 + exp((-37.6 - v) / 27.24))
        mtau = (3.53 - 1.85) / (1 + exp((v - 45) / 13.71)) + 1.85
UNITSON
}