UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX Ih
    NONSPECIFIC_CURRENT i
    RANGE gh
    GLOBAL eh
}

PARAMETER {
    v (mV)
    dt (ms)
    gh    = 0.001 (mho/cm2)
    
    celsius
}

STATE {
    f
}

ASSIGNED { 
    eh (mV)
    i (mA/cm2)
    finf  
    ftau (ms)
}

BREAKPOINT {
    SOLVE integrate METHOD cnexp
    i = gh*f*(v-eh)
}

UNITSOFF

INITIAL {
    setinf(v)
    f = finf
}

DERIVATIVE integrate {
    setinf(v)
    f' = (finf - f)/ftau
}

PROCEDURE setinf(v) {
    
    LOCAL vhalf, sfactor
    TABLE finf, ftau DEPEND celsius FROM -100 TO 100 WITH 400

    vhalf = -75     
    sfactor = 5.5   
    finf = 1.0/(1+exp((v-vhalf)/sfactor))
    ftau = (1.0/(exp(-14.59-0.086*v)+exp(-1.87+0.0701*v)))
}

UNITSON