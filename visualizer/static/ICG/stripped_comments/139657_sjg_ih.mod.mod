UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (nA) = (nanoamp)
}

NEURON {
    SUFFIX sjg_ih
    NONSPECIFIC_CURRENT i
    RANGE ghbar, gh
    GLOBAL uinf, utau, eh
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
    v (mV)
    celsius = 22 (degC)
    dt (ms)
    ghbar = 0.00318 (mho/cm2) <0,1e9>
    
}

STATE {
    u
}

ASSIGNED {
    eh (mV)
    gh (mho/cm2)
    i (mA/cm2)
    uinf
    utau (ms)
}

LOCAL uexp

BREAKPOINT {
    SOLVE states
    
    gh = ghbar*u
    i = gh*(v - eh)
}

UNITSOFF

INITIAL {
    trates(v)
    u = uinf
}

PROCEDURE states() {  
    trates(v)      
    u = u + uexp*(uinf-u)
VERBATIM
    return 0;
ENDVERBATIM
}

LOCAL q10
PROCEDURE rates(v) {  
                      
		      
    q10 = 3^((celsius - 22)/10)
    uinf = 1 / (1+exp((v + 101) / 11))
    utau = (10000 / (235.55*exp(0.0782*(v+23.76)) + 0.33*exp(-0.0614*(v+23.76)))) + 154.57

}

PROCEDURE trates(v) {  
                      
    LOCAL tinc
    TABLE uinf, uexp
    DEPEND dt, celsius FROM -200 TO 150 WITH 350

    rates(v)    
        
        

    tinc = -dt * q10
    uexp = 1 - exp(tinc/utau)
}

FUNCTION vtrap(x,y) {  
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/(exp(x/y) - 1)
    }
}

UNITSON