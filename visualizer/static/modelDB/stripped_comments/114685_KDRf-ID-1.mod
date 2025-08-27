UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX KDRf
    USEION k READ ki,ek WRITE ik
    RANGE gk
    GLOBAL rate_k,gmax_k
}

PARAMETER {
    v             (mV)
    dt            (ms)
    gk = 0.0038   (mho/cm2)  
    ek
    ki
    celsius
}

STATE {
    pfast 
}

ASSIGNED { 
    ik (mA/cm2)
    pinf
    ptau (ms)
    rate_k
    gmax_k
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = (gk*gmax_k)*pfast*(v-ek)
}

UNITSOFF

INITIAL {
    rate_k = 2.05
    gmax_k = 2.05
    settables(v)
    pfast = pinf
}

DERIVATIVE states {  
    settables(v)
    pfast' = (pinf-pfast)/ptau
}

PROCEDURE settables(v) {
	
    TABLE pinf, ptau DEPEND celsius FROM -100 TO 100 WITH 400
    
    pinf = 1/(1+exp(-(v+16.2)/8.6))                              
    ptau = 6.7/(exp(-(v+21.7)/21.2)+exp((v-11.7)/21.2))/rate_k   

}

UNITSON