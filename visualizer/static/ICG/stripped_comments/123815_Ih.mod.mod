UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Ih
        
        NONSPECIFIC_CURRENT i
        RANGE gkhbar
        GLOBAL rinf, rexp, tau_r, eh
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        p = 5 (degC)
        dt (ms)
        gkhbar = 0.001385 (mho/cm2)			
        
}
 
STATE {
        r
}
 
ASSIGNED {
        eh (mV)
        i (mA/cm2)
	rinf rexp
	tau_r
}
 
BREAKPOINT {
        SOLVE deriv METHOD derivimplicit
        i = gkhbar*r*(v - eh)
}
 
INITIAL {
	rates(v)
	r = rinf
}

DERIVATIVE deriv { 
	rates(v)
	r' = (rinf - r)/tau_r
}

PROCEDURE rates(v) {  
                      
        TABLE rinf, rexp, tau_r DEPEND dt, p FROM -200
TO 100 WITH 300
	rinf = 1/(1 + exp((v+84.1)/10.2))
	rexp = 1 - exp(-dt/(tau_r))
	tau_r = 100 + 1/(exp(-17.9-0.116*v)+exp(-1.84+0.09*v))
}
 
UNITSON