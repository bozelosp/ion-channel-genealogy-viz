UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Ih
        
        NONSPECIFIC_CURRENT i
        RANGE gkhbar
        GLOBAL rinf, tau_r,eh
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
	rinf 
	tau_r	(ms)
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        i = gkhbar*r*(v - eh)
}
 
INITIAL {
	rates(v)
	r = rinf
}

DERIVATIVE state { 
	rates(v)
	r' = (rinf - r)/tau_r
}

PROCEDURE rates(v(mV)) {  
                      

	rinf = 1/(1 + exp((v+84.1(mV))/10.2(mV)))
	tau_r = 100(ms) + 1(ms)/(exp(-17.9-0.116(/mV)*v)+exp(-1.84+0.09(/mV)*v))

}
 
UNITSON