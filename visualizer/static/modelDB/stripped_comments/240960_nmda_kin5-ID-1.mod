NEURON {
	POINT_PROCESS NMDA_KIN5
	RANGE Cmax, Cdur, kon, koff, CC, CO, Beta, Alpha, Erev, Kd, gamma, mg, gmax, g, B, kin_scale, Ro
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax        = 1.        (mM)        
	Cdur	    = 0.3       (ms)		
	kon         = 86.89     (/mM/ms)    
    koff        = 0.69      (/ms)       
    CC          = 9.64      (/ms)       
    CO          = 2.60      (/ms)       
    Beta	    = 0.68      (/ms)	    
    Alpha       = 0.079     (/ms)       
	Erev	    = 0.        (mV)		
	Kd          = 9.98      (mM)        
    gamma       = 0.101     (/mV)       
    mg          = 1.0       (mM)        
    gmax	    = 0.003026  (umho)	    
    kin_scale   = 1.83      (1)         
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	C       (mM)        
    B                   
    kin_factor          
    scale               
}

STATE {
    Ru                  
    Rb                  
    Rc                  
    Ro                  
}

INITIAL {
	C = 0.
    Ru = 1.
    Rb = 0.
    Rc = 0.
    Ro = 0.
    scale = 1.
    mgblock(v)
}

BREAKPOINT {
    mgblock(v)
    SOLVE kstates METHOD sparse
	g = scale * gmax * B * Ro
    i = g * (v - Erev)
}

KINETIC kstates {
    ~ Ru <-> Rb     (C * kon, koff)
    ~ Rb <-> Rc     (CC, CO * kin_factor)
    ~ Rc <-> Ro     (Beta, Alpha * kin_factor)
}

NET_RECEIVE(weight) {
    if (flag == 0) { 
        C = Cmax
        scale = weight
        net_send(Cdur, 1)
    } else {    
        C = 0.
    }
}

PROCEDURE mgblock(v(mV)) {
	
    B = 1. / (1. + exp(gamma * (-v)) * (mg / Kd))
    kin_factor = (1 - kin_scale) * B + kin_scale
}