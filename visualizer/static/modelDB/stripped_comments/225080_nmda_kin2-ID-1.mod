NEURON {
	POINT_PROCESS NMDA_KIN2
	RANGE Cmax, Cdur, kon, koff, CC, CO, Beta, Alpha, Erev, Kd, gamma, mg, gmax, g, sh, B, kin_slope, kin_offset
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
	Kd          = 9.888     (mM)        
    gamma       = 0.09137   (/mV)       
    mg          = 1.0       (mM)        
	sh 			= 2.222		(mV) 		
    gmax	    = 0.002223  (umho)	    
    kin_slope   = -0.00754  (1)         
	kin_offset  = 1.302     (1)         
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
	
    B = 1. / (1. + exp(gamma * (-v+sh)) * (mg / Kd))
	if (v <= -70) {
		kin_factor = kin_slope * (-70) + kin_offset
		}
	else {
		if (v >= 40) {
			kin_factor = 1
			}
		else {
    		kin_factor = kin_slope * (v) + kin_offset
			}
		}
}