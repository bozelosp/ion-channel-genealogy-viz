NEURON {
	POINT_PROCESS GABA_A_KIN
	RANGE Cmax, Cdur, kon, koff, CC, CO, Beta, Alpha, Erev, gmax, g, Ro
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax    = 1.        (mM)    	    
	Cdur	= 0.5       (ms)		    
	kon     = 5.397     (/ms/mM)        
    koff    = 4.433     (/ms)           
    CC      = 20.945    (/ms)           
    CO      = 1.233     (/ms)           
    Beta	= 283.090   (/ms)	        
    Alpha   = 254.520   (/ms)           
	Erev	= -73.0     (mV)		    
	gmax	= 0.000590  (umho)	        
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	C                   
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
}

BREAKPOINT {
    SOLVE kstates METHOD sparse
	g = scale * gmax * Ro
    i = g * (v - Erev)
}

KINETIC kstates {
    ~ Ru <-> Rb     (C * kon, koff)
    ~ Rb <-> Rc     (CC, CO)
    ~ Rc <-> Ro     (Beta, Alpha)
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