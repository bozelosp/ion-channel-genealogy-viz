NEURON {
	POINT_PROCESS AMPA_KIN
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
	Cdur	= 0.3       (ms)		    
	kon     = 12.88     (/ms/mM)        
    koff    = 6.47      (/ms)           
    CC      = 69.97     (/ms)           
    CO      = 6.16      (/ms)           
    Beta	= 100.63    (/ms)	        
    Alpha   = 173.04    (/ms)           
	Erev	= 0.        (mV)		    
	gmax	= 0.001     (umho)	        
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