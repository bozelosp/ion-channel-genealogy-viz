NEURON {
    POINT_PROCESS sEPSP
    RANGE A, onset, i
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
	A		= 0.0 (nA)	
	taur	= 0.3 (ms)	
	tauf	= 3.0 (ms)	
	onset	= 0.0 (ms)	
}

ASSIGNED {
    i     (nA)        
}


BREAKPOINT {

    if ((t < onset) || (t > onset+20)) {
    	i = 0
	} else { 
		i = -A*(1 - exp(-1*(t-onset)/taur)) * (exp(1 - (t-onset) / tauf))
    }
}