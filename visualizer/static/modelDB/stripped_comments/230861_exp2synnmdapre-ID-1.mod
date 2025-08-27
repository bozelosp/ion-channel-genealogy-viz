NEURON {

	POINT_PROCESS Exp2SynNMDApre
	RANGE e, i
	NONSPECIFIC_CURRENT i
	POINTER pre
	
	RANGE gmax
}

UNITS {

	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {

	e    = 0   (mV)
	gmax = 0   (uS)
	Smax = 125 (1)
}

ASSIGNED {

    pre (mV)
	v (mV)
	i (nA)
	mgblock
}

STATE { S }

INITIAL { S = 0 }

BREAKPOINT {

	mgblock = 1 / (1 + 0.28 * exp(-0.062(/mV) * v) )

	if (S>=Smax) {S=Smax} else {SOLVE state METHOD cnexp}

	i = gmax * S * (v - e) * mgblock
}

DERIVATIVE state { 

    S' = H(pre)/1(ms)-S/150(ms) 
}

FUNCTION H(x(mV)) {

    if (x>-50) {
    
        H = 1
        
    }else{
        
        H = 0
        
    }
}