NEURON {
	POINT_PROCESS Exp2SynAMPApre
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

	e    = 0 (mV)
	gmax = 1 (uS)
}

ASSIGNED {

    pre (mV)
	v   (mV)
	i   (nA)
}

STATE { W }

INITIAL { W = 0 }

BREAKPOINT {

	SOLVE state METHOD cnexp
	
	i = gmax * W * (v - e)
}

DERIVATIVE state {

    W' = H(pre)/1(ms)-W/2(ms)    
}

FUNCTION H(x(mV)) {

    if (x>-40) {
    
        H = 1
        
    }else{
        
        H = 0
        
    }
}