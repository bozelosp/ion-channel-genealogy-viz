NEURON {
	POINT_PROCESS Refrac_rel
	RANGE gr, e,  g, tau
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tau=1 (ms)   
	gr = 1 (umho) 
	e = -65 (mV) 
}

ASSIGNED {
	
	v (mV)
	i (nA)
}

STATE {
	g (umho)
}


INITIAL {
	g = 0 
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g*(v - e)
}

DERIVATIVE state {
	g' = -g/tau
}


NET_RECEIVE(w) { 
	if (flag == 0) { 
		g = g+gr   
	}
}