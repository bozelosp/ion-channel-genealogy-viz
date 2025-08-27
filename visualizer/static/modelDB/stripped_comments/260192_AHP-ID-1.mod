NEURON {
	POINT_PROCESS AHP
	USEION k READ ek 
	RANGE gr, e,  g, tau
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(uS)   = (microsiemens)
}

PARAMETER {
	tau=5 (ms)   
	gr = 300e-6 (uS) 
	e = -90 (mV) 
}

ASSIGNED {
	
	v (mV)
	i (nA)
	ek (mV)
}

STATE {
	g (uS)
}


INITIAL {
	g = 0 
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	i = g*(v - ek)
}

DERIVATIVE state {
	g' = -g/tau
}


NET_RECEIVE(weight (uS)) { 
	if (flag == 0) { 
		g = g+weight   
	}
}