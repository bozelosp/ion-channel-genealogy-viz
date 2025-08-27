NEURON {

	POINT_PROCESS Exp2SynNMDAperm
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i

	RANGE gmax
}

UNITS {

	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1 = 0 (ms)
	tau2 = 0 (ms)
	gmax = 0 (uS)
	e    = 0 (mV)
}

ASSIGNED {

	v       (mV)
	i       (nA)
	mgblock (1)
}

BREAKPOINT {
	    
	mgblock = 1.0 / (1.0 + 0.28 * exp(-0.062(/mV) * v) )
	i = gmax * (v - e) * mgblock
}

NET_RECEIVE(w) {}