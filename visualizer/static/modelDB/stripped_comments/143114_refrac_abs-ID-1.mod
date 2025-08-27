NEURON {
	POINT_PROCESS Refrac_abs
	RANGE gr, e, tr, g
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tr = 2 (ms) 
	gr = 1 (umho) 
	e = -65 (mV) 
}

ASSIGNED {
	g (umho)
	v (mV)
	i (nA)
}

INITIAL {
	g = 0 
}


BREAKPOINT {
	i = g*(v - e)
}

NET_RECEIVE(w) { 
	if (flag == 0) { 
		g = gr
		net_send(tr, 1)
	}else{ 
		g = 0
	}
}