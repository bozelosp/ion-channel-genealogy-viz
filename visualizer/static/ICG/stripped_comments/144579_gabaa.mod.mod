NEURON {
	POINT_PROCESS GABAA
	RANGE tau, e, i
	NONSPECIFIC_CURRENT i
	GLOBAL gfac




}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau = 0.1 (ms) <1e-9,1e9>
	e = 0	(mV)
	gfac = 1
}

ASSIGNED {
	v (mV)
	i (nA)






}

STATE {
	g (uS)
}

INITIAL {
	g=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = gfac*g*(v - e)


}

DERIVATIVE state {
	g' = -g/tau
}

NET_RECEIVE(weight (uS)) {
	g = g + weight
}