NEURON {
	POINT_PROCESS TrigSynCD
	RANGE tau, gmax, e, i, g
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau=.1 (ms)
	gmax=0 	(uS)
	e=0	(mV)
}

ASSIGNED {
	i (nA)
	bath (uS)
	g (uS)
	k (/ms)
	v	(mV)
}

CONSTANT { exp1 = 2.7182818284590452354} 

STATE { A (uS) G (uS) }

INITIAL {
	k = 1/tau
	A = 0
	G = 0
}


BREAKPOINT {
	SOLVE state METHOD sparse
	i = G*(v - e)
	g = G
}

KINETIC state {
	~ A <-> G	(k, 0)
	~ G -> 	(k)
}

NET_RECEIVE (weight (uS)) {
    state_discontinuity(A, A + gmax*exp1)
}