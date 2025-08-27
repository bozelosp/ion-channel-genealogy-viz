NEURON {
	POINT_PROCESS Exp2Syn_v2
	RANGE tau1, tau2, e, i, v1
	NONSPECIFIC_CURRENT i

	RANGE g
}

UNITS {
	(pA) = (picoamp)
	(mV) = (millivolt)
	(nS) = (nanosiemens)
}

PARAMETER {
	tau1 = 0.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	e=0	(mV)
}

ASSIGNED {
	v (mV)
	v1 (mV)
	i (pA)
	g (nS)
	factor
}

STATE {
	A (nS)
	B (nS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > 0.9999) {
		tau1 = 0.9999*tau2
	}
	if (tau1/tau2 < 1e-9) {
		tau1 = tau2*1e-9
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor

}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
	v1 = v
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (nS)) {
	A = A + weight*factor
	B = B + weight*factor
}