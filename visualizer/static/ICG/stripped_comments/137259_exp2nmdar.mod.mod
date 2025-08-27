NEURON {
	POINT_PROCESS Exp2NMDAR
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i

	RANGE g
	GLOBAL total

	RANGE mg
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=4 (ms) <1e-9,1e9>
	tau2=42 (ms) <1e-9,1e9>
	e=0	(mV)
	mg=1 (mM) 
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	total (uS)
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	total = 0
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor

	mgblock(v)
}

BREAKPOINT {
	SOLVE state METHOD cnexp

	g = B - A
	i = g*mgblock(v)*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	state_discontinuity(A, A + weight*factor)
	state_discontinuity(B, B + weight*factor)
	total = total+weight
}


FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}