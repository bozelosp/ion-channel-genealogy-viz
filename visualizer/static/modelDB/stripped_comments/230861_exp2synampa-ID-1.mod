NEURON {
	POINT_PROCESS Exp2SynAMPA
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i

	RANGE S, total
	GLOBAL c
}

UNITS {

	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {

	tau1  =   1 (ms) <1e-9,1e9>
	tau2  =   2 (ms) <1e-9,1e9>
	e     = 0	(mV)
	c     = 16.5
}

ASSIGNED {

	v      (mV)
	i      (nA)
	S      (uS)
	factor (1)
	total  (1)
}

STATE {

	A (uS)
	B (uS)
}

INITIAL {

	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
    tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	total = 0
}

BREAKPOINT {

	SOLVE state METHOD cnexp
	S = B - A
	
	i = S * (v - e)
}

DERIVATIVE state {

	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(w (uS) ) {

	
	
	
	state_discontinuity(A, A + factor*c*w)
	state_discontinuity(B, B + factor*c*w)
	total = total + 1
}