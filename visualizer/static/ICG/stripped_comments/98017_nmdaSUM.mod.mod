NEURON {
	POINT_PROCESS nmdaSUM
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i

	RANGE g, s
	GLOBAL total
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tau1 	= .1	(ms)
	tau2 	= 10	(ms)
	e	= 0	(mV)
	mag     = 1     (mM)
	eta	= 3.57  (mM)
	gamma   = 0.062 (/mV)
}

ASSIGNED {
	v (mV)
	i (nA)
	g (umho)
	s
	factor
	total (umho)
}

STATE {
	A (umho)
	B (umho)
}

INITIAL {
	LOCAL tp
	total = 0
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	s = B - A
	g = s * 1(umho) /(1 + mag * exp( - (gamma * v)) / eta )
	i = g * (v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (umho)) {
	state_discontinuity(A, A + weight*factor)
	state_discontinuity(B, B + weight*factor)
	total = total+weight
}