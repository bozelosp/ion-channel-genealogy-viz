NEURON {
	POINT_PROCESS EPSC
	RANGE tau1, tau2, i, inward, imax
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
}

PARAMETER {
	tau1 = 0.2   (ms)    <1e-9,1e9>
	tau2 = 5.0    (ms)   <1e-9,1e9>
    inward = -1.0                     
    imax = 1.0           <0,1>        
}

ASSIGNED {
	v (mV)
	i (nA)
	factor
}

STATE {
	A (nA)
	B (nA)
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
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = inward * imax * (B - A)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (nA)) {
	A = A + weight*factor
    B = B + weight*factor
}