NEURON {
	POINT_PROCESS Exp2EPSG_NMDA
	RANGE tau1, tau2, e, gmax, Kd, gamma, mg, M, i
	NONSPECIFIC_CURRENT i

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1    = .1    (ms) <1e-9,1e9>
	tau2    = 10    (ms) <1e-9,1e9>
    Kd      = 9.98  (mM)    
    gamma   = 0.101 (/mV)   
    mg      = 1.0   (mM)    
    e       = 0	    (mV)    
    gmax    = 1     (1)     
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
    M       
	factor
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
    mgblock(v)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
    mgblock(v)
	g = B - A
	i = g*gmax*M*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}

PROCEDURE mgblock(v(mV)) {
	
    M = 1. / (1. + exp(gamma * (-v)) * (mg / Kd))
}