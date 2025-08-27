NEURON {
	POINT_PROCESS GrCGABAexp
	RANGE tau1, tau2, e, i,tau3
	NONSPECIFIC_CURRENT i

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=  3  (ms) <1e-9,1e9> 
	tau2 = 4 (ms) <1e-9,1e9> 
	tau3 = 25 (ms)
	e=-80	(mV)

}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
}

STATE {
	A (uS)
	B (uS)
	C (uS)
}

INITIAL {


	LOCAL tp
     	if (tau1/tau2 > .9999) {
     		tau1 = .9999*tau2
     	}

     	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
     	factor = -exp(-tp/tau1) + exp(-tp/tau2)
     	factor = 1/factor


	A = 0
	B = 0
	C = 0

}

BREAKPOINT {
	SOLVE state METHOD cnexp

	g = ((B+C) - A)
        i = g*(v - e)

}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
	C' = -C/tau3

}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + 0.9*weight*factor
	C = C + 0.1*weight*factor

}