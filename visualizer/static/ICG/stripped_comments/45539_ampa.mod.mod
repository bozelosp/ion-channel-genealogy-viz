NEURON {
	POINT_PROCESS AMPA  
	RANGE tau, e, i 
			
			
	NONSPECIFIC_CURRENT i

	RANGE g
	GLOBAL near_unity, gfac





}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	near_unity = 0.999 (1) 
	tau = 10 (ms) <1e-9,1e9>
	e=0	(mV)
	gfac = 1
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	tau1 (ms)







}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	tau1 = near_unity * tau
	A = 0
	B = 0
	tp = (tau1*tau)/(tau - tau1) * log(tau/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau)
	factor = 1/factor




	factor = factor * tau * exp(-1)*1(/ms)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	g = gfac*g
	i = g*(v - e)


}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau
}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}