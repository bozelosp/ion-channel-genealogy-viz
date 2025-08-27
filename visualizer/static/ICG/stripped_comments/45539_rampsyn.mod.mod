NEURON {
	POINT_PROCESS RampSyn
	RANGE time_interval, e, i, weight, saturation_fact, k
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	time_interval = 5 (ms) <1e-9,1e9>
	e = 0	(mV)
	weight = 2.5e-5 (uS)	
			 	
	saturation_fact=1e10 (1) 
		
		
}

ASSIGNED {
	v (mV)
	i (nA)
	k (uS/ms)
}

STATE {
	g (uS)
}

INITIAL {
	g=0
	k=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	if (g > saturation_fact * weight) { g = saturation_fact * weight }
	i = g*(v - e)
}

DERIVATIVE state {
	g' = k
}

NET_RECEIVE(weight (uS)) {
	if (flag>=1) {
		
		k = k - weight/time_interval
		g = g - weight
	} else {
		
		net_send(time_interval, 1) 
		k = k + weight/time_interval
	}
}