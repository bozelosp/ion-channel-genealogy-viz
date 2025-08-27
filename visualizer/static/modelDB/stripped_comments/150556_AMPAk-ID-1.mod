NEURON {
	POINT_PROCESS AMPAk
	RANGE tauon, tauoff, gAmax, gA, Erev, i
	NONSPECIFIC_CURRENT i
	GLOBAL total
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(pS) = (picosiemens)
}

PARAMETER {
	Erev	= 0    (mV)	
	gAmax	= 30  (pS)	
	tauon	= 1.1  (ms)<1e-9,1e9>
	tauoff	= 5.75 (ms)<1e-9,1e9>
}

ASSIGNED {
	v (mV)
	i (nA)
	gA (uS)
	factor
	total (uS)
}

STATE {
	m (uS)
	h (uS)
}

INITIAL {
	LOCAL tp
	total = 0
	if (tauon/tauoff > .9999) {
		tauon = .9999*tauoff
	}
	m = 0
	h = 0
	tp = (tauon*tauoff)/(tauoff - tauon) * log(tauoff/tauon)
	factor = -exp(-tp/tauon) + exp(-tp/tauoff)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gA = (1e-6)*gAmax*(h-m) 	
        i = gA*(v - Erev)
}	

DERIVATIVE state {
	m' = -m/tauon
	h' = -h/tauoff
}

NET_RECEIVE(weight (uS)) {
	state_discontinuity(m, m + weight*factor)
	state_discontinuity(h, h + weight*factor)
	total = total+weight
}