NEURON {
	POINT_PROCESS NMDAk
	RANGE tauon, tauoff, gNmax, gN, Erev, i,mg
	NONSPECIFIC_CURRENT i
	GLOBAL total,vmin, vmax
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(pS) = (picosiemens)
	(mM) = (milli/liter)
}

PARAMETER {
	Erev	= 0    (mV)	
	gNmax	= 60  (pS)	
	tauon	= 2.23  (ms)<1e-9,1e9>	
	tauoff	= 75.68 (ms)<1e-9,1e9> 	
	mg	= 1.2    (mM)	
	vmin = -120	(mV)
	vmax = 100	(mV)
}

ASSIGNED {
	v (mV)
	i (nA)
	gN (uS)
	factor
	total (uS)
}

STATE {
	m (uS)
	h (uS)
	B		
}

INITIAL {
	LOCAL tp
	rates(v)
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
	rates(v)
	SOLVE state METHOD cnexp
	gN = (1e-6)*gNmax*(h-m) 	
        i = gN*B*(v - Erev)
	
}

DERIVATIVE state {
	m' = -m/tauon
	h' = -h/tauoff
}

PROCEDURE rates(v(mV)) {
	TABLE B
	DEPEND mg
	FROM vmin TO vmax WITH 200

	

	B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}


NET_RECEIVE(weight (uS)) {
	state_discontinuity(m, m + weight*factor)
	state_discontinuity(h, h + weight*factor)
	total = total+weight
}