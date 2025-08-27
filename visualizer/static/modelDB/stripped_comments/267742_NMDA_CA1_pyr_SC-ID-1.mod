NEURON {
	POINT_PROCESS NMDA_CA1_pyr_SC
	USEION ca READ eca WRITE ica
	RANGE tau1, tau2, e, i, mg, pf
	NONSPECIFIC_CURRENT i
	RANGE g, caf
}



UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1= 3.3 (ms) <1e-9,1e9>
	tau2 = 102.38 (ms) <1e-9,1e9>
	e=0	(mV)
	mg=1    (mM)		
	pf = 0.03  (1)      
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	eca (mV)
	ica (nA)
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
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = (B - A)*mgblock(v)
	i = g*(v - e)*(1-pf)
	ica = g*(v - eca)*pf
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	
	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))

}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}