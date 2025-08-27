NEURON {
	POINT_PROCESS NMDA 
	RANGE tau1, tau2, e, i, weight, A, B
	NONSPECIFIC_CURRENT i

	RANGE g, tp, factor, cond_scale, offset
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1= 0.5 (ms) 
	tau2 = 44 (ms) 
	e=0	(mV)
	mg=1    (mM)		
        weight = 0.5e-3 (uS) 
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor (1)
        rescale_epsps (1)
        cond_scale (1)
        tp (ms)
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
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
	g = B - A
	i = g*mgblock(v)*(v - e)
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















FUNCTION_TABLE t_table(g) (ms)  

NET_RECEIVE(weight (uS)) {
        LOCAL t_tmp
        t_tmp = t_table(g/weight)
        A = weight * factor * exp( -t_tmp/tau1)
        B = weight * factor * exp( -t_tmp/tau2)
}