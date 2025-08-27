NEURON {
	POINT_PROCESS Exp2Syn_depress
	RANGE tau1, tau2, e, i, tau_recover, attenuation
	NONSPECIFIC_CURRENT i

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 10 (ms) <1e-9,1e9>
	tau_recover = 100 (ms) <1e-9,1e9>
	attenuation = .8 <1e-9,1>
	e=0	(mV)
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
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS),tsyn (ms),weight_attenuate ) {

	INITIAL{
		weight_attenuate =1
		tsyn=-1000
	}
	weight_attenuate  = (weight_attenuate -1) * exp(-1/tau_recover*(t-tsyn)) +1
	A = A + weight*factor*weight_attenuate 
	B = B + weight*factor*weight_attenuate 
	tsyn=t

	
	weight_attenuate  = weight_attenuate *attenuation

}