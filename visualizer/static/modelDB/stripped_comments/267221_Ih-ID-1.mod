NEURON {
	SUFFIX Ih
	USEION h READ eh WRITE ih VALENCE 1
	RANGE gkhbar, ih
	RANGE slope, alpha, taumin, amp
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gkhbar = 0.0 (S/cm2)
	eh = -32.9 (mV)
	slope = 10.2 (mV)
	alpha = 84.1 (mV)
	taumin = 100 (ms)
	amp = 1 (ms)
}

ASSIGNED {
	v (mV)
	ih (mA/cm2)
	rinf (1)
	taur (ms)
}

STATE {
	r
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ih = gkhbar*r*(v - eh)
}

DERIVATIVE state {
	
	rates(v)
	r' = (rinf - r)/taur
}

INITIAL {
	rates(v)
	r = rinf
}

PROCEDURE rates(v (mV)) {
	
	
	LOCAL ar, br

	rinf = 1/(1 + exp((v + alpha)/slope))

	ar = exp(-(0.116(/mV)*v + 0.116))
	br = exp(0.09(/mV)*v - 1.84)
	taur = taumin + amp/(ar + br)
}