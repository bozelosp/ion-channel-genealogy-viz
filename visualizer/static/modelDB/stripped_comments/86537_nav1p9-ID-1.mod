NEURON {
	SUFFIX nav1p9
	NONSPECIFIC_CURRENT i
	RANGE gbar, ena, slow_inact, m, h, s, gate
	RANGE tau_m, tau_h, tau_s
	
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.0069005 (S/cm2)
	ena=62.94 (mV)

	A_am9 = 1.032 (/ms)
	B_am9 = 6.99 (mV) 
	C_am9 = -14.87115

	A_ah9 = 0.06435 (/ms) 
	B_ah9 = 73.26415 (mV) 
	C_ah9 = 3.71928 (mV) 

	A_as9 = 0.00000016 (/ms) 
	B_as9 = 0 (mV)
	gate = 0 (mV)
	C_as9 = 12 (mV)

	A_bs9 = 0.0005 (/ms)
	B_bs9 = 32 (mV)
	C_bs9 = 23 (mV)

	A_bm9 = 5.79 (/ms) 
	B_bm9 = 130.4 (mV) 
	C_bm9 = 22.9 (mV) 

	A_bh9 = 0.13496 
	B_bh9 = 10.27853 
	C_bh9 = -9.09334 
	
	slow_inact = 1 (1) 
}

ASSIGNED {
	v	(mV) 
	i	(mA/cm2)
	g	(S/cm2)
	tau_h	(ms)
	tau_m	(ms)
	tau_s	(ms)
	minf
	hinf
	sinf
}

STATE { m h s}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m * h * s
	i = g * (v-ena)
}

INITIAL {
	rates(v) 
	
	m = minf
	h = hinf
	s = sinf
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
	h' = (hinf - h)/tau_h
	s' = (sinf - s)/tau_s
}

FUNCTION alpham(Vm (mV)) (/ms) {
	alpham=A_am9/(1+exp((Vm+B_am9)/C_am9))
}

FUNCTION alphah(Vm (mV)) (/ms) {
	alphah=A_ah9/(1+exp((Vm+B_ah9)/C_ah9))
}

FUNCTION alphas(Vm (mV)) (/ms) {
	alphas=A_as9*exp(-(Vm+gate+B_as9)/C_as9)
}


FUNCTION betam(Vm (mV)) (/ms) {
	betam=A_bm9/(1+exp((Vm+B_bm9)/C_bm9))
}

FUNCTION betah(Vm (mV)) (/ms) {
	betah=A_bh9/(1+exp((Vm+B_bh9)/C_bh9))
}

FUNCTION betas(Vm (mV)) (/ms) {
	betas=A_bs9/(1+exp(-(Vm+gate+B_bs9)/C_bs9))
}

FUNCTION rates(Vm (mV)) (/ms) {
	tau_m = 1.0 / (alpham(Vm) + betam(Vm))
	minf = alpham(Vm) * tau_m

	tau_h = 1.0 / (alphah(Vm) + betah(Vm))
	hinf = alphah(Vm) * tau_h

	if (slow_inact) {
		tau_s = 1.0 / (alphas(Vm) + betas(Vm))
		sinf = alphas(Vm) * tau_s
	} else {
		tau_s = 0.1	
		sinf = 1.0	
	}
}