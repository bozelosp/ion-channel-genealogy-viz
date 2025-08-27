NEURON {
	SUFFIX bkfast
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar
	
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 37530e-6	(S/cm2) < 0, 1e9 >
	Erev = -82 (mV)
	
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	ninf
	tau_n (ms)
}

STATE {	n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * n*n
	i = g * (v - Erev)
}

INITIAL {
	
	n = alphan(v)/(alphan(v) + betan(v))
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/tau_n
}

FUNCTION alphan(Vm (mV)) (/ms) {
	UNITSOFF
	alphan = 5.82 /(1 + exp( -0.125 * (Vm + 3.3)))
	UNITSON
}

FUNCTION betan(Vm (mV)) (/ms) {
	UNITSOFF
	betan =  2.413 / (1 + exp( 0.0675 * (Vm + 46.35)))
	UNITSON
}

FUNCTION taun(Vm (mV)) (/ms) {
	UNITSOFF
	taun = 1.0 / (alphan(Vm) + betan(Vm))
	
	
	
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_n = taun(Vm)
	
	ninf = alphan(Vm) * tau_n	
}