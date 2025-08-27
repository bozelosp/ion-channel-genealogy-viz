NEURON {
	SUFFIX nmda
	NONSPECIFIC_CURRENT i
	RANGE i, Erev, gbar
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 2570e-6	(S/cm2) < 0, 1e9 >
	Erev = 0 (mV)
	mg_o = 2 (mM)	
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
	g = gbar * n
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
	alphan = 3.0 * exp(0.035 * (Vm + 10.0))
	UNITSON
}

FUNCTION betan(Vm (mV)) (/ms) {
	UNITSOFF
	betan = mg_o * exp(-0.035 * (Vm + 10.0))
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