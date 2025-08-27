NEURON {
	SUFFIX kca
	NONSPECIFIC_CURRENT i
	USEION ca READ cai
	POINTER gbar
	RANGE i, Erev
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(S)  = (siemens)
}

PARAMETER {
	gbar (S/cm2) 
	Erev = -80 (mV)




	cai (mM) 
	cao	= 3	(mM)  
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	minf
	tau_m (ms)
}

STATE {	m }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^4
	i = g * (v - Erev)
}

INITIAL {
	
	rates(v)
	m = minf
}
DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
}

FUNCTION taum(Vm (mV)) (ms) {
	UNITSOFF
	taum = 90.3-75.1/(1+exp(-(Vm+46)/22.7))
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_m = taum(Vm)
	UNITSOFF

	minf = (cai/(cai+3e-3)) * (1/(1+exp(-(Vm+28.3)/12.6)))
	UNITSON
}