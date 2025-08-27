NEURON {
	SUFFIX ka
	NONSPECIFIC_CURRENT i
	POINTER gbar
	RANGE i, Erev
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar (S/cm2) 
	Erev = -80 (mV)
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	minf
	hinf
	tau_h (ms)
	tau_m (ms)
}

STATE {	m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^3 * h
	i = g * (v - Erev)
}

INITIAL {
	
	rates(v)
	m = minf
	h = hinf
}
DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
	h' = (hinf - h)/tau_h
}

FUNCTION taum(Vm (mV)) (ms) {
	UNITSOFF
	taum = 11.6 - 10.4/(1+exp(-(Vm+32.9)/15.2))
	UNITSON
}

FUNCTION tauh(Vm (mV)) (ms) {
	UNITSOFF
	tauh = 38.6 - 29.2/(1+exp(-(Vm+38.9)/26.5))
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_h = tauh(Vm)
	tau_m = taum(Vm)
	UNITSOFF
	minf = 1/(1+exp(-(Vm+27.2)/8.7))
	hinf = 1/(1+exp((Vm+56.9)/4.9))
	UNITSON
}