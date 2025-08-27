NEURON {
	SUFFIX h
	NONSPECIFIC_CURRENT i
	
	RANGE i, gbar
	GLOBAL eh
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 1.0 (S/cm2) 
	
}

ASSIGNED {
        eh (mV)
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	minf
	tau_m (ms)
}

STATE {	m }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m
	i = g * (v - eh)
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
	taum = 272+1499/(1+exp(-(Vm+42.2)/8.73))
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	tau_m = taum(Vm)
	UNITSOFF
	minf = 1/(1+exp((Vm+70)/6))
	UNITSON
}