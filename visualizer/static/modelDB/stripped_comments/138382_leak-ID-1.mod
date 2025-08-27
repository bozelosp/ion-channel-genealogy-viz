NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE i, e, gbar
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
}

PARAMETER {
	gbar = 9e-5 (S/cm2)  
	e = -61 (mV)
}

ASSIGNED {
	i (mA/cm2)
	v (mV)
}

BREAKPOINT {
	i = gbar*(v - e)
}