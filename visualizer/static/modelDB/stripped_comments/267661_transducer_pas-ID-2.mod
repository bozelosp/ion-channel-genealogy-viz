UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	SUFFIX transducer_pas
	NONSPECIFIC_CURRENT i
	RANGE g, e
}

PARAMETER {
	g = .0025	(S/cm2)	<0,1e9>
	e = -60	(mV)
}

ASSIGNED {v (mV)  i (mA/cm2)}

BREAKPOINT {
	i = g *(v - e)
}