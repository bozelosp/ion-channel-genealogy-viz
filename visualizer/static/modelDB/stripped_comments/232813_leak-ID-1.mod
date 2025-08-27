UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE g, erev
}

PARAMETER {
	g = .001	(mho/cm2)
	erev = -85	(mV)
}

ASSIGNED {
i(mA/cm2)
v(mV)
}

BREAKPOINT {
	i = g*(v - erev)
}