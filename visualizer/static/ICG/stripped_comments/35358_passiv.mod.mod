UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT { v FROM -100 TO 50 WITH 50	(mV) }

NEURON {
	SUFFIX Pass
	NONSPECIFIC_CURRENT i
	RANGE g,erev
}

PARAMETER {
	g = .001	(mho/cm2)
	erev = -70	(mV)
}

ASSIGNED { i	(mA/cm2)}

BREAKPOINT {
	i = g*(v - erev)
}