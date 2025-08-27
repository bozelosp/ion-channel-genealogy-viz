UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	SUFFIX Pass
	NONSPECIFIC_CURRENT i
	RANGE g, erev
}

PARAMETER {
	g = .001	(mho/cm2)
	erev = -70	(mV)
}

ASSIGNED { 
	i	(mA/cm2)
	v 	(mV)
}

BREAKPOINT {
	i = g*(v - erev)
        VERBATIM
        in_passiv_breakpoint();
        ENDVERBATIM
}

VERBATIM
void in_passiv_breakpoint() {}
ENDVERBATIM