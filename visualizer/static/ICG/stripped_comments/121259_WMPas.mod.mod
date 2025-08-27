UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	SUFFIX WMPas
	NONSPECIFIC_CURRENT i
	RANGE g, eleak
}

PARAMETER {
	v 		(mV)
	g = .0004	(mho/cm2)	<0,1e9>
	eleak = -65	(mV)
}

ASSIGNED { 
	i	(mA/cm2)
}

BREAKPOINT {
	i = g*(v - eleak)
}