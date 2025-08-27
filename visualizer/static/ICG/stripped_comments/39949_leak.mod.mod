NEURON {
	SUFFIX leak
	NONSPECIFIC_CURRENT i
	RANGE g, i, iL
	RANGE e
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	g = 0.008 (mS/cm2)	<0,1e9>
	e = -90	(mV)		
}

ASSIGNED {
	v	(mV)
	iL	(uA/cm2)	
	i	(mA/cm2)
}

BREAKPOINT {
	iL = g*(v - e)
	i = (0.001)*iL
}