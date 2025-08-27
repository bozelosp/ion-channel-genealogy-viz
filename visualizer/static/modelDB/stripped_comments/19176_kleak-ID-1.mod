INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS kleak
	RANGE gmax
	GLOBAL Erev
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	gmax	= 0.004	(umho)		
	Erev	= -100	(mV)		
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
}

INITIAL {
}

BREAKPOINT {
	i = gmax * (v - Erev)
}