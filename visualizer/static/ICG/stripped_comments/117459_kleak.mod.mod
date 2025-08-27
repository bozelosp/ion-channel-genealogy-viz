INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	
	SUFFIX kleak
	USEION k READ ek WRITE ik
	RANGE gmax
	
	
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	gmax	= 0.004	(umho)		
	
}


ASSIGNED {
	ek (mV)
	v		(mV)		
	ik 		(nA)		
}

INITIAL {
}

BREAKPOINT {
	ik = gmax * (v - ek)
}