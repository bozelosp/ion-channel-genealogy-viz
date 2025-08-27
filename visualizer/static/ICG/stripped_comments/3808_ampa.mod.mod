INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS AmpaSynapse
	RANGE onset, gmaxEPSP, e, i, g, w
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
}

PARAMETER {
	onset		(ms)
	gmaxEPSP= 0	(nS)
	w= 1				
	e= 0.0		(mV)
	v		(mV)
	celsius		(degC)
}

ASSIGNED { 
	i 		(nA)  
	g 		(nS)
	tadj
}

UNITSOFF
INITIAL {
	tadj = 3^((celsius-23.5)/10)
}

BREAKPOINT { LOCAL tt
	tt= (t-onset)*tadj
	if (t>onset && tt<1630) { 
	  g = w*gmaxEPSP * exp(-tt/18) * (1-exp(-tt/2.2))/0.68
	}
	else {g = 0}
	i = g * (v - e)
}
UNITSON