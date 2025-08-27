INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDA
	POINTER pre
	RANGE C, R, R0, R1, g, gmax, lastrelease
	NONSPECIFIC_CURRENT i
	GLOBAL Cmax, Cdur, Alpha, Beta, Erev, Prethresh, Deadtime, Rinf, Rtau
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax	= 1	(mM)		
	Cdur	= 1.1	(ms)		
	Alpha	= 10	(/ms mM)	
	Beta	= 0.0125 (/ms)		
	Erev	= 0	(mV)		
	Prethresh = 0 			
	Deadtime = 0	(ms)		
	gmax		(umho)		
	eta     = 0.33  (/mM)
	mag     = 1     (mM)
	gamma   = 0.06  (/mV)
}

ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	C		(mM)		
	R				
	R0				
	R1				
	Rinf				
	Rtau		(ms)		
	pre 				
	lastrelease	(ms)		
}

INITIAL {
	R = 0
	C = 0
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / ((Alpha * Cmax) + Beta)
	lastrelease = -9e9
}

BREAKPOINT {
	SOLVE release
	g = (gmax * R)/(1 + eta * mag * exp( - (gamma * v)))
	i = g*(v - Erev)
}

PROCEDURE release() { LOCAL q
	

	q = ((t - lastrelease) - Cdur)		

						
	if (q > Deadtime) {
		if (pre > Prethresh) {		
			C = Cmax			
			R0 = R
			lastrelease = t
		}
						
	} else if (q < 0) {			
	
		
	
	} else if (C == Cmax) {			
		R1 = R
		C = 0.
	}



	if (C > 0) {				

	   R = Rinf + (R0 - Rinf) * exptable (- (t - lastrelease) / Rtau)
				
	} else {				

  	   R = R1 * exptable (- Beta * (t - (lastrelease + Cdur)))
	}

	VERBATIM
	return 0;
	ENDVERBATIM
}

FUNCTION exptable(x) { 
	TABLE  FROM -10 TO 10 WITH 2000

	if ((x > -10) && (x < 10)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}