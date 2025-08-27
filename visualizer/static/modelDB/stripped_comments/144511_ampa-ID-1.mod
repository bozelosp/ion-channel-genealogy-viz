INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS AMPA
	POINTER pre
	RANGE C, R, R0, R1, g, gmax, lastrelease, TimeCount
	NONSPECIFIC_CURRENT i
	RANGE Cmax, Cdur, Alpha, Beta, Erev, Prethresh, Deadtime, Rinf, Rtau
    THREADSAFE
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	dt		(ms)
	Cmax	= 1	(mM)		
	Cdur	= 1	(ms)		
	Alpha	= 1.1	(/ms mM)	
	Beta	= 0.19	(/ms)		
	Erev	= 0	(mV)		
	Prethresh = 0 			
	Deadtime = 1	(ms)		
	gmax		(umho)		
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
	TimeCount	(ms)		
}

INITIAL {
	R = 0
	C = 0
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / ((Alpha * Cmax) + Beta)
	lastrelease = -1000
	R1=0
	TimeCount=-1
}

BREAKPOINT {
	SOLVE release
	g = gmax * R
	i = g*(v - Erev)
}

PROCEDURE release() {
	

	TimeCount=TimeCount-dt			

						
	if (TimeCount < -Deadtime) {
		if (pre > Prethresh) {		
			C = Cmax			
			R0 = R
			lastrelease = t
			TimeCount=Cdur
		}
						
	} else if (TimeCount > 0) {		
	
		
	
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