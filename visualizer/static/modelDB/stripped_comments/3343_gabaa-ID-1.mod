INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAa_S
	POINTER pre
	RANGE C, R, R0, R1, g, gmax, TimeCount
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
	dt		(ms)
	Cmax	= 0.5	(mM)		
	Cdur	= 0.3	(ms)		
	Alpha	= 10.5	(/ms mM)	
	Beta	= 0.166	(/ms)		
	Erev	= -80	(mV)		
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
	lastrelease = -9e9
	R1=0
	TimeCount=-1
}

BREAKPOINT {
	SOLVE release
	g = gmax * R
	i = g*(v - Erev)
}

PROCEDURE release() { 
	

	TimeCount = TimeCount - dt		

						
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