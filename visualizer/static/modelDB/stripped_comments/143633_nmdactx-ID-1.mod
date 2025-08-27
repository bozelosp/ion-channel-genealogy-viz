INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDActx
	POINTER pre
	RANGE B, C, R, R0, R1, g, gmax, lastrelease, F, F0, F1
	NONSPECIFIC_CURRENT i
	GLOBAL Cmax, Cdur, Alpha, Beta, Erev, mg
	GLOBAL Prethresh, Deadtime, Rinf, Rtau
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax	= 1	(mM)		
	Cdur	= 1	(ms)		
	Alpha	= 0.072	(/ms mM)	
	Beta	= 0.0066 (/ms)		
	Erev	= 0	(mV)		
	Prethresh = 0 			
	Deadtime = 1	(ms)	
	gmax		(umho)		
	mg	= 1    (mM)			
	f = 0.0 				
	Ftau = 50 (ms)			
	Finf = 1				
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
	B				
	F				
	F0				
	F1				
}

INITIAL {
	R = 0
	C = 0
	F = 1
	F0 = 1
	F1 = 1
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / ((Alpha * Cmax) + Beta)
	lastrelease = -9e9
	
}

BREAKPOINT {
	SOLVE release

	B = mgblock(v)		

	g = gmax * R * B * F 

	i = g*(v - Erev)
}

PROCEDURE release() { LOCAL q
	


	q = ((t - lastrelease) - Cdur)		

						
	if (q > Deadtime) {
		if (pre > Prethresh) {		
			C = Cmax			
			R0 = R
			F = F + f
			F0 = F
			lastrelease = t
		}
						
	} else if (q < 0) {			
	
		
	
	} else if (C == Cmax) {			
		R1 = R
		C = 0
		F1 = F

	}

	F = Finf + (F0 - Finf) * exp(- (t - lastrelease) / Ftau)

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


FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}