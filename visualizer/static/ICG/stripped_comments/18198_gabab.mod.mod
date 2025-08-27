INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAb
	POINTER pre
	RANGE C, R, G, g, gmax, lastrelease, TimeCount
	NONSPECIFIC_CURRENT i
	GLOBAL Cmax, Cdur, Prethresh, Deadtime
	GLOBAL K1, K2, K3, K4, KD, Erev
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
	Prethresh = 0 			
	Deadtime = 1	(ms)		



	K1	= 0.52	(/ms mM)	
	K2	= 0.0013 (/ms)		
	K3	= 0.098 (/ms)		
	K4	= 0.033 (/ms)		
	KD	= 100			
	n	= 4			
	Erev	= -95	(mV)		
	gmax		(umho)		
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	C		(mM)		
	Gn
	pre 				
	lastrelease	(ms)		
	TimeCount	(ms)		
}


STATE {
	R				
	G				
}


INITIAL {
	C = 0
	lastrelease = -1000

	R = 0
	G = 0
	TimeCount=-1
}

BREAKPOINT {
	SOLVE bindkin METHOD cnexp
	Gn = G^n
	g = gmax * Gn / (Gn+KD)
	i = g*(v - Erev)
}


DERIVATIVE bindkin {

	release()		

	R' = K1 * C * (1-R) - K2 * R
	G' = K3 * R - K4 * G

}


PROCEDURE release() {
	

	TimeCount=TimeCount-dt			

						
	if (TimeCount < -Deadtime) {
		if (pre > Prethresh) {		
			C = Cmax			
			lastrelease = t
			TimeCount=Cdur
		}
						
	} else if (TimeCount > 0) {		
	
		
	
	} else if (C == Cmax) {			
		C = 0.
	}

}