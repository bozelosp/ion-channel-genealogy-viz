INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAb
	POINTER pre
	RANGE C, R, G, g, gmax, lastrelease
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

	Cmax	= 1	(mM)		
	Cdur	= 1	(ms)		
	Prethresh = 0 			
	Deadtime = 1	(ms)		





	K1	= 0.2	(/ms mM)	
	K2	= 0.0012 (/ms)		
	K3	= 0.18 (/ms)		
	K4	= 0.034 (/ms)		
	KD	= 100			
	n	= 4			
	Erev	= -90	(mV)		
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
}

BREAKPOINT {
	SOLVE bindkin METHOD euler
	Gn = G^n
	g = gmax * Gn / (Gn+KD)
	i = g*(v - Erev)
}


DERIVATIVE bindkin {

	release()		

	R' = K1 * C * (1-R) - K2 * R
	G' = K3 * R - K4 * G

}


PROCEDURE release() { LOCAL q
	

	q = ((t - lastrelease) - Cdur)		

						
	if (q > Deadtime) {
		if (pre > Prethresh) {		
			C = Cmax			
			lastrelease = t
		}
						
	} else if (q < 0) {			
	
		
	
	} else if (C == Cmax) {			
		C = 0.
	}

}