INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAb3
	POINTER C
	RANGE R, D, G, g, gmax
	NONSPECIFIC_CURRENT i
	GLOBAL K1, K2, K3, K4, KD, d1, d2, Erev
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {





	K1	= 0.66	(/ms mM)	
	K2	= 0.020 (/ms)		
	K3	= 0.083 (/ms)		
	K4	= 0.0079 (/ms)		
	d1	= 0.017 (/ms)		
	d2	= 0.0053 (/ms)		
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
}


STATE {
	R				
	D				
	G				
}


INITIAL {
	R = 0
	D = 0
	G = 0
}

BREAKPOINT {
	SOLVE bindkin METHOD cnexp
	Gn = G^n
	g = gmax * Gn / (Gn+KD)
	i = g*(v - Erev)
}


DERIVATIVE bindkin {

	R' = K1 * C * (1-R-D) - K2 * R + d2 * D
	D' = d1 * R - d2 * D
	G' = K3 * R - K4 * G

}