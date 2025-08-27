INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAbKG
	POINTER pre, vext, pmodyn
	USEION k READ ek WRITE ik
	RANGE C, R, S, G, g, gmax, lastrelease, TimeCount, ek, K1, K2, K5, K6, nsm, i
	GLOBAL Cmax, Cdur, Prethresh, Deadtime
	GLOBAL K3, K4, KD
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
	stimon = 0
	nsm = 1



	K1	= 0.52   (/ms mM) 
	K2	= 0.0013 (/ms)	
	K3	= 0.098 (/ms)		
	K4	= 0.033 (/ms)		
	K5	= 0.52   (/ms) 
	K6	= 0.00013 (/ms)	
	KD	= 100			
	n	= 4			
	gmax		(umho)		
}

ASSIGNED {
	ek 		(mV)
	v		(mV)		
	ik 		(nA)		
	g 		(umho)		
	C		(mM)		
	Gn
	pre 				
	lastrelease	(ms)		
	TimeCount	(ms)		
	vext				
	pmodyn
	i 		(nA)
}


STATE {
	R				
	S				
	G				
}


INITIAL {
	C = 0
	lastrelease = -1000

	R = 0
	S = 0
	G = 0
	TimeCount=-1
}

BREAKPOINT {
	SOLVE bindkin METHOD euler
	Gn = G^n
	g = gmax * Gn / (Gn+KD)
	i = g*(v - ek) 
	ik = i
}


DERIVATIVE bindkin {

	release()		

	R' = K1 * C * (1-R) - K2 * R
	S' = K5 * pmodyn * (1-S) - K6 * S
	G' = K3 * (R + 2 * nsm * S) - K4 * G
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
	
	if(vext==0) {      
		stimon=0
	}else{
		stimon=1
	}


}