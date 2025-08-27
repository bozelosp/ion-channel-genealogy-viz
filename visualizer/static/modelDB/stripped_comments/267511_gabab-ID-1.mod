INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAB
	RANGE R, G, g
	NONSPECIFIC_CURRENT i
	GLOBAL Cmax, Cdur
	GLOBAL K1, K2, K3, K4, KD, Erev, cutoff
	RANGE warn
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax	= 0.5	(mM)		
	Cdur	= 0.3	(ms)		



	K1	= 0.52	(/ms mM)	
	K2	= 0.0013 (/ms)		
	K3	= 0.098 (/ms)		
	K4	= 0.033 (/ms)		
	KD	= 100			
	n	= 4			
	Erev	= -95	(mV)		
	warn	= 0			
        cutoff = 1e12
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	Gn
	R				
	edc
	synon
	Rinf
	Rtau (ms)
	Beta (/ms)
}

STATE {
	Ron Roff
	G				
}


INITIAL {
	R = 0
	G = 0
	Ron = 0
	Roff = 0
	synon = 0
	Rinf = K1*Cmax/(K1*Cmax + K2)
	Rtau = 1/(K1*Cmax + K2)
	Beta = K2

}

BREAKPOINT {
	SOLVE bindkin METHOD derivimplicit
	if (G < cutoff) {
		Gn = G*G*G*G 
		g = Gn / (Gn+KD)
	} else {
		if(!warn){
			printf("gabab.mod WARN
			warn = 1
		}
		g = 1
	}
	i = g*(v - Erev)
}


DERIVATIVE bindkin {
	Ron' = synon*K1*Cmax - (K1*Cmax + K2)*Ron
	Roff' = -K2*Roff
	R = Ron + Roff
	G' = K3 * R - K4 * G
}





NET_RECEIVE(weight,  r0, t0 (ms)) {
	if (flag == 1) { 
		r0 = weight*(Rinf + (r0 - Rinf)*exp(-(t - t0)/Rtau))
		t0 = t
		synon = synon - weight
		Ron = Ron - r0
		Roff = Roff + r0
        }else{ 
		r0 = weight*r0*exp(-Beta*(t - t0))
		t0 = t
		synon = synon + weight
		Ron = Ron + r0
		Roff = Roff - r0
		
		net_send(Cdur, 1)
        }
}