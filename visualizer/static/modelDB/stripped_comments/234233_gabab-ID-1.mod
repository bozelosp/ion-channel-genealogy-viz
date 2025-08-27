INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAb_S
	RANGE R, G, Gn, g, gmax, synon, Ron, Roff
	NONSPECIFIC_CURRENT i
	GLOBAL Cmax, Cdur
	GLOBAL K1, K2, K3, K4, KD, Erev
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



	K1	= 0.09	(/ms mM)	
	K2	= 0.0012 (/ms)		
	K3	= 0.18 (/ms)		
	K4	= 0.034 (/ms)		
	KD	= 100			
	n	= 4			
	Erev	= -95	(mV)		
	gmax		(umho)	
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
	synon = 0
	Rinf = K1*Cmax/(K1*Cmax + K2)
	Rtau = 1/(K1*Cmax + K2)
	Beta = K2
}

BREAKPOINT {
	SOLVE bindkin METHOD cnexp
	Gn = G*G*G*G 
	g = gmax * Gn / (Gn+KD)
	i = g*(v - Erev)
}


DERIVATIVE bindkin {
	Ron' = synon*K1*Cmax - (K1*Cmax + K2)*Ron
	Roff' = -K2*Roff
	R = Ron + Roff
	G' = K3 * R - K4 * G
}





NET_RECEIVE(weight (umho),  r0, t0 (ms)) {
	if (flag == 1) { 
		r0 = weight*(Rinf + (r0 - Rinf)*exp(-(t - t0)/Rtau))
		t0 = t
		synon = synon - weight
		state_discontinuity(Ron, Ron - r0)
		state_discontinuity(Roff, Roff + r0)
        }else{ 
		r0 = weight*r0*exp(-Beta*(t - t0))
		t0 = t
		synon = synon + weight
		state_discontinuity(Ron, Ron + r0)
		state_discontinuity(Roff, Roff - r0)
		
		net_send(Cdur, 1)
        }
}