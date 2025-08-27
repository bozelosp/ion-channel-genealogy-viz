NEURON {
	POINT_PROCESS GABABsyn
	RANGE C, R, G, B, g, gmax, tauD
	NONSPECIFIC_CURRENT i
	RANGE vgat,sst,npy,pv,xEff
	RANGE isOn
	GLOBAL K1, K2, K3, K4, KD, k1, k_1, k2, e, Bm
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(uS) = (microsiemens)
}

PARAMETER {

	tauD = 10	(ms)		
	K1	= 0.066	(/ms mM)	
	K2	= 0.008 (/ms)		
	K3	= 0.27 (/ms)		
	K4	= 0.044 (/ms)		
	KD	= 0.5				
	n	= 2			
	e	= -95	(mV)		
	gmax		(uS)		
    f   = 0.1              
	k1	= 30	(/ms mM)	
	k_1	= 0.1 (/ms)		
	k2	= 0.02 (/ms)		
	Bm = 1 (mM)			
	vgat=0
	sst=0
	npy=0
	pv=0
	xEff=-1
	isOn=0
}


ASSIGNED {
	v		(mV)		
	i		(nA)		
	g		(uS)		
	Gn
}


STATE {
	C	(mM)		
	R				
	G				
	B	(mM)		
}


INITIAL {
	C = 0
	R = 0
	G = 0
	B = 0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	Gn = G^n
	g = isOn * gmax * Gn / (Gn+KD)
	i = g *(v - e)
}


DERIVATIVE state {

	C' = (-C/tauD -k1 * C * (Bm - B) + k_1 * B) 
	R' = (K1 * C * (1-R) - K2 * R) 
	G' = (K3 * R * (1-G) - K4 * G) * f
	B' = (k1 * C * (Bm - B) - (k_1 + k2) * B) 

}


NET_RECEIVE(weight (mM)) {
	C = C + weight
}