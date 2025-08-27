INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAa5
	POINTER C
	RANGE C0, C1, C2, O1, O2
	RANGE g, gmax, f1, f2
	GLOBAL Erev, kf1, kf2, kb1, kb2, a1, b1, a2, b2
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

PARAMETER {

	Erev	= -80    (mV)	
	gmax	= 500  (pS)	
	


	

	kf1	= 0.02   (/uM /ms)	
	kf2	= 0.01   (/uM /ms)	
	kb1	= 4.6	(/ms)	
	kb2	= 9.2	(/ms)	
	a1	= 3.3	(/ms)	
	b1	= 9.8	(/ms)	
	a2	= 10.6	(/ms)	
	b2	= 0.41  (/ms)	
}



ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(pS)		
	C 		(mM)		

	f1		(/ms)    
	f2		(/ms)    
}

STATE {
	
	C0		
	C1		
	C2		
	O1		
	O2		
}

INITIAL {
	C0 = 1
	C1 = 0
	C2 = 0
	O1 = 0
	O2 = 0
}

BREAKPOINT {
	SOLVE kstates METHOD sparse

	g = gmax * (O1+O2)
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {
	
	f1 = kf1 * (1e3) * C 
	f2 = kf2 * (1e3) * C 

	~ C0 <-> C1	(f1,kb1)
	~ C1 <-> C2	(f2,kb2)
	~ C1 <-> O1	(a1,b1)
	~ C2 <-> O2	(a2,b2)

	CONSERVE C0+C1+C2+O1+O2 = 1
}