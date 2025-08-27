INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS AMPA5
	POINTER C
	RANGE C0, C1, C2, D1, D2, O
	RANGE g, gmax, rb
	GLOBAL Erev
	GLOBAL Rb, Ru1, Ru2, Rd, Rr, Ro, Rc
	GLOBAL vmin, vmax
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

	Erev	= 0    (mV)	
	gmax	= 500  (pS)	
	vmin = -120	(mV)
	vmax = 100	(mV)
	


	Rb	= 13   (/mM /ms)
				
	Ru1	= 0.0059  (/ms)	
	Ru2	= 86  (/ms)	
	Rd	= 0.9   (/ms)	
	Rr	= 0.064 (/ms)	
	Ro	= 2.7    (/ms)	
	Rc	= 0.2    (/ms)	
}

ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(pS)		
	C 		(mM)		

	rb		(/ms)    
}

STATE {
	
	C0		
	C1		
	C2		
 	D1		
 	D2		
	O		
}

INITIAL {
	C0=1
	C1=0
	C2=0
	D1=0
	D2=0
	O=0
}

BREAKPOINT {
	SOLVE kstates METHOD sparse

	g = gmax * O
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {
	
	rb = Rb * C 

	~ C0  <-> C1	(rb,Ru1)
	~ C1 <-> C2	(rb,Ru2)
	~ C1 <-> D1	(Rd,Rr)
	~ C2 <-> D2	(Rd,Rr)
	~ C2 <-> O	(Ro,Rc)

	CONSERVE C0+C1+C2+D1+D2+O = 1
}