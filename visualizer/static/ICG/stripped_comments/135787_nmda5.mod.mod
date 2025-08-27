INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDA5
	POINTER C
	RANGE C0, C1, C2, D, O, B
	RANGE g, gmax, rb
	GLOBAL Erev, mg, Rb, Ru, Rd, Rr, Ro, Rc
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
	mg	= 0    (mM)	
	vmin = -120	(mV)
	vmax = 100	(mV)
	


	
	Rb	= 5e-3    (/uM /ms)	
	Ru	= 12.9e-3  (/ms)	
	Rd	= 8.4e-3   (/ms)	
	Rr	= 6.8e-3   (/ms)	
	Ro	= 46.5e-3   (/ms)	
	Rc	= 73.8e-3   (/ms)	
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
	D		
	O		

	B		
}

INITIAL {
	rates(v)
	C0 = 1
}

BREAKPOINT {
	rates(v)
	SOLVE kstates METHOD sparse

	g = gmax * O * B
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {
	
	rb = Rb * (1e3) * C 

	~ C0 <-> C1	(rb,Ru)
	~ C1 <-> C2	(rb,Ru)
	~ C2 <-> D	(Rd,Rr)
	~ C2 <-> O	(Ro,Rc)

	CONSERVE C0+C1+C2+D+O = 1
}

PROCEDURE rates(v(mV)) {
	TABLE B
	DEPEND mg
	FROM vmin TO vmax WITH 200

	

	B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}