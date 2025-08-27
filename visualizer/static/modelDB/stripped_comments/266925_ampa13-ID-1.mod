INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS AMPA13
	POINTER C
	RANGE C0, C1, C2, C3, C4, D1, D2, D3, D4, O1, O2, O3, O4
	RANGE g, gmax, rb1, rb2, rb3, rb4, Q10_binding, Q10_desensitization, Q10_opening, Q10_unbinding
	GLOBAL Erev
	GLOBAL Rb1, Rb2, Rb3, Rb4, Ru1, Ru2, Ru3, Ru4, Rd1, Rd2, Rd3, Rd4, Rr1, Rr2, Rr3, Rr4, Ro1, Ro2, Ro3, Ro4, Rc1, Rc2, Rc3, Rc4
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
    Q10_binding = 2.4
    Q10_unbinding = 2.4
    Q10_desensitization = 2.4
    Q10_opening = 2.4
    
    celsius       (degC)
	Erev	= 0    (mV)	
	gmax	= 20  (pS)	
	vmin = -120	(mV)
	vmax = 100	(mV)
	


	Rb1	= 800   (/mM /ms) 
	Rb2	= 600   (/mM /ms) 
	Rb3	= 400   (/mM /ms) 
	Rb4	= 200   (/mM /ms) 

	Ru1	= 30	(/ms)	
	Ru2	= 40  (/ms)
	Ru3	= 60	(/ms)	
	Ru4	= 80	(/ms)	
	
	Rd1	= .25		(/ms)	
	Rd2	= .25		(/ms)	
	Rd3	= 1		(/ms)	
	Rd4	= 1		(/ms)	

	Rr1	= 0.05 (/ms)	
	Rr2	= 0.05 (/ms)	
	Rr3	= 0.022 (/ms)	
	Rr4	= 0.022 (/ms)	

	Ro1	= 3	(/ms)	
	Ro2	= 4	(/ms)	
	Ro3	= 4	(/ms)	
	Ro4	= 4	(/ms)	
    
	Rc1	= 1.5	(/ms)	
   	Rc2	= 1	(/ms)	
    Rc3	= 1	(/ms)	
  	Rc4	= 1.5	(/ms)	
}

ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(pS)		
	C 		(mM)		

	rb1		(/ms)    
	rb2		(/ms)    
	rb3		(/ms)    
	rb4		(/ms)    
    
    Q10b (1)
    Q10u (1)
    Q10dr (1)
    Q10oc (1)
}

STATE {
	
	C0		
	C1		
	C2		
	C3		
	C4		
 	D1		
 	D2		
	D3		
	D4		
	O1		
    O2		
    O3		
    O4		
}

INITIAL {
	C0=1
	C1=0
	C2=0
	C3=0
	C4=0
	D1=0
	D2=0
	D3=0
	D4=0
	O1=0
    O2=0
    O3=0
    O4=0
    
    Q10b = Q10_binding^((celsius-22)/10)
    Q10u = Q10_unbinding^((celsius-22)/10)
    Q10dr = Q10_desensitization^((celsius-22)/10)
    Q10oc = Q10_opening^((celsius-22)/10)
}

BREAKPOINT {
	SOLVE kstates METHOD sparse

	g = gmax * (O4 + 0.75*O3 + 0.5*O2 + 0.25*O1)
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {
	
	rb1 = Rb1 * C 
	rb2 = Rb2 * C
    rb3 = Rb3 * C
	rb4 = Rb4 * C
	~ C0 <-> C1	(rb1*Q10b,Ru1*Q10u)
	~ C1 <-> C2	(rb2*Q10b,Ru2*Q10u)
	~ C2 <-> C3	(rb3*Q10b,Ru3*Q10u)
	~ C3 <-> C4	(rb4*Q10b,Ru4*Q10u)
	~ C1 <-> D1	(Rd1*Q10dr,Rr1*Q10dr)
	~ C2 <-> D2	(Rd2*Q10dr,Rr2*Q10dr)
	~ C3 <-> D3	(Rd3*Q10dr,Rr3*Q10dr)
	~ C4 <-> D4	(Rd4*Q10dr,Rr4*Q10dr)
	~ C1 <-> O1	(Ro1*Q10oc,Rc1*Q10oc)
	~ C2 <-> O2	(Ro2*Q10oc,Rc2*Q10oc)
	~ C3 <-> O3	(Ro3*Q10oc,Rc3*Q10oc)
	~ C4 <-> O4	(Ro4*Q10oc,Rc4*Q10oc)

	CONSERVE C0+C1+C2+C3+C4+D1+D2+D3+D4+O1+O2+O3+O4 = 1
}