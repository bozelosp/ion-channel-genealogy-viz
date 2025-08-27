INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDA10_2_2
	RANGE C0, C1, C2, D, O, CB0, CB1, CB2, DB, OB
	RANGE g, gmax, rb, RMgB, RMgU, Rd, Rr, Ro, Rc, Rb, Ru
	RANGE T_max, T, tau, tRel, Erev, synon
	GLOBAL mg
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(uS) = (microsiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

PARAMETER {

	Erev = -0.7	(mV)	
	gmax = 50	(pS)	
	mg	 = 1	(mM)	
	

	Rb	= 2.83		(/mM /ms)	
	Ru	= 38.1e-3  	(/ms)		
	Rd	= 4.7161   	(/ms)		
	Rr	= 0.16116  	(/ms)		
	Ro	= 0.099631	(/ms)		
	Rc	= 0.056999	(/ms)		
	
	
	tau  = .3 (ms) <1e-9,1e9>
	T_max = 1.5 (mM)			
}

ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(uS)		

	rb		(/ms)		
	RMgB	(/ms)
	RMgU	(/ms)
	
	T		(mM)		
	tRel	(ms)		
	synon				
	w					
}

STATE {
	
	C0		
	C1		
	C2		
	D		
	O		
	CB0		
	CB1		
	CB2		
	DB		
	OB		
}

INITIAL {
	T = 0
	synon = 0
	tRel = 0
	
	rates(v)
	C0	= 1
	C1	= 0
	C2	= 0
	D	= 0
	O	= 0
	CB0	= 0
	CB1	= 0
	CB2	= 0
	DB	= 0
	OB	= 0
	
	net_send(590, 1)
}

BREAKPOINT {
	SOLVE kstates METHOD sparse

	g = w * gmax * O
	i = g * (v - Erev)
}

KINETIC kstates {	
	release(t)
	
	rb = Rb * T 

	rates(v)

	~ C0  <-> C1	((2 * rb),Ru)
	~ C1  <-> C2	(rb,(2 * Ru))
	~ C2  <-> D		(Rd,Rr)
	~ C2  <-> O		(Ro,Rc)
	~ O   <-> OB	(RMgB,RMgU)
	~ OB  <-> CB2	((3*Rc),Ro)
	~ CB2 <-> DB	(Rd,Rr)
	~ CB2 <-> CB1	((2 * Ru),rb)
	~ CB1 <-> CB0	(Ru,(2 * rb))

	CONSERVE C0+C1+C2+D+O+CB0+CB1+CB2+DB+OB = 1
}

NET_RECEIVE(weight) {
	if (flag == 0) {
		tRel = t	
		synon = 1	
					
		w = weight
	}
	if (flag == 1) {
	
		C0	= 1
		C1	= 0
		C2	= 0
		D	= 0
		O	= 0
		CB0	= 0
		CB1	= 0
		CB2	= 0
		DB	= 0
		OB	= 0
	}
}
PROCEDURE release(t(ms)) {
	T = T_max * (t - tRel) / tau * exp(1 - (t - tRel) / tau) * synon
	
}

PROCEDURE rates(v(mV)) {
	RMgB = 610e-3 * exp(1 (/mV) * -v / 17) * (mg / 1 (mM)) * 1 (/ms)	
	RMgU = 5400e-3 * exp(1 (/mV) * v / 47) * 1 (/ms)					
}