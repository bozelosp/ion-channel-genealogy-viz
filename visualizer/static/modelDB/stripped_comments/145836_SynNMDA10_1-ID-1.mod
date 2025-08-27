INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDA10_1
	RANGE T_max, T, tau, tRel, Erev, synon
	RANGE U, Cl, D1, D2, O, UMg, ClMg, D1Mg, D2Mg, OMg
	RANGE g, gmax, rb, rmb, rmu, rbMg,rmc1b,rmc1u,rmc2b,rmc2u
	GLOBAL mg, Rb, Ru, Rd1, Rr1, Rd2, Rr2, Ro, Rc, Rmb, Rmu
	GLOBAL RbMg, RuMg, Rd1Mg, Rr1Mg, Rd2Mg, Rr2Mg, RoMg, RcMg
	GLOBAL Rmd1b,Rmd1u,Rmd2b,Rmd2u,rmd1b,rmd1u,rmd2b,rmd2u
	GLOBAL Rmc1b,Rmc1u,Rmc2b,Rmc2u
	GLOBAL valence, memb_fraction
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

	Erev	= 0    	(mV)	
	gmax	= 50  	(pS)	
	mg		= 1  	(mM)	
	valence = -2			
	memb_fraction = 0.8
	
	
	tau  = .3 (ms) <1e-9,1e9>
	T_max = 1.5 (mM)		



	Rb		= 10	   		(/mM /ms)	
	Ru		= 5.6e-3  		(/ms)		
	Ro		= 10e-3   		(/ms)		
	Rc		= 273e-3   		(/ms)		
	Rd1		= 2.2e-3   		(/ms)		
	Rr1		= 1.6e-3   		(/ms)		
	Rd2 	= 0.43e-3 		(/ms)		
	Rr2 	= 0.5e-3		(/ms)		
	Rmb		= 0.05			(/mM /ms)	
	Rmu		= 12800e-3		(/ms)		
	Rmc1b	= 0.00005		(/mM /ms)	
	Rmc1u	= 2.438312e-3	(/ms)		
	Rmc2b	= 0.00005		(/mM /ms)	
	Rmc2u	= 5.041915e-3	(/ms)		
	Rmd1b	= 0.00005		(/mM /ms)	
	Rmd1u	= 2.98874e-3	(/ms)		
	Rmd2b	= 0.00005		(/mM /ms)	
	Rmd2u	= 2.953408e-3	(/ms)		
	RbMg	= 10			(/mM /ms)	
	RuMg	= 17.1e-3		(/ms)		
	RoMg	= 10e-3			(/ms)		
	RcMg	= 548e-3		(/ms)		
	Rd1Mg	= 2.1e-3		(/ms)		
	Rr1Mg	= 0.87e-3		(/ms)		
	Rd2Mg	= 0.26e-3		(/ms)		
	Rr2Mg	= 0.42e-3		(/ms)		
}

ASSIGNED {
	v		(mV)	
	i 		(nA)	
	g 		(uS)	
	
	T		(mM)	
	tRel	(ms)	
	synon			
	w				

	rb		(/ms)   
	rmb		(/ms)	
	rmu		(/ms)	
	rbMg	(/ms)	
	rmc1b	(/ms)	
	rmc1u	(/ms)	
	rmc2b	(/ms)	
	rmc2u	(/ms)	
	rmd1b	(/ms)	
	rmd1u	(/ms)	
	rmd2b	(/ms)	
	rmd2u	(/ms)	
}

STATE {
	
	U		
	Cl		
	D1		
	D2		
	O		
	UMg		
	ClMg	
	D1Mg	
	D2Mg	
	OMg		
}

INITIAL {
	T = 0
	synon = 0
	tRel = 0
	
	U = 1
}

BREAKPOINT {
	SOLVE kstates METHOD sparse

	g = w * gmax * O
	i = g * (v - Erev)
}

KINETIC kstates {
	release(t)
	
	rb 		= Rb   * T
	rbMg 	= RbMg * T
	
	rmb 	= Rmb 	* mg * exp((1 (/mV) * v-40) * valence * memb_fraction /25)
	rmu 	= Rmu 	* exp((-1 (/mV))*(v-40) * valence * (1-memb_fraction) /25)
	rmc1b 	= Rmc1b * mg * exp((1 (/mV) * v-40) * valence * memb_fraction /25)
	rmc1u 	= Rmc1u * exp((-1 (/mV))*(v-40) * valence * (1-memb_fraction) /25)
	rmc2b 	= Rmc2b * mg * exp((1 (/mV) * v-40) * valence * memb_fraction /25)
	rmc2u 	= Rmc2u * exp((-1 (/mV))*(v-40) * valence * (1-memb_fraction) /25)
	rmd1b 	= Rmd1b * mg * exp((1 (/mV) * v-40) * valence * memb_fraction /25)
	rmd1u 	= Rmd1u * exp((-1 (/mV))*(v-40) * valence * (1-memb_fraction) /25)
	rmd2b 	= Rmd2b * mg * exp((1 (/mV) * v-40) * valence * memb_fraction /25)
	rmd2u 	= Rmd2u * exp((-1 (/mV))*(v-40) * valence * (1-memb_fraction) /25)

	~ U <-> Cl	(rb,Ru)
	~ Cl <-> O	(Ro,Rc)
	~ Cl <-> D1	(Rd1,Rr1)
	~ D1 <-> D2	(Rd2,Rr2)
	~ O <-> OMg	(rmb,rmu)
	~ UMg <-> ClMg 	(rbMg,RuMg)
	~ ClMg <-> OMg 	(RoMg,RcMg)
	~ ClMg <-> D1Mg (Rd1Mg,Rr1Mg)
	~ D1Mg <-> D2Mg (Rd2Mg,Rr2Mg)
	~ U <-> UMg     (rmc1b,rmc1u)
	~ Cl <-> ClMg	(rmc2b,rmc2u)
	~ D1 <-> D1Mg	(rmd1b,rmd1u)
	~ D2 <-> D2Mg	(rmd2b,rmd2u)

	CONSERVE U+Cl+D1+D2+O+UMg+ClMg+D1Mg+D2Mg+OMg = 1
}

NET_RECEIVE(weight) {
	if (flag == 0) {
		tRel = t	
		synon = 1	
					
		w = weight
	}
}
PROCEDURE release(t(ms)) {
	T = T_max * (t - tRel) / tau * exp(1 - (t - tRel) / tau) * synon
}