NEURON {
	POINT_PROCESS NMDA_Mg_T
	POINTER C
	RANGE U, Cl, D1, D2, O, UMg, ClMg, D1Mg, D2Mg, OMg
	RANGE g, gmax, rb, rmb, rmu, rbMg,rmc1b,rmc1u,rmc2b,rmc2u
	GLOBAL Erev, mg, Rb, Ru, Rd1, Rr1, Rd2, Rr2, Ro, Rc, Rmb, Rmu
	GLOBAL RbMg, RuMg, Rd1Mg, Rr1Mg, Rd2Mg, Rr2Mg, RoMg, RcMg
	GLOBAL Rmd1b,Rmd1u,Rmd2b,Rmd2u,rmd1b,rmd1u,rmd2b,rmd2u
	GLOBAL Rmc1b,Rmc1u,Rmc2b,Rmc2u
	GLOBAL vmin, vmax, valence, memb_fraction
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

	Erev	= 5    	(mV)	
	gmax	= 500  	(pS)	
	mg	= 1  	(mM)	
	vmin 	= -120	(mV)
	vmax 	= 100	(mV)
	valence = -2		
	memb_fraction = 0.8



	Rb		= 10e-3    	(/uM /ms)	
	Ru		= 0.02016 	(/ms)	
	Ro		= 46.5e-3   	(/ms)	
	Rc		= 91.6e-3   	(/ms)	
	Rd1		= 0.02266  	(/ms)	
	Rr1		= 0.00736  	(/ms)	
	Rd2 		= 0.004429	(/ms)	
	Rr2 		= 0.0023	(/ms)	
	Rmb		= 0.05e-3	(/uM /ms)	
	Rmu		= 12800e-3	(/ms)	
	Rmc1b		= 0.00005e-3	(/uM /ms)	
	Rmc1u		= 0.06	(/ms)	
	Rmc2b		= 0.00005e-3	(/uM /ms)	
	Rmc2u		= 0.06	(/ms)	
	Rmd1b		= 0.00005e-3	(/uM /ms)	
	Rmd1u		= 0.06	(/ms)	
	Rmd2b		= 0.00005e-3	(/uM /ms)	
	Rmd2u		= 0.06	(/ms)	
	RbMg		= 10e-3		(/uM /ms)	
	RuMg		= 0.06156	(/ms)	
	RoMg		= 46.5e-3		(/ms)	
	RcMg		= 91.6e-3	(/ms)	
	Rd1Mg		= 0.02163	(/ms)	
	Rr1Mg		= 0.004002	(/ms)	
	Rd2Mg		= 0.002678	(/ms)	
	Rr2Mg		= 0.001932	(/ms)	
}

ASSIGNED {
	v		(mV)	
	i 		(nA)	
	g 		(pS)	
	C 		(mM)	

	rb		(/ms)   
	rmb		(/ms)	
	rmu		(/ms)	
	rbMg		(/ms)	
	rmc1b		(/ms)	
	rmc1u		(/ms)	
	rmc2b		(/ms)	
	rmc2u		(/ms)	
	rmd1b		(/ms)	
	rmd1u		(/ms)	
	rmd2b		(/ms)	
	rmd2u		(/ms)	
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
	U = 1
  
  
  SOLVE kstates STEADYSTATE sparse 
}

BREAKPOINT {
	SOLVE kstates METHOD sparse

	g = gmax * O
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {

	rb 	= Rb 	* (1e3) * C
	rbMg 	= RbMg 	* (1e3) * C
	rmb 	= Rmb 	* mg * (1e3) * exp((v-40) * valence * memb_fraction /25)
	rmu 	= Rmu 	* exp((-1)*(v-40) * valence * (1-memb_fraction) /25)
	rmc1b 	= Rmc1b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25)
	rmc1u 	= Rmc1u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25)
	rmc2b 	= Rmc2b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25)
	rmc2u 	= Rmc2u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25)
	rmd1b 	= Rmd1b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25)
	rmd1u 	= Rmd1u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25)
	rmd2b 	= Rmd2b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25)
	rmd2u 	= Rmd2u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25)

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