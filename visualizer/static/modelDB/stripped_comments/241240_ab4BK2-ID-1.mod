NEURON {
	SUFFIX ab4BK2
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gbar,  g, i
	GLOBAL	inf_n, tau_n
}
UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(pS)	= (picosiemens)
	(um)	= (micrometer)
	(S) = (siemens)

}
PARAMETER {
	gbar 	= 3	(pS/um2)
	
	
	vhalf_n 	= -56.5 	(mV)
	slope_n 	= 11.8	(mV)

	
	tauA_n	= 14.5		(ms)
	tauG_n	= 0.48	
	tau0_n	= 1		(ms)  

	
	bv 	= 	0.003		(mM)
	nv 	= 2
	ah 	= 120		(mV)

	
	bt	=	0.0022		(mM)
	nt	=	4
	th	=	140		(ms)

	
	scale	= 1	(1)	
	cai0 	= .00005	(mM)

}
ASSIGNED {
	v		(mV)
	celsius	(degC)
	ek		(mV)
	ik		(mA/cm2)
	ica		(mA/cm2)
	i		(mA/cm2)
	area		(microm2)
	g		(pS/um2)
	inf_n		
	tau_n		(ms)
	vh		(mV)
	ta		(ms)
 	cai		(mM)

}

STATE { n }

BREAKPOINT {
	SOLVE state METHOD cnexp
	g=0
	if( n >=0 ) {  
	   g = gbar * n
	}
	i = g * (v - ek) * (1e-4)
	ik = i
 
}

INITIAL {
	rates(v, cai)
	n = inf_n
}

DERIVATIVE state {

	rates (v, cai0+scale*(cai-cai0))
	n' = (inf_n - n) / tau_n
}


PROCEDURE rates (v (mV), cai (mM)) {

	vh 	= vhalf_n + ah /( 1+(cai/bv)^nv ) 

	inf_n	= 1 / (1 + exp((v-vh)/-slope_n))

	ta	= tauA_n + th / ( 1 + (cai/bt)^nt) 

	tau_n	= tau0_n + 4*ta*sqrt(tauG_n*(1-tauG_n))/(exp(tauG_n*(v - vh)/slope_n)+exp(-(1-tauG_n)*(v - vh)/slope_n))

}