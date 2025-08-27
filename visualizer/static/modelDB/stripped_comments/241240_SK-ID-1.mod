NEURON {
	SUFFIX 	SK
	USEION 	ca READ cai, cao, ica
	USEION	k READ ek WRITE ik
	RANGE 	gbar, g, i, inf_c, tau_c

}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
	(pS)	=	(picosiemens)
	(um)	=	(micrometer)


	
	
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	gbar		= 	5 	(pS/um2)  
	
	hill_c	=	4.64	
	hill_t	=	2
	K_c		=	.00056	(mM)	
	tauA_c	=	45  	(ms)  
	tau0_c	=	5	(ms)


	diff		=	1	(1)	
	cai0		=	50e-6	(mM)
	scale		=	1	(1)	
}

ASSIGNED { 

	v		(mV)
	i		(mA/cm2)	
	ik		(mA/cm2)
	g		(pS/um2)
	ek		(mV)
	ica		(mA/cm2)
	cai		(mM)
	cao		(mM)
	inf_c
	tau_c		(ms)

}

STATE {
	c	
}		

BREAKPOINT {
	SOLVE states METHOD cnexp
	g 	= 0
	g 	= gbar  * c
	i 	= g * (v - ek)*(1e-4)
	ik 	= i
}

INITIAL {
	cai=cai0
	rates( cai0)
	c = inf_c
}

DERIVATIVE states {
	rates(cai0+scale*(cai-cai0)) 
	c' = ( inf_c - c) / tau_c
}

PROCEDURE rates ( cai ( mM)) {
	inf_c = 1/(1 + (K_c/cai)^hill_c)
 	tau_c = tau0_c + tauA_c/(1 + (cai/K_c)^hill_t)

}