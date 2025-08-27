NEURON {
	SUFFIX Cabuffer
	USEION ca READ ica,cai,cao WRITE cai, cao
	USEION nca READ inca WRITE ncai VALENCE 0
	USEION lca READ ilca WRITE lcai VALENCE 0
	USEION tca READ itca WRITE tcai VALENCE 0
	GLOBAL depth,cao0,cai0
	RANGE cai,cao,ncai,lcai,brat, tau
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	(um) = (micrometer)
}

PARAMETER {
	tau = 9				(ms)
	depth = .2 		(um)
	cai0  				(mM)	
	cao0   				(mM)	
	Fa = 96485.3365 (coulomb)
	brat = 1  
}

ASSIGNED {
	ica		(mA/cm2)
	diam	(um)
	VSR (um)
	ncai		(mM)
	inca		(mA/cm2) 
	lcai		(mM)
	ilca		(mA/cm2) 
	tcai		(mM)
	itca		(mA/cm2) 
	B 			(mM*cm2/mA)
}

STATE { 
cai (mM) 		<1e-5> 
cao (mM)}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {	
	ncai = - inca * B  
	lcai = - ilca * B 
	tcai = - itca * B  
	cai' = -ica * B / brat -(cai-cai0)/tau	
	cao' = 0
}

INITIAL {
	if (2*depth >= diam) {
		VSR = 0.25*diam 
	}else{
		VSR = depth*(1-depth/diam)
	}
	B = (1e4)/(2*Fa*VSR)
	cao0 	= 		cao
	cai0 	= 		cai
}