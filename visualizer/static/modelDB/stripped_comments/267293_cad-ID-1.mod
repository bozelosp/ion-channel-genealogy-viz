INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	RANGE ca,depth,cainf,taur,gamma
}

UNITS {
	(molar) = (1/liter)			
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
	FARADAY = (faraday) (coulomb)
}


PARAMETER {
	depth	= .1	(um)		
	taur	= 80	(ms)		
	cainf	= 100e-6(mM)
	cai		(mM)
	gamma = 1.0
}

STATE {
	ca		(mM) 
}

INITIAL {
	ca = cainf
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ica * gamma / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0. }	

	ca' = drive_channel + (cainf-ca)/taur
	cai = ca
}