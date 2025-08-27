INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cadyn
	USEION ca READ ica, cai WRITE cai
	RANGE depth,kt,kd,cainf,taur
}

UNITS {
	(molar) = (1/liter)			
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
}

CONSTANT {
	FARADAY = 96489		(coul)		
}

PARAMETER {
	depth	= .1	(um)		
	taur	= 1e10	(ms)		
	cainf	= 2.4e-4 (mM)
	kt	= 1e-4	(mM/ms)
	kd	= 1e-4	(mM)
}

STATE {
	cai		(mM) 
}

INITIAL {
	cai = kd
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
	drive_pump	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)

	if (drive_channel <= 0.) { drive_channel = 0. }	

	drive_pump = -kt * cai / (cai + kd )		

	cai' = drive_channel + drive_pump + (cainf-cai)/taur
}