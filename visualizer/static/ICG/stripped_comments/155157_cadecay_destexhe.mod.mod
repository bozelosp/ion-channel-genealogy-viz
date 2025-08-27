INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cadecay_destexhe
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
	taur	= 1000	(ms)		
	cainf	= 2e-4	(mM)
	kt	= 0	(mM/ms)		
	kd	= 0	(mM)		
}

STATE {
	cai		(mM) 
}

INITIAL {
	cai = cainf
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD euler
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)

	if (drive_channel <= 0.) { drive_channel = 0. }	

	cai' = drive_channel + (cainf-cai)/taur
}