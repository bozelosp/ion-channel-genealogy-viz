INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cas
	USEION ca READ ica, cai WRITE cai
	RANGE kd,cainf,taur
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
	diam (micron)
	taur	= 1.5	(ms)		
	cainf	= 0.0001 (mM)
	kd	= 0.0011	(mM)
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
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
					
}

DERIVATIVE state { 

	drive_channel = - ( (2 * ica) / (FARADAY * diam) ) 

	if (drive_channel <= 0.) { drive_channel = 0. }	

	cai' = drive_channel - ( (cai - cainf) /  taur )
}