NEURON {
	SUFFIX CAD
	USEION ca READ ica,cai WRITE cai  VALENCE 2 
	RANGE depth,kd,cainf,taur, mt
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
	depth	= 0.41	(um)  
	taur	= 1800	(ms)		
	cainf	= 0.00015 (mM)     
	kd	= 0.0003	(mM)
	mt = -100
}

STATE {
	cai		(mM) 
}

INITIAL {
	cai = kd
}

ASSIGNED {
        ica     (mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

	drive_channel =  - (1000) * ica / (2 * FARADAY * depth)

	if (drive_channel <= 0.) { drive_channel = 0. }	

	cai' = drive_channel - mt *(cainf-cai)/taur
}