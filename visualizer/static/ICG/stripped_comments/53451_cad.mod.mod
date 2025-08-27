INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cad
	USEION ca READ ica WRITE cai
	GLOBAL depth,cainf,taur
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
	taur	= 200	(ms)		
	cainf	= 100e-6(mM)
}

STATE {
	cai		(mM) <1e-6>
}

INITIAL {

	cai = cainf	
			
			
			
			

}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit

}

DERIVATIVE state { 

	
	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0. }	

	cai' = drive_channel + (cainf-cai)/taur

}