NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai	
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
	taur	= 2	(ms)		

	cainf	= 2.4e-4 (mM)

}

STATE {
	cai		(mM) 
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
	if (drive_channel <= 0.) { drive_channel = 0.  }   
         
	cai' = drive_channel + (cainf-cai)/taur
}