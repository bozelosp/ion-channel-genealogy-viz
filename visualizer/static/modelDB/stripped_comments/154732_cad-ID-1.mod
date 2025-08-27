NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai	
        RANGE ca
	GLOBAL depth,cainf,taur,bcap
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

	taur	= 28.6	(ms)		
	bcap	= 17	(1)		
	cainf	= 100e-6(mM)
	cai		(mM)
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
        SOLVE state METHOD derivimplicit 
    
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0.  }   
         
	
        ca' = drive_channel/(1+bcap) + (cainf-ca)/taur
	cai = ca
}