INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cadspine
	USEION ca READ ica, cai WRITE cai
	RANGE ca
	GLOBAL depth,cainf,taur,cadend,D,lneck,sneck,vspine
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
	depth	= 0.125	(um)       
	taur	= 15	(ms)       
	cainf	= 5e-7  (mM)
	cai		(mM)
	D = 0.22	(um2/ms) 
	cadend		(mM)       
	lneck		(um)       
	sneck		(um2)    
	vspine		(um3) 
}

STATE {
	ca		(mM)
}

INITIAL {
	ca = cainf
        cai=ca
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)

}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit 
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth * 18.0)
	if (drive_channel <= 0.) { drive_channel = 0. }	

	ca' = drive_channel + (cainf-ca)/taur-D*(ca-cadend)*sneck/(lneck*vspine)
	cai = ca
}