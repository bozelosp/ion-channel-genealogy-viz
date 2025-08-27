INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	RANGE ca
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
	depth	= 0.35	(um)		
	taur	= 70	(ms)		
	cainf	= 5e-7(mM)
	cai		(mM)

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

	ca' = drive_channel + (cainf-ca)/taur
	cai = ca
}