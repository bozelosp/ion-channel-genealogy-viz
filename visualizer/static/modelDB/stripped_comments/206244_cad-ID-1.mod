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
	depth	= .1	(um)		
	depth_corr 	(um)
	taur	= 200	(ms)		
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
	diam		(um)
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 
	if (diam>2*depth) {
		depth_corr = depth*(diam-depth)/diam 
	} else {
		depth_corr = diam/4
	}
	drive_channel =  - (10000) * ica / (2 * FARADAY * depth_corr)
	if (drive_channel <= 0.) { drive_channel = 0.  }   
         
	
        ca' = drive_channel/18 + (cainf -ca)/taur*7
	cai = ca
}