NEURON {
	SUFFIX cadyn_new
	USEION ca READ cai,ica WRITE cai 
	RANGE  ca 
	GLOBAL depth,cainf,taur 
     
}

UNITS {
	(molar) = (1/liter)		
	(mM)    = (milli/liter)
	(um)	= (micron) 
	(mA)    = (milliamp)
	(msM)	= (ms mM)  
	FARADAY = (faraday) (coul)
}

PARAMETER {
	depth	= .1	(um)		
	taur    =  200  (ms)	
	cainf	= 50e-6 (mM)	
	cai		        (mM)
}

ASSIGNED {
	ica		      (mA/cm2)
	diam          (um)
	VSR           (um)
	B             (mM*cm2/mA) 
	drive_channel (mM/ms)
}

STATE {
	ca (mM) 
}

BREAKPOINT {
	SOLVE state METHOD euler
}

DERIVATIVE state {
	drive_channel =  - ica * B

	if (drive_channel <= 0.) {
		
		drive_channel = 0.
	}
    
    ca' = drive_channel/18 + (cainf -ca)/taur*11
    
	cai = ca
}

INITIAL {



	if (2*depth >= diam) {
		VSR = 0.25 * diam 
	} else {
		VSR = depth * ( 1 - depth/diam)
	}

	B  = (1e4) / (2*FARADAY*VSR)
	ca = cainf
}