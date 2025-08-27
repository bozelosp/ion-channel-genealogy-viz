INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Cadyn
	USEION Ca READ iCa, Cai WRITE Cai VALENCE 2
	RANGE cainf, taur, drive
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
	pump = 0.02					
}

PARAMETER {
    drive   = 10000  (1)
	depth	= 0.1	(um)		
	cainf	= 1e-5	(mM)		
	taur	= 43	(ms)		
	kt	= 1e-4	(mM/ms)			
	kd	= 1e-4	(mM)
}

STATE {
	Cai		(mM) 
}

INITIAL {
	Cai = cainf
}

ASSIGNED {
	iCa		(mA/cm2)
	drive_channel	(mM/ms)
	drive_pump	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit 
}

DERIVATIVE state { 
	drive_channel =  - drive * iCa / (2 * FARADAY * depth)
	    
	    

	if (drive_channel <= 0.) { drive_channel = 0. }	

	drive_pump = -kt * Cai / (Cai + kd )		
	    
	    
	    
	
	Cai' = ( drive_channel + pump*drive_pump + (cainf-Cai)/taur )
	    
	    
}