INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX caldyn
	USEION cal READ ical, cali WRITE cali VALENCE 2
	RANGE pump, cainf, taur, drive
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
    drive   = 10000  (1)
	depth	= 0.1	(um)		
	cainf	= 1e-5	(mM)		
	taur	= 43	(ms)		
	kt	= 1e-4	(mM/ms)			
	kd	= 1e-4	(mM)
	
	pump = 0.05897					
}

STATE {
	cali		(mM) 
}

INITIAL {
	cali = cainf
}

ASSIGNED {
	ical		(mA/cm2)
	drive_channel	(mM/ms)
	drive_pump	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 
	drive_channel =  - drive * ical / (2 * FARADAY * depth)
	    
	    

	if (drive_channel <= 0.) { drive_channel = 0. }	

	drive_pump = -kt * cali / (cali + kd )		
	    
	    
	    
	
		cali' = ( drive_channel + pump*drive_pump + (cainf-cali)/taur )
	    
	    
}