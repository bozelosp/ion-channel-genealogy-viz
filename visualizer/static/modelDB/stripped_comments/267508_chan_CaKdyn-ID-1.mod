INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX CaKdyn
	USEION ca READ ica, cai WRITE cai
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
	
	pump = 0.02					
}

STATE {
	cai		(mM) 
}

INITIAL {
	cai = cainf
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
	drive_pump	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD euler
}

DERIVATIVE state { 
	drive_channel =  - drive * ica / (2 * FARADAY * depth)
	    
	    

	if (drive_channel <= 0.) { drive_channel = 0. }	

	drive_pump = -kt * cai / (cai + kd )		
	    
	    
	    
	
	cai' = ( drive_channel + pump*drive_pump + (cainf-cai)/taur )
	    
	    
}