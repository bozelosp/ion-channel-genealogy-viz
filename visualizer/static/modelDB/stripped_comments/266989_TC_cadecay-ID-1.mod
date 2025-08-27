NEURON {
	SUFFIX TC_cad
        
	USEION cal1 READ ical1 WRITE cal1i VALENCE 2
        USEION cal2 READ ical2 WRITE cal2i VALENCE 2
	USEION cat READ icat  WRITE cati  VALENCE 2
        
        RANGE delta_L1, taur_L1, cal1_inf,  delta_L2, taur_L2, cal2_inf,  delta_T, taur_T, cat_inf
}

UNITS {

	(mM)	= (milli/liter)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
	FARADAY = (faraday) (coulombs)
}

PARAMETER {
	
	
        
        delta_L1   = 0.5    (/um)
        taur_L1    = 5      (ms)  
        cal1_inf   = 5e-5   (mM)  
		
        delta_L2   = 0.5    (/um)
        taur_L2    = 5      (ms)  
        cal2_inf   = 5e-5   (mM)  
        
        delta_T    = 0.5    (/um)
        taur_T     = 5      (ms)  
        cat_inf    = 5e-5   (mM)  

}

STATE {
	cal1i		(mM) 
	cal2i		(mM) 
	cati		(mM) 
}

INITIAL {
	cal1i = cal1_inf
	cal2i = cal2_inf
	cati  = cat_inf
}

ASSIGNED {
	ical1		(mA/cm2)
	ical2		(mA/cm2)
	icat		(mA/cm2)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
	
}

DERIVATIVE state {
        
        cal1i' = -10000*(ical1 * delta_L1/(2*FARADAY)) - (cal1i-cal1_inf)/taur_L1
        cal2i' = -10000*(ical2 * delta_L2/(2*FARADAY)) - (cal2i-cal2_inf)/taur_L2
        cati'  = -10000*(icat  * delta_T/(2*FARADAY))  - (cati-cat_inf)/taur_T
}