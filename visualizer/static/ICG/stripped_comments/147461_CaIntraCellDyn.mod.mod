INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX CaIntraCellDyn
	USEION ca READ ica, cai WRITE cai	
  RANGE cai_new, depth, cai_inf, cai_tau
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
	depth	= 0.1	  (um)		
	cai_tau	= 2.0     (ms)		
	cai_inf	= 50.0e-6 (mM)		
	cai		  (mM)
}

STATE {
	cai_new		(mM) 
}

INITIAL {

	cai_new = cai_inf
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)
	if (drive_channel <= 0.) { drive_channel = 0.  }   
         
	cai_new' = drive_channel + (cai_inf-cai_new)/cai_tau
	cai = cai_new
}