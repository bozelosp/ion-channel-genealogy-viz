INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Cad
	USEION Ca READ iCa, Cai WRITE Cai VALENCE 2
	RANGE Cainf,taur,k
}

UNITS {
	(molar) = (1/liter)			
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
}


PARAMETER {
	depth	= .1(um)		
	taur	= 50	(ms)		
	Cainf	= 5e-5	(mM)  
	Cainit  = 5e-5 (mM)	
      k       = 0.0155458135   (mmol/C cm)  
}


STATE {
	Cai		(mM) <1e-8> 
}


INITIAL {
	Cai = Cainit
}


ASSIGNED {
	iCa		(mA/cm2)
	drive_channel	(mM/ms)
	drive_pump	(mM/ms)
}

	
BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state { 
	drive_channel =  - k * iCa
	if (drive_channel<=0.) { drive_channel = 0. }
	Cai' = drive_channel +(Cainf-Cai)/taur
}