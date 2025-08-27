INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Cad
	USEION Ca READ iCa, Cai WRITE Cai VALENCE 2
	RANGE depth,kt,kd,Cainf,taur,k,Cainf2,taur2
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
	depth	= .1	(um)		
	taur	= 700	(ms)		
	taur2	= 70	(ms)		
	Cainf	= 1e-8	(mM)
	Cainf2	= 5e-5	(mM)
	Cainit  = 5e-5
	kt	= 1	(mM/ms)		
	kd	= 5e-4	(mM)		
        k       = 1
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
	drive_channel =  - (k*10000) * iCa / (2 * FARADAY * depth)
	if (drive_channel<=0.) { drive_channel = 0. }

	drive_pump = -kt * Cai / (Cai + kd )		
	Cai' = drive_channel+drive_pump+(Cainf-Cai)/taur+(Cainf2-Cai)/taur2
}