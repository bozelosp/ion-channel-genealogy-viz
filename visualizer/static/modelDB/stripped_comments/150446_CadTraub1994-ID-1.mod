NEURON {
	SUFFIX Cad
	USEION ca READ ica WRITE cai
	RANGE  phi, beta
}

UNITS {
	(mA)	= (milliamp)
 	(molar) = (1/liter)			
 	(mM)	= (millimolar)
}

PARAMETER {
	phi	 	
	beta 	(/ms)
}

STATE {	cai (mM) }

INITIAL { 
	cai = 0
}

ASSIGNED { 
	ica		(mA/cm2) 
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	if( cai < 0 ){ cai = 0 }
}

DERIVATIVE state { 
	cai' = - phi * ica - beta * (cai)
}