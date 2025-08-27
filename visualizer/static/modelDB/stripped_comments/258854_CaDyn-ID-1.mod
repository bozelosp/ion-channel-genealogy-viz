NEURON {
	SUFFIX cad
	USEION ca READ ica WRITE cai
	RANGE  phi, tau
	GLOBAL ceiling
}

UNITS {
	(mA)	= (milliamp)
	(mM) 	= (milli/liter)
}

PARAMETER {
	phi		= 1e-3 			(100/coulomb meter)
	tau		= 50 			(ms)
	ceiling	= 1e6			(mM)
}

STATE {	
	cai (mM)
	}

INITIAL { 
	cai = 0 
}

ASSIGNED { 
	ica		(mA/cm2) 
}
	
BREAKPOINT {
	SOLVE state METHOD cnexp
	if( cai < 0 ){ cai = 0 }
	if( cai > ceiling ){ cai = ceiling }
}

DERIVATIVE state { 

	cai' = - phi * ica - (cai) / tau
}