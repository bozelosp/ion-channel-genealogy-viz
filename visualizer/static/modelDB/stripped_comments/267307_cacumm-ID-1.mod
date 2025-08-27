NEURON {
	SUFFIX cacum
	USEION ca READ ica WRITE cai
	NONSPECIFIC_CURRENT i
	RANGE depth, tau, cai0, cmax, gamma
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)
}

PARAMETER {
	depth = 0.1 (um)	
	irest = 0  (mA/cm2)		
	tau = 100 (ms)
	cai0 = 50e-6 (mM)	
			
			
			
	gamma = 1 
}

ASSIGNED {
	ica (mA/cm2)
	cmax
	i  	 (mA/cm2)
}

STATE {
	cai (mM)
}

INITIAL {
	cai = cai0
	
	cmax=cai
}

BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
	if (cai>cmax) {cmax=cai}
	i=0
}

DERIVATIVE integrate {
        
	cai' = (1e4) * (gamma*(irest-ica)/(2*F*depth)) + (cai0 - cai)/tau
}