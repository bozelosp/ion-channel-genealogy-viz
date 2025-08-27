NEURON {
	SUFFIX cacum
	USEION ca READ ica WRITE cai
	NONSPECIFIC_CURRENT i
	RANGE depth, tau, cai0, cmax
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
	irest = ica
	cmax=cai
}

BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
	if (cai>cmax) {cmax=cai}
	i=0
}

DERIVATIVE integrate {
	cai' = (irest-ica)/depth/F/2 * (1e4) + (cai0 - cai)/tau
}