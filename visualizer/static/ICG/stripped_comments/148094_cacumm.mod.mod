NEURON {
	SUFFIX cacum
	USEION ca READ ica WRITE cai
	NONSPECIFIC_CURRENT i
	RANGE depth, tau, cai0, cmax, irest, ca_tmax
}

UNITS {
        (um) = (micron)
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
        ica      (mA/cm2)
        cmax     (milli/liter)
        ca_tmax  (ms)
        i        (mA/cm2)
}

STATE {
	cai (mM)
}

INITIAL {
	cai = cai0

	cmax=cai
	ca_tmax=0
}

BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
	if (cai>cmax) {cmax=cai ca_tmax=t}
	if (cai<0) {cai=0}
	i=0
}

DERIVATIVE integrate {
	cai' = -(ica+irest)/depth/F/2 * (1e4) + (cai0 - cai)/tau
}