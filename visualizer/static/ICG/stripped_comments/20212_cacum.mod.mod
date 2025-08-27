NEURON {
	SUFFIX cacum
	USEION ca READ ica WRITE cai
	RANGE depth, tau, cai0
}

UNITS {
	(mM) = (milli/liter)
	(mA) = (milliamp)
	F = (faraday) (coulombs)
}

PARAMETER {
	depth = 1 (nm)	
	tau = 10 (ms)
	cai0 = 50e-6 (mM)	
			
			
}

ASSIGNED {
	ica (mA/cm2)
}

STATE {
	cai (mM)
}

BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
}

DERIVATIVE integrate {
	cai' = -ica/depth/F * (1e7) + (cai0 - cai)/tau
}