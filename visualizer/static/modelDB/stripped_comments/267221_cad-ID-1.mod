NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	RANGE ca, ica, cai
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulomb)
}

PARAMETER {
	depth = .1 (um)		
	taur = 200 (ms)		
	cainf = 100e-6 (mM)
}

ASSIGNED {
	cai	(mM)
	ica	(mA/cm2)
	drive_channel (mM/ms)
}

STATE { ca (mM) }
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 
	drive_channel =  - (1e4) * ica / (2 * FARADAY * depth)
	
	if (drive_channel <= 0.) { drive_channel = 0. }		
	
	ca' = drive_channel/18 + (cainf - ca)/taur
	cai = ca
}

INITIAL { ca = cainf }