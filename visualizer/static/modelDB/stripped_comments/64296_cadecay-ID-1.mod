INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms) }

NEURON{
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	RANGE ica, channel_flow, depth, B
	GLOBAL tau, cainf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
}

CONSTANT {
      FARADAY = 96154 (coul)
	
						
						
}

PARAMETER {
	dt (ms)
	depth = 1 	(um)		
					
					
	tau = 1e-5 	(ms)		
	cainf = 1e-5	(mM)		
	ica		(mA/cm2)
}

STATE {
	cai		(mM)
}

INITIAL {
	cai = cainf
}

ASSIGNED {
	channel_flow	(mM/ms)
	B		(mM cm2/ms/mA)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	B = -(1e4)/(2*FARADAY*depth)
	channel_flow = B*ica
	if (channel_flow <= 0.0 ) { channel_flow = 0.0 }	
	cai' = channel_flow  - (cai - cainf)/tau
}