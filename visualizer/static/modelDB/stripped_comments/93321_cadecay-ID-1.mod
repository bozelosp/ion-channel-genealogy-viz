INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms) }

NEURON{
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	RANGE ica, channel_flow, depth, B
	GLOBAL tau, cainf, setB
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um) = (micron)
}

CONSTANT {
        FARADAY=96485.309   
}

PARAMETER {
	dt (ms)
        depth = 1 (um)       
	tau = 20 	(ms) 
	cainf = 5e-5	(mM) 
	ica		(mA/cm2)
        setB = -4.7e-2 (cm2 mM/mA/ms)
}

STATE {
	cai		(mM)
}

INITIAL {
	cai = cainf
        B = -4.7e-2 
}

ASSIGNED {
	channel_flow	(mM/ms)
	B		(mM cm2/ms/mA)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {


        B = setB
	channel_flow = B*ica
	if (channel_flow <= 0.0 ) { channel_flow = 0.0 }	
	cai' = channel_flow  - (cai - cainf)/tau


}