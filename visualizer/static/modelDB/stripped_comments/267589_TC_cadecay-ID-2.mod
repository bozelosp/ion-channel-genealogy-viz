NEURON {
	SUFFIX TC_cad
	USEION ca READ ica WRITE cai
	RANGE depth,kt,kd,cainf,taur, cai_rec, gamma
}

UNITS {

	(mM)	= (milli/liter)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
	FARADAY = (faraday) (coulombs)
}

PARAMETER {
	depth	= .1	 (um)		
	gamma   = 0.05 	 (1)		
	taur	= 5	 (ms)		
	cainf	= 5e-5 (mM)		
}

STATE {
	cai		(mM) 
}

INITIAL {
	cai = cainf
}

ASSIGNED {
	ica		(mA/cm2)
	cai_rec		(mM)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
	
}

DERIVATIVE state { 

	cai' = -(10000)*(ica*gamma/(2*FARADAY*depth)) - (cai - cainf)/taur
	cai_rec = cai

}