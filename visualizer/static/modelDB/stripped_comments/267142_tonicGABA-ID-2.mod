NEURON{
	SUFFIX tonicGABA
        NONSPECIFIC_CURRENT i
	RANGE  i, e, g
}

UNITS {
	(uS)  = (micromho)
	(nA)  = (nanoamp)
	(mV)  = (millivolt)
}


PARAMETER {
	g    = 0.0001 (siemens/cm2) <0, 1e9>
	e    = 0 (mV) 	
}

ASSIGNED{
	i (milliamp/cm2)
	v	(mV)		
}


BREAKPOINT {
	i = g * (v-e)
}